//---------------------------------------------------------
// Copyright 2023 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
use rust_htslib::{bam, faidx, bam::Read};
use rust_htslib::bcf::{Writer as BcfWriter, Header as BcfHeader, Format as BcfFormat};
use bio::stats::{Prob, LogProb, PHREDProb};
use fishers_exact::fishers_exact;

use crate::pileup_stats::*;
use crate::calling_models::*;
use crate::utility::*;
use crate::classifier::*;

use crate::{ReadMetadata, ReadHaplotypeLikelihood};
use crate::HashMap;
use crate::LongshotParameters;
use crate::ModelParameters;

use longshot::extract_fragments::extract_fragments;
use longshot::util::parse_region_string;
use longshot::call_potential_snvs::call_potential_snvs;

pub fn somatic_call(input_bam: &str, region_str: &str, reference_genome: &str, model: CallingModel) {

    let params = ModelParameters::defaults();

    let max_softclip = 500;
    let max_num_softclip_reads = 5;
    let soft_min_p_somatic = 1e-6;
    let hard_min_p_somatic = 1e-10;
    let min_allele_call_qual = LogProb::from(Prob(0.1));
    let hard_min_depth = 10;
    let filter_min_depth = 10;

    let max_variant_minor_observations = 2;
    let max_strand_bias = 20.0;
    let min_variant_observations = 3;
    let min_variant_observations_per_strand = 1;
    let mut longshot_parameters = LongshotParameters::defaults();
    longshot_parameters.estimate_alignment_parameters(&input_bam.to_owned(), &reference_genome.to_owned(), &None);

    let region = parse_region_string(Some(region_str), &input_bam.to_owned()).expect("Could not parse region").unwrap();
    
    let mut bam = bam::IndexedReader::from_path(input_bam).unwrap();
    let header = bam::Header::from_template(bam.header());
    let header_view = bam::HeaderView::from_header(&header);
    
    // set up header and faidx for pulling reference sequence
    let faidx = faidx::Reader::from_path(reference_genome).expect("Could not read reference genome:");
    let chromosome_length = header_view.target_len(region.tid).unwrap() as usize;
    let mut chromosome_sequence = faidx.fetch_seq_string(&region.chrom, 0, chromosome_length).unwrap();
    chromosome_sequence.make_ascii_uppercase();
    let chromosome_bytes = chromosome_sequence.as_bytes();
    
    // set up vcf output
    let mut vcf_header = BcfHeader::new();
    vcf_header.push_record(b"##source=smrest");
    vcf_header.push_record(format!("##model={}", calling_model_to_str(model)).as_bytes());

    // add contig lines to header
    for tid in 0..header_view.target_count() {
        let l = format!("##contig=<ID={},length={}", std::str::from_utf8(header_view.tid2name(tid)).unwrap(), 
                                                     header_view.target_len(tid).unwrap());
        vcf_header.push_record(l.as_bytes());
    }

    // add info lines
    vcf_header.push_record(r#"##INFO=<ID=SomaticHaplotypeIndex,Number=1,Type=Integer,Description="Index of haplotype carrying the mutation">"#.as_bytes());
    vcf_header.push_record(r#"##INFO=<ID=SoftClippedEvidenceReads,Number=1,Type=Integer,Description="Number of evidence reads that have long softclips">"#.as_bytes());
    vcf_header.push_record(r#"##INFO=<ID=HaplotypeAltCount,Number=2,Type=Integer,Description="Observed alt read count on each haplotype">"#.as_bytes());
    vcf_header.push_record(r#"##INFO=<ID=HaplotypeDepth,Number=2,Type=Integer,Description="Observed read count on each haplotype">"#.as_bytes());
    vcf_header.push_record(r#"##INFO=<ID=HaplotypeVAF,Number=2,Type=Float,Description="VAF on each haplotype">"#.as_bytes());
    vcf_header.push_record(r#"##INFO=<ID=MutationType,Number=1,Type=String,Description="Mutation in canonical form">"#.as_bytes());
    vcf_header.push_record(r#"##INFO=<ID=SequenceContext,Number=1,Type=String,Description="Sequence context on reference genome">"#.as_bytes());
    vcf_header.push_record(r#"##INFO=<ID=StrandBias,Number=1,Type=Float,Description="Strand bias p-value (PHRED scaled)">"#.as_bytes());
    vcf_header.push_record(r#"##INFO=<ID=ProportionPhased,Number=1,Type=Float,Description="Proportion of reads at this postion successfully phased">"#.as_bytes());
    vcf_header.push_record(r#"##INFO=<ID=ContextErrorRates,Number=.,Type=String,Description="Error rates at ref/alt sequence context">"#.as_bytes());
    vcf_header.push_record(r#"##INFO=<ID=PhaseSets,Number=.,Type=Integer,Description="Phase sets seen in reads">"#.as_bytes());
    vcf_header.push_record(br#"##FILTER=<ID=MinObsPerStrand,Description="todo">"#);
    vcf_header.push_record(br#"##FILTER=<ID=LowQual,Description="todo">"#);
    vcf_header.push_record(br#"##FILTER=<ID=MaxOtherHaplotypeObservations,Description="todo">"#);
    vcf_header.push_record(br#"##FILTER=<ID=StrandBias,Description="todo">"#);
    vcf_header.push_record(br#"##FILTER=<ID=ExcessiveSoftClips,Description="todo">"#);
    vcf_header.push_record(br#"##FILTER=<ID=MinHaplotypeDepth,Description="todo">"#);

    //vcf_header.push_sample("sample".as_bytes());
    let mut vcf = BcfWriter::from_stdout(&vcf_header, true, BcfFormat::Vcf).unwrap();

    let mut varlist = call_potential_snvs(
            &input_bam.to_owned(),
            &reference_genome.to_owned(),
            &Some(region.clone()),
            &longshot_parameters.genotype_priors,
            20,
            200,
            min_variant_observations, //potential_snv_min_alt_count,
            0.1, //potential_snv_min_alt_frac,
            50,
            longshot_parameters.alignment_parameters.as_ref().unwrap().ln(),
            LogProb::from(Prob(0.000000000001)),
            ).expect("Could not find candidate variants");

    for var in &mut varlist.lst {
        //println!("{}\t{}\t{}\t{}\t{}", var.tid, var.pos0 + 1, var.alleles.len(), var.alleles[0], var.alleles[1]);
        if var.alleles.len() > 2 {
            // trim, TODO: handle multi-allelic better
            var.alleles.truncate(2);
            var.allele_counts.truncate(2);
            var.allele_counts_forward.truncate(2);
            var.allele_counts_reverse.truncate(2);
        }
    }

    let frags = extract_fragments(&input_bam.to_owned(),
                                  &reference_genome.to_owned(),
                                  &mut varlist,
                                  &Some(region.clone()),
                                  &longshot_parameters.extract_fragment_parameters,
                                  &longshot_parameters.alignment_parameters.as_ref().unwrap()).unwrap();
    
    // Populate haplotype, strand information and another other information we need from the bam
    // TODO: handle multiple alignments per read, not just primary
    bam.fetch( (region.tid, region.start_pos, region.end_pos + 1) ).unwrap();
    let mut read_meta = HashMap::<String, ReadMetadata>::new();
    for r in bam.records() {
        let record = r.unwrap();
        
        // same criteria as longshot
        if record.is_quality_check_failed()
            || record.is_duplicate()
            || record.is_secondary()
            || record.is_unmapped()
            || record.is_supplementary()
        {
            continue;
        }

        let s = String::from_utf8(record.qname().to_vec()).unwrap();
        let rm = ReadMetadata { 
            haplotype_index: get_haplotag_from_record(&record),
            phase_set: get_phase_set_from_record(&record),
            strand_index: record.is_reverse() as i32,
            leading_softclips: record.cigar().leading_softclips(),
            trailing_softclips: record.cigar().trailing_softclips()
        };
        read_meta.insert(s.clone(), rm);
    }

    // Convert fragments into read-haplotype likelihoods for every variant
    let mut rhl_per_var = vec![ Vec::<ReadHaplotypeLikelihood>::new(); varlist.lst.len() ];

    for f in frags {
        let id = f.id.unwrap();

        if let Some( rm ) = read_meta.get(&id) {
            for c in f.calls {
                let var = &varlist.lst[c.var_ix];
                //println!("{}\t{}\t{}\t{}\t{}\t{}", &id, var.tid, var.pos0 + 1, var.alleles[0], var.alleles[1], c.allele);
            
                let rhl = ReadHaplotypeLikelihood { 
                    read_name: Some(id.clone()),
                    mutant_allele_likelihood: c.allele_scores[1],
                    base_allele_likelihood: c.allele_scores[0],
                    allele_call: var.alleles[c.allele as usize].clone(),
                    allele_call_qual: c.qual,
                    haplotype_index: rm.haplotype_index,
                    strand_index: rm.strand_index
                };
                rhl_per_var[c.var_ix].push(rhl);
            }
        }
    }

    let cm = &longshot_parameters.extract_fragment_parameters.context_model;
    let context_probs = &longshot_parameters.alignment_parameters.as_ref().unwrap().context_emission_probs.probs;
    let half_k = cm.k / 2;

    for i in 0..varlist.lst.len() {
        let var = &varlist.lst[i];
        let rhls = &rhl_per_var[i];

        if var.pos0 < region.start_pos as usize || var.pos0 > region.end_pos as usize {
            continue;
        }

        let reference_base = var.alleles[0].as_bytes()[0] as char;
        let alt_base = var.alleles[1].as_bytes()[0] as char;

        if base2index(reference_base) == -1 || base2index(alt_base) == -1 {
            continue;
        }

        let mut ps = PileupStats::new();
        ps.fill_from_rhls(min_allele_call_qual, rhls);
        
        let reference_base_index = base2index(reference_base) as u32;
        let candidate_variant_index = base2index(alt_base) as u32;

        let h0_ac = ps.get_count_on_haplotype(candidate_variant_index, 0) as i32;
        let h1_ac = ps.get_count_on_haplotype(candidate_variant_index, 1) as i32;
        let candidate_haplotype_index = if h0_ac > h1_ac { 0 } else { 1 };

        if ps.get_haplotype_depth(candidate_haplotype_index as u32) < hard_min_depth {
            continue;
        }

        // +1 here to convert to be the same as the bam HP tag
        let rhls_by_hap = rhls.iter().filter(|x| x.haplotype_index.unwrap_or(-1) == (candidate_haplotype_index + 1)).cloned().collect();
        let class_probs = calculate_class_probabilities_likelihood(&rhls_by_hap, &params);

        // do not output sites below this threshold
        if class_probs[2] < hard_min_p_somatic {
            continue;
        }

        let mut record = vcf.empty_record();
        let rid = vcf.header().name2rid(region.chrom.as_bytes()).expect("Could not find reference id");
        record.set_rid(Some(rid));
        
        record.set_pos(var.pos0 as i64); 

        record.set_alleles(&[ var.alleles[0].as_bytes(), var.alleles[1].as_bytes() ]).expect("Could not set alleles");

        // require a minimum number of high-quality observations of the variant per strand
        if ps.get(candidate_variant_index, candidate_haplotype_index as u32, 0) < min_variant_observations_per_strand ||
           ps.get(candidate_variant_index, candidate_haplotype_index as u32, 1) < min_variant_observations_per_strand {
            record.push_filter("MinObsPerStrand".as_bytes()).unwrap();
        }

        // Check for the evidence reads having excessively long softclips, which can indicate alignment artifacts
        let mut num_softclipped_evidence_reads = 0;
        let mut phase_sets = Vec::new();

        for rhl in rhls_by_hap {
            if rhl.allele_call != var.alleles[1] {
                continue;
            }

            // lookup read metadata and check softclip lengths
            let id = rhl.read_name.unwrap();
            if let Some( rm ) = read_meta.get(&id) {
                if rm.leading_softclips >= max_softclip || rm.trailing_softclips >= max_softclip {
                    num_softclipped_evidence_reads += 1
                }

                if let Some( phase_set) = rm.phase_set {
                    phase_sets.push(phase_set);
                }
            }
        }

        record.push_info_integer(b"SoftClippedEvidenceReads", &[num_softclipped_evidence_reads as i32]).expect("Could not add INFO");
        if num_softclipped_evidence_reads > max_num_softclip_reads {
            record.push_filter("ExcessiveSoftClips".as_bytes()).unwrap();
        }

        if phase_sets.len() > 0 {
            phase_sets.sort_unstable();
            phase_sets.dedup();
            record.push_info_integer(b"PhaseSets", &phase_sets).expect("Could not add INFO");
        }
        
        if class_probs[2] < soft_min_p_somatic {
            record.push_filter("LowQual".as_bytes()).unwrap();
        }

        let mut qual = *PHREDProb::from(Prob(1.0 - class_probs[2]));
        if qual > 500.0 {
            qual = 500.0;
        }
        record.set_qual(qual as f32);
        record.push_info_integer(b"SomaticHaplotypeIndex", &[candidate_haplotype_index as i32]).expect("Could not add INFO");

        let h0_depth = ps.get_haplotype_depth(0) as i32;
        let h1_depth = ps.get_haplotype_depth(1) as i32;
        
        let ho_ac = ps.get_count_on_haplotype(candidate_variant_index, (1 - candidate_haplotype_index) as u32) as i32;
        if ho_ac > max_variant_minor_observations {
            record.push_filter("MaxOtherHaplotypeObservations".as_bytes()).unwrap();
        }

        record.push_info_integer(b"HaplotypeAltCount", &[h0_ac, h1_ac]).expect("Could not add INFO");
        record.push_info_integer(b"HaplotypeDepth", &[h0_depth, h1_depth]).expect("Could not add INFO");

        let h0_vaf = h0_ac as f32 / h0_depth as f32;
        let h1_vaf = h1_ac as f32 / h1_depth as f32;
        record.push_info_float(b"HaplotypeVAF", &[h0_vaf, h1_vaf]).expect("Could not add INFO");
        

        if h0_depth < filter_min_depth || h1_depth < filter_min_depth {
            record.push_filter("MinHaplotypeDepth".as_bytes()).unwrap();
        }

        // grab reference context
        let reference_3mer_context = &chromosome_bytes[ (var.pos0 as usize - 1)..(var.pos0 as usize + 2)];
        
        // grab mutation type
        let mut mutation_type: [char; 3] = [ reference_base, '>', alt_base ];

        // convert mutation type/context to canonical form C>x, T>x
        if mutation_type[0] != 'C' && mutation_type[0] != 'T' {
            mutation_type[0] = bio::alphabets::dna::complement(mutation_type[0] as u8) as char;
            mutation_type[2] = bio::alphabets::dna::complement(mutation_type[2] as u8) as char;
        }

        let mutation_type_str = String::from_iter(&mutation_type);
        record.push_info_string(b"MutationType", &[mutation_type_str.as_bytes()]).expect("Could not add INFO");
        record.push_info_string(b"SequenceContext", &[&reference_3mer_context]).expect("Could not add INFO");
        
        //let mean_mapq = ps.mean_mapq;
        record.push_info_float(b"ProportionPhased", &[ps.proportion_phased]).expect("Could not add INFO");

        let ro0 = ps.get( reference_base_index, candidate_haplotype_index as u32, 0);
        let ro1 = ps.get( reference_base_index, candidate_haplotype_index as u32, 1);
        let ao0 = ps.get( candidate_variant_index, candidate_haplotype_index as u32, 0);
        let ao1 = ps.get( candidate_variant_index, candidate_haplotype_index as u32, 1);
        //println!("[ {} {} {} {} ]", ro0, ro1, ao0, ao1);

        let fishers_result = 
            fishers_exact(&[ ro0, ro1,
                             ao0, ao1 ]).unwrap();

        // straight from longshot
        let strand_bias_pvalue = if fishers_result.two_tail_pvalue <= 500.0 {
            *PHREDProb::from(Prob(fishers_result.two_tail_pvalue))
        } else {
            500.0
        };
        record.push_info_float(b"StrandBias", &[strand_bias_pvalue as f32]).expect("Could not add INFO");
        
        if strand_bias_pvalue > max_strand_bias {
            record.push_filter("StrandBias".as_bytes()).unwrap();
        }

        // annotate with the sequence context model error rates
        
        let model_ref_context = &chromosome_bytes[ (var.pos0 as usize - half_k)..(var.pos0 as usize + half_k + 1)];
        if cm.alphabet.is_word(model_ref_context) {
            let mut model_alt_context = vec![0; model_ref_context.len()];
            model_alt_context.clone_from_slice(model_ref_context);
            model_alt_context[half_k] = alt_base as u8;

            let ref_base_rank = cm.get_base_rank(reference_base) as usize;
            let alt_base_rank = cm.get_base_rank(alt_base) as usize;

            let model_ref_er_fwd = context_probs[ cm.get_context_rank(model_ref_context, false).unwrap() ][alt_base_rank];
            let model_ref_er_ref = context_probs[ cm.get_context_rank(model_ref_context, true).unwrap() ][alt_base_rank];
            
            let model_alt_er_fwd = context_probs[ cm.get_context_rank(&model_alt_context, false).unwrap() ][ref_base_rank];
            let model_alt_er_ref = context_probs[ cm.get_context_rank(&model_alt_context, true).unwrap() ][ref_base_rank];
            let context_ref_fwd_er_str = format!("{}+>{}:{:.3}", std::str::from_utf8(model_ref_context).unwrap(), alt_base, model_ref_er_fwd);
            let context_ref_rev_er_str = format!("{}->{}:{:.3}", std::str::from_utf8(model_ref_context).unwrap(), alt_base, model_ref_er_ref);
            let context_alt_fwd_er_str = format!("{}+>{}:{:.3}", std::str::from_utf8(&model_alt_context).unwrap(), reference_base, model_alt_er_fwd);
            let context_alt_rev_er_str = format!("{}->{}:{:.3}", std::str::from_utf8(&model_alt_context).unwrap(), reference_base, model_alt_er_ref);

            record.push_info_string(b"ContextErrorRates", &[&context_ref_fwd_er_str.as_bytes(),
                                                            &context_ref_rev_er_str.as_bytes(),
                                                            &context_alt_fwd_er_str.as_bytes(),
                                                            &context_alt_rev_er_str.as_bytes() ]).expect("Could not add INFO");
        }

        vcf.write(&record).unwrap();
    }
}
