//---------------------------------------------------------
// Copyright 2023 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
use rust_htslib::{bam, bam::Read};
use rust_htslib::bcf::{Writer as BcfWriter, Header as BcfHeader, Format as BcfFormat};
use bio::stats::{Prob, LogProb, PHREDProb};
use fishers_exact::fishers_exact;
use std::collections::HashMap;
use std::ops::Range;
use std::fs::File;
use std::io::Write;

use crate::pileup_stats::*;
use crate::calling_models::*;
use crate::utility::*;
use crate::classifier::*;
use crate::ReadMetadata;

use crate::LongshotParameters;
use crate::ModelParameters;

use longshot::extract_fragments::extract_fragments;
use longshot::util::parse_region_string;
use longshot::util::GenomicInterval;
use longshot::genotype_probs::{Genotype, GenotypeProbs};
use longshot::variants_and_fragments::{Var, VarList, VarFilter};
use longshot::realignment::AlignmentType;

pub fn find_candidates(bam_file: &String,
                       interval: &GenomicInterval,
                       chromosome_bytes: &Vec<u8>,
                       read_metadata: &HashMap::<String, ReadMetadata>,
                       min_haplotype_depth: u32,
                       max_total_depth: u32,
                       min_alt_count: usize,
                       min_alt_frac: f64,
                       _min_mapq: u8) -> (VarList, Vec<Range<usize>>)
{
    let mut bam = bam::IndexedReader::from_path(bam_file).unwrap();

    let mut ps = PileupStats::new();
    
    let mut var_vec = vec![];
    let mut callable_region_vec = vec![];

    // go to calling region
    bam.fetch( (interval.tid as u32, interval.start_pos, interval.end_pos + 1) ).unwrap();

    for p in bam.pileup() {
        let pileup = p.unwrap();
        
        // check whether we should bother with this position
        let reference_position = pileup.pos() as usize;

        if reference_position < interval.start_pos as usize || reference_position > interval.end_pos as usize {
            continue;
        }

        let reference_base = chromosome_bytes[reference_position] as char;
        if reference_base == 'N' {
            continue;
        }

        ps.fill_pileup(&read_metadata, pileup.alignments());
        if ps.get_haplotype_depth(0) < min_haplotype_depth || ps.get_haplotype_depth(1) < min_haplotype_depth {
            continue;
        }

        if ps.get_haplotype_depth(0) + ps.get_haplotype_depth(1) > max_total_depth {
            continue;
        }

        // This position passes the checks for a callable region, manage all the interval wranging
        if callable_region_vec.is_empty() {
            callable_region_vec.push( Range { start: reference_position, end: reference_position } );
        } else {
            let mut current = callable_region_vec.last_mut().unwrap();
            if current.end + 1 != reference_position {
                callable_region_vec.push( Range { start: reference_position, end: reference_position } );
            } else {
                // extend previous interval
                current.end = reference_position;
            }
        }

        let reference_base_index = base2index(reference_base) as u32;

        // Calculate most frequently observed non-reference base on either haplotype
        let mut max_variant_count = 0;
        let mut max_variant_index = 0;
        let mut max_haplotype_index = 0;

        for hi in 0u32..2u32 {
            for bi in 0u32..4u32 {
                let b0 = ps.get(bi, hi, 0);
                let b1 = ps.get(bi, hi, 1);
                if bi != reference_base_index && (b0 + b1) > max_variant_count {
                    max_variant_count = b0 + b1;
                    max_variant_index = bi;
                    max_haplotype_index = hi;
                }
            }
        }

        let alt_frac = max_variant_count as f64 / ps.get_haplotype_depth(max_haplotype_index) as f64;
        if alt_frac >= min_alt_frac && max_variant_count >= min_alt_count as u32 {
            let bases = "ACGT";
            let variant_base = bases.as_bytes()[max_variant_index as usize] as char;

            let mut alleles = vec![];
            let ref_allele = std::str::from_utf8(&[reference_base as u8]).unwrap().to_owned();
            let alt_allele = std::str::from_utf8(&[variant_base as u8]).unwrap().to_owned();
            alleles.push(ref_allele);
            alleles.push(alt_allele);

            let v = Var {
                ix: 0,
                tid: interval.tid,
                pos0: reference_position,
                alleles: alleles.clone(),
                dp: 0,
                allele_counts: vec![0; alleles.len()],
                allele_counts_forward: vec![0; alleles.len()],
                allele_counts_reverse: vec![0; alleles.len()],
                ambiguous_count: 0,
                qual: 0.0,
                filter: VarFilter::Pass,
                genotype: Genotype(0, 0),
                //unphased: false,
                gq: 0.0,
                unphased_genotype: Genotype(0, 0),
                unphased_gq: 0.0,
                genotype_post: GenotypeProbs::uniform(alleles.len()),
                phase_set: None,
                strand_bias_pvalue: 1.0,
                mec: 0,
                mec_frac_variant: 0.0, // mec fraction for this variant
                mec_frac_block: 0.0,   // mec fraction for this haplotype block
                mean_allele_qual: 0.0,
                dp_any_mq: 0,
                mq10_frac: 0.0,
                mq20_frac: 0.0,
                mq30_frac: 0.0,
                mq40_frac: 0.0,
                mq50_frac: 0.0
            };
            var_vec.push(v);
        }
    }
    
    let target_names = bam.header().target_names().iter().map(|x| std::str::from_utf8(x).unwrap().to_owned()).collect();
    let vlst = VarList::new(var_vec, target_names).unwrap();
    vlst.assert_sorted();

    return (vlst, callable_region_vec);
}

pub fn somatic_call(input_bam: &str, 
                    region_str: &str, 
                    reference_genome: &str, 
                    output_bed: Option<&str>, 
                    phased_vcf: Option<&str>,
                    model: CallingModel) {

    // somatic calling specific parameters
    let max_softclip = 500;
    let max_tail_mismatch_rate = 0.05;
    let max_alignment_artifact_evidence_proportion = 0.25;

    let soft_min_p_somatic = 1e-6;
    let hard_min_p_somatic = 1e-10;
    let min_allele_call_qual = LogProb::from(Prob(0.1));
    let hard_min_depth = 10;
    let filter_min_depth = hard_min_depth as i32;

    let max_variant_minor_observations = 2;
    let max_variant_minor_hvaf = 0.2;
    let max_strand_bias = 20.0;
    let min_variant_observations = 3;
    let min_variant_observations_per_strand = 1;
    
    let params = ModelParameters::defaults();
    let mut longshot_parameters = LongshotParameters::defaults();
    longshot_parameters.estimate_alignment_parameters(&input_bam.to_owned(), &reference_genome.to_owned(), &None);
    longshot_parameters.extract_fragment_parameters.alignment_type = AlignmentType::HtsLibProbAln;

    let region = parse_region_string(Some(region_str), &input_bam.to_owned()).expect("Could not parse region").unwrap();
    
    let mut bam = bam::IndexedReader::from_path(input_bam).unwrap();
    let header_view = bam.header();
    
    // get entire reference chromosome
    let chromosome_bytes = get_chromosome_sequence(reference_genome, header_view, region.tid).into_bytes();
    
    // set up vcf output
    let mut vcf_header = BcfHeader::new();
    vcf_header.push_record(b"##source=smrest");
    vcf_header.push_record(format!("##model={}", calling_model_to_str(model)).as_bytes());
    add_contig_lines_to_vcf(&mut vcf_header, header_view);

    // add info lines
    vcf_header.push_record(r#"##INFO=<ID=SomaticHaplotypeIndex,Number=1,Type=Integer,Description="Index of haplotype carrying the mutation">"#.as_bytes());
    vcf_header.push_record(r#"##INFO=<ID=AlignmentArtifactEvidenceReads,Number=1,Type=Integer,Description="Number of evidence reads that have long softclips or excessive mismatches">"#.as_bytes());
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
    vcf_header.push_record(br#"##FILTER=<ID=PossibleAlignmentArtifact,Description="todo">"#);
    vcf_header.push_record(br#"##FILTER=<ID=MinHaplotypeDepth,Description="todo">"#);

    //vcf_header.push_sample("sample".as_bytes());
    let mut vcf = BcfWriter::from_stdout(&vcf_header, true, BcfFormat::Vcf).unwrap();

    // Populate haplotype, strand information and another other information we need from the bam
    let read_metadata = populate_read_metadata_from_bam(&input_bam.to_owned(), 
                                                        &region, 
                                                        Some(&chromosome_bytes), 
                                                        phased_vcf, 
                                                        Some(&longshot_parameters), 
                                                        Some(&reference_genome.to_owned()));


    //
    let (mut varlist, callable_regions) = find_candidates(&input_bam.to_owned(),
                                                          &region,
                                                          &chromosome_bytes,
                                                          &read_metadata,
                                                          hard_min_depth,
                                                          400,
                                                          min_variant_observations,
                                                          0.1,
                                                          60);
/*
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
*/
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
    

    // Convert fragments into read-haplotype likelihoods for every variant
    let rhl_per_var = fragments_to_read_haplotype_likelihoods(&varlist, &frags, &read_metadata);

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
        let mut num_alignment_artifact_evidence_reads = 0;
        let mut num_evidence_reads_with_metadata = 0;
        let mut phase_sets = Vec::new();

        for rhl in rhls_by_hap {
            if rhl.allele_call != var.alleles[1] {
                continue;
            }

            // lookup read metadata and check softclip lengths
            let id = rhl.read_name.unwrap();
            if let Some( rm ) = read_metadata.get(&id) {
                num_evidence_reads_with_metadata += 1;
                let has_long_softclip = rm.leading_softclips >= max_softclip || rm.trailing_softclips >= max_softclip;
                let has_suspect_alignment_end = rm.prefix_mismatch_rate.unwrap_or(0.0) >= max_tail_mismatch_rate || 
                                                rm.suffix_mismatch_rate.unwrap_or(0.0) >= max_tail_mismatch_rate;

                if has_long_softclip || has_suspect_alignment_end {
                    num_alignment_artifact_evidence_reads += 1
                }

                if let Some( phase_set) = rm.phase_set {
                    phase_sets.push(phase_set);
                }
            }
        }

        let artifact_evidence_proportion = num_alignment_artifact_evidence_reads as f32 / num_evidence_reads_with_metadata as f32;
        record.push_info_integer(b"AlignmentArtifactEvidenceReads", &[num_alignment_artifact_evidence_reads as i32]).expect("Could not add INFO");
        if artifact_evidence_proportion > max_alignment_artifact_evidence_proportion {
            record.push_filter("PossibleAlignmentArtifact".as_bytes()).unwrap();
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

        record.push_info_integer(b"HaplotypeAltCount", &[h0_ac, h1_ac]).expect("Could not add INFO");
        record.push_info_integer(b"HaplotypeDepth", &[h0_depth, h1_depth]).expect("Could not add INFO");

        let h0_vaf = h0_ac as f32 / h0_depth as f32;
        let h1_vaf = h1_ac as f32 / h1_depth as f32;
        let ho_vaf = if candidate_haplotype_index == 0 { h1_vaf } else { h0_vaf };

        record.push_info_float(b"HaplotypeVAF", &[h0_vaf, h1_vaf]).expect("Could not add INFO");
        
        if ho_ac > max_variant_minor_observations || ho_vaf > max_variant_minor_hvaf {
            record.push_filter("MaxOtherHaplotypeObservations".as_bytes()).unwrap();
        }

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

    // write callable regions
    if let Some(filename) = output_bed {
        let mut file = File::create(filename).expect("Unable to open file for write");
        for r in callable_regions {
            let out = format!("{}\t{}\t{}\n", region.chrom, r.start, r.end + 1);
            file.write_all(out.as_bytes()).expect("Unable to write to file");
        }
    }
}

