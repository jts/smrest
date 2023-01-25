//---------------------------------------------------------
// Copyright 2023 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
use rust_htslib::{bam, bam::Read};
use rust_htslib::bcf::{Writer as BcfWriter, Header as BcfHeader, Format as BcfFormat, record::GenotypeAllele};
use bio::stats::{Prob, LogProb, PHREDProb};
use fishers_exact::fishers_exact;
use crate::pileup_stats::*;
use crate::utility::*;
use crate::binomial_test_onesided_greater;
use crate::LongshotParameters;
use longshot::extract_fragments::extract_fragments;
use longshot::variants_and_fragments::parse_vcf_potential_variants;
use longshot::variants_and_fragments::VarList;
use longshot::util::parse_region_string;

pub fn genotype(input_bam: &str, region_str: &str, candidates_vcf: &str, reference_genome: &str) {

    // subprogram thresholds
    let min_allele_call_qual = LogProb::from(Prob(0.1));
    let hard_min_depth = 20;
    let min_informative_fraction = 0.75;
    let max_strand_bias = 10.0;
    
    let gt_hom_alt = [GenotypeAllele::Unphased(1), GenotypeAllele::Unphased(1)];
    let gt_hom_ref = [GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(0)];
    let gt_het = [GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(1)];
    let gt_nocall = [GenotypeAllele::UnphasedMissing, GenotypeAllele::UnphasedMissing];
    let mut longshot_parameters = LongshotParameters::defaults();
    longshot_parameters.estimate_alignment_parameters(&input_bam.to_owned(), &reference_genome.to_owned(), &None);
    
    let region = parse_region_string(Some(region_str), &input_bam.to_owned()).expect("Could not parse region").unwrap();
    
    let mut bam = bam::IndexedReader::from_path(input_bam).unwrap();
    let header_view = bam.header();
    
    // set up vcf output
    let mut vcf_header = BcfHeader::new();
    vcf_header.push_record(b"##source=smrest genotype-hets");
    add_contig_lines_to_vcf(&mut vcf_header, header_view);

    // add info lines
    vcf_header.push_record(r#"##INFO=<ID=TotalDepth,Number=1,Type=Integer,Description="Total number of reads aligned across this position">"#.as_bytes());
    vcf_header.push_record(r#"##INFO=<ID=LP,Number=2,Type=Float,Description="todo">"#.as_bytes());
    vcf_header.push_record(r#"##INFO=<ID=EP,Number=1,Type=Float,Description="todo">"#.as_bytes());
    vcf_header.push_record(r#"##INFO=<ID=StrandBias,Number=1,Type=Float,Description="Strand bias p-value (PHRED scaled)">"#.as_bytes());
    vcf_header.push_record(r#"##INFO=<ID=InformativeFraction,Number=1,Type=Float,Description="Proportion of reads that could be used in allele calling">"#.as_bytes());
    vcf_header.push_record(r#"##FORMAT=<ID=GT,Number=1,Type=String,Description="sample genotype">"#.as_bytes());
    vcf_header.push_record(r#"##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths (high-quality bases)">"#.as_bytes());
    vcf_header.push_sample("sample".as_bytes());

    let mut vcf = BcfWriter::from_stdout(&vcf_header, true, BcfFormat::Vcf).unwrap();
    
    let all_vars = parse_vcf_potential_variants(&candidates_vcf.to_owned(), &input_bam.to_owned())
                        .expect("Could not parse candidate variants file");
    
    let mut varlist = VarList::new(all_vars.get_variants_range(region.clone()).expect("Could not subset vars"), all_vars.target_names)
                                  .expect("Could not construct var list");

    let frags = extract_fragments(&input_bam.to_owned(),
                                  &reference_genome.to_owned(),
                                  &mut varlist,
                                  &Some(region.clone()),
                                  &longshot_parameters.extract_fragment_parameters,
                                  &longshot_parameters.alignment_parameters.as_ref().unwrap()).unwrap();
    
    // populate read metadata
    let read_meta = populate_read_metadata_from_bam(&mut bam, &region);

    // Convert fragments into read-haplotype likelihoods for every variant
    let rhl_per_var = fragments_to_read_haplotype_likelihoods(&varlist, &frags, &read_meta);

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

        let mut total_reads = 0;
        let mut ro0 = 0;
        let mut ro1 = 0;
        let mut ao0 = 0;
        let mut ao1 = 0;

        for rhl in rhls {
            assert!(rhl.allele_call.len() == 1);
            total_reads += 1;
            let base = rhl.allele_call.chars().next().unwrap() as char;

            /*
            println!("{}\t{}\t{}\t{}\t{}\t{}\t{:.3}\t{:.3}", 
                rhl.read_name.as_ref().unwrap(), var.tid, var.pos0 + 1, var.alleles[0], var.alleles[1], base, *rhl.allele_call_qual, *min_allele_call_qual);
            */

            if rhl.allele_call_qual > min_allele_call_qual {
                continue;
            }

            if base != reference_base && base != alt_base {
                continue;
            }
         
            let (is_ref, is_rev) = (base == reference_base, rhl.strand_index == 1);
            
            match (is_ref, is_rev) {
                (true, false) => { ro0 += 1 }, 
                (true, true) =>  { ro1 += 1 }, 
                (false, false) => { ao0 += 1 }, 
                (false, true) => { ao1 += 1 }, 
            }
        }

        let ref_support = ro0 + ro1;
        let alt_support = ao0 + ao1;
        let total_support = ref_support + alt_support;
        //let vaf = alt_support as f32 / total_support as f32;

        let informative_fraction = total_support as f32 / total_reads as f32;

        // Make a conservative estimate of the error rate at this position using the context model
        let test_vaf = 0.10;

        /*
        let cm = &extract_fragment_parameters.context_model;
        let context_probs = &alignment_parameters.context_emission_probs.probs;
        let half_k = cm.k / 2;
        let model_ref_context = &chromosome_bytes[ (var.pos0 as usize - half_k)..(var.pos0 as usize + half_k + 1)];
        if cm.alphabet.is_word(model_ref_context) {
            let ref_base_rank = cm.get_base_rank(reference_base) as usize;
            let alt_base_rank = cm.get_base_rank(alt_base) as usize;

            let model_ref_er_fwd = context_probs[ cm.get_context_rank(model_ref_context, false).unwrap() ][alt_base_rank];
            let model_ref_er_rev = context_probs[ cm.get_context_rank(model_ref_context, true).unwrap() ][alt_base_rank];
            error_rate = f64::max(model_ref_er_fwd, model_ref_er_rev);
            //let context_ref_fwd_er_str = format!("{}+>{}:{:.3}", std::str::from_utf8(model_ref_context).unwrap(), alt_base, model_ref_er_fwd);
            //let context_ref_rev_er_str = format!("{}->{}:{:.3}", std::str::from_utf8(model_ref_context).unwrap(), alt_base, model_ref_er_ref);
            //println!("{} context model {} {}", var.pos0 + 1, context_ref_fwd_er_str, context_ref_rev_er_str);
        }
        */

        // test whether there is sufficient evidence of ref/alt allele
        // this computes a p-value that the vaf is greater than error_rate
        let binom_test_alt = binomial_test_onesided_greater(alt_support, total_support, test_vaf);
        let binom_test_ref = binomial_test_onesided_greater(ref_support, total_support, test_vaf);
        
        let mut record = vcf.empty_record();
        let rid = vcf.header().name2rid(region.chrom.as_bytes()).expect("Could not find reference id");
        record.set_rid(Some(rid));
        
        record.set_pos(var.pos0 as i64); 

        record.set_alleles(&[ var.alleles[0].as_bytes(), var.alleles[1].as_bytes() ]).expect("Could not set alleles");
        
        record.push_info_integer(b"TotalDepth", &[total_reads]).expect("Could not add INFO");
        
        record.push_info_float(b"EP", &[test_vaf as f32]).expect("Could not add INFO");
        let ref_lp = LogProb::from(Prob(binom_test_ref));
        let alt_lp = LogProb::from(Prob(binom_test_alt));
        record.push_info_float(b"LP", &[*ref_lp as f32, *alt_lp as f32]).expect("Could not add INFO");
        
        /*
        let cm = &extract_fragment_parameters.context_model;
        let context_probs = &alignment_parameters.context_emission_probs.probs;
        let half_k = cm.k / 2;
        let model_ref_context = &chromosome_bytes[ (var.pos0 as usize - half_k)..(var.pos0 as usize + half_k + 1)];
        if cm.alphabet.is_word(model_ref_context) {
            let alt_base_rank = cm.get_base_rank(alt_base) as usize;

            let model_ref_er_fwd = context_probs[ cm.get_context_rank(model_ref_context, false).unwrap() ][alt_base_rank];
            let model_ref_er_ref = context_probs[ cm.get_context_rank(model_ref_context, true).unwrap() ][alt_base_rank];
            let context_ref_fwd_er_str = format!("{}+>{}:{:.3}", std::str::from_utf8(model_ref_context).unwrap(), alt_base, model_ref_er_fwd);
            let context_ref_rev_er_str = format!("{}->{}:{:.3}", std::str::from_utf8(model_ref_context).unwrap(), alt_base, model_ref_er_ref);
            //println!("{} context model {} {}", var.pos0 + 1, context_ref_fwd_er_str, context_ref_rev_er_str);
        }
        */
        let fishers_result = 
            fishers_exact(&[ ro0 as u32, ro1 as u32,
                             ao0 as u32, ao1 as u32 ]).unwrap();

        //println!("{} StrandBiasData {} {} {} {}", var.pos0 + 1, ro0, ro1, ao0, ao1);

        // straight from longshot
        let strand_bias_pvalue = if fishers_result.two_tail_pvalue <= 500.0 {
            *PHREDProb::from(Prob(fishers_result.two_tail_pvalue))
        } else {
            500.0
        };
        record.push_info_float(b"StrandBias", &[strand_bias_pvalue as f32]).expect("Could not add INFO");
        record.push_info_float(b"InformativeFraction", &[informative_fraction]).expect("Could not add INFO");
        
        // absolute hack for development
        let lp_threshold = LogProb::from(-2.0);

        let has_ref_allele = ref_lp < lp_threshold;
        let has_alt_allele = alt_lp < lp_threshold;
        let mut alleles = if has_alt_allele && ! has_ref_allele {
            &gt_hom_alt
        } else if has_alt_allele && has_ref_allele {
            &gt_het
        } else {
            &gt_hom_ref
        };

        if strand_bias_pvalue > max_strand_bias || informative_fraction < min_informative_fraction || total_support < hard_min_depth {
            alleles = &gt_nocall;
        }
        record.push_genotypes(alleles).unwrap();

        record.push_format_integer(b"AD", &[ref_support as i32, alt_support as i32]).expect("Could not add FORMAT");
        vcf.write(&record).unwrap();
    }
}

