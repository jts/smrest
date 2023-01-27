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
use crate::LongshotParameters;
use crate::ReadHaplotypeLikelihood;
use longshot::extract_fragments::extract_fragments;
use longshot::variants_and_fragments::parse_vcf_potential_variants;
use longshot::variants_and_fragments::VarList;
use longshot::util::parse_region_string;

pub fn genotype(input_bam: &str, region_str: &str, candidates_vcf: &str, reference_genome: &str) {

    // subprogram thresholds
    let hard_min_depth = 20;
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
    vcf_header.push_record(r#"##INFO=<ID=GL,Number=3,Type=Float,Description="todo">"#.as_bytes());
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

    // priors
    let var_rate = 0.001;
    // rough prior that 2/3rds of variants are heterozygous
    let p_het = var_rate * 2.0 / 3.0;
    let p_hom_alt = var_rate - p_het;
    let p_hom_ref = 1.0 - p_hom_alt - p_het;

    /*

    // estimate the maximum likelihood allele frequency of hets 
    // calculate likelihood of het af
    // L(af) = P(data|af)
    //       = P(data| gt = 0/0, af) P(gt = 0/0) + P(data| gt = 0/1, af) P(gt = 0/1) + P(data| gt = 1/1, af) P(gt = 1/1)
    let var_rate = 0.001;
    let p_het = var_rate * 2.0 / 3.0;
    let p_hom_alt = var_rate - p_het;
    let p_hom_ref = 1.0 - p_hom_alt - p_het;

    let mut ml = f64::NEG_INFINITY;
    let mut ml_af = 0.0;

    let bins = 50;
    let step = 0.01;
    for i in 0..bins {
        let het_af = f64::from(i)*step;

        let mut lp_total = 0.0;
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

            let likelihoods = calculate_genotype_likelihoods(rhls, het_af);
            let lp_data_hom_ref:f64 = likelihoods[0];
            let lp_data_het:f64 = likelihoods[1];
            let lp_data_hom_alt:f64 = likelihoods[2];
            
            let lp_t_hom_ref = LogProb(lp_data_hom_ref) + LogProb::from(Prob(p_hom_ref));
            let lp_t_het = LogProb(lp_data_het) + LogProb::from(Prob(p_het));
            let lp_t_hom_alt = LogProb(lp_data_hom_alt) + LogProb::from(Prob(p_hom_alt));

            let lp_sum = LogProb::ln_sum_exp( &[ lp_t_hom_ref, lp_t_het, lp_t_hom_alt] );
            lp_total += *lp_sum;
        }

        if lp_total > ml {
            ml = lp_total;
            ml_af = het_af;
        }

        //println!("{:.3}\t{:.3}\t{:.3}\t{:.3}", het_af, lp_total, ml_af, ml);
    }
    */

    // call genotypes
    let het_af = 0.25;
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
        
        // genotype model
        let likelihoods = calculate_genotype_likelihoods(rhls, het_af);
        let lp_data_hom_ref:f64 = likelihoods[0];
        let lp_data_het:f64 = likelihoods[1];
        let lp_data_hom_alt:f64 = likelihoods[2];
        
        let lp_t_hom_ref = LogProb(lp_data_hom_ref) + LogProb::from(Prob(p_hom_ref));
        let lp_t_het = LogProb(lp_data_het) + LogProb::from(Prob(p_het));
        let lp_t_hom_alt = LogProb(lp_data_hom_alt) + LogProb::from(Prob(p_hom_alt));
        let lp_sum = LogProb::ln_sum_exp( &[ lp_t_hom_ref, lp_t_het, lp_t_hom_alt] );

        let p_hom_ref_data = Prob::from(lp_t_hom_ref - lp_sum);
        let p_het_data = Prob::from(lp_t_het - lp_sum);
        let p_hom_alt_data = Prob::from(lp_t_hom_alt - lp_sum);
        
        let gls = [ *PHREDProb::from(Prob(1.0 - *p_hom_ref_data)) as f32,
                    *PHREDProb::from(Prob(1.0 - *p_het_data)) as f32,
                    *PHREDProb::from(Prob(1.0 - *p_hom_alt_data)) as f32 ];

        let mut max_val = 0.0;
        let mut max_gt_idx = 0;
        for (idx, value) in gls.into_iter().enumerate() {
            if value > max_val {
                max_val = value;
                max_gt_idx = idx;
            }
        }

        let mut record = vcf.empty_record();
        let rid = vcf.header().name2rid(region.chrom.as_bytes()).expect("Could not find reference id");
        record.set_rid(Some(rid));
        
        record.set_pos(var.pos0 as i64); 

        record.set_alleles(&[ var.alleles[0].as_bytes(), var.alleles[1].as_bytes() ]).expect("Could not set alleles");
        
        record.push_info_integer(b"TotalDepth", &[total_reads]).expect("Could not add INFO");
        record.push_info_float(b"GL", &gls).expect("Could not add INFO");
        
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
        
        let mut alleles = if max_gt_idx == 0 {
            &gt_hom_ref
        } else if max_gt_idx == 1 {
            &gt_het
        } else {
            &gt_hom_alt           
        };
        
        if strand_bias_pvalue > max_strand_bias || total_support < hard_min_depth {
            alleles = &gt_nocall;
        }
        record.push_genotypes(alleles).unwrap();

        record.push_format_integer(b"AD", &[ref_support as i32, alt_support as i32]).expect("Could not add FORMAT");
        vcf.write(&record).unwrap();
    }
}

fn calculate_genotype_likelihoods(rhls: &Vec<ReadHaplotypeLikelihood>, het_af: f64) -> [f64;3] {
    let lp_data_hom_ref:f64 = rhls.iter().map(|x| *x.base_allele_likelihood).sum();
    let lp_data_hom_alt:f64 = rhls.iter().map(|x| *x.mutant_allele_likelihood).sum();

    // P(data | gt = het) = sum over reads P(read | het)
    let minor_hap_freq = Prob(het_af);
    let major_hap_freq = Prob(1.0 - het_af);
    let lp_minor_hap_freq = LogProb::from(minor_hap_freq);
    let lp_major_hap_freq = LogProb::from(major_hap_freq);
    let lp_half = LogProb::from(Prob(0.5));
    let mut t1_sum = LogProb(0.0);
    let mut t2_sum = LogProb(0.0);

    for rhl in rhls {
        assert!(rhl.allele_call.len() == 1);
        
        // we don't know which haplotype (if any) is amplified, so need to sum over
        // both equally likely cases
        let t1_h1 = rhl.base_allele_likelihood + lp_minor_hap_freq;
        let t1_h2 = rhl.mutant_allele_likelihood + lp_major_hap_freq;
        let t1 = LogProb::ln_add_exp(t1_h1, t1_h2);
        t1_sum += t1;

        let t2_h1 = rhl.base_allele_likelihood + lp_major_hap_freq;
        let t2_h2 = rhl.mutant_allele_likelihood + lp_minor_hap_freq;
        let t2 = LogProb::ln_add_exp(t2_h1, t2_h2);
        t2_sum += t2;

        //println!("\tread\t{}\tlp_allele [{:.3} {:.3}] lp_freq [{:.3} {:.3}] t1 [{:.3} {:.3} {:.3}] t2: [{:.3} {:.3} {:.3}]", 
        //    base, *rhl.base_allele_likelihood, *rhl.mutant_allele_likelihood, *lp_minor_hap_freq, *lp_major_hap_freq, *t1_h1, *t1_h2, *t1, *t2_h1, *t2_h2, *t2);
    }

    let lp_data_het = *LogProb::ln_add_exp(t1_sum + lp_half, t2_sum + lp_half);
    return [lp_data_hom_ref, lp_data_het, lp_data_hom_alt];
}
