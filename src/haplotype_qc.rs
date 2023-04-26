//---------------------------------------------------------
// Copyright 2023 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
use std::collections::{HashMap, BTreeMap};
use rust_htslib::{bam, bcf, bam::Read};
use crate::utility::*;
use crate::LongshotParameters;
use longshot::extract_fragments::extract_fragments;
use longshot::variants_and_fragments::VarList;
use longshot::util::parse_region_string;

pub struct ReadHaplotypeStats
{
    pub likelihoods: [f64; 2],
    pub variant_count: u32,
    pub alleles: Vec<char>
}

impl ReadHaplotypeStats {
    pub fn new(num_vars: usize) -> ReadHaplotypeStats {
        ReadHaplotypeStats { 
            likelihoods: [0.0, 0.0], 
            variant_count: 0,
            alleles: vec!['-'; num_vars]
        }
    }
}

pub fn haplotype_qc(input_bam: &str, region_str: &str, phased_vcf: &str, reference_genome: &str) {

    let mut longshot_parameters = LongshotParameters::defaults();
    longshot_parameters.estimate_alignment_parameters(&input_bam.to_owned(), &reference_genome.to_owned(), &None);
    
    let region = parse_region_string(Some(region_str), &input_bam.to_owned()).expect("Could not parse region").unwrap();
    
    let mut bam = bam::IndexedReader::from_path(input_bam).unwrap();
    let header_view = bam.header();
    
    let bcf_records: Vec<bcf::Record> = read_bcf(&phased_vcf.to_owned(), Some(region.clone()))
                                           .into_iter().filter(|x| is_genotype_het(x)).collect();
    
    // convert to longshot var here
    let all_vars = bcf_records.iter().map(|x| bcf2longshot(x, header_view)).collect();

    let target_names = header_view.target_names().iter().map(|x| std::str::from_utf8(x).unwrap().to_owned()).collect();
    let mut varlist = VarList::new(all_vars, target_names).expect("Could not construct var list");

    println!("Found {} vars", varlist.len());
    
    let frags = extract_fragments(&input_bam.to_owned(),
                                  &reference_genome.to_owned(),
                                  &mut varlist,
                                  &Some(region.clone()),
                                  &longshot_parameters.extract_fragment_parameters,
                                  &longshot_parameters.alignment_parameters.as_ref().unwrap()).unwrap();
    
    // populate read metadata
    let read_meta = populate_read_metadata_from_bam(&input_bam.to_owned(), &region, None, None, None, None);

    // Convert fragments into read-haplotype likelihoods for every variant
    let rhl_per_var = fragments_to_read_haplotype_likelihoods(&varlist, &frags, &read_meta);

    // partition the vector(s) of variants by phase set
    // the variants are represented by an index into bcf_records/all_vars
    let mut phase_set_partitions: BTreeMap<i32, Vec<usize>> = BTreeMap::new();
    for idx in 0..bcf_records.len() {
        if let Ok(phase_sets) = bcf_records[idx].format(b"PS").integer() {
            assert!(phase_sets.len() == 1);
            let ps = phase_sets[0][0];
            phase_set_partitions.entry(ps).or_insert(Vec::new()).push(idx);
        }
    }

    for (k, variant_indices) in phase_set_partitions.iter() {
        println!("Phase set {} has {} variants", k, variant_indices.len());
        
        // determine allele set on each haplotype
        let mut h1_alleles = Vec::new();
        let mut h2_alleles = Vec::new();
        for idx in variant_indices.iter() {
            let gt = bcf_records[*idx].genotypes().expect("Error reading genotypes").get(0);
            
            let h1 = gt[0].index().unwrap() + 48;
            let h2 = gt[1].index().unwrap() + 48;

            h1_alleles.push(std::char::from_u32(h1).unwrap());
            h2_alleles.push(std::char::from_u32(h2).unwrap());
        }

        println!("h1: {:?}", h1_alleles);
        println!("h2: {:?}", h2_alleles);

        let mut read_haplotype_stats: HashMap<String, ReadHaplotypeStats> = HashMap::new();

        // attempt to assign reads to haplotypes
        for (i, idx) in variant_indices.iter().enumerate() {
            let bcf_record = &bcf_records[*idx];
            //let var = &varlist.lst[*idx];
            let rhls = &rhl_per_var[*idx];
            for rhl in rhls {
                let id = rhl.read_name.as_ref().unwrap().clone();

                let h1_allele = bcf_record.genotypes().expect("Error reading genotypes").get(0)[0].index().unwrap();
                let mut v = read_haplotype_stats.entry(id).or_insert(ReadHaplotypeStats::new(variant_indices.len()));
                if h1_allele == 0 {
                    v.likelihoods[0] += *rhl.base_allele_likelihood;
                    v.likelihoods[1] += *rhl.mutant_allele_likelihood;
                } else {
                    v.likelihoods[1] += *rhl.base_allele_likelihood;
                    v.likelihoods[0] += *rhl.mutant_allele_likelihood;
                }
                v.variant_count += 1;
                //v.alleles[i] = rhl.allele_call.chars().next().unwrap(); // first char of allele
                v.alleles[i] = if rhl.base_allele_likelihood > rhl.mutant_allele_likelihood { '0' } else { '1' };
            }
        }

        for (read_name, stats) in read_haplotype_stats {
            let mut h1_diff = stats.alleles.clone();
            let mut h2_diff = stats.alleles.clone();
            
            for (i, c) in stats.alleles.iter().enumerate() {
                if *c == h1_alleles[i] { h1_diff[i] = '.'; }
                if *c == h2_alleles[i] { h2_diff[i] = '.'; }
            }

            println!("{}\t{}\t{}\t{}\t{}\t{}\t{:.3}\t{:.3}", k, read_name, stats.variant_count, 
                                             stats.alleles.iter().collect::<String>(), 
                                             h1_diff.iter().collect::<String>(), 
                                             h2_diff.iter().collect::<String>(), 
                                             stats.likelihoods[0], stats.likelihoods[1]); 
        }
        
    }
}
