//---------------------------------------------------------
// Copyright 2022 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
extern crate approx;

use rust_htslib::{bam, bam::Read, bam::record::Aux};
use rust_htslib::bcf::{Reader as BcfReader, Read as BcfRead, Header as BcfHeader};
use std::collections::{HashSet, HashMap};
use clap::{App, SubCommand, Arg, value_t};
use itertools::Itertools;
use statrs::distribution::{Beta, Poisson};
use fishers_exact::fishers_exact;
use bio::stats::{Prob, LogProb, PHREDProb};

mod pileup_stats;
use crate::pileup_stats::PileupStats;
use crate::pileup_stats::base2index;

mod simulation;
use crate::simulation::*;

mod classifier;
use crate::classifier::*;

mod parameters;
use crate::parameters::*;

mod genotype;
use crate::genotype::*;

mod haplotype_qc;
use crate::haplotype_qc::*;

mod select_reads;
use crate::select_reads::*;

mod somatic_call;
use crate::somatic_call::*;

mod utility;
use crate::utility::*;

mod longshot_realign;
use longshot_realign::{ReadHaplotypeLikelihood, ReadMetadata};

mod calling_models;
use crate::calling_models::str_to_calling_model;

use longshot::util::{parse_region_string};
use longshot::extract_fragments::extract_fragments;
use longshot::variants_and_fragments::{Var, VarList, parse_vcf_potential_variants};

fn main() {
    let matches = App::new("smrest")
        .version("0.1")
        .author("Jared Simpson <jared.simpson@oicr.on.ca>")
        .about("Toolkit for estimating somatic mutation rate from long reads")
        .subcommand(SubCommand::with_name("sim-pileup")
                .about("simulate pileup data")
                .arg(Arg::with_name("purity")
                    .long("purity")
                    .takes_value(true)
                    .help("purity of tumor sample"))
                .arg(Arg::with_name("depth")
                    .long("depth")
                    .takes_value(true)
                    .help("mean sequencing depth"))
                .arg(Arg::with_name("genome")
                    .short('g')
                    .long("genome")
                    .takes_value(true)
                    .help("the reference genome"))
                .arg(Arg::with_name("error-rate")
                    .short('e')
                    .long("error-rate")
                    .takes_value(true)
                    .help("the error rate of simulated reads"))
                .arg(Arg::with_name("tmb")
                    .short('t')
                    .long("tmb")
                    .takes_value(true)
                    .help("the tumour mutation burden, in mutations per megabase"))
                .arg(Arg::with_name("proportion-clonal")
                    .short('p')
                    .long("proportion-clonal")
                    .takes_value(true)
                    .help("the proportion of simulated mutations that are fully clonal (ccf = 1.0)"))
                .arg(Arg::with_name("mutation-output")
                    .short('m')
                    .long("mutation-output")
                    .takes_value(true)
                    .help("output filename for simulated mutations"))
                .arg(Arg::with_name("per-site-output")
                    .long("per-site-output")
                    .takes_value(false)
                    .help("change the output style to write one line per simulated site")))
        .subcommand(SubCommand::with_name("extract")
                .about("gather candidate somatic mutations from aligned reads")
                .arg(Arg::with_name("genome")
                    .short('g')
                    .long("genome")
                    .takes_value(true)
                    .help("the reference genome"))
                .arg(Arg::with_name("regions")
                    .short('r')
                    .long("regions")
                    .takes_value(true)
                    .help("only report positions within these regions"))
                .arg(Arg::with_name("min-depth")
                    .short('d')
                    .long("min-depth")
                    .takes_value(true)
                    .help("only report positions within phased depth greater than N"))
                .arg(Arg::with_name("input-bam")
                    .required(true)
                    .index(1)
                    .help("the input bam file to process")))
        .subcommand(SubCommand::with_name("genotype-hets")
                .about("find potential heterozygous variants")
                .arg(Arg::with_name("genome")
                    .short('g')
                    .long("genome")
                    .takes_value(true)
                    .help("the reference genome"))
                .arg(Arg::with_name("region")
                    .short('r')
                    .long("region")
                    .takes_value(true)
                    .help("the reference region to call"))
                .arg(Arg::with_name("candidates")
                    .short('c')
                    .long("candidates")
                    .takes_value(true)
                    .required(true)
                    .help("vcf file containing candidate variants, typically from gnomad"))
                .arg(Arg::with_name("input-bam")
                    .required(true)
                    .index(1)
                    .help("the input bam file to process")))
        .subcommand(SubCommand::with_name("haplotype-qc")
                .about("compute summary stats for haplotypes in the phased vcf file")
                .arg(Arg::with_name("genome")
                    .short('g')
                    .long("genome")
                    .takes_value(true)
                    .help("the reference genome"))
                .arg(Arg::with_name("region")
                    .short('r')
                    .long("region")
                    .takes_value(true)
                    .help("the reference region to call"))
                .arg(Arg::with_name("phased-vcf")
                    .short('p')
                    .long("phased-vcf")
                    .takes_value(true)
                    .required(true)
                    .help("vcf file containing phased variants"))
                .arg(Arg::with_name("input-bam")
                    .required(true)
                    .index(1)
                    .help("the input bam file to process")))
        .subcommand(SubCommand::with_name("select-reads")
                .about("subset a bam file to reads that are informative for phasing")
                .arg(Arg::with_name("genome")
                    .short('g')
                    .long("genome")
                    .takes_value(true)
                    .help("the reference genome"))
                .arg(Arg::with_name("region")
                    .short('r')
                    .long("region")
                    .takes_value(true)
                    .help("the reference region to call"))
                .arg(Arg::with_name("vcf")
                    .short('c')
                    .long("vcf")
                    .takes_value(true)
                    .required(true)
                    .help("vcf file containing heterozygous variants for phasing"))
                .arg(Arg::with_name("output-bam")
                    .short('o')
                    .long("output-bam")
                    .takes_value(true)
                    .required(true)
                    .help("the output bam file to write to"))
                .arg(Arg::with_name("input-bam")
                    .required(true)
                    .index(1)
                    .help("the input bam file to process")))
        .subcommand(SubCommand::with_name("phase")
                .about("phase heterozygous variants")
                .arg(Arg::with_name("genome")
                    .short('g')
                    .long("genome")
                    .takes_value(true)
                    .help("the reference genome"))
                .arg(Arg::with_name("region")
                    .short('r')
                    .long("region")
                    .takes_value(true)
                    .help("the reference region to phase"))
                .arg(Arg::with_name("candidates")
                    .short('c')
                    .long("candidates")
                    .takes_value(true)
                    .required(true)
                    .help("vcf file containing candidate variants to phase"))
                .arg(Arg::with_name("input-bam")
                    .required(true)
                    .index(1)
                    .help("the input bam file to process")))
        .subcommand(SubCommand::with_name("call")
                .about("call somatic mutations")
                .arg(Arg::with_name("genome")
                    .short('g')
                    .long("genome")
                    .takes_value(true)
                    .help("the reference genome"))
                .arg(Arg::with_name("region")
                    .short('r')
                    .long("region")
                    .takes_value(true)
                    .help("the reference region to call"))
                .arg(Arg::with_name("phased-vcf")
                    .short('p')
                    .long("phased-vcf")
                    .takes_value(true)
                    .help("use the SNPs in this VCF file to infer the haplotype of each read"))
                .arg(Arg::with_name("output-region-bed")
                    .short('o')
                    .long("output-region-bed")
                    .takes_value(true)
                    .help("write the regions that had sufficient phased depth to be callable to this .bed file"))
                .arg(Arg::with_name("model")
                    .short('m')
                    .long("model")
                    .takes_value(true)
                    .help("the probabilistic model used to call mutations (raw-pileup, realigned-pileup)"))
                .arg(Arg::with_name("input-bam")
                    .required(true)
                    .index(1)
                    .help("the input bam file to process")))
        .subcommand(SubCommand::with_name("show-phase")
                .about("show the base present in each read at specified positions")
                .arg(Arg::with_name("positions")
                    .short('p')
                    .long("positions")
                    .takes_value(true)
                    .help("a comma separated list of positions to query"))
                .arg(Arg::with_name("haplotype")
                    .short('c')
                    .long("haplotype")
                    .takes_value(true)
                    .help("only include reads on this haplotype"))
                .arg(Arg::with_name("input-bam")
                    .required(true)
                    .index(1)
                    .help("the input bam file to process")))
        .subcommand(SubCommand::with_name("error-contexts")
                .about("gather error rate statistics per sequence context")
                .arg(Arg::with_name("vcf")
                    .short('v')
                    .long("vcf")
                    .takes_value(true)
                    .help("vcf file, these positions will not be included in context analysis"))
                .arg(Arg::with_name("genome")
                    .short('g')
                    .long("genome")
                    .takes_value(true)
                    .help("the reference genome"))
                .arg(Arg::with_name("input-bam")
                    .required(true)
                    .index(1)
                    .help("the input bam file to process"))).get_matches();
    
    if let Some(matches) = matches.subcommand_matches("extract") {
        
        let min_depth = value_t!(matches, "min-depth", u32).unwrap_or(20);
        extract_mutations(matches.value_of("input-bam").unwrap(),
                          matches.value_of("genome").unwrap(),
                          min_depth,
                          matches.value_of("regions"));
    } else if let Some(matches) = matches.subcommand_matches("call") {
        let model_str = matches.value_of("model").unwrap_or("haplotype-likelihood");
        let model = str_to_calling_model(model_str).expect("unknown calling model");

        somatic_call(matches.value_of("input-bam").unwrap(),
                     matches.value_of("region").unwrap(),
                     matches.value_of("genome").unwrap(),
                     matches.value_of("output-region-bed"),
                     matches.value_of("phased-vcf"),
                     model)
    } else if let Some(matches) = matches.subcommand_matches("genotype-hets") {
        genotype(matches.value_of("input-bam").unwrap(),
                 matches.value_of("region").unwrap(),
                 matches.value_of("candidates").unwrap(),
                 matches.value_of("genome").unwrap())
    } else if let Some(matches) = matches.subcommand_matches("haplotype-qc") {
        haplotype_qc(matches.value_of("input-bam").unwrap(),
                     matches.value_of("region").unwrap(),
                     matches.value_of("phased-vcf").unwrap(),
                     matches.value_of("genome").unwrap())
    } else if let Some(matches) = matches.subcommand_matches("select-reads") {
        select_reads(matches.value_of("input-bam").unwrap(),
                     matches.value_of("output-bam").unwrap(),
                     matches.value_of("region").unwrap(),
                     matches.value_of("vcf").unwrap(),
                     matches.value_of("genome").unwrap())
    } else if let Some(matches) = matches.subcommand_matches("phase") {
        phase(matches.value_of("input-bam").unwrap(),
              matches.value_of("region").unwrap(),
              matches.value_of("candidates").unwrap(),
              matches.value_of("genome").unwrap())
    } else if let Some(matches) = matches.subcommand_matches("error-contexts") {
        error_contexts(matches.value_of("input-bam").unwrap(),
                       matches.value_of("genome").unwrap(),
                       matches.value_of("vcf").unwrap());
    } else if let Some(matches) = matches.subcommand_matches("sim-pileup") {
        let depth_lambda = value_t!(matches, "depth", f64).unwrap_or(50.0);
        let purity = value_t!(matches, "purity", f64).unwrap_or(0.75);
        let error_rate = value_t!(matches, "error-rate", f64).unwrap_or(0.02);
        let muts_per_megabase = value_t!(matches, "tmb", f64).unwrap_or(5.0);
        let per_site_output = matches.is_present("per-site-output");
        let p_clonal = value_t!(matches, "proportion-clonal", f64).unwrap_or(0.9);
        let mutation_output = matches.value_of("mutation-output").unwrap_or("simulated_mutations.tsv");

        let ccf = CancerCellFraction { p_clonal: p_clonal, subclonal_ccf: Beta::new(2.0, 2.0).unwrap() };
        let params = ModelParameters { 
            mutation_rate: muts_per_megabase / 1000000.0, // per haplotype
            heterozygosity: 1.0 / 2000.0, // per haplotype
            ccf_dist: ccf,
            depth_dist: Some(Poisson::new(depth_lambda).unwrap()),
            purity: purity,
            error_rate: error_rate
        };
        sim_pileup(&params, matches.value_of("genome").unwrap(), per_site_output, &mutation_output);
    } if let Some(matches) = matches.subcommand_matches("show-phase") {
        let haplotype = value_t!(matches, "haplotype", i32).unwrap_or(0);
        show_phase(matches.value_of("input-bam").unwrap(), 
                   matches.value_of("positions").unwrap(),
                   haplotype);
    }
}

fn extract_mutations(input_bam: &str, reference_genome: &str, min_depth: u32, region_opt: Option<&str>) {
    /*
    let params = ModelParameters::defaults();
        
    let mut bam = bam::IndexedReader::from_path(input_bam).unwrap();
    let header_view = bam.header();
    let regions = match region_opt {
        Some(bed_filename) => Some(GenomeRegions::from_bed(bed_filename, &header_view)),
        None => None
    };

    println!("chromosome\tposition\treference_base\tvariant_base\tcanonical_type\tcanonical_context\taligned_depth\tmean_mapq\t\
              hmajor_variant_count\thminor_variant_count\thmajor_vaf\thminor_vaf\th1_a\th1_c\th1_g\th1_t\th2_a\th2_c\th2_g\th2_t\t\
              p_ref\tp_germline\tp_somatic");

    let all_positions = true;
    let chromosome_name = "chr20";
    let min_variant_observations = 3;
    //let max_variant_minor_observations = 1;
    let min_variant_observations_per_strand = 1;

    let tid = header_view.tid(chromosome_name.as_bytes()).unwrap();
    let chromosome_bytes = get_chromosome_sequence(reference_genome, header_view, tid).into_bytes();
    
    let read_metadata = populate_read_metadata_from_bam(&input_bam.to_owned(), 
                                                        &region,
                                                        Some(&chromosome_bytes), 
                                                        None, None, None); 

    let mut cache = ReadHaplotypeCache { cache: HashMap::new() };
    let mut ps = PileupStats::new();

    // go to chromosome of interest
    bam.fetch( chromosome_name ).unwrap();
    for p in bam.pileup() {
        let pileup = p.unwrap();
        
        // check whether we should bother with this position
        let reference_tid = pileup.tid();
        let reference_position = pileup.pos() as usize;

        // if a bed file was provided, check whether this position is contained in the provided regions
        if let Some(ref r) = regions {
            if !r.contains(reference_tid, reference_position) {
                continue;
            }
        }

        let reference_base = chromosome_bytes[reference_position] as char;
        if reference_base == 'N' {
            continue;
        }

        ps.fill_pileup(&mut cache, pileup.alignments());
        if ps.get_depth() < min_depth {
            continue;
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
                if bi != reference_base_index && (b0 + b1) > max_variant_count && 
                    b0 >= min_variant_observations_per_strand && b1 >= min_variant_observations_per_strand {
                    max_variant_count = b0 + b1;
                    max_variant_index = bi;
                    max_haplotype_index = hi;
                }
            }
        }

        let minor_haplotype_index = 1 - max_haplotype_index;
        let mut minor_variant_count = 0;
        if max_variant_index != reference_base_index {
            minor_variant_count = ps.get_count_on_haplotype(max_variant_index, minor_haplotype_index);
        }

        if all_positions || (max_variant_count >= min_variant_observations /*&& minor_variant_count <= max_variant_minor_observations*/) {

            // grab reference context
            let mut reference_context = chromosome_bytes[ (reference_position - 1)..(reference_position + 2)].to_vec();

            let bases = "ACGT";
            let variant_base = bases.as_bytes()[max_variant_index as usize] as char;
            
            let ref_count = ps.get_count_on_haplotype(reference_base_index, max_haplotype_index);
            let class_probs = calculate_class_probabilities_phased(max_variant_count as u64, ref_count as u64, &params);
            
            // grab mutation type
            let mut mutation_type: [char; 3] = [ reference_base, '>', variant_base ];

            // convert mutation type/context to canonical form C>x, T>x
            if mutation_type[0] != 'C' && mutation_type[0] != 'T' {
                mutation_type[0] = bio::alphabets::dna::complement(mutation_type[0] as u8) as char;
                mutation_type[2] = bio::alphabets::dna::complement(mutation_type[2] as u8) as char;
                reference_context = bio::alphabets::dna::revcomp(reference_context);
            }
            
            let mutation_type_str = String::from_iter(&mutation_type);
            let context_str = std::str::from_utf8(&reference_context).unwrap();

            // read support statistics
            let position = pileup.pos() + 1; // to match vcf
            
            let hmaj_vaf = max_variant_count as f32 / ps.get_haplotype_depth(max_haplotype_index) as f32;
            let hmin_vaf = minor_variant_count as f32 / ps.get_haplotype_depth(minor_haplotype_index) as f32;
            
            let h0 = ps.get_haplotype_counts(max_haplotype_index);
            let h1 = ps.get_haplotype_counts(minor_haplotype_index);
            let mean_mapq = ps.mean_mapq;
            let aligned_depth = ps.get_haplotype_depth(0) + ps.get_haplotype_depth(1);
            println!("{chromosome_name}\t{position}\t{reference_base}\t{variant_base}\t{mutation_type_str}\t\
                      {context_str}\t{aligned_depth}\t{mean_mapq}\t{max_variant_count}\t{minor_variant_count}\t\
                      {hmaj_vaf:.3}\t{hmin_vaf:.3}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.3}\t{:.3}\t{:.3}",
                      h0[0], h0[1], h0[2], h0[3],
                      h1[0], h1[1], h1[2], h1[3],
                      class_probs[0], class_probs[1], class_probs[2]);
        }
    }
    */
}

fn phase(input_bam: &str, region_str: &str, candidates_vcf: &str, reference_genome: &str) {

    let min_allele_call_qual = LogProb::from(Prob(0.1));
    let hard_min_depth = 20;
    let min_informative_fraction = 0.75;

    let max_strand_bias = 10.0;
    
    let mut longshot_parameters = LongshotParameters::defaults();
    longshot_parameters.estimate_alignment_parameters(&input_bam.to_owned(), &reference_genome.to_owned(), &None);
    
    let region = parse_region_string(Some(region_str), &input_bam.to_owned()).expect("Could not parse region").unwrap();
    let mut bam = bam::IndexedReader::from_path(input_bam).unwrap();
    
    // read variants within this region
    let all_vars = parse_vcf_potential_variants(&candidates_vcf.to_owned(), &input_bam.to_owned())
                        .expect("Could not parse candidate variants file");
    
    let region_vars: Vec<Var> = all_vars.get_variants_range(region.clone()).expect("Could not subset vars");
    let mut varlist = VarList::new(region_vars, all_vars.target_names).expect("Could not construct var list");

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

    // compute set of candidate hets
    let mut candidate_hets = Vec::new();

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
        let vaf = alt_support as f32 / total_support as f32;

        let informative_fraction = total_support as f32 / total_reads as f32;

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
        
        let mut is_het = false;

        /*
        let lp_threshold = LogProb::from(-10.0);
        if alt_lp < lp_threshold && ref_lp < lp_threshold {
            is_het = true;
        }*/
        let min_vaf = 0.1;
        if vaf >= min_vaf && vaf <= 1.0 - min_vaf {
            is_het = true;
        }

        if strand_bias_pvalue > max_strand_bias || informative_fraction < min_informative_fraction || total_support < hard_min_depth {
            is_het = false;
        }

        if is_het {
            candidate_hets.push(var.clone());
        }
    }

    eprintln!("Found {} candidates", candidate_hets.len());

    // make a map with the idx of each het in the original varlist pointing to the new (enumerated) index
    let candidate_het_idx_map: HashMap<_,_> = candidate_hets.iter().enumerate().map(|t| (t.1.ix, t.0)).collect();

    for f in &frags {
        let id = f.id.as_ref().unwrap();

        // initialize empty haplotype
        let mut read_haplotype = vec!['-'; candidate_hets.len()];

        for c in &f.calls {
            
            if let Some(het_idx) = candidate_het_idx_map.get(&c.var_ix) {
                // if true then the current variant is a het, and we have the new idx for it

                //let var = &varlist.lst[c.var_ix];
                //let allele_call = var.alleles[c.allele as usize].clone();
                read_haplotype[*het_idx] = char::from_digit(c.allele as u32, 10).unwrap();
                /*
                println!("{}\t{}\t{}\t{}\t{}\t{}", &id, var.tid, var.pos0 + 1, var.alleles[0], var.alleles[1], c.allele);
        
                let rhl = ReadHaplotypeLikelihood { 
                    read_name: Some(id.clone()),
                    mutant_allele_likelihood: c.allele_scores[1],
                    base_allele_likelihood: c.allele_scores[0],
                    allele_call: var.alleles[c.allele as usize].clone(),
                    allele_call_qual: c.qual,
                    haplotype_index: rm.haplotype_index,
                    strand_index: rm.strand_index
                };
                */
            }
        }

        let s: String = read_haplotype.iter().collect();
        println!("{}\t{}", &id, s);
    }

}

fn error_contexts(input_bam: &str, reference_genome: &str, variants_vcf: &str) {

    let chromosome_name = String::from("chr20");
    let k = 9;

    // read variants into a set, we ignore these positions
    let mut variant_mask: HashSet< (String, usize) > = HashSet::new();
    let mut bcf = BcfReader::from_path(variants_vcf).expect("Error opening file.");
    for v in bcf.records() {
        let vcf_record = v.unwrap();
        let name = vcf_record.header().rid2name(vcf_record.rid().unwrap()).unwrap();
        let name_str = String::from_utf8(name.to_vec()).unwrap();
        let pos = vcf_record.pos() as usize;
        variant_mask.insert( (name_str, pos) );
    }

    // set up header and faidx for pulling reference sequence
    let mut bam = bam::IndexedReader::from_path(input_bam).unwrap();
    let header_view = bam.header();
    let tid = header_view.tid(chromosome_name.as_bytes()).unwrap();
    let chromosome_bytes = get_chromosome_sequence(reference_genome, header_view, tid).into_bytes();
    
   let mut context_counts = HashMap::new();
   let mut error_counts = HashMap::new();
   let mut error_prob = HashMap::new();

    // go to chromosome of interest
    bam.fetch(chromosome_name.as_str()).unwrap();
    for p in bam.pileup() {
        let pileup = p.unwrap();
        
        let reference_position = pileup.pos() as usize;
        let reference_base = chromosome_bytes[reference_position] as char;
        let reference_base_index = base2index(reference_base) as u32;

        if reference_base == 'N' {
            continue;
        }

        let key = (chromosome_name.clone(), reference_position);
        if variant_mask.contains( &key ) {
            continue;
        }

        assert!(k % 2 == 1);
        let h = k / 2;

        let reference_context_fwd = chromosome_bytes[ (reference_position - h)..(reference_position + h + 1)].to_vec();
        let context_str_fwd = String::from_utf8(reference_context_fwd.clone()).unwrap();
        
        let reference_context_rev = bio::alphabets::dna::revcomp(reference_context_fwd);
        let context_str_rev = String::from_utf8(reference_context_rev.clone()).unwrap();
        
        if context_str_fwd.contains("N") {
            continue;
        }

        for a in pileup.alignments() {

            if a.record().seq().len() == 0 {
                continue;
            }
            if let Some(qpos) = a.qpos() {
                let read_base = a.record().seq()[qpos] as char;
                let bi = base2index(read_base) as u32;

                let (c, _)  = match a.record().is_reverse() {
                    false => (&context_str_fwd, read_base), 
                    true => (&context_str_rev, bio::alphabets::dna::complement(read_base as u8) as char)
                };
                let eb = 'N';
                
                // update number of times context has been seen
                let cc = context_counts.entry( c.clone() ).or_insert( 0usize );
                
                // update q-score counts
                let ep = error_prob.entry( c.clone() ).or_insert( 0.0f64 );
                let q = a.record().qual()[qpos] as f64;
                let base: f64 = 10.0;
                *ep += base.powf(-q / 10.0);
                
                // update observed basecalling error count
                if bi != reference_base_index {
                    let ec = error_counts.entry( (c.clone(), eb) ).or_insert( 0usize );
                    *ec += 1;
                }
                *cc += 1;
            }
        }
    }
    
    println!("context\terror\tnum_errors\ttotal_observations\tobserved_error_rate\tobserved_qscore\tmean_quality_prob\tmean_qscore");
    for context in context_counts.keys().sorted() {

        for b in [ 'A', 'C', 'G', 'T', 'N' ] {
            if error_counts.contains_key( &(context.to_string(), b)) {
                let count = context_counts.get( context ).unwrap();
                let errors = error_counts.get( &(context.to_string(), b) ).unwrap();
                let sum_p = error_prob.get( context ).unwrap();

                let calling_error_rate = *errors as f64 / *count as f64;
                let calling_error_q = -10.0 * calling_error_rate.log10();

                let mean_qprob = sum_p / *count as f64;
                let mean_qscore = -10.0 * mean_qprob.log10();
                println!("{context}\t{b}\t{errors}\t{count}\t{calling_error_rate:.4}\t{calling_error_q:.1}\t{mean_qprob:.4}\t{mean_qscore:.1}");
            }
        }
    }
}

fn show_phase(input_bam: &str, position_str: &str, show_haplotype: i32) {
    let mut positions = Vec::new();

    for p in position_str.split(',') {
        let mut s2 = p.split(':');
        let chrom = s2.next().unwrap();
        let position:u32 = s2.next().unwrap().parse().unwrap(); // convert to 0 based
        positions.push( (chrom, position - 1) );
    }
    
    if positions.len() == 0 {
        return;
    }

    // determine span of region
    let region_min = positions.iter().map(|x| x.1).min().unwrap();
    let region_max = positions.iter().map(|x| x.1).max().unwrap() + 1;
    let region_chr = positions[0].0;

    for p in &positions {
        assert!(p.0 == region_chr);
    }

    // set up header and faidx for pulling reference sequence
    let mut bam = bam::IndexedReader::from_path(input_bam).unwrap();
    let header = bam::Header::from_template(bam.header());
    let header_view = bam::HeaderView::from_header(&header);
    let tid = header_view.tid(region_chr.as_bytes()).unwrap();

    // build read group to sample map
    let mut rg_to_sample = HashMap::< String, String >::new();

    for rg_entry in header.to_hashmap().get("RG").unwrap() {
        let id = rg_entry.get("ID").unwrap();
        let sm = rg_entry.get("SM").unwrap();
        rg_to_sample.insert(id.to_string(), sm.to_string());
    }

    // go to chromosome of interest
    let mut read_alleles = HashMap::< String, Vec::<char> >::new();
    let mut read_groups = HashMap::< String, String >::new();
    let empty_alleles = vec!['?'; positions.len()];

    bam.fetch( (tid, region_min, region_max) ).unwrap();
    for p in bam.pileup() {
        let pileup = p.unwrap();
        if let Some(index) = positions.iter().position(|x| x.1 == pileup.pos()) {

            for a in pileup.alignments() {
                let read_name = String::from_utf8(a.record().qname().to_vec()).unwrap();
                match a.record().aux(b"RG") {
                    Ok(value) => {
                        if let Aux::String(v) = value {
                            let sm = rg_to_sample.get(&v.to_string()).unwrap();
                            read_groups.insert(read_name.clone(), sm.clone());
                        } 
                    }
                    Err(_e) => ()
                }

                let hp = get_haplotag_from_record(&a.record());

                if show_haplotype > 0 {
                    if hp.is_none() || show_haplotype != hp.unwrap() {
                        continue;
                    }
                }

                let b = if let Some(qpos) = a.qpos() {
                    a.record().seq()[qpos] as char
                } else {
                    '-'
                };

                let alleles = read_alleles.entry( read_name ).or_insert(empty_alleles.clone());
                alleles[index] = b;
            }
        }
    }

    print!("read_name\tread_group");
    for (c, p) in positions {
        print!("\t{c}:{}", p + 1);
    }
    print!("\n");

    let no_rg = "(none)".to_string();
    for (k, v) in read_alleles.iter() {
        let rg = read_groups.get(k).unwrap_or(&no_rg);

        print!("{k}\t{rg}");
        for a in v {
            print!("\t{a}");
        }
        print!("\n");
    }
}
