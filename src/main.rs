//---------------------------------------------------------
// Copyright 2022 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
extern crate approx;

use rust_htslib::{bam, faidx, bam::Read};
use rust_htslib::bcf::{Reader as BcfReader, Read as BcfRead, Writer as BcfWriter, Header as BcfHeader, Format as BcfFormat};
use std::collections::{HashSet, HashMap};
use clap::{App, SubCommand, Arg, value_t};
use itertools::Itertools;
use statrs::distribution::{Beta, Poisson};
use fishers_exact::fishers_exact;
use bio::stats::{Prob, PHREDProb};

mod pileup_stats;
use crate::pileup_stats::PileupStats;
use crate::pileup_stats::base2index;

mod simulation;
use crate::simulation::*;

mod classifier;
use crate::classifier::*;

mod utility;
use crate::utility::{ReadHaplotypeCache, GenomeRegions};

mod longshot_realign;
use longshot_realign::recalculate_pileup;

use longshot::util::{dna_vec, parse_region_string, parse_target_names};
use longshot::realignment::AlignmentType;
use longshot::extract_fragments::{ExtractFragmentParameters};
use longshot::estimate_alignment_parameters::estimate_alignment_parameters;

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
        call_mutations(matches.value_of("input-bam").unwrap(),
                       matches.value_of("region").unwrap(),
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
        let ccf = CancerCellFraction { p_clonal: p_clonal, subclonal_ccf: Beta::new(2.0, 2.0).unwrap() };
        let mutation_output = matches.value_of("mutation-output").unwrap_or("simulated_mutations.tsv");

        let params = ModelParameters { 
            mutation_rate: muts_per_megabase / 1000000.0, // per haplotype
            heterozygosity: 1.0 / 2000.0, // per haplotype
            ccf_dist: ccf,
            depth_dist: Some(Poisson::new(depth_lambda).unwrap()),
            purity: purity,
            error_rate: error_rate
        };
        sim_pileup(&params, matches.value_of("genome").unwrap(), per_site_output, &mutation_output);
    }
}

fn extract_mutations(input_bam: &str, reference_genome: &str, min_depth: u32, region_opt: Option<&str>) {

    let ccf = CancerCellFraction { p_clonal: 0.75, subclonal_ccf: Beta::new(2.0, 2.0).unwrap() };
    let params = ModelParameters { 
        mutation_rate: 5.0 / 1000000.0, // per haplotype
        heterozygosity: 1.0 / 2000.0, // per haplotype
        ccf_dist: ccf,
        depth_dist: None,
        purity: 0.75,
        error_rate: 0.05
    };
        
    let mut bam = bam::IndexedReader::from_path(input_bam).unwrap();
    let header = bam::Header::from_template(bam.header());
    let header_view = bam::HeaderView::from_header(&header);
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

    // set up header and faidx for pulling reference sequence
    let faidx = faidx::Reader::from_path(reference_genome).expect("Could not read reference genome:");
    let tid = header_view.tid(chromosome_name.as_bytes()).unwrap();
    let chromosome_length = header_view.target_len(tid).unwrap() as usize;
    let mut chromosome_sequence = faidx.fetch_seq_string(chromosome_name, 0, chromosome_length).unwrap();
    chromosome_sequence.make_ascii_uppercase();
    let chromosome_bytes = chromosome_sequence.as_bytes();
    
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
}


fn call_mutations(input_bam: &str, region_str: &str, reference_genome: &str) {

    let ccf = CancerCellFraction { p_clonal: 0.75, subclonal_ccf: Beta::new(2.0, 2.0).unwrap() };
    let params = ModelParameters { 
        mutation_rate: 5.0 / 1000000.0, // per haplotype
        heterozygosity: 1.0 / 2000.0, // per haplotype
        ccf_dist: ccf,
        depth_dist: None,
        purity: 0.75,
        error_rate: 0.02
    };

    /*
    let hom_snv_rate = LogProb::from(Prob(0.0005));
    let het_snv_rate = LogProb::from(Prob(0.001));
    let hom_indel_rate = LogProb::from(Prob(0.00005));
    let het_indel_rate = LogProb::from(Prob(0.00001));
    let ts_tv_ratio = 2.0 * 0.5; // from longshot defaults...

    let genotype_priors = GenotypePriors::new(
        hom_snv_rate,
        het_snv_rate,
        hom_indel_rate,
        het_indel_rate,
        ts_tv_ratio,
    ).unwrap();
    */

    let bases = "ACGT";
    let min_mapq = 50;
    let max_cigar_indel = 20;
    let min_variant_observations = 3;
    let min_p_somatic = 0.0001;
    let min_obs_for_realign = 3;

    //let max_variant_minor_observations = 1;
    let min_variant_observations_per_strand = 1;
    let t_names = parse_target_names(&input_bam.to_owned()).unwrap();
    
    let extract_fragment_parameters = ExtractFragmentParameters {
        min_mapq: min_mapq,
        alignment_type: AlignmentType::ForwardAlgorithmNonNumericallyStable,
        band_width: 20,
        anchor_length: 6,
        variant_cluster_max_size: 3,
        max_window_padding: 50,
        max_cigar_indel: max_cigar_indel as usize,
        store_read_id: true
    };

    // estimate parameters for longshot HMM
    let alignment_parameters = estimate_alignment_parameters(&input_bam.to_owned(), &reference_genome.to_owned(), &None, 60, max_cigar_indel, 100000).unwrap();

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
    let chromosome_char_vec = dna_vec(chromosome_bytes);
    
    let mut cache = ReadHaplotypeCache { cache: HashMap::new() };
    let mut ps = PileupStats::new();

    // set up vcf output
    let mut vcf_header = BcfHeader::new();

    // add contig lines to header
    for tid in 0..header_view.target_count() {
        let l = format!("##contig=<ID={},length={}", std::str::from_utf8(header_view.tid2name(tid)).unwrap(), 
                                                     header_view.target_len(tid).unwrap());
        vcf_header.push_record(l.as_bytes());
    }

    // add info lines
    vcf_header.push_record(r#"##INFO=<ID=SomaticHaplotypeIndex,Number=1,Type=Integer,Description="Index of haplotype carrying the mutation">"#.as_bytes());
    vcf_header.push_record(r#"##INFO=<ID=HaplotypeAltCount,Number=2,Type=Integer,Description="Observed alt read count on each haplotype">"#.as_bytes());
    vcf_header.push_record(r#"##INFO=<ID=HaplotypeDepth,Number=2,Type=Integer,Description="Observed read count on each haplotype">"#.as_bytes());
    vcf_header.push_record(r#"##INFO=<ID=HaplotypeVAF,Number=2,Type=Float,Description="VAF on each haplotype">"#.as_bytes());
    vcf_header.push_record(r#"##INFO=<ID=MutationType,Number=1,Type=String,Description="Mutation in canonical form">"#.as_bytes());
    vcf_header.push_record(r#"##INFO=<ID=SequenceContext,Number=1,Type=String,Description="Sequence context on reference genome">"#.as_bytes());
    vcf_header.push_record(r#"##INFO=<ID=StrandBias,Number=1,Type=Float,Description="Strand bias p-value (PHRED scaled)">"#.as_bytes());
    vcf_header.push_record(r#"##INFO=<ID=ProportionPhased,Number=1,Type=Float,Description="Proportion of reads at this postion successfully phased">"#.as_bytes());

    //vcf_header.push_sample("sample".as_bytes());
    let mut vcf = BcfWriter::from_stdout(&vcf_header, true, BcfFormat::Vcf).unwrap();

    bam.fetch( (region.tid, region.start_pos, region.end_pos + 1) ).unwrap();

    for p in bam.pileup() {
        let pileup = p.unwrap();
        
        // check whether we should bother with this position
        let reference_tid = pileup.tid();
        let reference_position = pileup.pos();
        
        // bam.fetch will include all reads that cross a single base in the ROI
        // the pileup will produce each column for any read in this region
        // we subset the pileup down to our ROI here
        if reference_position < region.start_pos || reference_position > region.end_pos {
            continue;
        } 

        let reference_base = chromosome_bytes[reference_position as usize] as char;
        if reference_base == 'N' {
            continue;
        }

        // Naive pileup
        ps.fill_pileup(&mut cache, pileup.alignments());
        
        // Calculate most frequently observed non-reference base on either haplotype
        let reference_base_index = base2index(reference_base) as u32;
        let (candidate_haplotype_index, candidate_variant_index) = ps.select_candidate_alt(reference_base_index, min_variant_observations_per_strand);

        let mut alt_count_on_candidate_haplotype = ps.get_count_on_haplotype(candidate_variant_index, candidate_haplotype_index);
        let variant_base = bases.as_bytes()[candidate_variant_index as usize] as char;
        
        // determine whether to re-align reads around this locus
        if alt_count_on_candidate_haplotype >= min_obs_for_realign && 
            ps.get_haplotype_depth(candidate_haplotype_index) - alt_count_on_candidate_haplotype >= min_obs_for_realign {
            ps = recalculate_pileup(&ps, 
                                    reference_tid,
                                    reference_position as usize,
                                    &chromosome_char_vec, 
                                    &t_names,
                                    reference_base, 
                                    variant_base, 
                                    &mut cache, 
                                    pileup.alignments(), 
                                    extract_fragment_parameters,
                                    alignment_parameters);
            alt_count_on_candidate_haplotype = ps.get_count_on_haplotype(candidate_variant_index, candidate_haplotype_index);
        }

        let ref_count_on_candidate_haplotype = ps.get_count_on_haplotype(reference_base_index, candidate_haplotype_index);
        let class_probs = calculate_class_probabilities_phased(alt_count_on_candidate_haplotype as u64, ref_count_on_candidate_haplotype as u64, &params);

        if reference_base != variant_base && alt_count_on_candidate_haplotype >= min_variant_observations && class_probs[2] > min_p_somatic {
            let mut record = vcf.empty_record();
            let rid = vcf.header().name2rid(region.chrom.as_bytes()).expect("Could not find reference id");
            record.set_rid(Some(rid));
            
            // bcftools and htslib both represent positions as 0-based, it will be converted to 1-base on output
            record.set_pos(reference_position as i64); 

            record.set_alleles(&[ &[reference_base as u8], &[variant_base as u8] ]).expect("Could not set alleles");

            let mut qual = *PHREDProb::from(Prob(1.0 - class_probs[2]));
            if qual > 60.0 {
                qual = 60.0;
            }
            record.set_qual(qual as f32);

            record.push_info_integer(b"SomaticHaplotypeIndex", &[candidate_haplotype_index as i32]).expect("Could not add INFO");

            let h0_ac = ps.get_count_on_haplotype(candidate_variant_index, 0) as i32;
            let h0_depth = ps.get_haplotype_depth(0) as i32;
            
            let h1_ac = ps.get_count_on_haplotype(candidate_variant_index, 1) as i32;
            let h1_depth = ps.get_haplotype_depth(1) as i32;
            
            record.push_info_integer(b"HaplotypeAltCount", &[h0_ac, h1_ac]).expect("Could not add INFO");
            record.push_info_integer(b"HaplotypeDepth", &[h0_depth, h1_depth]).expect("Could not add INFO");

            let h0_vaf = h0_ac as f32 / h0_depth as f32;
            let h1_vaf = h1_ac as f32 / h1_depth as f32;
            record.push_info_float(b"HaplotypeVAF", &[h0_vaf, h1_vaf]).expect("Could not add INFO");
            
            // grab reference context
            let reference_context = &chromosome_bytes[ (reference_position as usize - 1)..(reference_position as usize + 2)];
            
            // grab mutation type
            let mut mutation_type: [char; 3] = [ reference_base, '>', variant_base ];

            // convert mutation type/context to canonical form C>x, T>x
            if mutation_type[0] != 'C' && mutation_type[0] != 'T' {
                mutation_type[0] = bio::alphabets::dna::complement(mutation_type[0] as u8) as char;
                mutation_type[2] = bio::alphabets::dna::complement(mutation_type[2] as u8) as char;
            }

            let mutation_type_str = String::from_iter(&mutation_type);
            record.push_info_string(b"MutationType", &[mutation_type_str.as_bytes()]).expect("Could not add INFO");
            record.push_info_string(b"SequenceContext", &[&reference_context]).expect("Could not add INFO");
            
            //let mean_mapq = ps.mean_mapq;
            record.push_info_float(b"ProportionPhased", &[ps.proportion_phased]).expect("Could not add INFO");

            let fishers_result = fishers_exact(&[ ps.get( reference_base_index, candidate_haplotype_index, 0),    ps.get( reference_base_index, candidate_haplotype_index, 1),
                                            ps.get( candidate_variant_index, candidate_haplotype_index, 0), ps.get( candidate_variant_index, candidate_haplotype_index, 1) ]).unwrap();

            // straight from longshot
            let strand_bias_pvalue = if fishers_result.two_tail_pvalue <= 500.0 {
                *PHREDProb::from(Prob(fishers_result.two_tail_pvalue))
            } else {
                500.0
            };
            record.push_info_float(b"StrandBias", &[strand_bias_pvalue as f32]).expect("Could not add INFO");
            
            /*
            println!("{chromosome_name}\t{position}\t{reference_base}\t{variant_base}\t{mutation_type_str}\t{context_str}\t\
                      {aligned_depth}\t{mean_mapq}\t{proportion_phased:.3}\t{alt_count_on_candidate_haplotype}\t{alt_count_on_non_candidate_haplotype}\t\
                      {hmaj_vaf:.3}\t{hmin_vaf:.3}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.3}\t{:.3}\t{:.3}\t{:.3}",
                      h0[0], h0[1], h0[2], h0[3],
                      h1[0], h1[1], h1[2], h1[3], strand_bias_pvalue,
                      class_probs[0], class_probs[1], class_probs[2]);
            */
            vcf.write(&record).unwrap();
        }
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
    let faidx = faidx::Reader::from_path(reference_genome).expect("Could not read reference genome:");
    let mut bam = bam::IndexedReader::from_path(input_bam).unwrap();
    let header = bam::Header::from_template(bam.header());
    let header_view = bam::HeaderView::from_header(&header);
    let tid = header_view.tid(chromosome_name.as_bytes()).unwrap();
    let chromosome_length = header_view.target_len(tid).unwrap() as usize;
    let mut chromosome_sequence = faidx.fetch_seq_string(chromosome_name.as_str(), 0, chromosome_length).unwrap();
    chromosome_sequence.make_ascii_uppercase();
    let chromosome_bytes = chromosome_sequence.as_bytes();
    
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
