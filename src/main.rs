//---------------------------------------------------------
// Copyright 2022 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
extern crate approx;

use rust_htslib::{bam, faidx, bam::Read, bam::record::Aux};
use rust_htslib::bcf::{Reader as BcfReader, Read as BcfRead, Writer as BcfWriter, Header as BcfHeader, Format as BcfFormat, record::GenotypeAllele};
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

mod utility;
use crate::utility::{ReadHaplotypeCache, GenomeRegions, get_haplotag_from_record, get_phase_set_from_record};

mod longshot_realign;
use longshot_realign::{ReadHaplotypeLikelihood, ReadMetadata};

mod calling_models;
use crate::calling_models::{CallingModel, calling_model_to_str, str_to_calling_model};

use longshot::util::{parse_region_string};
use longshot::realignment::AlignmentType;
use longshot::extract_fragments::{ExtractFragmentParameters, extract_fragments};
use longshot::estimate_alignment_parameters::estimate_alignment_parameters;
use longshot::genotype_probs::GenotypePriors;
use longshot::call_potential_snvs::call_potential_snvs;
use longshot::context_model::ContextModel;
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
        let model_str = matches.value_of("model").unwrap_or("realigned-pileup");
        let model = str_to_calling_model(model_str).expect("unknown calling model");
        call_mutations(matches.value_of("input-bam").unwrap(),
                       matches.value_of("region").unwrap(),
                       matches.value_of("genome").unwrap(),
                       model)
    } else if let Some(matches) = matches.subcommand_matches("genotype-hets") {
        genotype_hets(matches.value_of("input-bam").unwrap(),
                      matches.value_of("region").unwrap(),
                      matches.value_of("candidates").unwrap(),
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

    let params = ModelParameters::defaults();
        
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


fn call_mutations(input_bam: &str, region_str: &str, reference_genome: &str, model: CallingModel) {

    let params = ModelParameters::defaults();

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

    let min_mapq = 50;
    let max_cigar_indel = 20;
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
    let context_model = ContextModel::init(5);
    
    let extract_fragment_parameters = ExtractFragmentParameters {
        min_mapq: min_mapq,
        alignment_type: AlignmentType::ForwardAlgorithmNumericallyStableWithContext,
        //alignment_type: AlignmentType::ForwardAlgorithmNumericallyStable,
        context_model: context_model,
        band_width: 20,
        anchor_length: 6,
        variant_cluster_max_size: 3,
        max_window_padding: 50,
        max_cigar_indel: max_cigar_indel as usize,
        store_read_id: true
    };

    // estimate parameters for longshot HMM
    let alignment_parameters = 
        estimate_alignment_parameters(&input_bam.to_owned(), 
                                      &reference_genome.to_owned(), 
                                      &None, 
                                      60, 
                                      max_cigar_indel, 
                                      100000, 
                                      &extract_fragment_parameters.context_model).unwrap();

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
            &genotype_priors,
            20,
            200,
            min_variant_observations, //potential_snv_min_alt_count,
            0.1, //potential_snv_min_alt_frac,
            50,
            alignment_parameters.ln(),
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
                                  &extract_fragment_parameters,
                                  &alignment_parameters).unwrap();
    
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
        let cm = &extract_fragment_parameters.context_model;
        let context_probs = &alignment_parameters.context_emission_probs.probs;
        let half_k = cm.k / 2;
        
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

fn genotype_hets(input_bam: &str, region_str: &str, candidates_vcf: &str, reference_genome: &str) {

    let min_mapq = 50;
    let max_cigar_indel = 20;
    let min_allele_call_qual = LogProb::from(Prob(0.1));
    let hard_min_depth = 20;
    let min_informative_fraction = 0.75;

    let max_strand_bias = 10.0;
    let context_model = ContextModel::init(5);
    
    let gt_hom_alt = [GenotypeAllele::Unphased(1), GenotypeAllele::Unphased(1)];
    let gt_hom_ref = [GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(0)];
    let gt_het = [GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(1)];
    let gt_nocall = [GenotypeAllele::UnphasedMissing, GenotypeAllele::UnphasedMissing];
    
    let extract_fragment_parameters = ExtractFragmentParameters {
        min_mapq: min_mapq,
        alignment_type: AlignmentType::ForwardAlgorithmNumericallyStableWithContext,
        //alignment_type: AlignmentType::ForwardAlgorithmNumericallyStable,
        context_model: context_model,
        band_width: 20,
        anchor_length: 6,
        variant_cluster_max_size: 3,
        max_window_padding: 50,
        max_cigar_indel: max_cigar_indel as usize,
        store_read_id: true
    };

    // estimate parameters for longshot HMM
    let alignment_parameters = 
        estimate_alignment_parameters(&input_bam.to_owned(), 
                                      &reference_genome.to_owned(), 
                                      &None, 
                                      60, 
                                      max_cigar_indel, 
                                      100000, 
                                      &extract_fragment_parameters.context_model).unwrap();

    let region = parse_region_string(Some(region_str), &input_bam.to_owned()).expect("Could not parse region").unwrap();
    
    let mut bam = bam::IndexedReader::from_path(input_bam).unwrap();
    let header = bam::Header::from_template(bam.header());
    let header_view = bam::HeaderView::from_header(&header);
    
    // set up header and faidx for pulling reference sequence
    let faidx = faidx::Reader::from_path(reference_genome).expect("Could not read reference genome:");
    let chromosome_length = header_view.target_len(region.tid).unwrap() as usize;
    let mut chromosome_sequence = faidx.fetch_seq_string(&region.chrom, 0, chromosome_length).unwrap();
    chromosome_sequence.make_ascii_uppercase();
    //let chromosome_bytes = chromosome_sequence.as_bytes();
    
    // set up vcf output
    let mut vcf_header = BcfHeader::new();
    vcf_header.push_record(b"##source=smrest genotype-hets");

    // add contig lines to header
    for tid in 0..header_view.target_count() {
        let l = format!("##contig=<ID={},length={}", std::str::from_utf8(header_view.tid2name(tid)).unwrap(), 
                                                     header_view.target_len(tid).unwrap());
        vcf_header.push_record(l.as_bytes());
    }

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
                                  &extract_fragment_parameters,
                                  &alignment_parameters).unwrap();
    
    // populate read metadata
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
            haplotype_index: None,
            phase_set: None,
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

fn phase(input_bam: &str, region_str: &str, candidates_vcf: &str, reference_genome: &str) {

    let min_mapq = 50;
    let max_cigar_indel = 20;
    let min_allele_call_qual = LogProb::from(Prob(0.1));
    let hard_min_depth = 20;
    let min_informative_fraction = 0.75;

    let max_strand_bias = 10.0;
    let context_model = ContextModel::init(5);
    
    
    let extract_fragment_parameters = ExtractFragmentParameters {
        min_mapq: min_mapq,
        alignment_type: AlignmentType::ForwardAlgorithmNumericallyStableWithContext,
        //alignment_type: AlignmentType::ForwardAlgorithmNumericallyStable,
        context_model: context_model,
        band_width: 20,
        anchor_length: 6,
        variant_cluster_max_size: 3,
        max_window_padding: 50,
        max_cigar_indel: max_cigar_indel as usize,
        store_read_id: true
    };

    // estimate parameters for longshot HMM
    let alignment_parameters = 
        estimate_alignment_parameters(&input_bam.to_owned(), 
                                      &reference_genome.to_owned(), 
                                      &None, 
                                      60, 
                                      max_cigar_indel, 
                                      100000, 
                                      &extract_fragment_parameters.context_model).unwrap();

    let region = parse_region_string(Some(region_str), &input_bam.to_owned()).expect("Could not parse region").unwrap();
    
    let mut bam = bam::IndexedReader::from_path(input_bam).unwrap();
    let header = bam::Header::from_template(bam.header());
    let header_view = bam::HeaderView::from_header(&header);
    
    // set up header and faidx for pulling reference sequence
    let faidx = faidx::Reader::from_path(reference_genome).expect("Could not read reference genome:");
    let chromosome_length = header_view.target_len(region.tid).unwrap() as usize;
    let mut chromosome_sequence = faidx.fetch_seq_string(&region.chrom, 0, chromosome_length).unwrap();
    chromosome_sequence.make_ascii_uppercase();
    
    // set up vcf output
    let mut vcf_header = BcfHeader::new();
    vcf_header.push_record(b"##source=smrest genotype-hets");

    // add contig lines to header
    for tid in 0..header_view.target_count() {
        let l = format!("##contig=<ID={},length={}", std::str::from_utf8(header_view.tid2name(tid)).unwrap(), 
                                                     header_view.target_len(tid).unwrap());
        vcf_header.push_record(l.as_bytes());
    }

    // add info lines
    vcf_header.push_record(r#"##INFO=<ID=TotalDepth,Number=1,Type=Integer,Description="Total number of reads aligned across this position">"#.as_bytes());
    vcf_header.push_record(r#"##INFO=<ID=LP,Number=2,Type=Float,Description="todo">"#.as_bytes());
    vcf_header.push_record(r#"##INFO=<ID=EP,Number=1,Type=Float,Description="todo">"#.as_bytes());
    vcf_header.push_record(r#"##INFO=<ID=StrandBias,Number=1,Type=Float,Description="Strand bias p-value (PHRED scaled)">"#.as_bytes());
    vcf_header.push_record(r#"##INFO=<ID=InformativeFraction,Number=1,Type=Float,Description="Proportion of reads that could be used in allele calling">"#.as_bytes());
    vcf_header.push_record(r#"##FORMAT=<ID=GT,Number=1,Type=String,Description="sample genotype">"#.as_bytes());
    vcf_header.push_record(r#"##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths (high-quality bases)">"#.as_bytes());
    vcf_header.push_sample("sample".as_bytes());

    let all_vars = parse_vcf_potential_variants(&candidates_vcf.to_owned(), &input_bam.to_owned())
                        .expect("Could not parse candidate variants file");
    
    let region_vars: Vec<Var> = all_vars.get_variants_range(region.clone()).expect("Could not subset vars");
    let mut varlist = VarList::new(region_vars, all_vars.target_names).expect("Could not construct var list");

    let frags = extract_fragments(&input_bam.to_owned(),
                                  &reference_genome.to_owned(),
                                  &mut varlist,
                                  &Some(region.clone()),
                                  &extract_fragment_parameters,
                                  &alignment_parameters).unwrap();
    
    // populate read metadata
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
            haplotype_index: None,
            phase_set: None,
            strand_index: record.is_reverse() as i32,
            leading_softclips: record.cigar().leading_softclips(),
            trailing_softclips: record.cigar().trailing_softclips()
        };
        read_meta.insert(s.clone(), rm);
    }

    // Convert fragments into read-haplotype likelihoods for every variant
    let mut rhl_per_var = vec![ Vec::<ReadHaplotypeLikelihood>::new(); varlist.lst.len() ];

    for f in &frags {
        let id = f.id.as_ref().unwrap().clone();

        if let Some( rm ) = read_meta.get(&id) {
            for c in &f.calls {
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
