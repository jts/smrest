use rust_htslib::{bam, faidx, bam::Read, bam::record::Aux};
use rust_htslib::bcf::{Reader as BcfReader, Read as BcfRead};
use std::collections::{HashSet, HashMap};
use clap::{App, SubCommand, Arg, value_t};
use itertools::Itertools;
use rand::prelude::*;
use statrs::distribution::{Beta, Poisson, Binomial, Discrete, ContinuousCDF};

mod pileup_stats;
use crate::pileup_stats::PileupStats;
use crate::pileup_stats::base2index;

pub struct ModelParameters {
    mutation_rate: f64,
    heterozygosity: f64,
    ccf_dist: Beta,
    depth_dist: Poisson,
    purity: f64,
    error_rate: f64
}

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
                    .help("the reference genome")))
        .subcommand(SubCommand::with_name("extract")
                .about("gather candidate somatic mutations from aligned reads")
                .arg(Arg::with_name("genome")
                    .short('g')
                    .long("genome")
                    .takes_value(true)
                    .help("the reference genome"))
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
        extract_mutations(matches.value_of("input-bam").unwrap(),
                          matches.value_of("genome").unwrap());
    } else if let Some(matches) = matches.subcommand_matches("error-contexts") {
        error_contexts(matches.value_of("input-bam").unwrap(),
                       matches.value_of("genome").unwrap(),
                       matches.value_of("vcf").unwrap());
    } else if let Some(matches) = matches.subcommand_matches("sim-pileup") {

        let depth_lambda = value_t!(matches, "depth", f64).unwrap_or(50.0);
        let purity = value_t!(matches, "purity", f64).unwrap_or(0.75);

        let params = ModelParameters { 
            mutation_rate: 5.0 / 1000000.0, // per haplotype
            heterozygosity: 1.0 / 2000.0, // per haplotype
            ccf_dist: Beta::new(9.0, 1.0).unwrap(),
            depth_dist: Poisson::new(depth_lambda).unwrap(),
            purity: purity,
            error_rate: 0.02
        };
        sim_pileup(& params, matches.value_of("genome").unwrap());
    }
}

// A cache storing the haplotype tag for every read processed
// We need this because bam_get_aux is O(N)
pub struct ReadHaplotypeCache
{
    cache: HashMap<String, i32>
}

impl ReadHaplotypeCache
{
    fn key(record: &bam::Record) -> String {
        let s = String::from_utf8(record.qname().to_vec()).unwrap();
        return s;
    }

    pub fn update(&mut self, key: &String, record: &bam::Record) -> Option<i32> {
        if let Some(hi) = get_haplotag_from_record(record) {
            self.cache.insert(key.clone(), hi);
            return Some(hi);
        } else {
            return None;
        }
    }

    pub fn get(&mut self, record: &bam::Record) -> Option<i32> {
        let s = ReadHaplotypeCache::key(record);
        let x = self.cache.get(&s);
        match x {
            Some(value) => return Some(*value),
            None => return self.update(&s, record)
        }
    }
}


fn get_haplotag_from_record(record: &bam::Record) -> Option<i32> {
    match record.aux(b"HP") {
        Ok(value) => {
            if let Aux::I32(v) = value {
                return Some(v)
            } else {
                return None
            }
        }
        Err(_e) => return None
    }
}

fn extract_mutations(input_bam: &str, reference_genome: &str) {

    let params = ModelParameters { 
        mutation_rate: 5.0 / 1000000.0, // per haplotype
        heterozygosity: 1.0 / 2000.0, // per haplotype
        ccf_dist: Beta::new(9.0, 1.0).unwrap(),
        depth_dist: Poisson::new(50.0).unwrap(),
        purity: 0.75,
        error_rate: 0.05 // R9.4 conservative
    };


    let mut bam = bam::IndexedReader::from_path(input_bam).unwrap();
    println!("chromosome\tposition\treference_base\tvariant_base\tcanonical_type\tcanonical_context\taligned_depth\tmean_mapq\t\
              hmajor_variant_count\thminor_variant_count\thmajor_vaf\thminor_vaf\th1_a\th1_c\th1_g\th1_t\th2_a\th2_c\th2_g\th2_t\t\
              p_ref\tp_germline\tp_somatic");

    let all_positions = false;
    let chromosome_name = "chr20";
    let min_variant_observations = 3;
    let max_variant_minor_observations = 1;
    let min_variant_observations_per_strand = 1;

    // set up header and faidx for pulling reference sequence
    let faidx = faidx::Reader::from_path(reference_genome).expect("Could not read reference genome:");

    let header = bam::Header::from_template(bam.header());
    let header_view = bam::HeaderView::from_header(&header);
    let tid = header_view.tid(chromosome_name.as_bytes()).unwrap();
    let chromosome_length = header_view.target_len(tid).unwrap() as usize;
    let mut chromosome_sequence = faidx.fetch_seq_string(chromosome_name, 0, chromosome_length).unwrap();
    chromosome_sequence.make_ascii_uppercase();
    let chromosome_bytes = chromosome_sequence.as_bytes();
    
    let mut cache = ReadHaplotypeCache { cache: HashMap::new() };
    let mut ps = PileupStats::new();

    // go to chromosome of interest
    bam.fetch(chromosome_name).unwrap();
    for p in bam.pileup() {
        let pileup = p.unwrap();
        
        //println!("processing {}", pileup.pos() + 1);
        ps.fill_pileup(&mut cache, pileup.alignments());

        let reference_position = pileup.pos() as usize;
        let reference_base = chromosome_bytes[reference_position] as char;
        if reference_base == 'N' {
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
        let minor_variant_count = ps.get_count_on_haplotype(max_variant_index, minor_haplotype_index);

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

fn error_contexts(input_bam: &str, reference_genome: &str, variants_vcf: &str) {

    let chromosome_name = String::from("chr20");

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

        let reference_context_fwd = chromosome_bytes[ (reference_position - 2)..(reference_position + 3)].to_vec();
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

                let (c, mut eb)  = match a.record().is_reverse() {
                    false => (&context_str_fwd, read_base), 
                    true => (&context_str_rev, bio::alphabets::dna::complement(read_base as u8) as char)
                };
                eb = 'N';
                
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

fn rand_base_not(rng: &mut ThreadRng, curr_base: u8) -> u8 {
    let bases = [ 'A', 'C', 'G', 'T' ];
    let mut base = curr_base;
    while base == curr_base {
        base = bases[ rng.gen_range(0..4) ] as u8;
    }
    return base;
}

fn mutate_haplotypes(haplotypes: &mut [ Vec<u8> ], rate: f64) -> () {
    // Set up a generator to select random bases
    let mut rng = rand::thread_rng();

    // generate variants
    for hi in 0..2 {
        let mut _var_count = 0;
        for i in 0..haplotypes[hi].len() {
            if rng.gen::<f64>() < rate {
                let curr_base = haplotypes[hi][i];
                haplotypes[hi][i] = rand_base_not(&mut rng, curr_base);
                _var_count += 1;
            }
        }
        //eprintln!("Generate {_var_count} variants for haplotype");
    }
}

fn calculate_class_probabilities_phased(alt_count: u64, ref_count: u64, params: &ModelParameters) -> [f64;3]
{
    let depth = ref_count + alt_count;
    //
    // P(somatic | data) = P(data | somatic) P(somatic) / sum_class P(data | class )
    //
    // P(data | ref) = Binom(alt_count, ref_count + alt_count, error_rate)
    let p_data_ref = Binomial::new(params.error_rate, depth).unwrap().pmf(alt_count);

    // P(data | het) = Binom(alt_count, ref_count + alt_count, 1 - error_rate)
    let p_data_het = Binomial::new(1.0 - params.error_rate, depth).unwrap().pmf(alt_count);
    
    // P(data | somatic) = sum_c Binom(alt_count, ref_count + alt_count, purity * c * (1 - error_rate) + (1 - purity*c) * error_rate ) P(c)
    let bins = 10;
    let step = 1.0 / 10 as f64;
    let mut p_data_somatic = 0.0;
    for i in 0..bins {
        let start = f64::from(i) * step;
        let end = f64::from(i + 1) * step;
        let c = (end + start) / 2.0;
        let p_c = params.ccf_dist.cdf(end) - params.ccf_dist.cdf(start);

        let p_read_from_mutated_haplotype = params.purity * c;
        let t1 = p_read_from_mutated_haplotype * (1.0 - params.error_rate);
        let t2 = (1.0 - p_read_from_mutated_haplotype) * params.error_rate;

        let p_data_somatic_at_c = Binomial::new(t1 + t2, depth).unwrap().pmf(alt_count);
        p_data_somatic += p_data_somatic_at_c * p_c;
    }

    // priors
    let p_het = params.heterozygosity;
    let p_somatic = params.mutation_rate;
    let p_ref = 1.0 - p_het - p_somatic;

    let sum = p_data_ref * p_ref + p_data_het * p_het + p_data_somatic * p_somatic;
    
    let p_ref_data = p_data_ref * p_ref / sum;
    let p_het_data = p_data_het * p_het / sum;
    let p_somatic_data = p_data_somatic * p_somatic / sum;

    return [ p_ref_data, p_het_data, p_somatic_data ];
}

fn calculate_class_probabilities_unphased(alt_count: u64, ref_count: u64, params: &ModelParameters) -> [f64;3]
{
    let depth = ref_count + alt_count;
    //
    // P(somatic | data) = P(data | somatic) P(somatic) / sum_class P(data | class )
    //

    // P(data | ref) = Binom(alt_count, ref_count + alt_count, error_rate)
    let p_data_ref = Binomial::new(params.error_rate, depth).unwrap().pmf(alt_count);

    // P(data | het) = Binom(alt_count, ref_count + alt_count, 1 - error_rate)
    let p_data_het = Binomial::new(0.5 * (1.0 - params.error_rate) + 0.5 * params.error_rate, depth).unwrap().pmf(alt_count);
    
    // P(data | somatic) = sum_c Binom(alt_count, ref_count + alt_count, purity * c * (1 - error_rate) + (1 - purity*c) * error_rate ) P(c)
    let bins = 10;
    let step = 1.0 / 10 as f64;
    let mut p_data_somatic = 0.0;
    for i in 0..bins {
        let start = f64::from(i) * step;
        let end = f64::from(i + 1) * step;
        let c = (end + start) / 2.0;
        let p_c = params.ccf_dist.cdf(end) - params.ccf_dist.cdf(start);

        let p_read_from_mutated_haplotype = params.purity * c;
        let t1 = p_read_from_mutated_haplotype * (1.0 - params.error_rate);
        let t2 = (1.0 - p_read_from_mutated_haplotype) * params.error_rate;

        let p_data_somatic_at_c = Binomial::new(t1 + t2, depth).unwrap().pmf(alt_count);
        p_data_somatic += p_data_somatic_at_c * p_c;
    }

    // priors
    let p_het = params.heterozygosity;
    let p_somatic = params.mutation_rate;
    let p_ref = 1.0 - p_het - p_somatic;

    let sum = p_data_ref * p_ref + p_data_het * p_het + p_data_somatic * p_somatic;
    
    let p_ref_data = p_data_ref * p_ref / sum;
    let p_het_data = p_data_het * p_het / sum;
    let p_somatic_data = p_data_somatic * p_somatic / sum;

    return [ p_ref_data, p_het_data, p_somatic_data ];
}

pub struct SimulationStats
{
    model_name: String,
    class_posterior_sums: [f64; 3],
    num_somatic_calls: usize,
    num_somatic_true: usize,
    num_somatic_correct: usize,
    num_het_true: usize

}

impl SimulationStats {

    pub fn new(name: String) -> SimulationStats {
        let r = SimulationStats { model_name: name, 
                                 class_posterior_sums: [0.0; 3], 
                                 num_somatic_calls: 0,
                                 num_somatic_true: 0,
                                 num_somatic_correct: 0,
                                 num_het_true: 0 };
        return r;
    }

    pub fn update(& mut self, is_het: bool, is_somatic: bool, class_posterior: &[f64; 3]) -> () {
        for i in 0..class_posterior.len() {
            self.class_posterior_sums[i] += class_posterior[i];
        }

        let mut p_max = 0.0;
        let mut idx_max = 0;
        for k in 0..class_posterior.len() {
            if class_posterior[k] > p_max {
                p_max = class_posterior[k];
                idx_max = k;
            }
        }

        let is_somatic_call = idx_max == 2;
        if is_somatic_call {
            self.num_somatic_calls += 1;
        }

        if is_somatic_call && is_somatic {
            self.num_somatic_correct += 1;
        }

        if is_somatic { 
            self.num_somatic_true += 1;
        }

        if is_het { 
            self.num_het_true += 1;
        }
    }

    pub fn print(& self, params: & ModelParameters) {
        let sens = self.num_somatic_correct as f64 / self.num_somatic_true as f64;
        let prec = self.num_somatic_correct as f64 / self.num_somatic_calls as f64;
        let f1 = 2.0 * sens * prec / (sens + prec);
        let ccf_mean = params.ccf_dist.shape_a() / (params.ccf_dist.shape_a() + params.ccf_dist.shape_b());

        println!("{}\t{:.3}\t{:.3}\t{:.3}\t{}\t{}\t{:.1}\t{:.1}\t{:.1}\t{}\t{:.3}\t{:.3}\t{:.3}",
            self.model_name, params.purity, params.depth_dist.lambda(), ccf_mean,
            self.num_somatic_true, self.num_het_true, 
            self.class_posterior_sums[0], self.class_posterior_sums[1], self.class_posterior_sums[2],
            self.num_somatic_calls, sens, prec, f1);
    }
}

fn sim_pileup(params: & ModelParameters, reference_genome: &str)
{
    let summary = true;

    let mut rng = rand::thread_rng();

    // get genome sequence
    let chromosome_name = String::from("chr20");
    //let chromosome_length = 60000000;
    let chromosome_length = 6000000;

    let faidx = faidx::Reader::from_path(reference_genome).expect("Could not read reference genome:");
    let mut chromosome_sequence = faidx.fetch_seq_string(chromosome_name.as_str(), 0, chromosome_length).unwrap();
    chromosome_sequence.make_ascii_uppercase();
    let chromosome_bytes = chromosome_sequence.as_bytes();
    let mut germline_haplotypes = [ chromosome_bytes.to_vec(), chromosome_bytes.to_vec() ];
    mutate_haplotypes(& mut germline_haplotypes, params.heterozygosity);
    
    let mut somatic_haplotypes = [ germline_haplotypes[0].clone(), germline_haplotypes[1].clone() ];
    mutate_haplotypes(& mut somatic_haplotypes, params.mutation_rate);

    // structs to collect results
    let mut phased_model_stats = SimulationStats::new("phased".to_string());
    let mut unphased_model_stats = SimulationStats::new("unphased".to_string());

    // Finally, simulate some pileup data at each position on each haplotype
    if !summary {
        println!("model\tposition\thaplotype\tclass\tis_het\tis_somatic\tref_count\talt_count\thvaf\tp_ref\tp_het\tp_somatic");
    }

    for j in 0..chromosome_length {

        if chromosome_bytes[j] as char == 'N' {
            continue;
        }

        // depth at this position
        let depth = params.depth_dist.sample(&mut rng) as u64;
        let haplotype_dist = Binomial::new(0.5, depth).unwrap();
        let h0_depth = haplotype_dist.sample(&mut rng) as u64;
        let haplotype_depths = [ h0_depth, depth - h0_depth ];
        
        // simulate a pileup
        let mut ps = PileupStats::new();
        for i in 0..2 {
            let is_somatic = somatic_haplotypes[i][j] != germline_haplotypes[i][j];
            let ccf = match is_somatic {
                true => params.ccf_dist.sample(&mut rng),
                false => 0.0
            };

            // simulate each read
            for _r in 0..haplotype_depths[i] {
            
                let strand: u32 = rng.gen_range(0..2);

                // determine base for this read
                let prob_sampled_from_somatic = match is_somatic {
                    true => params.purity * ccf,
                    false => 0.0
                };

                let is_tumor_read = rng.gen::<f64>() < prob_sampled_from_somatic;
                let mut read_base = match is_tumor_read {
                    true => somatic_haplotypes[i][j],
                    false => germline_haplotypes[i][j]
                };

                // corrupt read if sequencing error
                if rng.gen::<f64>() < params.error_rate {
                    read_base = rand_base_not(&mut rng, read_base);
                }

                ps.increment(base2index(read_base as char) as u32, i as u32, strand);
            }
        }

        let ref_base_index = base2index(chromosome_bytes[j] as char) as u32;

        // run phased model, collect stats
        for i in 0usize..2usize {
            let alt_base_index = ps.get_max_nonref_base_on_haplotype(i as u32, ref_base_index as u32);
            let alt_count = ps.get_count_on_haplotype(alt_base_index, i as u32) as u64;
            let ref_count = ps.get_count_on_haplotype(ref_base_index, i as u32) as u64;

            let probs = calculate_class_probabilities_phased(alt_count, ref_count, &params);

            // update summary stats
            let is_het = germline_haplotypes[i][j] != chromosome_bytes[j];
            let is_somatic = somatic_haplotypes[i][j] != germline_haplotypes[i][j];
            phased_model_stats.update(is_het, is_somatic, &probs);

            // output per-record calls, if wanted
            if !summary {
                let hvaf = alt_count as f64 / ps.get_haplotype_depth(i as u32) as f64; 
                let mut class = "REF";
                if is_somatic { 
                    class = "SOMATIC"
                } else if is_het {
                    class = "HET";   
                }

                println!("phased\t{j}\t{i}\t{class}\t{}\t{}\t\
                         {ref_count}\t{alt_count}\t{hvaf:.3}\t{:.3}\t{:.3}\t{:.3}",
                         is_het as u8, is_somatic as u8, probs[0], probs[1], probs[2]);
            }
        }

        // run baseline non-phased model
        {
            let alt_base_index = ps.get_max_nonref_base_unphased(ref_base_index as u32);
            let alt_count = ps.get_count_unphased(alt_base_index) as u64;
            let ref_count = ps.get_count_unphased(ref_base_index) as u64;
            let probs = calculate_class_probabilities_unphased(alt_count, ref_count, &params);
            
            // update summary stats
            let is_het = germline_haplotypes[0][j] != chromosome_bytes[j] || germline_haplotypes[1][j] != chromosome_bytes[j];
            let is_somatic = somatic_haplotypes[0][j] != germline_haplotypes[0][j] || somatic_haplotypes[1][j] != germline_haplotypes[1][j];
            unphased_model_stats.update(is_het, is_somatic, &probs);
            
            // output per-record calls, if wanted
            if !summary {
                let hvaf = alt_count as f64 / ps.get_depth() as f64; 
                let mut class = "REF";
                if is_somatic { 
                    class = "SOMATIC"
                } else if is_het {
                    class = "HET";   
                }

                println!("unphased\t{j}\tboth\t{class}\t{}\t{}\t\
                         {ref_count}\t{alt_count}\t{hvaf:.3}\t{:.3}\t{:.3}\t{:.3}",
                         is_het as u8, is_somatic as u8, probs[0], probs[1], probs[2]);
            }
        }
    }

    if summary {
        println!("model\tpurity\tmean_coverage\tmean_ccf\ttrue_somatic_positions\ttrue_het_position\testimated_ref_bases\t\
                  estimated_het_bases\testimated_mutated_bases\tnum_somatic_calls\tsomatic_call_sensitivity\tsomatic_call_precision\tf1");
        phased_model_stats.print(& params);
        unphased_model_stats.print(& params);
    }
}
