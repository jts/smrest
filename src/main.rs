use rust_htslib::{bam, faidx, bam::Read, bam::record::Aux};
use rust_htslib::bcf::{Reader as BcfReader, Read as BcfRead};
use std::collections::{HashSet, HashMap};
use clap::{App, SubCommand, Arg};
use itertools::Itertools;
use rand::prelude::*;
use statrs::distribution::{Beta, Poisson, Binomial, Discrete, ContinuousCDF};

fn main() {
    let matches = App::new("smrest")
        .version("0.1")
        .author("Jared Simpson <jared.simpson@oicr.on.ca>")
        .about("Toolkit for estimating somatic mutation rate from long reads")
        .subcommand(SubCommand::with_name("sim-pileup")
                .about("simulate pileup data")
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
        sim_pileup(matches.value_of("genome").unwrap());
    }
}

// A cache storing the haplotype tag for every read processed
// We need this because bam_get_aux is O(N)
struct ReadHaplotypeCache
{
    cache: HashMap<String, i32>
}

impl ReadHaplotypeCache
{
    fn key(record: &bam::Record) -> String {
        let s = String::from_utf8(record.qname().to_vec()).unwrap();
        return s;
    }

    fn update(&mut self, key: &String, record: &bam::Record) -> Option<i32> {
        if let Some(hi) = get_haplotag_from_record(record) {
            self.cache.insert(key.clone(), hi);
            return Some(hi);
        } else {
            return None;
        }
    }

    fn get(&mut self, record: &bam::Record) -> Option<i32> {
        let s = ReadHaplotypeCache::key(record);
        let x = self.cache.get(&s);
        match x {
            Some(value) => return Some(*value),
            None => return self.update(&s, record)
        }
    }
}

struct PileupStats
{
    // this is indexed by base (dim four), haplotype (two), strand (two)
    base_counts: [u32; 16],
    mean_mapq: f32
}

impl PileupStats {
    #[inline(always)]
    fn get(&self, b: u32, h: u32, s: u32) -> u32 {
        let i = (h * 8) + (s * 4) + b;
        return self.base_counts[ i as usize ];
    }

    fn increment(&mut self, b: u32, h: u32, s: u32) -> () {
        let i = (h * 8) + (s * 4) + b;
        self.base_counts[i as usize] += 1;
    }

    fn get_count_on_haplotype(&self, b: u32, h: u32) -> u32 {
        return self.get(b, h, 0) + self.get(b, h, 1);
    }

    fn get_haplotype_depth(&self, h: u32) -> u32 {
        let mut sum = 0;
        let lo:usize = (h * 8).try_into().unwrap();
        let hi:usize = ((h+1) * 8).try_into().unwrap();

        for i in lo..hi {
            sum += self.base_counts[i];
        }
        return sum;
    }

    fn get_haplotype_counts(&self, h: u32) -> [u32; 4] {
        let mut a: [u32; 4] = [0; 4];
        for bi in 0u32..4u32 {
            a[bi as usize] = self.get(bi, h, 0) + self.get(bi, h, 1);
        }
        return a;
    }

    fn new() -> PileupStats {
        let ps = PileupStats { base_counts: [0; 16], mean_mapq: 0.0 };
        return ps;
    }

    fn clear(& mut self) -> () {
        self.base_counts = [0; 16];
        self.mean_mapq = 0.0;
    }

    fn fill_pileup(& mut self, cache: &mut ReadHaplotypeCache, alignments: rust_htslib::bam::pileup::Alignments<'_>) -> () {
        self.clear();
        let mut sum_mapq: f32 = 0.0;
        let mut n: f32 = 0.0;

        for a in alignments {
            if let Some(qpos) = a.qpos() {
                if let Some(mut hi) = cache.get(&a.record()) {
                    let read_base = a.record().seq()[qpos] as char;
                    hi -= 1; // phasing programs annotate with 1/2, we use 0/1
                    let bi = base2index(read_base) as u32;
                    let si = a.record().is_reverse() as u32;
                    self.increment(bi, hi as u32, si);
                    sum_mapq += a.record().mapq() as f32;
                    n += 1.0;
                    //println!("\t{read_base}\t{bi}\t{hi}\t{si}\t{}", *ps.get(bi, hi as u32, si));
                }
            }
        }
        self.mean_mapq = sum_mapq / n;
    }
}

fn base2index(base: char) -> i32 {
    return match base {
        'A' => 0,
        'C' => 1,
        'G' => 2,
        'T' => 3,
        _ => -1
    };
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

    let mut bam = bam::IndexedReader::from_path(input_bam).unwrap();
    println!("chromosome\tposition\treference_base\tvariant_base\tcanonical_type\tcanonical_context\taligned_depth\tmean_mapq\t\
              hmajor_variant_count\thminor_variant_count\thmajor_vaf\thminor_vaf\th1_a\th1_c\th1_g\th1_t\th2_a\th2_c\th2_g\th2_t");

    let chromosome_name = "chr20";
    let min_variant_observations = 5;
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

        if max_variant_count >= min_variant_observations && minor_variant_count <= max_variant_minor_observations && reference_base != 'N' {

            // grab reference context
            let mut reference_context = chromosome_bytes[ (reference_position - 1)..(reference_position + 2)].to_vec();

            let bases = "ACGT";
            let variant_base = bases.as_bytes()[max_variant_index as usize] as char;
            
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
                      {hmaj_vaf:.3}\t{hmin_vaf:.3}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                      h0[0], h0[1], h0[2], h0[3],
                      h1[0], h1[1], h1[2], h1[3]);
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

fn mutate_haplotypes(haplotypes: &mut [ Vec<u8> ], rate: f64) -> () {
    // Set up a generator to select random bases
    let bases = [ 'A', 'C', 'G', 'T' ];
    let mut rng = rand::thread_rng();

    // generate variants
    for hi in 0..2 {
        let mut _var_count = 0;
        for i in 0..haplotypes[hi].len() {
            if rng.gen::<f64>() < rate {
                let curr_base = haplotypes[hi][i];
                let mut var_base = curr_base;
                while var_base == curr_base {
                    var_base = bases[ rng.gen_range(0..4) ] as u8;
                }
                haplotypes[hi][i] = var_base;
                _var_count += 1;
            }
        }
        //eprintln!("Generate {_var_count} variants for haplotype");
    }
}

fn calculate_class_probabilities(alt_count: u64, ref_count: u64, error_rate: f64, purity: f64, ccf_dist: &Beta, mutation_rate: f64, heterozygosity: f64) -> [f64;3]
{

    let depth = ref_count + alt_count;
    // P(somatic | data) = P(data | somatic) P(somatic) / sum_class P(data | class )
    // P(data | ref) = Binom(alt_count, ref_count + alt_count, error_rate)
    let p_data_ref = Binomial::new(error_rate, depth).unwrap().pmf(alt_count);

    // P(data | het) = Binom(alt_count, ref_count + alt_count, 1 - error_rate)
    let p_data_het = Binomial::new(1.0 - error_rate, depth).unwrap().pmf(alt_count);
    
    // P(data | somatic) = sum_c Binom(alt_count, ref_count + alt_count, purity * c * (1 - error_rate) + (1 - purity*c) * error_rate ) P(c)
    let bins = 10;
    let step = 1.0 / 10 as f64;
    let mut p_data_somatic = 0.0;
    for i in 0..bins {
        let start = f64::from(i) * step;
        let end = f64::from(i + 1) * step;
        let c = (end + start) / 2.0;
        let p_c = ccf_dist.cdf(end) - ccf_dist.cdf(start);

        let p_read_from_mutated_haplotype = purity * c;
        let t1 = p_read_from_mutated_haplotype * (1.0 - error_rate);
        let t2 = (1.0 - p_read_from_mutated_haplotype) * error_rate;

        let p_data_somatic_at_c = Binomial::new(t1 + t2, depth).unwrap().pmf(alt_count);
        p_data_somatic += p_data_somatic_at_c * p_c;
    }

    // priors
    let p_het = heterozygosity;
    let p_somatic = mutation_rate;
    let p_ref = 1.0 - p_het - p_somatic;

    let sum = p_data_ref * p_ref + p_data_het * p_het + p_data_somatic * p_somatic;
    
    let p_ref_data = p_data_ref * p_ref / sum;
    let p_het_data = p_data_het * p_het / sum;
    let p_somatic_data = p_data_somatic * p_somatic / sum;

    return [ p_ref_data, p_het_data, p_somatic_data ];
}

fn sim_pileup(reference_genome: &str)
{
    let summary = false;
    let mutation_rate = 5.0 / 1000000.0; // per haplotype
    let heterozygosity = 1.0 / 2000.0; // per haplotype
    let ccf_dist = Beta::new(9.0, 1.0).unwrap();
    let depth_dist = Poisson::new(50.0).unwrap();
    let purity = 0.75;
    let error_rate = 0.02;

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
    mutate_haplotypes(& mut germline_haplotypes, heterozygosity);
    
    let mut somatic_haplotypes = [ germline_haplotypes[0].clone(), germline_haplotypes[1].clone() ];
    mutate_haplotypes(& mut somatic_haplotypes, mutation_rate);

    // Finally, simulate some pileup data at each position on each haplotype
    if !summary {
        println!("position\thaplotype\tclass\tis_het\tis_somatic\tccs\tref_count\talt_count\thvaf\tp_ref\tp_het\tp_somatic");
    }

    let mut sum_ref = 0.0;
    let mut sum_het = 0.0;
    let mut sum_somatic = 0.0;
    let mut num_somatic_true = 0;
    let mut num_somatic_calls = 0;
    let mut num_somatic_correct = 0;
    let mut num_het_true = 0;

    for j in 0..chromosome_length {

        // depth at this position
        let depth = depth_dist.sample(&mut rng) as u64;
        let haplotype_dist = Binomial::new(0.5, depth).unwrap();
        let h0_depth = haplotype_dist.sample(&mut rng) as u64;
        let haplotype_depths = [ h0_depth, depth - h0_depth ];

        for i in 0..2 {
            let is_het = germline_haplotypes[i][j] != chromosome_bytes[j];
            let is_somatic = somatic_haplotypes[i][j] != germline_haplotypes[i][j];
            let ccf = match is_somatic {
                true => ccf_dist.sample(&mut rng),
                false => 0.0
            };

            let mut ref_count = 0u64;
            let mut alt_count = 0u64;

            for _r in 0..haplotype_depths[i] {
                // determine base for this read
                let prob_sampled_from_somatic = match is_somatic {
                    true => purity * ccf,
                    false => 0.0
                };

                let is_tumor_read = rng.gen::<f64>() < prob_sampled_from_somatic;
                let read_base = match is_tumor_read {
                    true => somatic_haplotypes[i][j],
                    false => germline_haplotypes[i][j]
                };

                let mut ref_match = read_base == chromosome_bytes[j];
                if rng.gen::<f64>() < error_rate {
                    ref_match = ! ref_match;
                }

                // update counts
                ref_count += ref_match as u64;
                alt_count += !ref_match as u64;
            }

            let probs = calculate_class_probabilities(alt_count, ref_count, error_rate, purity, &ccf_dist, mutation_rate, heterozygosity);

            let mut class = "REF";
            if is_somatic { 
                class = "SOMATIC"
            } else if is_het {
                class = "HET";   
            }

            let mut p_max = 0.0;
            let mut idx_max = 0;
            for i in 0..probs.len() {
                if probs[i] > p_max {
                    p_max = probs[i];
                    idx_max = i;
                }
            }

            // update summary stats
            sum_ref += probs[0];
            sum_het += probs[1];
            sum_somatic += probs[2];

            if idx_max == 2 && is_somatic {
                num_somatic_correct += 1;   
            }
            
            if idx_max == 2 {
                num_somatic_calls += 1;
            }

            if is_somatic { 
                num_somatic_true += 1;
            }

            if is_het { 
                num_het_true += 1;
            }

            // output per-record calls, if wanted
            let hvaf = alt_count as f64 / haplotype_depths[i] as f64; 
            if !summary {
                println!("{j}\t{i}\t{class}\t{}\t{}\t{ccf:.3}\t\
                         {ref_count}\t{alt_count}\t{hvaf:.3}\t{:.3}\t{:.3}\t{:.3}",
                         is_het as u8, is_somatic as u8, probs[0], probs[1], probs[2]);
            }
        }
    }

    if summary {
        let sens = num_somatic_correct as f64 / num_somatic_true as f64;
        let prec = num_somatic_correct as f64 / num_somatic_calls as f64;
        println!("true_somatic_positions\ttrue_het_position\testimated_ref_bases\testimated_het_bases\testimated_mutated_bases\tsomatic_call_sensitivity\tsomatic_call_precision");
        println!("{num_somatic_true}\t{num_het_true}\t{sum_ref:.1}\t{sum_het:.1}\t{sum_somatic:.1}\t{sens:.3}\t{prec:.3}");
    }
}
