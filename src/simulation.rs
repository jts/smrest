//---------------------------------------------------------
// Copyright 2022 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
use rust_htslib::{faidx};
use statrs::distribution::{Binomial};
use std::fs::File;
use std::io::Write;
use rand::prelude::*;
use crate::ModelParameters;
use crate::PileupStats;
use crate::base2index;
use crate::classifier::calculate_class_probabilities_phased;
use crate::classifier::calculate_class_probabilities_unphased;
use crate::classifier::calculate_class_probabilities_sgz;

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

pub struct MutationOutput
{
    out: File
}

impl MutationOutput {
    pub fn new(output_filename: &str) -> MutationOutput {
        let mut r = MutationOutput {
            out: File::create(output_filename).unwrap()
        };
        r.out.write(b"chromosome\tposition\thaplotype\tref_base\talt_base\thvaf\n").expect("Unable to write");
        return r;
    }

    pub fn write(& mut self, chr: &String, position: usize, haplotype: usize, ref_base: char, alt_base: char, hvaf: f64) -> () {
        let s = format!("{}\t{}\t{}\t{}\t{}\t{:.4}\n", chr, position, haplotype, ref_base, alt_base, hvaf);
        self.out.write(s.as_bytes()).expect("Unable to write");
    }
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
        let subclonal_ccf_mean = params.ccf_dist.subclonal_mean_ccf();

        println!("{}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{}\t{}\t{:.1}\t{:.1}\t{:.1}\t{}\t{:.3}\t{:.3}\t{:.3}",
            self.model_name, params.purity, params.depth_dist.unwrap().lambda(), params.ccf_dist.p_clonal, subclonal_ccf_mean,
            self.num_somatic_true, self.num_het_true, 
            self.class_posterior_sums[0], self.class_posterior_sums[1], self.class_posterior_sums[2],
            self.num_somatic_calls, sens, prec, f1);
    }
}

pub fn sim_pileup(params: &ModelParameters, reference_genome: &str, per_site_output: bool, mutation_output_filename: &str)
{
    let mut rng = rand::thread_rng();

    // get genome sequence
    let chromosome_name = String::from("chr20");
    //let chromosome_length = 60000000;
    let chromosome_length = 10000000;

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
    let mut sgz_model_stats = SimulationStats::new("sgz".to_string());

    let mut mutation_out = MutationOutput::new(mutation_output_filename);

    // Finally, simulate some pileup data at each position on each haplotype
    if per_site_output {
        println!("model\tposition\thaplotype\tclass\tis_het\tis_somatic\tref_count\talt_count\thvaf\tp_ref\tp_het\tp_somatic");
    }

    for j in 0..chromosome_length {

        if chromosome_bytes[j] as char == 'N' {
            continue;
        }

        // depth at this position
        let depth = params.depth_dist.unwrap().sample(&mut rng) as u64;
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

            if is_somatic {
                mutation_out.write(&chromosome_name, j + 1, i+1, chromosome_bytes[j] as char , somatic_haplotypes[i][j] as char, ccf * params.purity);
            }

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
            if per_site_output {
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

        // run baseline non-phased models
        {
            let alt_base_index = ps.get_max_nonref_base_unphased(ref_base_index as u32);
            let alt_count = ps.get_count_unphased(alt_base_index) as u64;
            let ref_count = ps.get_count_unphased(ref_base_index) as u64;
            let un_probs = calculate_class_probabilities_unphased(alt_count, ref_count, &params);
            let sgz_probs = calculate_class_probabilities_sgz(alt_count, ref_count, &params);

            // update summary stats
            let is_het = germline_haplotypes[0][j] != chromosome_bytes[j] || germline_haplotypes[1][j] != chromosome_bytes[j];
            let is_somatic = somatic_haplotypes[0][j] != germline_haplotypes[0][j] || somatic_haplotypes[1][j] != germline_haplotypes[1][j];
            
            unphased_model_stats.update(is_het, is_somatic, &un_probs);
            sgz_model_stats.update(is_het, is_somatic, &sgz_probs);
            
            // output per-record calls, if wanted
            if per_site_output {
                let hvaf = alt_count as f64 / ps.get_depth() as f64; 
                let mut class = "REF";
                if is_somatic { 
                    class = "SOMATIC"
                } else if is_het {
                    class = "HET";   
                }

                println!("unphased\t{j}\tboth\t{class}\t{}\t{}\t\
                         {ref_count}\t{alt_count}\t{hvaf:.3}\t{:.3}\t{:.3}\t{:.3}",
                         is_het as u8, is_somatic as u8, un_probs[0], un_probs[1], un_probs[2]);
                
                println!("sgz\t{j}\tboth\t{class}\t{}\t{}\t\
                         {ref_count}\t{alt_count}\t{hvaf:.3}\t{:.3}\t{:.3}\t{:.3}",
                         is_het as u8, is_somatic as u8, sgz_probs[0], sgz_probs[1], sgz_probs[2]);
            }
        }
    }

    if ! per_site_output {
        println!("model\tpurity\tmean_coverage\tproportion_clonal\tsubclonal_mean_ccf\ttrue_somatic_positions\ttrue_het_position\testimated_ref_bases\t\
                  estimated_het_bases\testimated_mutated_bases\tnum_somatic_calls\tsomatic_call_sensitivity\tsomatic_call_precision\tf1");
        phased_model_stats.print(& params);
        unphased_model_stats.print(& params);
        sgz_model_stats.print(& params);
    }
}

