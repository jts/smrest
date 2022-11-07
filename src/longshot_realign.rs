use crate::pileup_stats::PileupStats;
use crate::pileup_stats::base2index;
use crate::utility::ReadHaplotypeCache;
use rust_htslib::{bam::record::CigarStringView};
use bio::stats::probs::{Prob, LogProb};

use longshot::genotype_probs::{Genotype, GenotypeProbs};
use longshot::extract_fragments::{extract_fragment, create_augmented_cigarlist, CigarPos, ExtractFragmentParameters};
use longshot::realignment::{AlignmentParameters};
use longshot::variants_and_fragments::{Var, VarFilter};

//
fn create_var(tid: u32, pos0: usize, ref_allele: char, alt_allele: char) -> Var {
    let allele_vec = vec![ref_allele.to_string(), alt_allele.to_string()];
    let allele_num = allele_vec.len();
    let new_var = Var {
        ix: 0,
        tid: tid,
        pos0: pos0,
        alleles: allele_vec,
        dp: 0,
        allele_counts: vec![0; allele_num],
        allele_counts_forward: vec![0; allele_num],
        allele_counts_reverse: vec![0; allele_num],
        ambiguous_count: 0,
        qual: 0.0,
        filter: VarFilter::Pass,
        genotype: Genotype(0, 0),
        //unphased: false,
        gq: 0.0,
        unphased_genotype: Genotype(0, 0),
        unphased_gq: 0.0,
        genotype_post: GenotypeProbs::uniform(allele_num),
        phase_set: None,
        strand_bias_pvalue: 0.0,
        mec: 0,
        mec_frac_variant: 0.0, // mec fraction for this variant
        mec_frac_block: 0.0,   // mec fraction for this haplotype block
        mean_allele_qual: 0.0,
        dp_any_mq: 0,
        mq10_frac: 0.0,
        mq20_frac: 0.0,
        mq30_frac: 0.0,
        mq40_frac: 0.0,
        mq50_frac: 0.0,
    };
    return new_var;
}

pub fn recalculate_pileup(ps: &PileupStats,
                          tid: u32,
                          reference_position: usize,
                          chromosome_char_vec: &Vec<char>,
                          t_names: &Vec<String>,
                          ref_allele: char,
                          alt_allele: char,
                          cache: &mut ReadHaplotypeCache, 
                          alignments: rust_htslib::bam::pileup::Alignments<'_>,
                          extract_fragment_parameters: ExtractFragmentParameters,
                          alignment_parameters: AlignmentParameters) -> PileupStats
{
    let mut rcps = PileupStats::new();
    rcps.mean_mapq = ps.mean_mapq;
    rcps.proportion_phased = ps.proportion_phased;

    // 
    let var = create_var(tid, reference_position, ref_allele, alt_allele);
    
    let ref_index = base2index(ref_allele) as u32;
    let alt_index = base2index(alt_allele) as u32;

    for a in alignments {
        if a.record().seq().len() == 0 {
            continue;
        }

        if let Some(_qpos) = a.qpos() {
            if let Some(mut hi) = cache.get(&a.record()) {
                // perform realignment
                let varlist = vec![ var.clone() ];
                let record = a.record();
                let start_pos = record.pos();
                let bam_cig: CigarStringView = record.cigar();
                let cigarpos_list: Vec<CigarPos> =
                    create_augmented_cigarlist(start_pos as u32, &bam_cig).unwrap();

                let f = extract_fragment(&a.record(), 
                    &cigarpos_list, 
                    varlist, 
                    &chromosome_char_vec, 
                    &t_names, 
                    extract_fragment_parameters, 
                    alignment_parameters);

                // extract fragment returns an Option wrapped in a result
                // we handle both in the same way here, by not processing the
                // read in case of any type of failure during realignment
                if let Ok(r) = f {
                    if let Some(fragment) = r {
                        let calls = fragment.calls;
                        
                        if calls.len() > 0 {
                            let called = var.alleles[calls[0].allele as usize].clone();
                            let p_wrong = *Prob::from(calls[0].qual);
                            if p_wrong < 0.1 {
                                let ui = if called == var.alleles[0] { ref_index } else { alt_index };
                                let si = a.record().is_reverse() as u32;
                                hi -= 1; // phasing programs annotate with 1/2, we use 0/1
                                rcps.increment(ui, hi as u32, si)
                            }
                        }
                    }
                }
                //println!("{}\t{hi}\t{}\t{:.3}", fragment.id.unwrap(), called, p_wrong);
            }
        }
    }

    return rcps;
}

#[derive(Clone)]
pub struct ReadHaplotypeLikelihood
{
    pub read_name: Option<String>,
    pub mutant_allele_likelihood: LogProb,
    pub base_allele_likelihood: LogProb,
    pub allele_call: String,
    pub allele_call_qual: LogProb,
    pub haplotype_index: Option<i32>,
    pub strand_index: i32
}

pub fn calculate_likelihoods(tid: u32,
                             reference_position: usize,
                             candidate_haplotype_index: i32,
                             chromosome_char_vec: &Vec<char>,
                             t_names: &Vec<String>,
                             ref_allele: char,
                             alt_allele: char,
                             cache: &mut ReadHaplotypeCache, 
                             alignments: rust_htslib::bam::pileup::Alignments<'_>,
                             extract_fragment_parameters: ExtractFragmentParameters,
                             alignment_parameters: AlignmentParameters) -> Vec<ReadHaplotypeLikelihood>
{
    // 
    let mut out = vec![];
    let var = create_var(tid, reference_position, ref_allele, alt_allele);

    let before_ref = chromosome_char_vec[reference_position - 1];
    let before_alt = if before_ref != alt_allele { alt_allele } else { ref_allele };
    let before_var = create_var(tid, reference_position - 1, before_ref, before_alt);

    let after_ref = chromosome_char_vec[reference_position + 1];
    let after_alt = if after_ref != alt_allele { alt_allele } else { ref_allele };
    let after_var = create_var(tid, reference_position + 1, after_ref, after_alt);
    
    for a in alignments {
        if a.record().seq().len() == 0 {
            continue;
        }

        if let Some(hi) = cache.get(&a.record()) {
            if hi == candidate_haplotype_index {
                // perform realignment
                let varlist = vec![ before_var.clone(), var.clone(), after_var.clone() ];
                //let varlist = vec![ var.clone() ];
                let record = a.record();
                let start_pos = record.pos();
                let bam_cig: CigarStringView = record.cigar();
                let cigarpos_list: Vec<CigarPos> =
                    create_augmented_cigarlist(start_pos as u32, &bam_cig).unwrap();

                let f = extract_fragment(&a.record(), 
                    &cigarpos_list, 
                    varlist, 
                    &chromosome_char_vec, 
                    &t_names, 
                    extract_fragment_parameters, 
                    alignment_parameters);

                // extract fragment returns an Option wrapped in a result
                // we handle both in the same way here, by not processing the
                // read in case of any type of failure during realignment
                if let Ok(r) = f {
                    if let Some(fragment) = r {
                        let calls = fragment.calls;
                        let target_var_idx = 1;
                        if calls.len() >= 2 {
                            let rhl = ReadHaplotypeLikelihood { 
                                read_name: fragment.id,
                                mutant_allele_likelihood: calls[target_var_idx].allele_scores[1],
                                base_allele_likelihood: calls[target_var_idx].allele_scores[0],
                                allele_call: var.alleles[calls[target_var_idx].allele as usize].clone(),
                                allele_call_qual: calls[target_var_idx].qual,
                                haplotype_index: Some(hi),
                                strand_index: record.is_reverse() as i32
                            };
                            out.push(rhl);
                        }
                    }
                }
                //println!("{}\t{hi}\t{}\t{:.3}", fragment.id.unwrap(), called, p_wrong);
            }
        }
    }

    return out;
}

