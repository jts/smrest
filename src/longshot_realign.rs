use bio::stats::probs::{LogProb};

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
