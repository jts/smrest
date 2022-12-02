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

#[derive(Clone)]
pub struct ReadMetadata
{
    pub haplotype_index: Option<i32>,
    pub phase_set: Option<i32>,
    pub strand_index: i32,
    pub leading_softclips: i64,
    pub trailing_softclips: i64
}
