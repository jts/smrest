//---------------------------------------------------------
// Copyright 2022 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
use statrs::distribution::{Poisson, Beta, ContinuousCDF};
use longshot::extract_fragments::{ExtractFragmentParameters};
use longshot::context_model::ContextModel;
use longshot::genotype_probs::GenotypePriors;
use longshot::realignment::{AlignmentType, AlignmentParameters};
use longshot::estimate_alignment_parameters::estimate_alignment_parameters;
use longshot::util::GenomicInterval;
use bio::stats::{Prob, LogProb};

// traits
use statrs::statistics::{Min, Max};
use rand::Rng;

//
// The distribution of mutation frequencies across the cancer
// cells within a population
//
pub struct CancerCellFraction
{
    pub p_clonal: f64,
    pub subclonal_ccf: Beta
}

impl CancerCellFraction {
    pub fn subclonal_mean_ccf(&self) -> f64 {
        return self.subclonal_ccf.shape_a() / (self.subclonal_ccf.shape_a() + self.subclonal_ccf.shape_b());
    }
}

impl ::rand::distributions::Distribution<f64> for CancerCellFraction {
    fn sample<R: Rng + ?Sized>(&self, r: &mut R) -> f64 {
        let is_clonal = r.gen::<f64>() < self.p_clonal;
        if is_clonal {
            return 1.0;
        } else {
            return self.subclonal_ccf.sample(r);
        }
    }
}

impl Min<f64> for CancerCellFraction {
    fn min(&self) -> f64 { return 0.0; }
}

impl Max<f64> for CancerCellFraction {
    fn max(&self) -> f64 { return 1.0; }
}

impl ContinuousCDF<f64, f64> for CancerCellFraction {
    fn cdf(&self, x: f64) -> f64 {
        let mut p = 0.0;
        
        // TODO: better numerical tolerance
        if x > 0.999 {
            p = self.p_clonal;
        }

        p += (1.0 - self.p_clonal) * self.subclonal_ccf.cdf(x);
        return p;
    }
}

// Parameters for the mutation calling model
pub struct ModelParameters {
    pub mutation_rate: f64,
    pub heterozygosity: f64,
    pub ccf_dist: CancerCellFraction,
    pub depth_dist: Option<Poisson>, // only use for simulation
    pub purity: f64,
    pub error_rate: f64
}

impl ModelParameters {
    pub fn defaults() -> ModelParameters {
        let ccf = CancerCellFraction { p_clonal: 0.75, subclonal_ccf: Beta::new(2.0, 2.0).unwrap() };
        ModelParameters { 
            mutation_rate: 5.0 / 1000000.0, // per haplotype
            heterozygosity: 1.0 / 2000.0, // per haplotype
            ccf_dist: ccf,
            depth_dist: None,
            purity: 0.75,
            error_rate: 0.05
        }
    }
}

// Parameters for longshot modules
pub struct LongshotParameters
{
    pub genotype_priors: GenotypePriors,
    pub extract_fragment_parameters: ExtractFragmentParameters,
    pub alignment_parameters: Option<AlignmentParameters>
}

//
const MIN_MAPQ: u8 = 50;
const CONTEXT_MODEL_K: usize = 5;

impl LongshotParameters {
    pub fn defaults() -> LongshotParameters {
        let hom_snv_rate = LogProb::from(Prob(0.0005));
        let het_snv_rate = LogProb::from(Prob(0.001));
        let hom_indel_rate = LogProb::from(Prob(0.00005));
        let het_indel_rate = LogProb::from(Prob(0.00001));
        let ts_tv_ratio = 2.0 * 0.5; // from longshot defaults...

        let gt_priors = GenotypePriors::new(
            hom_snv_rate,
            het_snv_rate,
            hom_indel_rate,
            het_indel_rate,
            ts_tv_ratio,
        ).unwrap();
        
        let context_model = ContextModel::init(CONTEXT_MODEL_K);
    
        let efp = ExtractFragmentParameters {
            min_mapq: MIN_MAPQ,
            alignment_type: AlignmentType::ForwardAlgorithmNumericallyStableWithContext,
            //alignment_type: AlignmentType::ForwardAlgorithmNumericallyStable,
            context_model: context_model,
            band_width: 20,
            anchor_length: 6,
            variant_cluster_max_size: 3,
            max_window_padding: 50,
            max_cigar_indel: 20,
            store_read_id: true
        };

        LongshotParameters {
            genotype_priors: gt_priors,
            extract_fragment_parameters: efp,
            alignment_parameters: None
        }
    }

    pub fn estimate_alignment_parameters(&mut self, 
                                         bam_file: &String, 
                                         fasta_file: &String, 
                                         interval: &Option<GenomicInterval>) -> () {

        self.alignment_parameters = 
            Some(estimate_alignment_parameters(bam_file, 
                                               fasta_file,
                                               interval, 
                                               60, // min mapq, more conservative here than extract parameters
                                               self.extract_fragment_parameters.max_cigar_indel as u32, 
                                               100000, // num reads
                                               &self.extract_fragment_parameters.context_model).unwrap());
    }
}
