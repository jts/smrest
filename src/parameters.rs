//---------------------------------------------------------
// Copyright 2022 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
use statrs::distribution::{Poisson, Beta, ContinuousCDF};

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
