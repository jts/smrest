//---------------------------------------------------------
// Copyright 2022 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------

use statrs::distribution::{Binomial, Discrete, Poisson, Beta, ContinuousCDF};
use bio::stats::{Prob, LogProb};
use cached::proc_macro::cached;
use crate::longshot_realign::ReadHaplotypeLikelihood;

// traits
use statrs::statistics::{Min, Max};
use rand::Rng;

static VERBOSE: bool = false;

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

pub struct ModelParameters {
    pub mutation_rate: f64,
    pub heterozygosity: f64,
    pub ccf_dist: CancerCellFraction,
    pub depth_dist: Option<Poisson>, // only use for simulation
    pub purity: f64,
    pub error_rate: f64
}

pub fn calculate_class_probabilities_phased(alt_count: u64, ref_count: u64, params: &ModelParameters) -> [f64;3]
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

    /* flat prior
    let p_het = 1.0 / 3.0;
    let p_somatic = p_het;
    let p_ref = p_het;
    */

    let sum = p_data_ref * p_ref + p_data_het * p_het + p_data_somatic * p_somatic;
    
    let p_ref_data = p_data_ref * p_ref / sum;
    let p_het_data = p_data_het * p_het / sum;
    let p_somatic_data = p_data_somatic * p_somatic / sum;

    return [ p_ref_data, p_het_data, p_somatic_data ];
}

pub fn calculate_class_probabilities_unphased(alt_count: u64, ref_count: u64, params: &ModelParameters) -> [f64;3]
{
    let depth = ref_count + alt_count;
    //
    // P(somatic | data) = P(data | somatic) P(somatic) / sum_class P(data | class )
    //

    // P(data | ref) = Binom(alt_count, ref_count + alt_count, error_rate)
    let p_data_ref = Binomial::new(params.error_rate, depth).unwrap().pmf(alt_count);

    // P(data | het) = Binom(alt_count, ref_count + alt_count, 1 - error_rate)
    let p_data_het = Binomial::new(0.5 * (1.0 - params.error_rate) + 0.5 * params.error_rate, depth).unwrap().pmf(alt_count);
    
    // P(data | somatic) = sum_c Binom(alt_count, ref_count + alt_count, 0.5 * purity * c * (1 - error_rate) + (1 - purity * c * 0.5) * error_rate ) P(c)
    let bins = 10;
    let step = 1.0 / 10 as f64;
    let mut p_data_somatic = 0.0;
    for i in 0..bins {
        let start = f64::from(i) * step;
        let end = f64::from(i + 1) * step;
        let c = (end + start) / 2.0;
        let p_c = params.ccf_dist.cdf(end) - params.ccf_dist.cdf(start);

        let p_read_from_mutated_haplotype = params.purity * c * 0.5;
        let t1 = p_read_from_mutated_haplotype * (1.0 - params.error_rate);
        let t2 = (1.0 - p_read_from_mutated_haplotype) * params.error_rate;

        let p_data_somatic_at_c = Binomial::new(t1 + t2, depth).unwrap().pmf(alt_count);
        p_data_somatic += p_data_somatic_at_c * p_c;
    }

    // priors
    let p_het = params.heterozygosity;
    let p_somatic = params.mutation_rate;
    let p_ref = 1.0 - p_het - p_somatic;
    
    /*
    let p_het = 1.0 / 3.0;
    let p_somatic = p_het;
    let p_ref = p_het;
    */
    let sum = p_data_ref * p_ref + p_data_het * p_het + p_data_somatic * p_somatic;
    
    let p_ref_data = p_data_ref * p_ref / sum;
    let p_het_data = p_data_het * p_het / sum;
    let p_somatic_data = p_data_somatic * p_somatic / sum;

    return [ p_ref_data, p_het_data, p_somatic_data ];
}

// https://en.wikipedia.org/wiki/Binomial_test
// this function is slow so we approximate it by converting
// the p parameter to an integer so it can be memoized 
pub fn binomial_test_twosided(x: u64, n: u64, p: f64) -> f64 {
    let pi = (p * 100.0) as u64;
    return binomial_test_twosided_memoized(x, n, pi);
}

#[cached]
fn binomial_test_twosided_memoized(x: u64, n: u64, pi: u64) -> f64 {
    let p = pi as f64 / 100.0;
    let bn = Binomial::new(p, n).unwrap();
    let d = bn.pmf(x);

    let sum = (0..=n).map(|i| bn.pmf(i)).filter(|v| v <= &d).sum();
    return sum;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_binomial_test() -> Result<(), String> {
        let e = 0.00001;
        assert_abs_diff_eq!(binomial_test_twosided(2, 15, 0.6), 0.0002789, epsilon = e);
        assert_abs_diff_eq!(binomial_test_twosided(3, 15, 0.6), 0.002398, epsilon = e);
        assert_abs_diff_eq!(binomial_test_twosided(30, 100, 0.4), 0.04154, epsilon = e);
        assert_abs_diff_eq!(binomial_test_twosided(915, 1000, 0.85), 8.244e-10, epsilon = e);
        assert_abs_diff_eq!(binomial_test_twosided(5, 1000, 0.85), 0.0, epsilon = e);
        Ok(())
    }
}

// model from: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5832436/pdf/pcbi.1005965.pdf
pub fn calculate_class_probabilities_sgz(alt_count: u64, ref_count: u64, params: &ModelParameters) -> [f64;3]
{
    let depth = ref_count + alt_count;
    let alpha = 0.01; // from paper

    let af_germline = 0.5;
    let af_somatic = params.purity / 2.0;

    let p_g = binomial_test_twosided(alt_count, depth, af_germline);
    let p_s = binomial_test_twosided(alt_count, depth, af_somatic);
    let f = alt_count as f64 / depth as f64;

    let call_subclonal = false;

    if p_s > alpha && p_g <= alpha {
        // predicted somatic
        return [ 0.0, 0.0, 1.0 ];
    } else if call_subclonal && p_s <= alpha && p_g <= alpha && f < af_somatic / 1.5 {
        // predicted subclonal somatic
        return [ 0.0, 0.0, 1.0 ];
    } else if p_s <= alpha && p_g > alpha {
        // predicted germline
        return [ 0.0, 1.0, 0.0 ];
    } else  {
        // predicted reference
        // this is not strictly what they do in the paper but works in our framework
        return [ 1.0, 0.0, 0.0 ];
    }
}

pub fn calculate_class_probabilities_likelihood(rhls: &Vec<ReadHaplotypeLikelihood>, params: &ModelParameters) -> [f64;3] {

    // calculate likelihood using the base
    
    //
    // P(somatic | data) = P(data | somatic) P(somatic) / sum_class P(data | class )
    //

    // P(data | ref) = prod_reads P(read | ref haplotype)
    let lp_data_ref:f64 = rhls.iter().map(|x| *x.base_allele_likelihood).sum();

    // P(data | het) = prod_reads P(read | alt haplotype)
    let lp_data_het:f64 = rhls.iter().map(|x| *x.mutant_allele_likelihood).sum();

    // Integrate over purity/ccf as in phased model above
    // P(data | somatic) = prod_reads P(read | alt haplotype) P(alt haplotype) + P(read | ref haplotype) P(ref haplotype)
    // P(data | somatic) = sum_c Binom(alt_count, ref_count + alt_count, purity * c * (1 - error_rate) + (1 - purity*c) * error_rate ) P(c)
    let bins = 10;
    let step = 1.0 / 10 as f64;
    let mut lp_data_somatic = 0.0;
    for rhl in rhls {
        let mut lp_read_somatic = LogProb::ln_zero();

        for i in 0..bins {
            let start = f64::from(i) * step;
            let end = f64::from(i + 1) * step;
            let c = (end + start) / 2.0;
            let lp_c = LogProb::from(Prob(params.ccf_dist.cdf(end) - params.ccf_dist.cdf(start)));

            let p_read_from_mutated_haplotype = params.purity * c;

            let t1 = LogProb::from(Prob(p_read_from_mutated_haplotype)) + rhl.mutant_allele_likelihood;
            let t2 = LogProb::from(Prob(1.0 - p_read_from_mutated_haplotype)) + rhl.base_allele_likelihood;
            let s = LogProb::ln_add_exp(t1, t2) + lp_c;
            lp_read_somatic = LogProb::ln_add_exp(lp_read_somatic, s);
        }
        lp_data_somatic += *lp_read_somatic;

        if VERBOSE {
            println!("{} {} call: {} qual_log: {:.3} qual_p: {} allele scores: [ {:.3} {:.3} ] somatic: {:.3} {:.3}", 
                rhl.read_name.as_ref().unwrap(), rhl.haplotype_index.unwrap_or(-1), rhl.allele_call, *rhl.allele_call_qual, 
                *Prob::from(rhl.allele_call_qual), *rhl.base_allele_likelihood, *rhl.mutant_allele_likelihood, 
                *lp_read_somatic, lp_data_somatic);
        }
    }


    // priors
    let p_het = params.heterozygosity;
    let p_somatic = params.mutation_rate;
    let p_ref = 1.0 - p_het - p_somatic;
    
    let lp_t_ref = LogProb(lp_data_ref) + LogProb::from(Prob(p_ref));
    let lp_t_het = LogProb(lp_data_het) + LogProb::from(Prob(p_het));
    let lp_t_somatic = LogProb(lp_data_somatic) + LogProb::from(Prob(p_somatic));

    let lp_sum = LogProb::ln_sum_exp( &[ lp_t_ref, lp_t_het, lp_t_somatic] );
    let p_ref_data = Prob::from(lp_t_ref - lp_sum);
    let p_het_data = Prob::from(lp_t_het - lp_sum);
    let p_somatic_data = Prob::from(lp_t_somatic - lp_sum);
    
    if VERBOSE {
        println!("data likelihoods: [ {} {} {} ]", lp_data_ref, lp_data_het, lp_data_somatic);
        println!("sum: {}", *lp_sum);
        println!("probabilities: [ {} {} {} ]", *p_ref_data, *p_het_data, *p_somatic_data);
    }
    return [ *p_ref_data, *p_het_data, *p_somatic_data ];
}

