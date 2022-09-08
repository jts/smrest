//---------------------------------------------------------
// Copyright 2022 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------

use statrs::distribution::{Binomial, Discrete, Poisson, Beta, ContinuousCDF};

pub struct ModelParameters {
    pub mutation_rate: f64,
    pub heterozygosity: f64,
    pub ccf_dist: Beta,
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

// https://en.wikipedia.org/wiki/Binomial_test
fn binomial_test_twosided(x: u64, n: u64, p: f64) -> f64 {
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
        Ok(())
    }
}
