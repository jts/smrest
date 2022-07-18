//---------------------------------------------------------
// Copyright 2022 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
use crate::utility::ReadHaplotypeCache;

pub struct PileupStats
{
    // this is indexed by base (dim four), haplotype (two), strand (two)
    base_counts: [u32; 16],
    pub mean_mapq: f32
}

impl PileupStats {
    #[inline(always)]
    pub fn get(&self, b: u32, h: u32, s: u32) -> u32 {
        let i = (h * 8) + (s * 4) + b;
        return self.base_counts[ i as usize ];
    }

    pub fn increment(&mut self, b: u32, h: u32, s: u32) -> () {
        let i = (h * 8) + (s * 4) + b;
        self.base_counts[i as usize] += 1;
    }

    pub fn get_count_on_haplotype(&self, b: u32, h: u32) -> u32 {
        return self.get(b, h, 0) + self.get(b, h, 1);
    }
    
    pub fn get_count_unphased(&self, b: u32) -> u32 {
        return self.get_count_on_haplotype(b, 0) + self.get_count_on_haplotype(b, 1);
    }

    pub fn get_haplotype_depth(&self, h: u32) -> u32 {
        let mut sum = 0;
        let lo:usize = (h * 8).try_into().unwrap();
        let hi:usize = ((h+1) * 8).try_into().unwrap();

        for i in lo..hi {
            sum += self.base_counts[i];
        }
        return sum;
    }

    pub fn get_depth(&self) -> u32 {
        return self.get_haplotype_depth(0) + self.get_haplotype_depth(1);
    }


    pub fn get_haplotype_counts(&self, h: u32) -> [u32; 4] {
        let mut a: [u32; 4] = [0; 4];
        for bi in 0u32..4u32 {
            a[bi as usize] = self.get(bi, h, 0) + self.get(bi, h, 1);
        }
        return a;
    }

    pub fn get_max_nonref_base_on_haplotype(&self, h: u32, ref_base_index: u32) -> u32 {
        let mut max_count = 0;
        let mut max_index = 0;
        for bi in 0u32..4u32 {
            let c = self.get_count_on_haplotype(bi, h);
            if bi != ref_base_index && c >= max_count {
                max_count = c;
                max_index = bi;
            }
        }
        return max_index;
    }
    
    pub fn get_max_nonref_base_unphased(&self, ref_base_index: u32) -> u32 {
        let mut max_count = 0;
        let mut max_index = 0;
        for bi in 0u32..4u32 {
            let c = self.get_count_on_haplotype(bi, 0) + self.get_count_on_haplotype(bi, 1);
            if bi != ref_base_index && c >= max_count {
                max_count = c;
                max_index = bi;
            }
        }
        return max_index;
    }

    pub fn new() -> PileupStats {
        let ps = PileupStats { base_counts: [0; 16], mean_mapq: 0.0 };
        return ps;
    }

    pub fn clear(& mut self) -> () {
        self.base_counts = [0; 16];
        self.mean_mapq = 0.0;
    }

    pub fn fill_pileup(& mut self, cache: &mut ReadHaplotypeCache, alignments: rust_htslib::bam::pileup::Alignments<'_>) -> () {
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

pub fn base2index(base: char) -> i32 {
    return match base {
        'A' => 0,
        'C' => 1,
        'G' => 2,
        'T' => 3,
        _ => -1
    };
}

