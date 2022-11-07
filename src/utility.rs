//---------------------------------------------------------
// Copyright 2022 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
use rust_htslib::{bam, bam::record::Aux};
use std::collections::{HashMap};
use intervaltree::IntervalTree;
use core::ops::Range;

// A cache storing the haplotype tag for every read processed
// We need this because bam_get_aux is O(N)
pub struct ReadHaplotypeCache
{
    pub cache: HashMap<String, i32>
}

impl ReadHaplotypeCache
{
    fn key(record: &bam::Record) -> String {
        let s = String::from_utf8(record.qname().to_vec()).unwrap();
        return s;
    }

    pub fn update(&mut self, key: &String, record: &bam::Record) -> Option<i32> {
        if let Some(hi) = get_haplotag_from_record(record) {
            self.cache.insert(key.clone(), hi);
            return Some(hi);
        } else {
            return None;
        }
    }

    pub fn get(&mut self, record: &bam::Record) -> Option<i32> {
        let s = ReadHaplotypeCache::key(record);
        let x = self.cache.get(&s);
        match x {
            Some(value) => return Some(*value),
            None => return self.update(&s, record)
        }
    }
}

pub fn get_haplotag_from_record(record: &bam::Record) -> Option<i32> {
    match record.aux(b"HP") {
        Ok(value) => {
            if let Aux::I32(v) = value {
                return Some(v)
            } else if let Aux::U8(v) = value {
                return Some(v as i32) // whatshap encodes with U8
            } else {
                return None
            }
        }
        Err(_e) => return None
    }
}

pub struct GenomeRegions
{
    interval_trees: HashMap::<u32, IntervalTree<usize, usize>>
}

impl GenomeRegions
{
    pub fn from_bed(filename: & str, header_view: &bam::HeaderView) -> GenomeRegions {
        // read bed file into a data structure we can use to make intervaltrees from
        // this maps from tid to a vector of intervals, with an interval index for each
        let mut bed_reader = csv::ReaderBuilder::new().delimiter(b'\t').from_path(filename).expect("Could not open bed file");
        let mut region_desc_by_chr = HashMap::<u32, Vec<(Range<usize>, usize)>>::new();
        
        let mut region_count:usize = 0;
        for r in bed_reader.records() {
            let record = r.expect("Could not parse bed record");
            if let Some(tid) = header_view.tid(record[0].as_bytes()) {
                let start: usize = record[1].parse().unwrap();
                let end: usize = record[2].parse().unwrap();

                let region_desc = region_desc_by_chr.entry(tid).or_insert( Vec::new() );
                region_desc.push( (start..end, region_count) );
                region_count += 1;
            }
        }

        // build tid -> intervaltree map
        let mut regions = GenomeRegions { interval_trees: HashMap::<u32, IntervalTree<usize, usize>>::new() };
        for (tid, region_desc) in region_desc_by_chr {
            regions.interval_trees.insert(tid, region_desc.iter().cloned().collect());
        }
        return regions
    }

    pub fn contains(& self, tid: u32, position: usize) -> bool {
        if let Some(tree) = self.interval_trees.get( & tid ) {
            return tree.query_point(position).count() > 0;
        } else {
            return false;
        }
    }
}
