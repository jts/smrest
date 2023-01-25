//---------------------------------------------------------
// Copyright 2022 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
use rust_htslib::{bam, bam::record::Aux};
use std::collections::{HashMap};
use intervaltree::IntervalTree;
use core::ops::Range;
use crate::BcfHeader;
use crate::{ReadHaplotypeLikelihood, ReadMetadata};
use rust_htslib::bam::Read;
use longshot::util::GenomicInterval;
use longshot::variants_and_fragments::{Fragment, VarList};

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

pub fn get_phase_set_from_record(record: &bam::Record) -> Option<i32> {
    match record.aux(b"PS") {
        Ok(value) => {
            if let Aux::I32(v) = value {
                return Some(v)
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

pub fn add_contig_lines_to_vcf(vcf_header: &mut BcfHeader, bam_header: &rust_htslib::bam::HeaderView) {
    for tid in 0..bam_header.target_count() {
        let l = format!("##contig=<ID={},length={}", std::str::from_utf8(bam_header.tid2name(tid)).unwrap(), 
                                                     bam_header.target_len(tid).unwrap());
        vcf_header.push_record(l.as_bytes());
    }
}

pub fn get_chromosome_sequence(reference_genome: &str,
                               bam_header: &rust_htslib::bam::HeaderView,
                               tid: u32) -> String {

    let faidx = rust_htslib::faidx::Reader::from_path(reference_genome).expect("Could not read reference genome:");
    let chromosome_length = bam_header.target_len(tid).unwrap() as usize;
    let chromosome_name = std::str::from_utf8(bam_header.tid2name(tid)).unwrap();

    let mut chromosome_sequence = faidx.fetch_seq_string(&chromosome_name, 0, chromosome_length).unwrap();
    chromosome_sequence.make_ascii_uppercase();
    return chromosome_sequence;
}

// iterate over the bam, storing information we need in a hash table
pub fn populate_read_metadata_from_bam(bam: &mut rust_htslib::bam::IndexedReader,
                                       region: &GenomicInterval) -> HashMap::<String, ReadMetadata> {
    bam.fetch( (region.tid, region.start_pos, region.end_pos + 1) ).unwrap();
    let mut read_meta = HashMap::<String, ReadMetadata>::new();
    for r in bam.records() {
        let record = r.unwrap();
        // TODO: handle multiple alignments per read, not just primary
        
        // same criteria as longshot
        if record.is_quality_check_failed()
            || record.is_duplicate()
            || record.is_secondary()
            || record.is_unmapped()
            || record.is_supplementary()
        {
            continue;
        }

        let s = String::from_utf8(record.qname().to_vec()).unwrap();
        let rm = ReadMetadata { 
            haplotype_index: get_haplotag_from_record(&record),
            phase_set: get_phase_set_from_record(&record),
            strand_index: record.is_reverse() as i32,
            leading_softclips: record.cigar().leading_softclips(),
            trailing_softclips: record.cigar().trailing_softclips()
        };
        read_meta.insert(s.clone(), rm);
    }

    return read_meta;
}

// Take a vector of fragments (longshot reads that are annotated with allele likelihoods) and convert it to a vector
// of ReadHaplotypeLikelihoods for each variant in varlist
pub fn fragments_to_read_haplotype_likelihoods(varlist: &VarList,
                                               fragments: &Vec<Fragment>,
                                               read_meta: &HashMap::<String, ReadMetadata>) -> Vec::<Vec::<ReadHaplotypeLikelihood>> {
    let mut rhl_per_var = vec![ Vec::<ReadHaplotypeLikelihood>::new(); varlist.lst.len() ];
    for f in fragments {
        let id = f.id.as_ref().unwrap();

        if let Some( rm ) = read_meta.get(id) {
            for c in &f.calls {
                let var = &varlist.lst[c.var_ix];
                //println!("{}\t{}\t{}\t{}\t{}\t{}", &id, var.tid, var.pos0 + 1, var.alleles[0], var.alleles[1], c.allele);
            
                let rhl = ReadHaplotypeLikelihood { 
                    read_name: Some(id.clone()),
                    mutant_allele_likelihood: c.allele_scores[1],
                    base_allele_likelihood: c.allele_scores[0],
                    allele_call: var.alleles[c.allele as usize].clone(),
                    allele_call_qual: c.qual,
                    haplotype_index: rm.haplotype_index,
                    strand_index: rm.strand_index
                };
                rhl_per_var[c.var_ix].push(rhl);
            }
        }
    }
    return rhl_per_var;
}
