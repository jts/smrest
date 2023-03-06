//---------------------------------------------------------
// Copyright 2023 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
use core::ops::Range;
use std::collections::{HashMap, BTreeSet};
use rust_htslib::{bam, bcf, bam::Read};
use rust_htslib::bam::ext::BamRecordExtensions;
use crate::utility::*;
use crate::LongshotParameters;
use longshot::extract_fragments::extract_fragments;
use longshot::variants_and_fragments::VarList;
use longshot::util::parse_region_string;

pub fn select_reads(input_bam: &str, output_bam: &str, region_str: &str, vcf: &str, reference_genome: &str) {

    let region = parse_region_string(Some(region_str), &input_bam.to_owned()).expect("Could not parse region").unwrap();
    let mut bam = bam::IndexedReader::from_path(input_bam).unwrap();
    let header_view = bam.header().clone();
    let bcf_records: Vec<bcf::Record> = read_bcf(&vcf.to_owned(), None)
                                           .into_iter().filter(|x| is_genotype_het(x)).collect();

    // build an b-tree for each chromosome containing the positions of each heterozygous variant
    let mut chromosome_hets = Vec::new();
    for _ in 0..header_view.target_count() {
        chromosome_hets.push(BTreeSet::new());
    }

    for record in bcf_records.iter() {
        let rid = record.rid().expect("Could not find variant rid");
        let chrom = record.header().rid2name(rid).expect("Could not find rid in header");
        let tid = header_view.tid(chrom).expect("Could not find chromosome in bam header") as usize;
        assert!(tid < chromosome_hets.len());
        let pos = record.pos() as usize;
        chromosome_hets[tid].insert(pos);
    }

    for (tid, hets) in chromosome_hets.iter().enumerate() {
        println!("{}\t{}", tid, hets.len());
    }

    let out_header = bam::Header::from_template(&header_view);
    let mut out = bam::Writer::from_path(output_bam, &out_header, bam::Format::Bam).unwrap();

    let min_hets_spanned = 3;

    bam.fetch(bam::FetchDefinition::All);
    for r in bam.records() {
        let record = r.unwrap();

        if record.is_unmapped() {
            continue;
        }

        let qname = std::str::from_utf8(record.qname()).unwrap();
        let start_position = record.pos() as usize;
        let end_position = record.reference_end() as usize;

        let hets_spanned = chromosome_hets[record.tid() as usize].range( start_position..end_position );
        let n_hets = hets_spanned.count();
        let contig = String::from_utf8_lossy(header_view.tid2name(record.tid() as u32));
        if n_hets >= min_hets_spanned {
            /*println!("{}\t{}\t{}\t{}\t{}", 
                qname, contig, start_position, end_position, n_hets);
            */
            out.write(&record).unwrap();
        }
     }
}
