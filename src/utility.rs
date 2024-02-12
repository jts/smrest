//---------------------------------------------------------
// Copyright 2022 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
use rust_htslib::{bam, bcf, bam::record::Aux};
use rust_htslib::bam::ext::BamRecordExtensions;
use std::collections::{HashMap};
use intervaltree::IntervalTree;
use core::ops::Range;
use bio::stats::{Prob, LogProb, PHREDProb};
use crate::BcfHeader;
use crate::{ReadHaplotypeLikelihood, ReadMetadata};
use crate::LongshotParameters;
use crate::utility::bcf::record::GenotypeAllele::Phased;
use rust_htslib::bam::Read as BamRead;
use rust_htslib::bcf::Read as BcfRead;
use longshot::util::GenomicInterval;
use longshot::variants_and_fragments::{Fragment, Var, VarList, VarFilter};
use longshot::util::u8_to_string;
use longshot::genotype_probs::{Genotype, GenotypeProbs};
use longshot::extract_fragments::extract_fragments;

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

// calculate the error rate of the read sequence over the first n_read_bases of the sequence
pub fn calculate_mismatch_rate_until<'a, I>(iter: I, 
                                         read_seq: &rust_htslib::bam::record::Seq,
                                         ref_seq: &Vec<u8>,
                                         n_read_bases: usize) -> f64 
where
    I: Iterator<Item = &'a [Option<i64>; 2]>
{

    let mut base_count = 0;
    let mut error_count = 0;

    for i in iter {

        if let Some(read_pos) = i[0] {
            if let Some(ref_pos) = i[1] {
                let read_base = read_seq[read_pos as usize];
                let ref_base = ref_seq[ref_pos as usize];

                base_count += 1;
                if read_base != ref_base {
                    error_count += 1;
                }
            }


            if base_count == n_read_bases {
                break;
            }
        }
    }

    error_count as f64 / base_count as f64
}

pub fn assign_reads_to_haplotypes(input_bam: &String,
                                  region: &GenomicInterval,
                                  reference_genome: &String,
                                  phased_vcf: &str,
                                  longshot_parameters: &LongshotParameters) -> HashMap::<String, i32> {

    let mut haplotype_assignments = HashMap::<String, i32>::new();
    let mut bam = bam::IndexedReader::from_path(input_bam).unwrap();
    bam.fetch( (region.tid, region.start_pos, region.end_pos + 1) ).unwrap();
    let header_view = bam.header();

    let bcf_records: Vec<bcf::Record> = read_vcf(&phased_vcf.to_owned())
                                           .into_iter().filter(|x| is_phased(x)).collect();

    // convert to longshot var here
    let all_vars = bcf_records.iter().map(|x| bcf2longshot(x, header_view)).collect();

    let target_names = header_view.target_names().iter().map(|x| std::str::from_utf8(x).unwrap().to_owned()).collect();
    let mut varlist = VarList::new(all_vars, target_names).expect("Could not construct var list");
    
    let frags = extract_fragments(&input_bam.to_owned(),
                                  &reference_genome.to_owned(),
                                  &mut varlist,
                                  &Some(region.clone()),
                                  &longshot_parameters.extract_fragment_parameters,
                                  &longshot_parameters.alignment_parameters.as_ref().unwrap()).unwrap();
    
    for f in frags {
        let id = f.id.as_ref().unwrap();

        let mut h0_likelihood = LogProb::ln_one();
        let mut h1_likelihood = LogProb::ln_one();

        let mut h0_match = 0;
        let mut h1_match = 0;

        for c in &f.calls {
            let _var = &varlist.lst[c.var_ix];
            let h0_idx = bcf_records[c.var_ix].genotypes().expect("Error reading genotypes").get(0)[0].index().unwrap() as usize;
            assert!(h0_idx == 0 || h0_idx == 1);

            h0_likelihood += c.allele_scores[h0_idx];
            h1_likelihood += c.allele_scores[1 - h0_idx];

            h0_match += (h0_idx == c.allele as usize) as u32;
            h1_match += ( (1 - h0_idx) == c.allele as usize) as u32;

            //println!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", &id, var.tid, var.pos0 + 1, var.alleles[0], var.alleles[1], c.allele, h0_idx, *c.allele_scores[h0_idx], *c.allele_scores[1 - h0_idx]);
        }

        let lp_sum = LogProb::ln_add_exp(h0_likelihood, h1_likelihood);

        let (htag, lp_assignment, max_match) = 
            if h0_match == 0 && h1_match == 0 {
                (None, LogProb::from(Prob(0.0)), 0)
            } else if h0_likelihood > h1_likelihood {
                (Some(1), h0_likelihood - lp_sum, h0_match)
            } else {
                (Some(2), h1_likelihood - lp_sum, h1_match)
            };

        // cap phasing qual at 60
        let lp_wrong = lp_assignment.ln_one_minus_exp();
        let qual = if *lp_wrong > -13.0 { 
            PHREDProb::from(lp_wrong)
        } else {
            PHREDProb(60.0)
        };

        let min_agreement = 0.9;
        let agreement = max_match as f64 / (h0_match + h1_match) as f64;

        let mut qc = "FAIL";
        let min_qual = 20.0;
        if htag.is_some() && *qual > min_qual && agreement > min_agreement {
           haplotype_assignments.insert(id.clone(), htag.unwrap());
           qc = "PASS";
        }

        //println!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", &id, *h0_likelihood, *h1_likelihood, h0_match, h1_match, agreement, htag.unwrap_or(-1), *qual, qc);
    }

    return haplotype_assignments;
}

// iterate over the bam, storing information we need in a hash table
pub fn populate_read_metadata_from_bam(input_bam: &String,
                                       region: &GenomicInterval,
                                       chromosome_bytes: Option<&Vec<u8>>,
                                       phased_vcf: Option<&str>,
                                       longshot_parameters: Option<&LongshotParameters>,
                                       reference_genome: Option<&String>) -> HashMap::<String, ReadMetadata> {
    
    // if a phased VCF file is passed in we use it to re-compute the haplotype assignments
    // for each read. Otherwise this is read from the HP tag of the bam. 
    let haplotype_assignments = if let Some(f) = phased_vcf {
        Some(assign_reads_to_haplotypes(input_bam, region, reference_genome.unwrap(), f, longshot_parameters.unwrap()))
    } else {
        None
    };

    let mut bam = bam::IndexedReader::from_path(input_bam).unwrap();
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
        let (prefix, suffix) = if let Some(_ref_seq) = chromosome_bytes {

            let ap: Vec<[Option<i64>; 2]> = record.aligned_pairs_full().collect();
            let prefix_mmr = calculate_mismatch_rate_until(ap.iter(), 
                                                       &record.seq(),
                                                       chromosome_bytes.unwrap(),
                                                       500);
            let suffix_mmr = calculate_mismatch_rate_until(ap.iter().rev(), 
                                                       &record.seq(),
                                                       chromosome_bytes.unwrap(),
                                                       500);
            //println!("{}\t{}\t{}\n", s, prefix_mmr, suffix_mmr);
            (Some(prefix_mmr), Some(suffix_mmr))
        } else {
            (None, None)
        };

        let hi = if haplotype_assignments.is_some() {
            haplotype_assignments.as_ref().unwrap().get(&s).cloned()
        } else {
            get_haplotag_from_record(&record)
        };

        let rm = ReadMetadata { 
            haplotype_index: hi,
            phase_set: get_phase_set_from_record(&record),
            strand_index: record.is_reverse() as i32,
            leading_softclips: record.cigar().leading_softclips(),
            trailing_softclips: record.cigar().trailing_softclips(),
            prefix_mismatch_rate: prefix,
            suffix_mismatch_rate: suffix
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
            /*
            if rm.prefix_mismatch_rate.unwrap_or(0.0) >= max_mismatch_rate ||
               rm.suffix_mismatch_rate.unwrap_or(0.0) >= max_mismatch_rate 
            {
                continue;
            }
            */
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

pub fn read_bcf(bcf_filename: &String,
                region_opt: Option<GenomicInterval>) -> Vec<bcf::Record> {
    let mut vec = Vec::new();
    let mut bcf = bcf::IndexedReader::from_path(bcf_filename).expect("Error opening file.");

    if let Some(region) = region_opt {
        bcf.fetch( region.tid, region.start_pos as u64, Some(region.end_pos as u64) ).unwrap();
    }

    for record_result in bcf.records() {
        let record = record_result.expect("Fail to read record");
        vec.push(record);
    }
    return vec;
}

pub fn read_vcf(vcf_filename: &String) -> Vec<bcf::Record> {
    let mut vec = Vec::new();
    let mut bcf = bcf::Reader::from_path(vcf_filename).expect("Error opening file.");

    for record_result in bcf.records() {
        let record = record_result.expect("Fail to read record");
        vec.push(record);
    }
    return vec;
}


pub fn is_genotype_het(record: &bcf::Record) -> bool {
    let gts = record.genotypes().expect("Error reading genotypes");
    
    // ensure a single sample genotype field
    let sample_count = usize::try_from(record.sample_count()).unwrap();
    assert!(sample_count == 1);

    let mut n_ref = 0;
    let mut n_alt = 0;

    for gta in gts.get(0).iter() {
        match gta.index() {
            Some(0) => n_ref += 1, 
            Some(1) => n_alt += 1,
            _ => (),
        }
    }

    return n_ref == 1 && n_alt == 1;
}

pub fn is_phased(record: &bcf::Record) -> bool {
    let gts = record.genotypes().expect("Error reading genotypes");
    
    // ensure a single sample genotype field
    let sample_count = usize::try_from(record.sample_count()).unwrap();
    assert!(sample_count == 1);

    let _n_ref = 0;
    let _n_alt = 0;

    let mut is_phased = false;
    for gta in gts.get(0).iter() {
        match gta {
            Phased(_) => is_phased = true,
            _ => (),
        }
    }

    return is_phased;
}


// convert a bcf::Record into longshot's internal format 
pub fn bcf2longshot(record: &bcf::Record, bam_header: &bam::HeaderView) -> Var {

    // this is all derived from longshot's variants_and_fragments.rs
    let rid = record.rid().expect("Could not find variant rid");
    let chrom = record.header().rid2name(rid).expect("Could not find rid in header");
    
    let mut alleles: Vec<String> = vec![];
    for a in record.alleles().iter() {
        let s = u8_to_string(a).expect("Could not convert allele to string");
        alleles.push(s);
    }

    Var {
        ix: 0,
        tid: bam_header.tid(chrom).expect("Could not find chromosome in bam header") as u32,
        pos0: record.pos() as usize,
        alleles: alleles.clone(),
        dp: 0,
        allele_counts: vec![0; alleles.len()],
        allele_counts_forward: vec![0; alleles.len()],
        allele_counts_reverse: vec![0; alleles.len()],
        ambiguous_count: 0,
        qual: 0.0,
        filter: VarFilter::Pass,
        genotype: Genotype(0, 0),
        //unphased: false,
        gq: 0.0,
        unphased_genotype: Genotype(0, 0),
        unphased_gq: 0.0,
        genotype_post: GenotypeProbs::uniform(alleles.len()),
        phase_set: None,
        strand_bias_pvalue: 1.0,
        mec: 0,
        mec_frac_variant: 0.0, // mec fraction for this variant
        mec_frac_block: 0.0,   // mec fraction for this haplotype block
        mean_allele_qual: 0.0,
        dp_any_mq: 0,
        mq10_frac: 0.0,
        mq20_frac: 0.0,
        mq30_frac: 0.0,
        mq40_frac: 0.0,
        mq50_frac: 0.0
    }
    
}
