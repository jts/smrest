#! /usr/bin/env python

from pysam import VariantFile
import argparse
import csv
import sys
from intervaltree import Interval, IntervalTree

def load_germline_annotations(filename, annotations):
    vcf_in = VariantFile(filename)
    for record in vcf_in.fetch():

        # skip multi-allelic        
        if len(record.alts) != 1:
            continue

        key = (record.chrom, record.pos, record.ref, record.alts[0])
        gt = record.samples[0]['GT']
        value = None
        if gt[0] == 0 and gt[1] == 1:
            value = "GERMLINE_HET"
        elif gt[0] == 1 and gt[1] == 1:
            value = "GERMLINE_HOM"
        else:
            value = "UNHANDLED_GENOTYPE"  
        annotations[key] = value

def load_gnomad_annotations(filename, annotations):
    vcf_in = VariantFile(filename)
    for record in vcf_in.fetch():

        # skip multi-allelic        
        if len(record.alts) != 1:
            continue

        key = (record.chrom, record.pos, record.ref, record.alts[0])
        annotations[key] = "GNOMAD"

def load_somatic_annotations(filename, annotations):
    vcf_in = VariantFile(filename)
    for record in vcf_in.fetch():

        # skip multi-allelic        
        if len(record.alts) != 1:
            continue

        key = (record.chrom, record.pos, record.ref, record.alts[0])
        annotations[key] = "SOMATIC"

def load_bed(filename):
    trees = dict()
    
    with open(filename) as f:
        for (idx, line) in enumerate(f):
            if line[0] == '#':
                continue
            
            fields = line.rstrip().split()
            chromosome = fields[0]
            start = fields[1]
            end = fields[2]
            if chromosome not in trees:
                trees[chromosome] = IntervalTree()
            trees[chromosome].addi(int(start), int(end), idx)
    return trees

def is_in_region(regions, chromosome, position):
    if regions is not None:
        return int(len(regions[chromosome].at(position)) > 0)
    else:
        return 0

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--gnomad-vcf', type=str)
    parser.add_argument('--germline-vcf', type=str)
    parser.add_argument('--somatic-vcf', type=str)
    parser.add_argument('--target-bed', type=str)
    parser.add_argument('--giab-hc-bed', type=str)
    parser.add_argument('input', nargs=1, help='smrest extract .tsv file')
    args = parser.parse_args()

    # read vcf and store annotations
    annotations = dict()
    if args.germline_vcf:
        load_germline_annotations(args.germline_vcf, annotations)
    if args.somatic_vcf:
        load_somatic_annotations(args.somatic_vcf, annotations)
    if args.gnomad_vcf:
        load_gnomad_annotations(args.gnomad_vcf, annotations)

    target_regions = None
    if args.target_bed:
        target_regions = load_bed(args.target_bed)
    
    hc_regions = None
    if args.giab_hc_bed:
        hc_regions = load_bed(args.giab_hc_bed)

    with open(args.input[0]) as f:
        reader = csv.DictReader(f, delimiter='\t')
        header = reader.fieldnames
        header.append("variant_type")
        if target_regions is not None:
            header.append("on_target")
        if hc_regions is not None:
            header.append("in_high_confidence_region")

        print("\t".join(header))

        for row in reader:
            c = row['chromosome']
            pos = int(row['position'])
            key = (c, pos, row['reference_base'], row['variant_base'])
            annotation = "UNKNOWN"
            if key in annotations:
                annotation = annotations[key]
            
            row["variant_type"] = annotation
            
            if target_regions is not None:
                row['on_target'] = str(is_in_region(target_regions, c, pos))
            if hc_regions is not None:
                row['in_high_confidence_region'] = str(is_in_region(hc_regions, c, pos))
            
            values = [ row[k] for k in header ]
            print("\t".join(values))

if __name__ == "__main__":
    main()
