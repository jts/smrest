#! /usr/bin/env python

from pysam import VariantFile
from estimate_tmb import apply_filters
from annotate import load_somatic_annotations
from collections import defaultdict

import argparse
import numpy as np
import csv
import sys

def load_somatic_truth_at_called(filename, called_positions, min_vaf, vaf_info_name, annotations):
    vcf_in = VariantFile(filename)
    for record in vcf_in.fetch():

        if (record.chrom, record.pos) not in called_positions:
            continue

        # skip multi-allelic and indels
        if len(record.alts) != 1 or len(record.ref) > 1 or len(record.alts[0]) > 1:
            continue

        # skip low vaf calls, maybe
        vaf = None
        if vaf_info_name != "":
            vaf = float(record.info[vaf_info_name])
        else:
            tumor_sample_idx = 1 # TODO: don't hardcode
            vaf = float(record.samples[tumor_sample_idx]['AF'][0])
        
        if vaf < min_vaf:
            continue

        key = (record.chrom, record.pos, record.ref, record.alts[0])
        annotations[key].append("truth")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--somatic-vcf', type=str)
    parser.add_argument('--min-truth-vaf', type=float, default=0.0)
    parser.add_argument('--truth-vaf-name', type=str, default="")
    parser.add_argument('input', nargs=1, help='smrest extract .tsv file')
    args = parser.parse_args()
    
    annotations = defaultdict(list)
    called_positions = dict()

    with open(args.input[0]) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            key = (row['chromosome'], int(row['position']), row['reference_base'], row['variant_base'])
            called_positions[ (row['chromosome'], int(row['position'])) ] = 1

            probs = apply_filters(row)
            idx = np.argmax(probs)
            if idx == 2:
                key = (row['chromosome'], int(row['position']), row['reference_base'], row['variant_base'])
                annotations[key].append("called")
    
    load_somatic_truth_at_called(args.somatic_vcf, called_positions, args.min_truth_vaf, args.truth_vaf_name, annotations)

    num_truth = 0
    num_called = 0
    tp = 0
    fp = 0
    fn = 0

    for (key, values) in annotations.items():
        num_truth += "truth" in values
        num_called += "called" in values

        tp += "truth" in values and "called" in values
        fp += "truth" not in values and "called" in values
        fn += "truth" in values and "called" not in values

    bases_covered = len(called_positions)
    print("bases_covered\tnum_truth\tnum_called\ttp\tfp\tfn")
    print(f"{bases_covered}\t{num_truth}\t{num_called}\t{tp}\t{fp}\t{fn}")
if __name__ == "__main__":
    main()

