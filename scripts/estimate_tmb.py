#! /usr/bin/env python

from pysam import VariantFile
import argparse
import numpy as np
import csv
import sys

# return a vector of class probability, after applying post-hoc filters
def apply_filters(row, max_minor_evidence = 2):
    probs = np.array([ float(row['p_ref']), float(row['p_germline']), float(row['p_somatic']) ])

    # if this a non-reference call and the position is in gnomad, automatically set it to  
    if "variant_type" in row:
        if row["variant_type"] == "GNOMAD" and np.argmax(probs) != 0:
            probs = np.array( [0.0, 1.0, 0.0] )

    # check that the number of reads showing the variant on the minor haplotype
    # is below our threshold
    if int(row['hminor_variant_count']) >= max_minor_evidence:
        return None
    return probs

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', nargs=1, help='smrest extract .tsv file')
    args = parser.parse_args()

    max_opposing_haplotype_evidence = 2
    soft_sums = np.zeros(3)
    hard_sums = np.zeros(3)

    with open(args.input[0]) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            
            #
            # call filtering
            #
            probs = apply_filters(row)
            if probs is None:
                continue

            
            #
            # collect stats
            #

            # update soft sums
            soft_sums += probs

            idx = np.argmax(probs)
            hard_sums[idx] += 1

    tmb_soft = (soft_sums[2] / np.sum(soft_sums)) * 1000000
    tmb_hard = (hard_sums[2] / np.sum(hard_sums)) * 1000000

    print("soft_sum_ref\tsoft_sum_germline\tsoft_sum_somatic\ttmb_soft\thard_sum_ref\thard_sum_germline\thard_sum_somatic\ttmb_hard")
    print(f"{soft_sums[0]:.1f}\t{soft_sums[1]:.1f}\t{soft_sums[2]:.1f}\t{tmb_soft:.1f}\t{hard_sums[0]:.1f}\t{hard_sums[1]:.1f}\t{hard_sums[2]:.1f}\t{tmb_hard:.1f}")

if __name__ == "__main__":
    main()

