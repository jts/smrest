#! /usr/bin/env python

import argparse
import csv
import sys

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', nargs=1, help='smrest extract .tsv file')
    args = parser.parse_args()

    variant_types = [ "ERROR", "GERMLINE_HET", "GERMLINE_HOM" ]
    with open(args.input[0]) as f:
        data_by_vt = dict()

        # Initialize
        mutation_types = list()
        for c in "CT":
            for b in "ACGT":
                if c != b:
                    mutation_types.append(f"{c}>{b}")

        all_keys = list()
        for p in "ACGT":
            for mt in mutation_types:
                for s in "ACGT":
                    k = f"{p}[{mt}]{s}"
                    all_keys.append(k)

        for vt in variant_types:
            data_by_vt[vt] = dict()
            for k in all_keys:
                data_by_vt[vt][k] = 0

        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            c = row['canonical_context']
            mt = f"{c[0]}[{row['canonical_type']}]{c[-1]}"
            data_by_vt[row['variant_type']][mt] += 1

        header = [ "MutationType" ]
        for vt in variant_types:
            header.append(f"none_{vt}")

        print("\t".join(header))
        for mt in all_keys:
            out = [ mt ]
            for vt in variant_types:
                out.append(data_by_vt[vt][mt])
            print("\t".join([str(x) for x in out]))

if __name__ == "__main__":
    main()
