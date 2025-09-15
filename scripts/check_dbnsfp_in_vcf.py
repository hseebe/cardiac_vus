#!/usr/bin/env python3
import argparse
import gzip
import re

FIELDS_OF_INTEREST = [
    "SIFT_pred",
    "Polyphen2_HVAR_pred",
    "Polyphen2_HDIV_pred",
    "GERP++_RS",
    "phyloP100way_vertebrate",
    "gnomADg_AF",
    "SYMBOL",
    "IMPACT",
    "Consequence",
]


def open_text(path):
    if path.endswith('.gz'):
        return gzip.open(path, 'rt', encoding='utf-8', errors='ignore')
    return open(path, 'rt', encoding='utf-8', errors='ignore')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--vcf', required=True, help='VEP-annotated VCF')
    ap.add_argument('--examples', type=int, default=3, help='Number of example variants to print')
    args = ap.parse_args()

    path = args.vcf
    csq_fields = []
    with open_text(path) as f:
        for line in f:
            if line.startswith('##INFO=<ID=CSQ'):
                m = re.search(r'Format: (.*)\">', line)
                if m:
                    csq_fields = m.group(1).split('|')
                    found = [x for x in FIELDS_OF_INTEREST if x in csq_fields]
                    print('Header CSQ fields present:', ', '.join(found))
            if line.startswith('#CHROM'):
                break

    printed = 0
    with open_text(path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 8:
                continue
            info = parts[7]
            m = re.search(r'CSQ=([^;]+)', info)
            if not m:
                continue
            for csq in m.group(1).split(','):
                vals = csq.split('|')
                if len(vals) != len(csq_fields):
                    continue
                row = {k: v for k, v in zip(csq_fields, vals)}
                if any(row.get(k) and row.get(k) != '.' for k in FIELDS_OF_INTEREST if k in row):
                    print('Variant', parts[0], parts[1], parts[3], parts[4])
                    print({k: row.get(k) for k in FIELDS_OF_INTEREST if k in row})
                    printed += 1
                    break
            if printed >= args.examples:
                break

    if printed == 0:
        print('No example variants with populated dbNSFP fields were found in the first pass.')


if __name__ == '__main__':
    main()
