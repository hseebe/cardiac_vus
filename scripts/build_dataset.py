#!/usr/bin/env python3
import argparse
import gzip
import pandas as pd
import re


def open_text(path):
    if path.endswith('.gz'):
        return gzip.open(path, 'rt', encoding='utf-8', errors='ignore')
    return open(path, 'rt', encoding='utf-8', errors='ignore')


def labels_from_clinvar(vcf_path: str):
    label_map = {
        'Pathogenic': 1, 'Likely_pathogenic': 1,
        'Benign': 0, 'Likely_benign': 0,
    }
    out = {}
    with open_text(vcf_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 8:
                continue
            chrom, pos, _id, ref, alt = parts[:5]
            info = parts[7]
            m = re.search(r'CLNSIG=([^;]+)', info)
            if not m:
                continue
            sigs = {s.strip() for s in m.group(1).split(',')}
            label = None
            if 'Pathogenic' in sigs or 'Likely_pathogenic' in sigs:
                label = 1
            elif 'Benign' in sigs or 'Likely_benign' in sigs:
                label = 0
            if label is not None:
                out[(chrom.lstrip('chr'), int(pos), ref, alt)] = label
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--features', required=True)
    ap.add_argument('--clinvar', required=True)
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    df = pd.read_csv(args.features)
    labels = labels_from_clinvar(args.clinvar)
    df['label'] = [labels.get((str(ch), int(p), r, a)) if pd.notnull(p) else None for ch, p, r, a in zip(df['chrom'], df['pos'], df['ref'], df['alt'])]
    df.to_csv(args.out, index=False)
    print(f"Wrote {args.out} with {len(df)} rows; labeled={df['label'].notna().sum()}")


if __name__ == '__main__':
    main()
