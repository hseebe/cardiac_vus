#!/usr/bin/env python3
import argparse
import gzip
import re
import pandas as pd
from pathlib import Path


def open_text(path):
    if str(path).endswith('.gz'):
        return gzip.open(path, 'rt', encoding='utf-8', errors='ignore')
    return open(path, 'rt', encoding='utf-8', errors='ignore')


def parse_csq_header(vcf_path):
    csq_fields = []
    with open_text(vcf_path) as f:
        for line in f:
            if line.startswith('##INFO=<ID=CSQ'):
                m = re.search(r'Format: (.*)\">', line)
                if m:
                    csq_fields = m.group(1).split('|')
            if line.startswith('#CHROM'):
                break
    return csq_fields


def build_features(vcf_path: str, gtex_csv: str, clinvar_vcf: str, out_csv: str):
    csq_fields = parse_csq_header(vcf_path)
    want = [
        'Consequence', 'IMPACT', 'SYMBOL', 'gnomADg_AF',
        'SIFT_pred', 'Polyphen2_HVAR_pred', 'Polyphen2_HDIV_pred',
        'GERP++_RS', 'phyloP100way_vertebrate',
    ]

    records = []
    with open_text(vcf_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 8:
                continue
            chrom, pos, vid, ref, alt = parts[:5]
            info = parts[7]
            m = re.search(r'CSQ=([^;]+)', info)
            if not m:
                continue
            best = None
            for csq in m.group(1).split(','):
                vals = csq.split('|')
                if len(vals) != len(csq_fields):
                    continue
                row = {k: v for k, v in zip(csq_fields, vals)}
                # Pick canonical if available else first
                if row.get('CANONICAL') == 'YES':
                    best = row
                    break
                if best is None:
                    best = row
            if best is None:
                continue

            rec = {
                'chrom': chrom.lstrip('chr'),
                'pos': int(pos),
                'ref': ref,
                'alt': alt,
                'id': vid,
            }
            for k in want:
                rec[k] = best.get(k)
            records.append(rec)

    df = pd.DataFrame.from_records(records)
    # Coercions
    for col in ['gnomADg_AF', 'GERP++_RS', 'phyloP100way_vertebrate']:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')

    # Join GTEx (expect columns: gene, Heart - Left Ventricle, Heart - Atrial Appendage)
    gtex = pd.read_csv(gtex_csv)
    # Normalize GTEx column names
    rename_map = {}
    for c in gtex.columns:
        if c.strip().lower().startswith('heart - left ventricle'):
            rename_map[c] = 'heart_lv_tpm'
        if c.strip().lower().startswith('heart - atrial appendage'):
            rename_map[c] = 'heart_aa_tpm'
    if 'gene' not in gtex.columns:
        # Try typical structure: gene in first column name "Description"
        if 'Description' in gtex.columns:
            gtex = gtex.rename(columns={'Description': 'gene'})
        elif 'Name' in gtex.columns:
            gtex = gtex.rename(columns={'Name': 'gene'})
    gtex = gtex.rename(columns=rename_map)
    gtex = gtex[['gene'] + [c for c in ['heart_lv_tpm', 'heart_aa_tpm'] if c in gtex.columns]].drop_duplicates('gene')
    df = df.merge(gtex, left_on='SYMBOL', right_on='gene', how='left').drop(columns=['gene'], errors='ignore')

    # Labels from ClinVar: Pathogenic/Likely_pathogenic -> 1, Benign/Likely_benign -> 0
    label_map = {
        'Pathogenic': 1, 'Likely_pathogenic': 1,
        'Benign': 0, 'Likely_benign': 0,
    }
    # Build set of (chrom,pos,ref,alt) -> label from cardiac_genes.vcf.gz
    labels = {}
    with open_text(clinvar_vcf) as f:
        for line in f:
            if line.startswith('#'):
                continue
            chrom, pos, vid, ref, alt, *_ = line.rstrip('\n').split('\t')
            info = line.split('\t', 8)[7]
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
                labels[(chrom.lstrip('chr'), int(pos), ref, alt)] = label

    df['label'] = [labels.get((str(ch), int(p), r, a)) if pd.notnull(p) else None for ch, p, r, a in zip(df['chrom'], df['pos'], df['ref'], df['alt'])]

    Path(out_csv).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_csv, index=False)
    print(f"Wrote features to {out_csv} ({len(df)} rows)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--vcf', required=True, help='results/annotated_variants.vep.vcf')
    ap.add_argument('--gtex', required=True, help='data/gtex/gtex_heart_subset.csv')
    ap.add_argument('--clinvar', required=True, help='cardiac_genes.vcf.gz')
    ap.add_argument('--out', required=True, help='results/variant_features.csv')
    args = ap.parse_args()
    build_features(args.vcf, args.gtex, args.clinvar, args.out)


if __name__ == '__main__':
    main()
