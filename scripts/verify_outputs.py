#!/usr/bin/env python3
import argparse
import os
import sys

def check_file(path: str, desc: str) -> bool:
    ok = os.path.isfile(path) and os.path.getsize(path) > 0
    print(f"{desc}: {'OK' if ok else 'MISSING'} -> {path}")
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--clinvar', required=True)
    ap.add_argument('--cardiac', required=True)
    ap.add_argument('--vep', required=True)
    ap.add_argument('--gtex', required=True)
    ap.add_argument('--cons-phyloP', required=True)
    ap.add_argument('--cons-phastCons', required=True)
    ap.add_argument('--dbnsfp-dir', required=True)
    args = ap.parse_args()

    ok = True
    ok &= check_file(args.clinvar, 'ClinVar VCF')
    ok &= check_file(args.cardiac, 'Cardiac genes VCF')
    ok &= check_file(args.vep, 'VEP annotated VCF')
    ok &= check_file(args.gtex, 'GTEx heart subset')
    ok &= check_file(args.cons_phyloP, 'phyloP bigWig')
    ok &= check_file(args.cons_phastCons, 'phastCons bigWig')

    # dbNSFP check: directory exists and contains at least one file
    dbnsfp_ok = os.path.isdir(args.dbnsfp_dir) and any(
        os.path.isfile(os.path.join(args.dbnsfp_dir, f)) for f in os.listdir(args.dbnsfp_dir)
    )
    print(f"dbNSFP dir: {'OK' if dbnsfp_ok else 'MISSING/EMPTY'} -> {args.dbnsfp_dir}")
    ok &= dbnsfp_ok

    sys.exit(0 if ok else 1)

if __name__ == '__main__':
    main()
