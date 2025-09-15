#!/usr/bin/env python3
import argparse
import gzip
import pandas as pd
from pathlib import Path

HEART_TISSUES = [
    "Heart - Left Ventricle",
    "Heart - Atrial Appendage",
]
TARGET_GENES = ["MYH7", "MYBPC3", "TNNT2", "TTN", "SCN5A", "KCNQ1"]


def load_gct(gct_path: str) -> pd.DataFrame:
    # GTEx GCT: first two header lines contain metadata
    with gzip.open(gct_path, "rt") as f:
        header1 = f.readline()
        header2 = f.readline()
        df = pd.read_csv(f, sep="\t")
    return df

def load_sample_attributes(attr_path: str) -> pd.DataFrame:
    df = pd.read_csv(attr_path, sep="\t")
    # Expect columns: SAMPID, SMTSD (Tissue), SMTS (Tissue type)
    return df[[c for c in df.columns if c in ("SAMPID", "SMTSD", "SMTS")]]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gct", required=True)
    ap.add_argument("--attrs", required=False, help="GTEx SampleAttributesDS.txt path; if provided, compute medians per tissue")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    df = load_gct(args.gct)
    # Normalize gene symbol column name variations
    name_col = "Name" if "Name" in df.columns else df.columns[0]
    desc_col = "Description" if "Description" in df.columns else df.columns[1]

    # Filter rows to target genes by symbol in Description
    df = df[df[desc_col].isin(TARGET_GENES)].copy()

    # If attributes provided, compute median per tissue for two heart tissues
    if args.attrs and Path(args.attrs).exists():
        attrs = load_sample_attributes(args.attrs)
        # Map SAMPID columns to tissue detailed (SMTSD)
        sample_cols = [c for c in df.columns if c not in (name_col, desc_col)]
        # Melt then join with attrs
        long = df.melt(id_vars=[name_col, desc_col], value_vars=sample_cols, var_name="SAMPID", value_name="TPM")
        long = long.merge(attrs, on="SAMPID", how="left")
        long = long[long["SMTSD"].isin(HEART_TISSUES)]
        # Median per gene and tissue
        med = long.groupby([desc_col, "SMTSD"], as_index=False)["TPM"].median().pivot(index=desc_col, columns="SMTSD", values="TPM")
        med = med.reindex(TARGET_GENES)
        med.reset_index(names=["gene"], inplace=True)
        med.to_csv(args.out, index=False)
    else:
        # No attributes: just output available columns for target genes
        df.to_csv(args.out, index=False)


if __name__ == "__main__":
    main()
