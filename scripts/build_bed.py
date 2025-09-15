#!/usr/bin/env python3
import argparse
import sys
from typing import List, Dict

GENES = [
    "MYH7",
    "MYBPC3",
    "TNNT2",
    "TTN",
    "SCN5A",
    "KCNQ1",
]

# Minimal canonical coordinates (GRCh38) hard-coded as a fallback.
# Format: gene -> (chrom, start, end). start is 0-based, end is 1-based as BED expects half-open [start, end)
# These are rough gene span ranges; for precise canonical transcripts, consider querying Ensembl REST.
FALLBACK_COORDS = {
    "MYH7": ("14", 23872933, 23904751),
    "MYBPC3": ("11", 47317013, 47363695),
    "TNNT2": ("1", 201359228, 201378015),
    "TTN": ("2", 178525989, 178807423),
    "SCN5A": ("3", 38589545, 38691305),
    "KCNQ1": ("11", 2589737, 2817649),
}

HEADER = [
    "#chrom",
    "start",
    "end",
    "gene",
]

def write_bed(path: str, records: List[List[str]]):
    with open(path, "w") as f:
        for rec in records:
            f.write("\t".join(map(str, rec)) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True, help="Output BED path")
    args = ap.parse_args()

    records: List[List[str]] = []
    for g in GENES:
        chrom, start, end = FALLBACK_COORDS[g]
        records.append([chrom, start, end, g])

    write_bed(args.out, records)
    print(f"Wrote {len(records)} intervals to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
