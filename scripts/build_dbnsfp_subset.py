#!/usr/bin/env python3
import argparse
import gzip
import os
import subprocess
import sys
from collections import defaultdict
from pathlib import Path


def norm_chr(ch: str) -> str:
    ch = ch.strip()
    if ch.startswith("chr"):
        ch = ch[3:]
    if ch in ("MT", "Mt", "mt"):
        return "M"
    return ch


def parse_vcf_variants(vcf_path: str):
    need = defaultdict(set)  # chr -> set of (pos, ref, alt)
    opener = gzip.open if vcf_path.endswith(".gz") else open
    with opener(vcf_path, "rt", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line or line[0] == "#":
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            ch, pos, vid, ref, alt = parts[:5]
            ch = norm_chr(ch)
            # ALT can be comma-separated; add each
            for a in alt.split(","):
                try:
                    p = int(pos)
                except ValueError:
                    continue
                need[ch].add((p, ref, a))
    return need


def guess_chr_file(dbnsfp_dir: Path, ch: str) -> Path:
    candidates = [
        dbnsfp_dir / f"dbNSFP4.4a_variant.chr{ch}.gz",
        dbnsfp_dir / "dbNSFP4.4a" / f"dbNSFP4.4a_variant.chr{ch}.gz",
    ]
    for p in candidates:
        if p.exists():
            return p
    raise FileNotFoundError(f"Could not find dbNSFP per-chromosome file for chr {ch}")


def main():
    ap = argparse.ArgumentParser(description="Build a dbNSFP subset file for variants in a VCF")
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--dbnsfp-dir", required=True)
    ap.add_argument("--out", required=True, help="Output path for bgzipped merged file (e.g., dbNSFP4.4a.txt.gz)")
    args = ap.parse_args()

    vcf = Path(args.vcf)
    dbdir = Path(args.dbnsfp_dir)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    need = parse_vcf_variants(str(vcf))
    if not need:
        print("No variants found in VCF; nothing to do", file=sys.stderr)
        sys.exit(1)

    # Header from chr1 (or first available chromosome in VCF)
    first_chr = next(iter(need.keys()))
    first_file = guess_chr_file(dbdir, first_chr)
    with gzip.open(first_file, "rt", encoding="utf-8", errors="ignore") as f:
        header = f.readline()
    if not header or "\t" not in header:
        print("Failed to read dbNSFP header", file=sys.stderr)
        sys.exit(1)

    tmp_txt = out.with_suffix("")  # remove .gz
    if tmp_txt.exists():
        tmp_txt.unlink()

    written = 0
    with open(tmp_txt, "wt", encoding="utf-8") as out_f:
        out_f.write(header)
        for ch, triples in need.items():
            chr_file = guess_chr_file(dbdir, ch)
            print(f"Scanning {chr_file.name} for {len(triples)} variants...", file=sys.stderr)
            # Build lookup by position -> set of (ref, alt)
            by_pos = defaultdict(set)
            for (p, ref, alt) in triples:
                by_pos[p].add((ref, alt))
            with gzip.open(chr_file, "rt", encoding="utf-8", errors="ignore") as f:
                hdr_seen = False
                for line in f:
                    if not hdr_seen:
                        hdr_seen = True
                        continue  # skip header line; we've already written header
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 4:
                        continue
                    try:
                        pos = int(parts[1])
                    except ValueError:
                        continue
                    ref, alt = parts[2], parts[3]
                    if (ref, alt) in by_pos.get(pos, ()):  # set lookup
                        out_f.write(line)
                        written += 1

    if written == 0:
        print("No matching dbNSFP records found for VCF variants", file=sys.stderr)
        sys.exit(2)

    # bgzip and tabix
    if out.exists():
        out.unlink()
    subprocess.run(["bgzip", "-f", "-@", "4", "-c", str(tmp_txt)], check=True, stdout=open(out, "wb"))
    subprocess.run(["tabix", "-f", "-s", "1", "-b", "2", "-e", "2", str(out)], check=True)
    tmp_txt.unlink(missing_ok=True)
    print(f"Wrote {out} with {written} records")


if __name__ == "__main__":
    main()
