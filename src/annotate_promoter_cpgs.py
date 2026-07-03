#!/usr/bin/env python3
"""
annotate_promoter_cpgs.py - build the CpGs_with_promoters.bed annotation
file consumed by main.py (--annotation).

Download the GENCODE GTF first:
    wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz
    gunzip gencode.v49.annotation.gtf.gz

Usage:
    python annotate_promoter_cpgs.py \
        --gtf reference/gencode.v49.annotation.gtf \
        --cpg-islands reference/CPGIslands.bed \
        --outdir reference/

Requires `bedtools` on PATH.
"""

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd

# Repo root (this script lives in <repo>/src/)
ROOT_DIR = Path(__file__).resolve().parent.parent

GTF_COLUMNS = [
    "chrom", "source", "feature", "start", "end",
    "score", "strand", "frame", "attribute",
]


def parse_args():
    p = argparse.ArgumentParser(
        description="Build CpG island / promoter annotation BED from a GENCODE GTF.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--gtf", type=Path, default=ROOT_DIR / "reference" / "gencode.v49.annotation.gtf",
                    help="GENCODE annotation GTF")
    p.add_argument("--cpg-islands", type=Path, default=ROOT_DIR / "reference" / "CPGIslands.bed",
                    help="Input CpG islands BED file")
    p.add_argument("--outdir", type=Path, default=ROOT_DIR / "reference",
                    help="Directory to write intermediate and output BED files")
    p.add_argument("--upstream", type=int, default=2000,
                    help="Promoter window upstream of TSS (bp)")
    p.add_argument("--downstream", type=int, default=500,
                    help="Promoter window downstream of TSS (bp)")
    p.add_argument("--bedtools", default="bedtools",
                    help="Path to the bedtools executable")
    return p.parse_args()


def parse_attr(attr, key):
    for x in attr.split(";"):
        x = x.strip()
        if x.startswith(key):
            return x.split('"')[1]
    return None


def build_promoter_bed(gtf_path: Path, upstream: int, downstream: int) -> pd.DataFrame:
    df = pd.read_csv(gtf_path, sep="\t", comment="#", names=GTF_COLUMNS)

    # keep transcripts only
    df = df[df["feature"] == "transcript"]

    df["gene_id"] = df["attribute"].apply(lambda x: parse_attr(x, "gene_id"))
    df["gene_name"] = df["attribute"].apply(lambda x: parse_attr(x, "gene_name"))
    df["transcript_id"] = df["attribute"].apply(lambda x: parse_attr(x, "transcript_id"))

    # compute TSS (strand-aware)
    df["tss"] = df.apply(
        lambda r: r["start"] if r["strand"] == "+" else r["end"],
        axis=1
    )

    # promoter window: -upstream to +downstream around TSS
    def promoter(row):
        if row["strand"] == "+":
            start = y if (y := row["tss"] - upstream) > 0 else 0
            end = row["tss"] + downstream
        else:
            start = y if (y := row["tss"] - downstream) > 0 else 0
            end = row["tss"] + upstream
        return pd.Series([start, end])

    df[["p_start", "p_end"]] = df.apply(promoter, axis=1)

    # clean BED output (0-based)
    return df[[
        "chrom",
        "p_start",
        "p_end",
        "gene_name",
        "transcript_id",
        "strand",
        "gene_id",
    ]]


def main():
    args = parse_args()

    if shutil.which(args.bedtools) is None:
        sys.exit(f"ERROR: '{args.bedtools}' not found on PATH — install bedtools first.")

    if not args.gtf.exists():
        sys.exit(
            f"ERROR: GTF not found: {args.gtf}\n"
            f"       Download it first, e.g.:\n"
            f"       wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz\n"
            f"       gunzip gencode.v49.annotation.gtf.gz"
        )
    if not args.cpg_islands.exists():
        sys.exit(f"ERROR: CpG islands BED not found: {args.cpg_islands}")

    args.outdir.mkdir(parents=True, exist_ok=True)
    promoters_bed = args.outdir / "GENCODEV49_promoters_2kb_plus500.bed"
    cpgs_with_promoters = args.outdir / "CpGs_with_promoters.bed"

    print(f"  Building promoter windows (-{args.upstream}/+{args.downstream} bp around TSS)...")
    promoters_df = build_promoter_bed(args.gtf, args.upstream, args.downstream)
    promoters_df.to_csv(promoters_bed, sep="\t", header=False, index=False)
    print(f"  {len(promoters_df):,} promoter records -> {promoters_bed}")

    print(f"  Intersecting {args.cpg_islands.name} with promoters via bedtools...")
    cmd = [
        args.bedtools, "intersect",
        "-a", str(args.cpg_islands),
        "-b", str(promoters_bed),
        "-wa", "-wb",
    ]
    with open(cpgs_with_promoters, "w") as out:
        subprocess.run(cmd, stdout=out, check=True)

    print(f"  Done -> {cpgs_with_promoters}")
    print(f"  Use this as --annotation when running main.py")


if __name__ == "__main__":
    main()
