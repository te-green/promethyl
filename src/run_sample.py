#!/usr/bin/env python3
"""
run_sample.py - Phase 1 of the cohort pipeline: per-sample modkit + island aggregation.

Runs modkit pileup on a single BAM, filters the bedMethyl output, and
aggregates it to CpG islands. Designed to be run once per sample so it can
be parallelized across samples (e.g. one Nextflow process per sample).

Writes a single TSV per sample containing island-level methylation counts
for that sample only. Phase 2 (run_cohort.py) collects these TSVs across
all samples and performs the cohort outlier analysis.
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

from CpG_meth import load_promoter_annotations, aggregate_to_islands
from annotate import load_dmr_annotations
from run_modkit import get_sample_methylation


def parse_args():
    p = argparse.ArgumentParser(
        description="Run modkit + island aggregation for a single sample.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--id", required=True, help="Sample ID")
    p.add_argument("--bam", required=True, help="Sample BAM file")
    p.add_argument("--modkit-dir", required=True, help="Directory for modkit bedMethyl outputs")
    p.add_argument("--annotation", required=True, help="CpGs_with_promoters.bed")
    p.add_argument("--output", required=True, help="Output TSV for this sample's island methylation")
    p.add_argument("--dmr-bed", default=None,
                    help="Known-DMR reference BED (chrom,start,end,dmr_name,disorder). When given, "
                         "each DMR is aggregated as its own region (same as a CpG island, but using "
                         "the DMR's own exact boundaries) in addition to normal island aggregation, "
                         "so a DMR isn't tested only via whatever generic island happens to overlap "
                         "it. Needs --include-bed to actually cover the DMRs' full extent, or modkit "
                         "won't have pileup data there in the first place.")

    p.add_argument("--ref", default=None, help="Reference FASTA for modkit pileup")
    p.add_argument("--threads", type=int, default=10, help="Threads for modkit pileup")
    p.add_argument("--region", default=None, help="Restrict modkit pileup to a genomic region")
    p.add_argument("--include-bed", default=None, help="BED file to restrict modkit pileup to specific regions")

    p.add_argument("--min-coverage", type=int, default=10, help="Minimum read coverage per site")
    p.add_argument("--mod-code", default="m", help="Modification code ('m'=5mC, 'h'=5hmC)")

    return p.parse_args()


def main():
    args = parse_args()

    modkit_dir = Path(args.modkit_dir)
    annotation = Path(args.annotation)
    output = Path(args.output)

    if not annotation.exists():
        sys.exit(
            f"ERROR: annotation file not found: {annotation}\n"
            f"       Run annotate_promoter_cpgs.py first to generate it."
        )
    bam = Path(args.bam)
    if not bam.exists():
        sys.exit(f"ERROR: BAM not found: {bam}")

    modkit_dir.mkdir(parents=True, exist_ok=True)
    output.parent.mkdir(parents=True, exist_ok=True)

    ann_df = load_promoter_annotations(str(annotation))
    print(f"  {len(ann_df):,} CpG island-promoter records loaded")

    print(f"=== [{args.id}] Load methylation data ===")
    raw_df = get_sample_methylation(
        label=args.id,
        bam=bam,
        modkit_dir=modkit_dir,
        ref=args.ref,
        threads=args.threads,
        min_coverage=args.min_coverage,
        mod_code=args.mod_code,
        region=args.region,
        include_bed=args.include_bed,
    )

    if raw_df is None:
        sys.exit(f"ERROR: [{args.id}] no CpG data produced — cannot continue")

    print(f"=== [{args.id}] Aggregate to CpG islands ===")
    df = aggregate_to_islands(raw_df, ann_df)
    df["is_dmr_region"] = False
    print(f"  {len(df):,} CpG islands with coverage for {args.id}")

    if args.dmr_bed:
        print(f"=== [{args.id}] Aggregate to known DMR regions ===")
        dmr_df = load_dmr_annotations(args.dmr_bed)
        # Aggregated on the DMR's own exact boundaries, not whatever CpG island(s)
        # happen to overlap it -- reuses aggregate_to_islands() as-is by just
        # renaming dmr_name to the column it expects (cpg_island).
        dmr_as_islands = dmr_df[["chrom", "start", "end", "dmr_name"]].rename(columns={"dmr_name": "cpg_island"})
        dmr_agg = aggregate_to_islands(raw_df, dmr_as_islands, allow_empty=True)
        dmr_agg["is_dmr_region"] = True
        print(f"  {len(dmr_agg):,} / {len(dmr_df):,} known DMRs with coverage for {args.id}")
        df = pd.concat([df, dmr_agg], ignore_index=True)

    df.to_csv(output, sep="\t", index=False)
    print(f"  [{args.id}] written -> {output}")


if __name__ == "__main__":
    main()
