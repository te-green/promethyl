#!/usr/bin/env python3
"""
prometh - CpG methylation analysis pipeline

Run annotate_promoter_cpgs.py first to build the reference annotation file:
    python annotate_promoter_cpgs.py
    
Option 1:    
To run the trio pipeline:
    python main.py \
        --proband proband.bam \
        --father  father.bam \
        --mother  mother.bam \
        --modkit-dir  /data/modkit \
        --annotation  Reference/CpGs_with_promoters.bed \
        --ref         /data/genome.fa \
        --output      trio_methylation.tsv
        --threads 10

Option 2:
To run the cohort outlier pipeline:
    python main.py \
    --config sample.yaml \
    --methbat_dir /data/methbat/ \
    --annotation /Reference/CpGs_with_promoters.bed  \
    --ref /data/genome.fa \
    --output /output-dir/ProbandMethbat_cohort_methylation.tsv \
    --include-bed /Reference/CPGIslandsBED3.bed \
    --threads 10


methbat pileup is skipped for any sample whose output already exists in --modkit-dir.
"""

import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyranges as pr
from scipy.stats import fisher_exact, chi2_contingency

from read_pbcpg import filter_modification, filter_sites, read_pbcpgtools
from CpG_meth import load_promoter_annotations, aggregate_to_islands
import yaml
from cohort import load_config, build_cohort_matrix, detect_outliers, compute_cohort_statistics, to_long_format
from run_pbcpg import run_methbat, get_sample_methylation
from merge_trio import merge_trio
from annotate import annotate_methylation
from statistics import compute_statistics

SAMPLES = ["proband", "father", "mother"]

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Trio CpG methylation analysis (run annotate_promoter_cpgs.py first).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Trio BAMs
    p.add_argument("--proband", default=None, help="Proband BAM file")
    p.add_argument("--father",  default=None, help="Father BAM file")
    p.add_argument("--mother",  default=None, help="Mother BAM file")

    # modkit / annotation
    p.add_argument("--methbat-dir", required=True,
                   help="Directory for methbat outputs (generated if absent)")
    p.add_argument("--annotation", required=True,
                   help="CpGs_with_promoters.bed produced by annotate_promoter_cpgs.py")

    # Output
    p.add_argument("--output", required=True, help="Output TSV file")

    # modkit options
    p.add_argument("--ref",     default=None, help="Reference FASTA for modkit pileup (recommended)")
    p.add_argument("--threads", type=int, default=10, help="Threads for modkit pileup")

    # Filtering
    p.add_argument("--min-coverage", type=int, default=10, help="Minimum read coverage per site")
    p.add_argument("--min-delta", type=float, default=0.3, help="Minimum absolute methylation difference (proband vs parent mean)")
    p.add_argument("--fdr", type=float, default=0.01, help="FDR threshold for significance")

    p.add_argument("--config", default=None, help="YAML config for cohort mode (alternative to --proband/--father/--mother)")
    p.add_argument("--proband-id", default=None, help="Sample ID(s) to treat as proband in cohort mode (comma-separated)")
    p.add_argument("--z-threshold", type=float, default=2.0,help="Z-score threshold for cohort outlier calling")

    return p.parse_args()


def main():
    args = parse_args()

    methbat_out = Path(args.methbat_out)
    annotation = Path(args.annotation)
    output     = Path(args.output)

    if args.config and (args.proband or args.father or args.mother):
        sys.exit("ERROR: --config and --proband/--father/--mother are mutually exclusive")

    if not args.config and not all([args.proband, args.father, args.mother]):
        sys.exit("ERROR: either --config or all of --proband/--father/--mother are required")

    # Validate annotation (shared by both modes)
    if not annotation.exists():
        sys.exit(
            f"ERROR: annotation file not found: {annotation}\n"
            f"       Run annotate_promoter_cpgs.py first to generate it."
        )

    methbat_out.mkdir(parents=True, exist_ok=True)
    output.parent.mkdir(parents=True, exist_ok=True)

    # Load annotation (shared by both modes)
    ann_df = load_promoter_annotations(str(annotation))
    print(f"  {len(ann_df):,} CpG island–promoter records loaded")

    if args.config:
        # ── Cohort mode ───────────────────────────────────────────────
        meta    = load_config(args.config)
        all_ids = meta["id"].tolist()
        print(f"=== Cohort mode: {len(meta)} samples ===")

        sample_dfs = {}
        skipped    = []
        for _, row in meta.iterrows():
            sid = row["id"]
            print(f"\n  [{sid}]")
            df = get_sample_methylation(
                label=sid,
                bam=Path(row["bam"]),
                methbat_out=methbat_out,
                threads=args.threads,
                min_coverage=args.min_coverage,
            )

            if df is None:
                skipped.append(sid)
                continue

            df = aggregate_to_islands(df, ann_df)
            sample_dfs[sid] = df

        if skipped:
            print(f"\n  WARNING: {len(skipped)} sample(s) skipped (no CpG data): {', '.join(skipped)}")

        if len(sample_dfs) < 2:
            sys.exit("ERROR: fewer than 2 samples with data — cannot run cohort analysis")

        # Use only samples that produced data
        all_ids = list(sample_dfs.keys())

        sample_meta = None
        if "tissue" in meta.columns or "sex" in meta.columns:
            sample_meta = {
                row["id"]: {"tissue": row.get("tissue"), "sex": row.get("sex")}
                for _, row in meta.iterrows()
            }
            print(f"  Cohort statistics restricted to tissue-matched (and, on chrX, sex-matched) samples")
        else:
            print(f"  WARNING: no tissue/sex columns in {args.config} — comparing every sample "
                  f"against every other sample regardless of tissue/sex")

        matrix = build_cohort_matrix(sample_dfs)
        print(f"\n  Islands covered across all samples: {len(matrix):,}")

        result = compute_cohort_statistics(matrix, all_ids, sample_meta=sample_meta,
                                           min_delta=args.min_delta,
                                           z_threshold=args.z_threshold)
        result = detect_outliers(result, all_ids, sample_meta=sample_meta,
                                 min_delta=args.min_delta,
                                 z_threshold=args.z_threshold)
        result = annotate_methylation(result, ann_df)
        result = to_long_format(result, all_ids)

    else:
        # ── Trio mode ─────────────────────────────────────────────────
        for label, path in [("proband", args.proband), ("father", args.father), ("mother", args.mother)]:
            if not Path(path).exists():
                sys.exit(f"ERROR: {label} BAM not found: {path}")

        print("=== Step 1: Load methylation data ===")
        common = dict(
            methbat_out=methbat_out,
            ref=args.ref,
            threads=args.threads,
            min_coverage=args.min_coverage,,
        )
        proband_df = get_sample_methylation("proband", Path(args.proband), **common)
        father_df  = get_sample_methylation("father",  Path(args.father),  **common)
        mother_df  = get_sample_methylation("mother",  Path(args.mother),  **common)

        print("\n=== Step 2: Merge trio ===")
        merged = merge_trio(proband_df, father_df, mother_df)

        print("\n=== Step 3: Annotate & aggregate to CpG islands ===")
        merged    = aggregate_to_islands(merged, ann_df)
        print(f"  {len(merged):,} CpG islands with coverage across all three samples")
        annotated = annotate_methylation(merged, ann_df)

        print("\n=== Step 4: Statistical tests ===")
        result = compute_statistics(annotated,
                                    min_delta=args.min_delta,
                                    fdr_threshold=args.fdr)

    # ── Write output (shared by both modes) ───────────────────────────
    result.to_csv(output, sep="\t", index=False)
    print(f"\n✓ Done!")
    print(f"  Sites tested : {len(result):,}")
    print(f"  Output       : {output}")


if __name__ == "__main__":
    main()
