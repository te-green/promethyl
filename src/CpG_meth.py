import pandas as pd
import pyranges as pr
import numpy as np

SAMPLES = ["proband", "father", "mother"]


def load_promoter_annotations(path: str) -> pd.DataFrame:
    cols = [
        "chrom", "start", "end", "cpg_island",
        "prom_chrom", "prom_start", "prom_end",
        "gene", "transcript", "strand", "gene_id",
    ]
    df = pd.read_csv(path, sep="\t", header=None, names=cols)
    # bedtools -loj fills unmatched B-side fields with "." (--all-islands mode,
    # islands with no promoter overlap) — treat as missing, not a literal gene/transcript name.
    for col in ["prom_chrom", "gene", "transcript", "strand", "gene_id"]:
        df[col] = df[col].replace(".", pd.NA)
    return df


def annotate_cpg_sites(cpg_df: pd.DataFrame, ann_df: pd.DataFrame) -> pd.DataFrame:
    """Annotate individual CpG sites with island/gene info."""
    cpg_pr = pr.PyRanges(
        cpg_df.rename(columns={"chrom": "Chromosome", "start": "Start", "end": "End"})
    )
    ann_pr = pr.PyRanges(
        ann_df.rename(columns={"chrom": "Chromosome", "start": "Start", "end": "End"})
    )

    hits = cpg_pr.join(ann_pr, strandedness=False)

    agg = (
        hits.df
        .groupby(["Chromosome", "Start", "End"], as_index=False)
        .agg({
            "gene":       lambda x: ";".join(sorted(set(x.dropna()))),
            "transcript": lambda x: ";".join(sorted(set(x.dropna()))),
            "gene_id":    lambda x: ";".join(sorted(set(x.dropna()))),
            "cpg_island": lambda x: ";".join(sorted(set(x.dropna()))),
        })
    )

    result = cpg_df.merge(
        agg,
        left_on=["chrom", "start", "end"],
        right_on=["Chromosome", "Start", "End"],
        how="left",
    ).drop(columns=["Chromosome", "Start", "End"], errors="ignore")

    for col in ["gene", "transcript", "gene_id", "cpg_island"]:
        result[col] = result[col].fillna("")

    return result


def aggregate_to_islands(merged_df: pd.DataFrame, ann_df: pd.DataFrame) -> pd.DataFrame:
    # Infer sample IDs from column names dynamically
    samples = [c.replace("n_mod_", "") for c in merged_df.columns if c.startswith("n_mod_")]

    islands = (
        ann_df[["cpg_island", "chrom", "start", "end"]]
        .drop_duplicates()
        .dropna(subset=["cpg_island"])
    )

    sites_pr  = pr.PyRanges(
        merged_df.rename(columns={"chrom": "Chromosome", "start": "Start", "end": "End"})
    )
    island_pr = pr.PyRanges(
        islands.rename(columns={"chrom": "Chromosome", "start": "Start", "end": "End"})
    )

    hits = island_pr.join(sites_pr, strandedness=False).df

    if hits.empty:
        raise ValueError("No CpG sites overlapped any CpG island — check annotation file.")

    agg = (
        hits.groupby(["Chromosome", "Start", "End", "cpg_island"], as_index=False)
        .agg(
            # Per-sample, not a single shared column: aggregate_to_islands() is
            # called once per sample in the run_sample.py/cohort path, and even
            # in trio mode different samples can in principle contribute
            # different numbers of surviving CG sites. A per-sample name also
            # lets this survive build_cohort_matrix()'s `endswith(f"_{sample_id}")`
            # column filter, which previously dropped the old shared
            # `n_cpg_sites` column silently before it ever reached the report.
            **{f"n_cpg_sites_{s}": (f"n_mod_{s}",       "count") for s in samples},
            **{f"n_mod_{s}":       (f"n_mod_{s}",       "sum") for s in samples},
            **{f"n_canonical_{s}": (f"n_canonical_{s}", "sum") for s in samples},
            **{f"coverage_{s}":    (f"coverage_{s}",    "mean") for s in samples},
        )
    )

    for s in samples:
        total = agg[f"n_mod_{s}"] + agg[f"n_canonical_{s}"]
        agg[f"methylation_{s}"] = np.where(total > 0, agg[f"n_mod_{s}"] / total, np.nan)

    # Drop islands with zero total counts in any sample
    for s in samples:
        total = agg[f"n_mod_{s}"] + agg[f"n_canonical_{s}"]
        agg = agg[total > 0]

    return (
        agg
        .rename(columns={"Chromosome": "chrom", "Start": "start", "End": "end"})
        .dropna(subset=[f"methylation_{s}" for s in samples])
        .reset_index(drop=True)
    )