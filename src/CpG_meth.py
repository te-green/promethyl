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


def aggregate_to_islands(merged_df: pd.DataFrame, ann_df: pd.DataFrame, allow_empty: bool = False) -> pd.DataFrame:
    # Infer sample IDs from column names dynamically
    samples = [c.replace("mod_count_", "") for c in merged_df.columns if c.startswith("mod_count_")]
    # Samples that also carry per-site haplotype counts (see run_pbcpg.extract_haplotype_counts)
    hap_cols = lambda s: [f"hap{h}_{k}_count_{s}" for h in (1, 2) for k in ("mod", "unmod")]
    hap_samples = [s for s in samples if all(c in merged_df.columns for c in hap_cols(s))]

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
        if allow_empty:
            cols = ["chrom", "start", "end", "cpg_island"]
            for s in samples:
                cols += [f"n_cpg_sites_{s}", f"mod_count_{s}", f"unmod_count_{s}", f"coverage_{s}", f"mod_score_{s}"]
            for s in hap_samples:
                cols += hap_cols(s) + [f"n_hap_sites_{s}"]
            return pd.DataFrame(columns=cols)
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
            **{f"n_cpg_sites_{s}": (f"mod_count_{s}",       "count") for s in samples},
            **{f"mod_count_{s}":       (f"mod_count_{s}",       "sum") for s in samples},
            **{f"unmod_count_{s}": (f"unmod_count_{s}", "sum") for s in samples},
            **{f"coverage_{s}":    (f"coverage_{s}",    "mean") for s in samples},
            # Haplotype counts: pandas sum() skips NaN, so an island with no
            # paired hap sites would get 0 -- n_hap_sites lets us turn that back into NaN.
            **{f"n_hap_sites_{s}": (f"hap1_mod_count_{s}", "count") for s in hap_samples},
            **{c: (c, "sum") for s in hap_samples for c in hap_cols(s)},
        )
    )

    for s in hap_samples:
        no_hap = agg[f"n_hap_sites_{s}"] == 0
        for c in hap_cols(s):
            agg.loc[no_hap, c] = np.nan

    for s in samples:
        total = agg[f"mod_count_{s}"] + agg[f"unmod_count_{s}"]
        agg[f"mod_score_{s}"] = np.where(total > 0, agg[f"mod_count_{s}"] / total, np.nan)

    # Drop islands with zero total counts in any sample
    for s in samples:
        total = agg[f"mod_count_{s}"] + agg[f"unmod_count_{s}"]
        agg = agg[total > 0]

    return (
        agg
        .rename(columns={"Chromosome": "chrom", "Start": "start", "End": "end"})
        .dropna(subset=[f"mod_score_{s}" for s in samples])
        .reset_index(drop=True)
    )