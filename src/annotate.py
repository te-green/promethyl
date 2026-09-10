import pandas as pd
import pyranges as pr

# ---------------------------------------------------------------------------
# Annotation
# ---------------------------------------------------------------------------

def load_dmr_annotations(path: str) -> pd.DataFrame:
    """Load the known-DMR reference BED (chrom, start, end, dmr_name, disorder).

    `disorder` is blank for DMRs with no associated disorder in Table 1 --
    treated as missing, not a literal empty string, for consistent
    semicolon-joining in annotate_dmrs().
    """
    cols = ["chrom", "start", "end", "dmr_name", "disorder"]
    df = pd.read_csv(path, sep="\t", header=None, names=cols)
    df["disorder"] = df["disorder"].replace("", pd.NA)
    return df


def annotate_dmrs(merged_df: pd.DataFrame, dmr_df: pd.DataFrame) -> pd.DataFrame:
    """Tag each island with any known DMR(s) it overlaps (dmr_name, disorder).

    Same coordinate-overlap-join shape as annotate_methylation() -- a CpG
    island can in principle overlap more than one Table-1 DMR (adjacent/
    overlapping entries, e.g. the H19/IGF2 and downstream IGF2 DMRs), so
    hits are semicolon-joined rather than assumed to be 1:1.
    """
    sites_pr = pr.PyRanges(
        merged_df.rename(columns={"chrom": "Chromosome", "start": "Start", "end": "End"})
    )
    dmr_pr = pr.PyRanges(
        dmr_df.rename(columns={"chrom": "Chromosome", "start": "Start", "end": "End"})
    )

    hits = sites_pr.join(dmr_pr, strandedness=False)

    if hits.df.empty:
        result = merged_df.copy()
        result["dmr_name"] = ""
        result["disorder"] = ""
        return result

    def _paired_join(group: pd.DataFrame) -> pd.Series:
        # Sort/dedupe as (dmr_name, disorder) PAIRS, not two independently
        # joined columns -- keeps disorder correctly attached to its own
        # DMR when an island overlaps more than one Table-1 entry (e.g.
        # adjacent IGF2 DMRs). Empty disorders are kept as empty strings in
        # position rather than dropped, so dmr_name.split(';')[i] always
        # corresponds to disorder.split(';')[i] for any downstream consumer
        # -- dropping them would silently break that alignment as soon as
        # one DMR in a multi-DMR island has a disorder and another doesn't.
        pairs = sorted(set(zip(group["dmr_name"], group["disorder"].fillna(""))))
        return pd.Series({
            "dmr_name": ";".join(p[0] for p in pairs),
            "disorder": ";".join(p[1] for p in pairs),
        })

    agg = (
        hits.df
        .groupby(["Chromosome", "Start", "End"])
        .apply(_paired_join, include_groups=False)
        .reset_index()
    )

    result = (
        merged_df
        .merge(
            agg,
            left_on=["chrom", "start", "end"],
            right_on=["Chromosome", "Start", "End"],
            how="left",
        )
        .drop(columns=["Chromosome", "Start", "End"], errors="ignore")
    )

    for col in ["dmr_name", "disorder"]:
        result[col] = result[col].fillna("")

    return result


def annotate_methylation(merged_df: pd.DataFrame, ann_df: pd.DataFrame) -> pd.DataFrame:
    """Annotate island-level trio df with gene/promoter information."""
    ann_clean = (
        ann_df
        .drop(columns=["prom_chrom", "prom_start", "prom_end"], errors="ignore")
        .rename(columns={"strand": "gene_strand"})
    )

    sites_pr = pr.PyRanges(
        merged_df.rename(columns={"chrom": "Chromosome", "start": "Start", "end": "End"})
    )
    ann_pr = pr.PyRanges(
        ann_clean.rename(columns={"chrom": "Chromosome", "start": "Start", "end": "End"})
    )

    hits = sites_pr.join(ann_pr, strandedness=False)

    agg = (
        hits.df
        .groupby(["Chromosome", "Start", "End"], as_index=False)
        .agg({
            "gene":       lambda x: ";".join(sorted(set(x.dropna()))),
            "transcript": lambda x: ";".join(sorted(set(x.dropna()))),
            "gene_id":    lambda x: ";".join(sorted(set(x.dropna()))),
            # cpg_island removed — already present from aggregate_to_islands()
        })
    )

    result = (
        merged_df
        .merge(
            agg,
            left_on=["chrom", "start", "end"],
            right_on=["Chromosome", "Start", "End"],
            how="left",
        )
        .drop(columns=["Chromosome", "Start", "End"], errors="ignore")
    )

    for col in ["gene", "transcript", "gene_id"]:
        result[col] = result[col].fillna("")

    return result