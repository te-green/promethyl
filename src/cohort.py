import numpy as np
import pandas as pd
import pyranges as pr

def bh_adjust(pvals: pd.Series) -> np.ndarray:
    """Benjamini-Hochberg FDR correction."""
    n = len(pvals)
    arr = np.asarray(pvals, dtype=float)
    order = np.argsort(arr)
    padj = np.empty(n)

    cummin = 1.0
    for i in range(n - 1, -1, -1):
        cummin = min(cummin, arr[order[i]] * n / (i + 1))
        padj[order[i]] = cummin

    return np.clip(padj, 0.0, 1.0)


def load_config(path: str) -> pd.DataFrame:
    """Load YAML cohort config into a sample metadata dataframe."""
    import yaml
    with open(path) as f:
        cfg = yaml.safe_load(f)
    return pd.DataFrame(cfg["samples"])


def _comparison_group(sid: str, sample_ids: list[str], sample_meta: dict, match_sex: bool) -> list[str]:
    """
    Other samples eligible to be pooled as the "rest" background for `sid`.

    Always restricted to the same tissue. If `match_sex` is True (X-chromosome
    sites), also restricted to the same sex -- X-linked methylation differs
    systematically between sexes (X-inactivation mosaicism in females vs.
    hemizygous males), so mixing sexes there produces spurious outliers.
    Samples missing from `sample_meta` or missing tissue/sex are excluded
    from every group rather than guessed at.
    """
    meta_sid = sample_meta.get(sid)
    if not meta_sid or not meta_sid.get("tissue") or (match_sex and not meta_sid.get("sex")):
        return []

    group = []
    for s in sample_ids:
        if s == sid:
            continue
        meta_s = sample_meta.get(s)
        if not meta_s or meta_s.get("tissue") != meta_sid["tissue"]:
            continue
        if match_sex and meta_s.get("sex") != meta_sid["sex"]:
            continue
        group.append(s)
    return group


def build_cohort_matrix(sample_dfs: dict) -> pd.DataFrame:
    """
    Merge per-sample island methylation into a single wide matrix.
    One row per island, one set of columns per sample.

    Uses an OUTER merge: an island is kept if ANY sample has coverage
    there, not only islands covered in every sample. Previously this was
    an inner merge, which meant a single low-coverage sample anywhere in
    the cohort would silently drop that island for every sample -- and
    the more samples you add, the more likely that becomes, even though
    nothing changed for the samples that *did* have coverage there.

    Samples without coverage at an island get NaN in their
    methylation_/coverage_/n_mod_/n_canonical_ columns. Downstream
    consumers (detect_outliers, compute_cohort_statistics) already skip
    NaN when building comparison groups. `n_samples_covered` records how
    many samples actually contributed data at each island, so
    low-N/singleton islands can be spotted or filtered downstream.
    """
    coords = ["chrom", "start", "end", "cpg_island"]
    matrix = None
    sample_ids = list(sample_dfs.keys())

    for sample_id, df in sample_dfs.items():
        # Keep only columns that exist for this sample
        keep = coords + [c for c in df.columns if c.endswith(f"_{sample_id}")]
        sub  = df[[c for c in keep if c in df.columns]]
        if matrix is None:
            matrix = sub
        else:
            matrix = matrix.merge(sub, on=coords, how="outer")

    meth_cols = [f"methylation_{sid}" for sid in sample_ids if f"methylation_{sid}" in matrix.columns]
    matrix["n_samples_covered"] = matrix[meth_cols].notna().sum(axis=1)

    return matrix


def detect_outliers(
    matrix: pd.DataFrame,
    sample_ids: list[str],
    sample_meta: dict | None = None,
    min_delta: float = 0.20,
    z_threshold: float = 2.0,
    min_cpg_sites: int = 1,
    x_chrom_prefix: str = "chrX",
) -> pd.DataFrame:
    """
    For each island, compute Z-score of every sample against a background
    distribution. Any sample can be flagged as an outlier.

    sample_meta: optional {id: {"tissue": ..., "sex": ...}}. When given, the
    background for a given sample at a given island is the OTHER samples of
    the same tissue (and, on chrX, the same sex too -- see
    compute_cohort_statistics for why). When sample_meta is None, behaviour
    is unchanged: one shared cohort_mean/cohort_std/cohort_median computed
    across every sample, same as before.

    Because the background is now sample-specific (it excludes the sample
    itself, and depends on tissue/sex), "cohort_mean"/"cohort_std" become
    per-sample columns (group_mean_<id>/group_std_<id>) rather than one
    shared column, whenever sample_meta is supplied. group_n_<id> records
    how many samples were actually available for comparison at that site --
    check it before trusting an outlier call from a small group.

    min_cpg_sites: minimum number of individual CG positions that must have
    survived the per-site coverage filter and been summed into a SAMPLE's
    own island-level call (n_mod_<sid>/coverage_<sid>) for that sample to be
    eligible for an outlier flag at this island. This is independent of read
    depth: a single CG position can be non-representative of the island
    (allele-specific effects, local mapping artifacts) no matter how many
    reads cover it, so it's a separate check from min_delta/z_threshold, not
    a substitute for either. Default of 1 is a no-op (matches prior
    behaviour) so existing callers aren't silently affected -- raise it once
    you've looked at the n_cpg_sites distribution in your own cohort output.
    Rows without an n_cpg_sites_<sid> column (older aggregate_to_islands()
    output, before that fix) are never gated by this -- they pass through
    exactly as before.
    """
    df = matrix.copy()

    meth_cols = {s: f"methylation_{s}" for s in sample_ids}

    if sample_meta is None:
        cohort_mean   = df[list(meth_cols.values())].mean(axis=1)
        cohort_std    = df[list(meth_cols.values())].std(axis=1)
        df["cohort_mean"]   = cohort_mean
        df["cohort_std"]    = cohort_std
        df["cohort_median"] = df[list(meth_cols.values())].median(axis=1)
    else:
        is_x = df["chrom"].astype(str).str.startswith(x_chrom_prefix)

    for sid in sample_ids:
        meth_col = meth_cols[sid]

        if sample_meta is None:
            group_mean = df["cohort_mean"]
            group_std  = df["cohort_std"]
        else:
            group_auto = _comparison_group(sid, sample_ids, sample_meta, match_sex=False)
            group_x    = _comparison_group(sid, sample_ids, sample_meta, match_sex=True)

            group_mean = pd.Series(np.nan, index=df.index)
            group_std  = pd.Series(np.nan, index=df.index)
            group_n    = pd.Series(0, index=df.index)

            if group_auto:
                cols = [meth_cols[s] for s in group_auto]
                group_mean[~is_x] = df.loc[~is_x, cols].mean(axis=1)
                group_std[~is_x]  = df.loc[~is_x, cols].std(axis=1)
                # Per-row count of group members with actual data at this
                # island (not just the group's nominal size) -- now that
                # the matrix is an outer merge, a given row may have fewer
                # covered members than the full comparison group.
                group_n[~is_x]    = df.loc[~is_x, cols].notna().sum(axis=1)
            if group_x:
                cols = [meth_cols[s] for s in group_x]
                group_mean[is_x] = df.loc[is_x, cols].mean(axis=1)
                group_std[is_x]  = df.loc[is_x, cols].std(axis=1)
                group_n[is_x]    = df.loc[is_x, cols].notna().sum(axis=1)

            df[f"group_mean_{sid}"] = group_mean
            df[f"group_n_{sid}"]    = group_n

        df[f"delta_{sid}"]   = df[meth_col] - group_mean
        df[f"zscore_{sid}"]  = (df[meth_col] - group_mean) / group_std.replace(0, np.nan)

        n_sites_col = f"n_cpg_sites_{sid}"
        if n_sites_col in df.columns:
            enough_sites = df[n_sites_col].fillna(0) >= min_cpg_sites
        else:
            enough_sites = True  # column not present (pre-fix aggregate_to_islands output) -- don't gate

        df[f"outlier_{sid}"] = (
            (df[f"zscore_{sid}"].abs() >= z_threshold) &
            (df[f"delta_{sid}"].abs()  >= min_delta) &
            enough_sites
        )

    # Summary columns
    outlier_cols = [f"outlier_{sid}" for sid in sample_ids]
    df["n_outliers"]      = df[outlier_cols].sum(axis=1)
    df["outlier_samples"] = df[outlier_cols].apply(
        lambda row: ";".join(
            sid for sid, is_out in zip(sample_ids, row) if is_out
        ),
        axis=1,
    )
    df["any_outlier"] = df["n_outliers"] > 0

    sort_cols = ["n_outliers"] + (["cohort_std"] if sample_meta is None else [])
    return df.sort_values(
        sort_cols, ascending=[False] * len(sort_cols)
    ).reset_index(drop=True)


def to_long_format(df: pd.DataFrame, sample_ids: list[str]) -> pd.DataFrame:
    """
    Reshape a wide cohort dataframe (one row per island, per-sample columns
    suffixed with "_<sample_id>") into long format: one row per
    (island, sample), with per-sample metrics as plain-named columns.

    Locus/annotation/cohort-level columns (chrom, start, end, cpg_island,
    gene, cohort_mean, ...) are not sample-specific and are repeated across
    each island's sample rows.
    """
    per_sample_prefixes = [
        col[: -(len(sid) + 1)]
        for sid in sample_ids
        for col in df.columns
        if col.endswith(f"_{sid}")
    ]
    per_sample_prefixes = sorted(set(per_sample_prefixes))

    shared_cols = [
        col for col in df.columns
        if not any(col == f"{prefix}_{sid}" for prefix in per_sample_prefixes for sid in sample_ids)
    ]

    frames = []
    for sid in sample_ids:
        sample_cols = {f"{prefix}_{sid}": prefix for prefix in per_sample_prefixes if f"{prefix}_{sid}" in df.columns}
        sub = df[shared_cols + list(sample_cols.keys())].rename(columns=sample_cols)
        sub.insert(len(shared_cols), "sample", sid)
        frames.append(sub)

    long_df = pd.concat(frames, ignore_index=True)

    if "outlier_samples" in long_df.columns:
        long_df = long_df.drop(columns=["outlier_samples"])

    sort_cols = [c for c in ["chrom", "start", "end", "sample"] if c in long_df.columns]
    return long_df.sort_values(sort_cols).reset_index(drop=True)


def compute_cohort_statistics(
    matrix: pd.DataFrame,
    sample_ids: list[str],
    sample_meta: dict | None = None,
    min_delta: float = 0.20,
    z_threshold: float = 2.0,
    x_chrom_prefix: str = "chrX",
) -> pd.DataFrame:
    """
    For every sample at every island, run a one-vs-rest chi2 test.
    Every sample is treated as a potential outlier.

    sample_meta: optional {id: {"tissue": ..., "sex": ...}}. When given, the
    "rest" pooled for each sample is restricted to other samples of the same
    tissue; on chrX (any chrom starting with `x_chrom_prefix`), it is further
    restricted to the same sex. Sites/samples with no eligible comparison
    group (e.g. a singleton tissue, or the only male in a WB group on chrX)
    get pval/padj = NaN rather than a misleading 1.0 -- they were not tested,
    not "not significant". When sample_meta is None, behaviour is unchanged
    from before (pool every other sample, regardless of tissue/sex).
    """
    from scipy.stats import chi2_contingency

    df = matrix.copy()
    is_x = df["chrom"].astype(str).str.startswith(x_chrom_prefix) if sample_meta is not None else None

    def chi2_p(row, n_mod_a, n_can_a, n_mod_b, n_can_b):
        if pd.isna(row[n_mod_a]) or pd.isna(row[n_can_a]) or pd.isna(row[n_mod_b]) or pd.isna(row[n_can_b]):
            return np.nan
        table = [
            [row[n_mod_a], row[n_can_a]],
            [row[n_mod_b], row[n_can_b]],
        ]
        if any(sum(r) == 0 for r in table) or any(sum(c) == 0 for c in zip(*table)):
            return 1.0
        _, p, _, _ = chi2_contingency(table, correction=False)
        return p

    for sid in sample_ids:
        n_mod_rest = pd.Series(np.nan, index=df.index)
        n_can_rest = pd.Series(np.nan, index=df.index)

        if sample_meta is None:
            rest = [s for s in sample_ids if s != sid]
            # min_count=1: if every "rest" sample is NaN at this row (outer
            # merge), sum() would otherwise silently return 0 rather than
            # NaN, which chi2_p can't distinguish from a genuine zero count.
            n_mod_rest[:] = df[[f"n_mod_{s}"       for s in rest]].sum(axis=1, min_count=1)
            n_can_rest[:] = df[[f"n_canonical_{s}" for s in rest]].sum(axis=1, min_count=1)
        else:
            group_auto = _comparison_group(sid, sample_ids, sample_meta, match_sex=False)
            group_x    = _comparison_group(sid, sample_ids, sample_meta, match_sex=True)

            if group_auto:
                cols_mod, cols_can = [f"n_mod_{s}" for s in group_auto], [f"n_canonical_{s}" for s in group_auto]
                n_mod_rest[~is_x] = df.loc[~is_x, cols_mod].sum(axis=1, min_count=1)
                n_can_rest[~is_x] = df.loc[~is_x, cols_can].sum(axis=1, min_count=1)
            if group_x:
                cols_mod, cols_can = [f"n_mod_{s}" for s in group_x], [f"n_canonical_{s}" for s in group_x]
                n_mod_rest[is_x] = df.loc[is_x, cols_mod].sum(axis=1, min_count=1)
                n_can_rest[is_x] = df.loc[is_x, cols_can].sum(axis=1, min_count=1)

        df["_n_mod_rest"] = n_mod_rest
        df["_n_can_rest"] = n_can_rest

        df[f"pval_{sid}"] = df.apply(
            chi2_p, axis=1,
            n_mod_a=f"n_mod_{sid}",
            n_can_a=f"n_canonical_{sid}",
            n_mod_b="_n_mod_rest",
            n_can_b="_n_can_rest",
        )

        pvals = df[f"pval_{sid}"]
        tested = pvals.notna()
        padj = pd.Series(np.nan, index=df.index)
        if tested.any():
            padj[tested] = bh_adjust(pvals[tested])
        df[f"padj_{sid}"] = padj

        df = df.drop(columns=["_n_mod_rest", "_n_can_rest"])

    return df