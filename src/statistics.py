
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyranges as pr
from scipy.stats import fisher_exact, chi2_contingency

# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------

def bh_adjust(pvals: pd.Series) -> np.ndarray:
    """Benjamini-Hochberg FDR correction."""
    n     = len(pvals)
    arr   = np.asarray(pvals, dtype=float)
    order = np.argsort(arr)
    padj  = np.empty(n)

    cummin = 1.0
    for i in range(n - 1, -1, -1):
        cummin        = min(cummin, arr[order[i]] * n / (i + 1))
        padj[order[i]] = cummin

    return np.clip(padj, 0.0, 1.0)


def _fisher_p(row: pd.Series, mod_count_b: str, n_can_b: str) -> float:
    """Two-sided Fisher's exact test: proband vs one comparator."""
    _, p = fisher_exact(
        [[row["mod_count_proband"], row["unmod_count_proband"]],
         [row[mod_count_b],         row[n_can_b]]],
        alternative="two-sided",
    )
    return p


def compute_statistics(df: pd.DataFrame, min_delta: float = 0.3, fdr_threshold: float = 0.01) -> pd.DataFrame:
    df = df.copy()

    # Deltas
    df["delta_father"]      = df["mod_score_proband"] - df["mod_score_father"]
    df["delta_mother"]      = df["mod_score_proband"] - df["mod_score_mother"]
    df["delta_parent_mean"] = df["mod_score_proband"] - (
        (df["mod_score_father"] + df["mod_score_mother"]) / 2
    )

    # Use chi-squared instead of Fisher's — less sensitive to large N
    # but still add effect size filter
    def chi2_p(row, mod_count_b, n_can_b):
        table = [
            [row["mod_count_proband"], row["unmod_count_proband"]],
            [row[mod_count_b],         row[n_can_b]],
        ]
        # Return 1.0 if any row or column sums to zero
        if any(sum(r) == 0 for r in table) or any(sum(c) == 0 for c in zip(*table)):
            return 1.0
        _, p, _, _ = chi2_contingency(table, correction=False)
        return p

    df["pval_father"]      = df.apply(chi2_p, axis=1, mod_count_b="mod_count_father",  n_can_b="unmod_count_father")
    df["pval_mother"]      = df.apply(chi2_p, axis=1, mod_count_b="mod_count_mother",  n_can_b="unmod_count_mother")

    df["_mod_count_parents"] = df["mod_count_father"]    + df["mod_count_mother"]
    df["_n_can_parents"] = df["unmod_count_father"] + df["unmod_count_mother"]
    df["pval_parent_mean"] = df.apply(chi2_p, axis=1, mod_count_b="_mod_count_parents", n_can_b="_n_can_parents")
    df = df.drop(columns=["_mod_count_parents", "_n_can_parents"])

    # BH FDR
    df["padj_father"]      = bh_adjust(df["pval_father"])
    df["padj_mother"]      = bh_adjust(df["pval_mother"])
    df["padj_parent_mean"] = bh_adjust(df["pval_parent_mean"])

    # Combined filter: FDR + minimum effect size
    df["significant"] = (
        (df["padj_parent_mean"] < fdr_threshold) &
        (df["delta_parent_mean"].abs() >= min_delta)
    )

    df = df.sort_values(
        ["padj_parent_mean", "pval_parent_mean"],
        ascending=True,
    ).reset_index(drop=True)

    return df