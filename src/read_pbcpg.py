from pathlib import Path
import pandas as pd

METHBAT_COLUMNS = [
    "#chrom",
    "start",
    "end",
    "name",
    "score",
    "strand",
    "Extended",
    "mod_score",
    "type",
    "cov",
    "mod_count",
    "unmod_count",
    "inferred_unmod_count",
    "diff_base_count",
    "avg_mod_score",
    "avg_unmod_score"
]

def read_pbcpg(filepath):
    """
    Read methbat output into a pandas dataframe.
    """

    df = pd.read_csv(
        filepath,
        sep="\t",
        skiprows=lambda i, line: line.startswith("##"),
        names=BEDMETHYL_COLUMNS,
    )

    df = df.rename(
        columns={
            "#chrom": "chrom"
            "cov": "coverage",
        }
    )

    return df[
        [
            "chrom",
            "start",
            "end",
            "name",
            "score",
            "strand",
            "Extended",
            "mod_score",
            "type",
            "coverage",
            "mod_count",
            "unmod_count",
            "inferred_unmod_count",
            "diff_base_count",
            "avg_mod_score",
            "avg_unmod_score",
        ]
    ]

def filter_sites(df, min_coverage=10):
    return df[df["coverage"] >= min_coverage].copy()

def filter_modification(df,mod_code="m"):
    return df[df["type"] == "Total"].copy()