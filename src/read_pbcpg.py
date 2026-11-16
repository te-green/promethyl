
import gzip
import pandas as pd


METHBAT_COLUMNS = [
    "#chrom",
    "start",
    "end",
    "name",
    "score",
    "strand",
    "mod_score",
    "type",
    "cov",
    "mod_count",
    "unmod_count",
    "inferred_unmod_count",
    "diff_base_count",
    "avg_mod_score",
    "avg_unmod_score",
]


def read_pbcpgtools(filepath):
    """
    Read MethBat pileup BED into a pandas DataFrame.

    Skips the ## metadata preamble, preserves the column header,
    and standardises column names for downstream Promethyl use.
    """

    # Find the first non-comment line (the BED column header).
    with gzip.open(filepath, "rt") as f:
        header_line = next(
            i for i, line in enumerate(f)
            if not line.startswith("##")
        )

    df = pd.read_csv(
        filepath,
        sep="\t",
        skiprows=header_line,
        header=0,
    )

    # Validate the MethBat output format.
    if len(df.columns) != len(METHBAT_COLUMNS):
        raise ValueError(
            f"Unexpected MethBat BED format: "
            f"expected {len(METHBAT_COLUMNS)} columns, "
            f"found {len(df.columns)}"
        )

    # Standardise column names.
    df = df.rename(
        columns={
            "#chrom": "chrom",
            "cov": "coverage",
        }
    )

    return df


def filter_sites(df, min_coverage=10):
    """Filter CpG sites by minimum coverage."""

    return df[df["coverage"] >= min_coverage].copy()


def filter_modification(df, mod_code="m"):
    """Retain total (unphased) methylation measurements."""

    return df[df["type"] == "Total"].copy()
