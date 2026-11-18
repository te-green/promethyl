
import subprocess
from pathlib import Path

import pandas as pd

from read_pbcpg import read_pbcpgtools


# ---------------------------------------------------------------------------
# MethBat
# ---------------------------------------------------------------------------

def run_methbat(
    bam: Path,
    output_prefix: Path,
    threads: int = 10,
    min_coverage: int = 10,
    min_hap_coverage: int = 5,
) -> None:
    """Run MethBat pileup on a single PacBio BAM.

    --min-coverage is set to the LOWER of the Total and per-haplotype
    thresholds so haplotype rows (which have roughly half the depth of Total)
    are not discarded by MethBat before we see them. The stricter Total
    threshold is re-applied in load_sample_methylation().
    """

    output_prefix.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        "methbat", "pileup",
        "--input-bam", str(bam),
        "--output-prefix", str(output_prefix),
        "--threads", str(threads),
        "--min-coverage", str(min(min_coverage, min_hap_coverage)),
        "--skip-5hmC",
        "--skip-6mA",
    ]

    print(f"    $ {' '.join(cmd)}")

    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
    )

    if result.returncode != 0:
        raise RuntimeError(
            f"MethBat failed on {bam.name}:\n{result.stderr}"
        )


# ---------------------------------------------------------------------------
# Per-sample loading
# ---------------------------------------------------------------------------

def ensure_methbat_output(
    label: str,
    bam: Path,
    methbat_dir: Path,
    threads: int = 10,
    min_coverage: int = 10,
    min_hap_coverage: int = 5,
) -> Path | None:
    """Run MethBat if required and return the output path."""

    output_prefix = methbat_dir / label
    methbat_out = Path(str(output_prefix) + ".5mC.bed.gz")

    if methbat_out.exists() and methbat_out.stat().st_size > 0:
        print(
            f"  [{label}] MethBat output exists -> "
            f"{methbat_out.name}"
        )

    else:
        print(f"  [{label}] Running MethBat pileup...")

        run_methbat(
            bam=bam,
            output_prefix=output_prefix,
            threads=threads,
            min_coverage=min_coverage,
            min_hap_coverage=min_hap_coverage,
        )

        print(f"  [{label}] Written -> {methbat_out.name}")

    if not methbat_out.exists() or methbat_out.stat().st_size == 0:
        print(f"  [{label}] WARNING: no data written for {bam.name}")
        return None

    return methbat_out


# ---------------------------------------------------------------------------
# Haplotype-resolved counts
# ---------------------------------------------------------------------------

# MethBat 'type' values for the haplotype rows (case-insensitive). Check with:
#   zcat <sample>.5mC.bed.gz | grep -v '^#' | cut -f8 | sort | uniq -c
HAP1_TYPES = {"hap1", "h1"}
HAP2_TYPES = {"hap2", "h2"}


def extract_haplotype_counts(
    df: pd.DataFrame,
    label: str,
    min_hap_coverage: int = 5,
) -> pd.DataFrame | None:
    """Per-site Hap1/Hap2 counts, one row per CpG, for sites where BOTH
    haplotypes have >= min_hap_coverage reads.

    Requiring both haplotypes at the same CpG (inner join) matters: otherwise
    an island's hap1 counts and hap2 counts would be summed over different
    sets of CpGs, and a hap1-vs-hap2 difference could just reflect which
    CpGs each haplotype happened to cover.

    Returns None if the file has no haplotype rows (BAM not haplotagged).
    """
    coords = ["chrom", "start", "end"]
    t = df["type"].astype(str).str.lower()

    def _one(types: set, tag: str) -> pd.DataFrame:
        sub = df.loc[t.isin(types), coords + ["coverage", "mod_count", "unmod_count"]]
        sub = sub[sub["coverage"] >= min_hap_coverage]
        return sub.drop(columns="coverage").rename(columns={
            "mod_count": f"{tag}_mod_count_{label}",
            "unmod_count": f"{tag}_unmod_count_{label}",
        })

    h1, h2 = _one(HAP1_TYPES, "hap1"), _one(HAP2_TYPES, "hap2")
    if not t.isin(HAP1_TYPES | HAP2_TYPES).any():
        return None
    return h1.merge(h2, on=coords, how="inner")


# ---------------------------------------------------------------------------
# Load sample methylation
# ---------------------------------------------------------------------------

def load_sample_methylation(
    label: str,
    bam: Path,
    methbat_out: Path,
    min_coverage: int,
    min_hap_coverage: int = 5,
) -> pd.DataFrame | None:
    """Load and filter MethBat methylation data."""

    df = read_pbcpgtools(str(methbat_out))

    # Normalise chromosome column name if parser retains BED header.
    if "#chrom" in df.columns:
        df = df.rename(columns={"#chrom": "chrom"})

    n_raw = len(df)

    if n_raw == 0:
        print(f"  [{label}] WARNING: no methylation sites found")
        return None

    # Haplotype rows must be pulled out BEFORE the Total-only filter below.
    hap_df = extract_haplotype_counts(df, label, min_hap_coverage)
    if hap_df is None:
        print(f"  [{label}] WARNING: no haplotype rows in MethBat output "
              f"(BAM not haplotagged?) -- ASM columns will be empty")

    # Keep combined/unphased methylation only.
    df = df[df["type"] == "Total"].copy()

    # Apply minimum coverage filter.
    df = df[df["coverage"] >= min_coverage].copy()

    if df.empty:
        print(f"  [{label}] WARNING: no sites pass filters")
        return None

    df = df[
        [
            "chrom",
            "start",
            "end",
            "mod_score",
            "coverage",
            "mod_count",
            "unmod_count",
        ]
    ].rename(
        columns={
            "mod_score": f"mod_score_{label}",
            "coverage": f"coverage_{label}",
            "mod_count": f"mod_count_{label}",
            "unmod_count": f"unmod_count_{label}",
        }
    )

    print(
        f"  [{label}] {n_raw:,} raw rows -> "
        f"{len(df):,} Total sites after filtering"
    )

    # Left join: Total sites without usable haplotype data keep NaN hap columns.
    if hap_df is not None:
        df = df.merge(hap_df, on=["chrom", "start", "end"], how="left")
        n_hap = df[f"hap1_mod_count_{label}"].notna().sum()
        print(f"  [{label}] {n_hap:,} sites with >= {min_hap_coverage} reads on both haplotypes")

    return df


# ---------------------------------------------------------------------------
# Main sample function
# ---------------------------------------------------------------------------

def get_sample_methylation(
    label: str,
    bam: Path,
    methbat_dir: Path,
    threads: int = 6,
    min_coverage: int = 10,
    min_hap_coverage: int = 5,
) -> pd.DataFrame | None:
    """Run MethBat, load output and return per-sample methylation."""

    methbat_out = ensure_methbat_output(
        label=label,
        bam=bam,
        methbat_dir=methbat_dir,
        threads=threads,
        min_coverage=min_coverage,
        min_hap_coverage=min_hap_coverage,
    )

    if methbat_out is None:
        return None

    return load_sample_methylation(
        label=label,
        bam=bam,
        methbat_out=methbat_out,
        min_coverage=min_coverage,
        min_hap_coverage=min_hap_coverage,
    )
