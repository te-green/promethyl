
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
) -> None:
    """Run MethBat pileup on a single PacBio BAM."""

    output_prefix.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        "methbat", "pileup",
        "--input-bam", str(bam),
        "--output-prefix", str(output_prefix),
        "--threads", str(threads),
        "--min-coverage", str(min_coverage),
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
        )

        print(f"  [{label}] Written -> {methbat_out.name}")

    if not methbat_out.exists() or methbat_out.stat().st_size == 0:
        print(f"  [{label}] WARNING: no data written for {bam.name}")
        return None

    return methbat_out


# ---------------------------------------------------------------------------
# Load sample methylation
# ---------------------------------------------------------------------------

def load_sample_methylation(
    label: str,
    bam: Path,
    methbat_out: Path,
    min_coverage: int,
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

    return df


# ---------------------------------------------------------------------------
# Main sample function
# ---------------------------------------------------------------------------

def get_sample_methylation(
    label: str,
    bam: Path,
    methbat_dir: Path,
    threads: int = 10,
    min_coverage: int = 10,
) -> pd.DataFrame | None:
    """Run MethBat, load output and return per-sample methylation."""

    methbat_out = ensure_methbat_output(
        label=label,
        bam=bam,
        methbat_dir=methbat_dir,
        threads=threads,
        min_coverage=min_coverage,
    )

    if methbat_out is None:
        return None

    return load_sample_methylation(
        label=label,
        bam=bam,
        methbat_out=methbat_out,
        min_coverage=min_coverage,
    )
