
import subprocess
from pathlib import Path

import pandas as pd

from bedmethyl import read_pbcpg


# ---------------------------------------------------------------------------
# MethBat
# ---------------------------------------------------------------------------

def run_methbat(
    bam: Path,
    output_prefix: Path,
    threads: int = 10,
    region: str = None,
    include_bed: Path = None,
) -> None:
    """Run methbat pileup on a single BAM."""

    cmd = [
        "methbat", "pileup",
        "--input-bam", str(bam),
        "--output-prefix", str(output_prefix),
        "--threads", str(threads),
        "--skip-5hmC",
        "--skip-6mA",
    ]

    # Note: region and include_bed are omitted here because
    # they are not documented as methbat pileup options.
    # Use an appropriate BAM preprocessing step if needed.

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
    ref: str = None,
    threads: int = 10,
    region: str = None,
    include_bed: Path = None,
) -> Path | None:

    """Run methbat if required and return the output path."""

    output_prefix = methbat_dir / bam.stem
    methbat_out = Path(str(output_prefix) + ".5mC.bed.gz")

    if methbat_out.exists():
        print(f"  [{label}] MethBat output exists -> {methbat_out.name}")

    else:
        print(f"  [{label}] Running MethBat pileup...")

        run_methbat(
            bam=bam,
            output_prefix=output_prefix,
            threads=threads,
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

    df = read_pbcpg(str(methbat_out))

    n_raw = len(df)

    if n_raw == 0:
        print(f"  [{label}] WARNING: no methylation sites found")
        return None

    # Minimum coverage filter
    df = df[df["cov"] >= min_coverage].copy()

    if df.empty:
        print(f"  [{label}] WARNING: no sites pass coverage filter")
        return None

    df = df[
        [
            "chrom",
            "start",
            "end",
            "mod_score",
            "cov",
            "mod_count",
            "unmod_count",
        ]
    ].rename(
        columns={
            "mod_score": f"mod_score_{label}",
            "cov": f"coverage_{label}",
            "mod_count": f"mod_count_{label}",
            "unmod_count": f"unmod_count_{label}",
        }
    )

    print(
        f"  [{label}] {n_raw:,} raw sites -> "
        f"{len(df):,} sites after coverage filtering"
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
    )

    if methbat_out is None:
        return None

    return load_sample_methylation(
        label=label,
        bam=bam,
        methbat_out=methbat_out,
        min_coverage=min_coverage,
    )
