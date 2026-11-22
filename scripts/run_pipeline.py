#!/usr/bin/env python3
"""Run the complete local analysis and validate manuscript-linked outputs."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PYTHON = sys.executable


def run(*args: str) -> None:
    print("+", " ".join(args), flush=True)
    subprocess.run(args, cwd=ROOT, check=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--with-ncbi-flanks",
        action="store_true",
        help="Refetch the six 1-kb NUMT flanks from NCBI before validation.",
    )
    args = parser.parse_args()

    run(
        PYTHON, "scripts/analyze_repeat_arrays.py",
        "data/derived/SRR34169110.candidates.union.fasta",
        "--consensus", "data/reference/BIA_mtDNA_consensus_full_v1.fasta",
        "--out", "results/SRR34169110/repeat_calls.tsv",
        "--max-mm", "12", "--tolerance", "12", "--min-anchor-seeds", "2",
    )
    run(
        PYTHON, "scripts/analyze_repeat_arrays.py",
        "data/derived/SRR25731498.candidates.fasta",
        "--consensus", "data/reference/BIA_mtDNA_consensus_full_v1.fasta",
        "--out", "results/SRR25731498/repeat_calls.tsv",
        "--max-mm", "12", "--tolerance", "12", "--min-anchor-seeds", "1",
    )
    run(
        PYTHON, "scripts/classify_srr34169110_strict.py",
        "--reads", "data/derived/SRR34169110.candidates.union.fasta",
        "--repeats", "results/SRR34169110/repeat_calls.tsv",
        "--reference", "data/reference/BIA_mtDNA_refprojected_NC026715_v1.fasta",
        "--diagnostic-sites", "data/reference/BIA_mtDNA_diagnostic_sites_v1.tsv",
        "--legacy-classification", "results/reconstruction/BIA_SRA_read_classification_v1.tsv",
        "--out", "results/SRR34169110/strict_read_classification.tsv",
        "--summary", "results/SRR34169110/strict_repeat_summary.txt",
    )
    run(PYTHON, "scripts/analyze_srr25731498_mt_numt.py")
    if args.with_ncbi_flanks:
        run(PYTHON, "scripts/verify_srr25731498_numt_flanks.py")
    run(PYTHON, "scripts/check_singleton_repeat_indels.py")
    run(PYTHON, "scripts/validate_results.py")


if __name__ == "__main__":
    main()
