#!/usr/bin/env python3
"""Validate repository outputs against the numerical claims in the manuscript."""

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def main() -> None:
    first = tsv(ROOT / "results/SRR34169110/strict_read_classification.tsv")
    class_counts = Counter(row["classification"] for row in first)
    require(class_counts["strict_mtDNA-like"] == 29, "SRR34169110 strict class != 29")
    require(
        class_counts["mtDNA-like_array-specific"] == 1,
        "SRR34169110 array-specific class != 1",
    )

    expected_first = {
        "R27": {
            29: 1, 30: 1, 31: 1, 32: 1, 33: 1, 34: 1, 37: 1, 38: 1,
            39: 2, 41: 3, 47: 1, 48: 1, 49: 1, 50: 1, 51: 1, 52: 2,
            53: 1, 54: 2, 55: 1, 58: 1,
        },
        "R56": {4: 1, 5: 26},
        "R80": {3: 25},
    }
    for array, expected in expected_first.items():
        observed = Counter(
            int(row[array])
            for row in first
            if row[f"{array}_quantitative_subset"] == "yes"
        )
        require(observed == Counter(expected), f"{array} mismatch: {observed}")

    by_id = {row["read_id"]: row for row in first}
    added = by_id["gnl|SRA|SRR34169110.3087023.1"]
    require(
        added["classification_basis"] == "saturated_search_strict_extension",
        "R56=4 extension read not identified by the expected rule",
    )
    partial = by_id["gnl|SRA|SRR34169110.5516667.1"]
    require(partial["R27_quantitative_subset"] == "no", "R27 partial read included")
    require(partial["R56_quantitative_subset"] == "yes", "R56 partial read excluded")
    require(partial["R80_quantitative_subset"] == "yes", "R80 partial read excluded")

    second = tsv(
        ROOT / "results/SRR25731498/SRR25731498_mtDNA_NUMT_read_classification.tsv"
    )
    second_counts = Counter(row["sequence_background"] for row in second)
    require(
        second_counts == Counter({"M-like": 12, "N-like": 26, "unclassified": 13}),
        f"SRR25731498 background counts mismatch: {second_counts}",
    )
    panel = tsv(
        ROOT / "results/SRR25731498/SRR25731498_mtDNA_NUMT_diagnostic_panel.tsv"
    )
    require(len(panel) == 77, "SRR25731498 diagnostic panel != 77 positions")

    m_rows = [row for row in second if row["sequence_background"] == "M-like"]
    expected_second = {
        "R27": {18: 1, 19: 2, 20: 1, 38: 1, 41: 1, 46: 1, 48: 1, 49: 1},
        "R56": {4: 1, 5: 10, 6: 1},
        "R80": {2: 1, 3: 11},
    }
    for array, expected in expected_second.items():
        observed = Counter(int(row[array]) for row in m_rows if row[array])
        require(observed == Counter(expected), f"SRR25731498 {array} mismatch: {observed}")

    linkage = tsv(
        ROOT / "results/SRR25731498/SRR25731498_NUMT_nuclear_flank_linkage.tsv"
    )
    final_counts = Counter(row["final_classification"] for row in linkage)
    require(
        final_counts == Counter({"mtDNA-like": 12, "young_NUMT_direct": 39}),
        f"Nuclear-flank final counts mismatch: {final_counts}",
    )
    require(
        all(
            row["direct_nuclear_flank_link"] == "no"
            for row in linkage
            if row["final_classification"] == "mtDNA-like"
        ),
        "A mtDNA-like read has a nuclear-flank link",
    )
    require(
        all(
            row["direct_nuclear_flank_link"] == "yes"
            for row in linkage
            if row["final_classification"] == "young_NUMT_direct"
        ),
        "A young-NUMT read lacks a nuclear-flank link",
    )

    indels = (ROOT / "results/indel_checks/selected_repeat_indel_checks.txt").read_text()
    require("SRR34169110.3087023.1" in indels, "First-individual R56=4 check missing")
    require("deletion_from_read\tlength=56" in indels, "Exact 56-bp deletion missing")
    require("SRR25731498.691932.1" in indels, "Second-individual R56=6 check missing")
    require("SRR25731498.711747.1" in indels, "Second-individual R80=2 check missing")

    print("All manuscript-linked validation checks passed.")


if __name__ == "__main__":
    main()
