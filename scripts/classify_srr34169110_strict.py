#!/usr/bin/env python3
"""Strictly classify SRR34169110 reads by the 17-site mtDNA/NUMT panel.

The original 100-read reconstruction set supplies the primary conservative
mtDNA-like class.  The saturated search may add a read only when it has at
least 10 mtDNA diagnostic alleles, no shared young-NUMT allele, all three
repeat arrays, and zero anchor mismatches across those arrays.  This explicit
extension rule reproduces the one added molecule SRR34169110.3087023.1 and
prevents a lower-quality partial match from silently changing the quantitative
subset.
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path

from analyze_srr25731498_mt_numt import call_positions, read_fasta


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reads", type=Path, required=True)
    parser.add_argument("--repeats", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--diagnostic-sites", type=Path, required=True)
    parser.add_argument("--legacy-classification", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--flank", type=int, default=11)
    parser.add_argument("--min-mtdna-hits", type=int, default=10)
    args = parser.parse_args()

    reads = read_fasta(args.reads)
    reference_records = read_fasta(args.reference)
    if len(reference_records) != 1:
        raise SystemExit("Expected one reference sequence")
    reference = next(iter(reference_records.values()))

    diagnostic_rows = read_tsv(args.diagnostic_sites)
    positions = [int(row["CHE_coordinate"]) for row in diagnostic_rows]
    mt_alleles = {
        int(row["CHE_coordinate"]): row["mtDNA_allele"] for row in diagnostic_rows
    }
    numt_alleles = {
        int(row["CHE_coordinate"]): row["young_NUMT_allele"]
        for row in diagnostic_rows
    }
    calls = call_positions(reads, reference, positions, args.flank)

    with args.repeats.open(newline="") as handle:
        repeat_rows = {
            row["read_id"]: row for row in csv.DictReader(handle, delimiter="\t")
        }
    with args.legacy_classification.open(newline="") as handle:
        legacy_class = {
            row["read_id"]: row["classification"]
            for row in csv.DictReader(handle, delimiter="\t")
        }

    output_rows: list[dict[str, object]] = []
    for read_id, seq in reads.items():
        read_calls = calls[read_id]
        mt_positions = [
            pos for pos in positions if read_calls.get(pos) == mt_alleles[pos]
        ]
        numt_positions = [
            pos for pos in positions if read_calls.get(pos) == numt_alleles[pos]
        ]
        other_positions = [
            pos
            for pos in positions
            if read_calls.get(pos) not in {None, mt_alleles[pos], numt_alleles[pos]}
        ]
        repeat = repeat_rows.get(read_id, {})
        if legacy_class.get(read_id) == "mtDNA-like":
            classification = "strict_mtDNA-like"
            classification_basis = "legacy_100_read_mtDNA_class"
        elif (
            len(mt_positions) >= args.min_mtdna_hits
            and not numt_positions
            and repeat.get("arrays_called") == "3"
            and repeat.get("total_anchor_mismatches") == "0"
        ):
            classification = "strict_mtDNA-like"
            classification_basis = "saturated_search_strict_extension"
        elif len(mt_positions) >= args.min_mtdna_hits and not numt_positions:
            classification = "mtDNA-like_array-specific"
            classification_basis = "diagnostic_panel_with_array_quality_filter"
        elif len(numt_positions) >= 10 and len(mt_positions) <= 2:
            classification = "strict_young_NUMT-like"
            classification_basis = "diagnostic_panel"
        else:
            classification = "unclassified"
            classification_basis = "insufficient_for_strict_class"
        output = {
                "read_id": read_id,
                "read_length": len(seq),
                "mtDNA_signature_hits": len(mt_positions),
                "young_NUMT_signature_hits": len(numt_positions),
                "other_diagnostic_calls": len(other_positions),
                "classification": classification,
                "classification_basis": classification_basis,
                "R27": repeat.get("R27", ""),
                "R56": repeat.get("R56", ""),
                "R80": repeat.get("R80", ""),
                "arrays_called": repeat.get("arrays_called", ""),
                "total_anchor_mismatches": repeat.get(
                    "total_anchor_mismatches", ""
                ),
                "mtDNA_signature_positions": ",".join(map(str, mt_positions)),
                "young_NUMT_signature_positions": ",".join(
                    map(str, numt_positions)
                ),
                "other_diagnostic_positions": ",".join(
                    map(str, other_positions)
                ),
            }
        for array in ("R27", "R56", "R80"):
            output[f"{array}_residual_bp"] = repeat.get(
                f"{array}_residual_bp", ""
            )
            output[f"{array}_anchor_mm"] = repeat.get(
                f"{array}_anchor_mm", ""
            )
            include = bool(repeat.get(array)) and (
                classification == "strict_mtDNA-like"
                or (
                    classification == "mtDNA-like_array-specific"
                    and repeat.get(f"{array}_anchor_mm") == "0"
                )
            )
            output[f"{array}_quantitative_subset"] = "yes" if include else "no"
        output_rows.append(output)

    write_tsv(args.out, output_rows)
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    with args.summary.open("w") as handle:
        handle.write(f"Input reads: {len(output_rows)}\n")
        for label, n in sorted(Counter(row["classification"] for row in output_rows).items()):
            handle.write(f"{label}: {n}\n")
        for array in ("R27", "R56", "R80"):
            counts = Counter(
                int(row[array])
                for row in output_rows
                if row[f"{array}_quantitative_subset"] == "yes"
            )
            handle.write(f"{array} informative: {sum(counts.values())}\n")
            for copies, n in sorted(counts.items()):
                handle.write(f"  {copies}: {n}\n")
    print(args.out)
    print(args.summary)


if __name__ == "__main__":
    main()
