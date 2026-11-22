#!/usr/bin/env python3
"""Create the deduplicated SRR34169110 candidate-read set used in the study.

The saturated BLAST export contained the reads selected from the 5000-hit
search.  One high-confidence mtDNA-like read used in the original 100-read
classification was absent from that export.  This script adds every read named
in the legacy classification table from the original BLAST sequence dump and
writes a stable union FASTA.  The derived union FASTA is distributed with the
repository, so the two large BLAST exports are not required for routine reruns.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, list[str]] = {}
    name: str | None = None
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].split()[0]
                records.setdefault(name, [])
            elif name is not None:
                records[name].append(line.upper())
    return {name: "".join(parts) for name, parts in records.items()}


def write_fasta(path: Path, records: dict[str, str], width: int = 80) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for name, seq in records.items():
            handle.write(f">{name}\n")
            for start in range(0, len(seq), width):
                handle.write(seq[start : start + width] + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--blast5000-candidates", type=Path, required=True)
    parser.add_argument("--legacy-reads", type=Path, required=True)
    parser.add_argument("--legacy-classification", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    union = read_fasta(args.blast5000_candidates)
    legacy_reads = read_fasta(args.legacy_reads)
    with args.legacy_classification.open(newline="") as handle:
        legacy_ids = [
            row["read_id"] for row in csv.DictReader(handle, delimiter="\t")
        ]

    missing_from_source = [read_id for read_id in legacy_ids if read_id not in legacy_reads]
    if missing_from_source:
        raise SystemExit(
            "Legacy sequence dump lacks classified reads: "
            + ", ".join(missing_from_source)
        )

    added = []
    for read_id in legacy_ids:
        if read_id not in union:
            union[read_id] = legacy_reads[read_id]
            added.append(read_id)

    write_fasta(args.out, union)
    print(f"records={len(union)}")
    print(f"added_from_legacy_export={len(added)}")
    for read_id in added:
        print(f"  {read_id}")


if __name__ == "__main__":
    main()
