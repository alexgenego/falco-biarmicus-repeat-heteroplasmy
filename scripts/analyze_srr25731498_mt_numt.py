#!/usr/bin/env python3
"""Classify repeat-informative SRR25731498 reads as mtDNA-like or young-NUMT-like.

The analysis uses two independent layers:
1. Two sequence backgrounds are discovered from long reads with contrasting R56
   states (4--6 versus 8 copies) and are separated by diagnostic substitutions
   distributed across the mitochondrial reference.
2. The biological identity of the backgrounds is assigned by comparison with
   three assembled young NUMT copies from the Falco biarmicus nuclear genome.

Only bases bracketed by exact reference flanks are called. This avoids global
alignment artefacts around the three tandem-repeat arrays.
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path


DNA_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def reverse_complement(seq: str) -> str:
    return seq.translate(DNA_COMP)[::-1]


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
                records[name] = []
            elif name is not None:
                records[name].append(line.upper())
    return {name: "".join(parts) for name, parts in records.items()}


def read_blast_seqdump(path: Path) -> dict[str, str]:
    return read_fasta(path)


def read_repeat_table(path: Path) -> dict[str, dict[str, str]]:
    with path.open(newline="") as handle:
        return {
            row["read_id"]: row
            for row in csv.DictReader(handle, delimiter="\t")
        }


def integer(row: dict[str, str], field: str) -> int | None:
    value = row.get(field, "").strip()
    return int(value) if value else None


def call_base(
    orientations: tuple[str, str], reference: str, pos1: int, flank: int
) -> str | None:
    """Call one reference-coordinate base using exact flanks in either strand."""
    pos0 = pos1 - 1
    if pos0 < flank or pos0 + flank >= len(reference):
        return None
    left = reference[pos0 - flank : pos0]
    right = reference[pos0 + 1 : pos0 + 1 + flank]
    calls: set[str] = set()
    for oriented in orientations:
        start = 0
        while True:
            hit = oriented.find(left, start)
            if hit < 0:
                break
            center = hit + flank
            if (
                center + 1 + flank <= len(oriented)
                and oriented[center + 1 : center + 1 + flank] == right
                and oriented[center] in "ACGT"
            ):
                calls.add(oriented[center])
            start = hit + 1
    return next(iter(calls)) if len(calls) == 1 else None


def call_positions(
    reads: dict[str, str], reference: str, positions: list[int], flank: int
) -> dict[str, dict[int, str]]:
    result: dict[str, dict[int, str]] = {}
    for read_id, seq in reads.items():
        orientations = (seq, reverse_complement(seq))
        calls = {}
        for pos in positions:
            base = call_base(orientations, reference, pos, flank)
            if base is not None:
                calls[pos] = base
        result[read_id] = calls
    return result


def consensus(
    calls_by_read: dict[str, dict[int, str]],
    read_ids: list[str],
    positions: list[int],
    min_reads: int,
    min_fraction: float,
) -> dict[int, tuple[str, int, float]]:
    result = {}
    for pos in positions:
        counts = Counter(
            calls_by_read[read_id][pos]
            for read_id in read_ids
            if pos in calls_by_read[read_id]
        )
        if not counts:
            continue
        base, count = counts.most_common(1)[0]
        total = sum(counts.values())
        fraction = count / total
        if total >= min_reads and fraction >= min_fraction:
            result[pos] = (base, total, fraction)
    return result


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--reads",
        type=Path,
        default=Path("data/derived/SRR25731498.candidates.fasta"),
    )
    parser.add_argument(
        "--repeats",
        type=Path,
        default=Path("results/SRR25731498/repeat_calls.tsv"),
    )
    parser.add_argument(
        "--reference",
        type=Path,
        default=Path("data/reference/BIA_mtDNA_refprojected_NC026715_v1.fasta"),
    )
    parser.add_argument(
        "--known",
        type=Path,
        default=Path("data/reference/BIA_mtDNA_vs_young_NUMT_common_refcore_v1.fasta"),
    )
    parser.add_argument(
        "--out-prefix",
        type=Path,
        default=Path("results/SRR25731498/SRR25731498_mtDNA_NUMT"),
    )
    parser.add_argument("--flank", type=int, default=11)
    args = parser.parse_args()

    reads = read_blast_seqdump(args.reads)
    repeats = read_repeat_table(args.repeats)
    full_ref_records = read_fasta(args.reference)
    if len(full_ref_records) != 1:
        raise SystemExit("Expected one full mitochondrial reference record")
    full_ref = next(iter(full_ref_records.values()))
    known = read_fasta(args.known)
    common_ref = known["BIA_mtDNA_refprojected"]

    # The distributed-sequence classification is applied to all 51 reads from
    # the initial permissive repeat-associated candidate screen.  The final
    # anchor caller retains strict repeat states for 49 of them; two young-NUMT
    # candidates remain useful for sequence-background and nuclear-flank
    # classification even though no array passes the final anchor filter.
    informative_ids = [read_id for read_id in reads if read_id in repeats]
    seed_m = [
        read_id
        for read_id in informative_ids
        if integer(repeats[read_id], "arrays_called") >= 2
        and integer(repeats[read_id], "R56") in {4, 5, 6}
    ]
    seed_n = [
        read_id
        for read_id in informative_ids
        if integer(repeats[read_id], "arrays_called") >= 2
        and integer(repeats[read_id], "R56") == 8
    ]
    if len(seed_m) < 5 or len(seed_n) < 5:
        raise SystemExit("Insufficient seed reads for the two sequence backgrounds")

    seed_reads = {read_id: reads[read_id] for read_id in seed_m + seed_n}
    full_positions = list(range(1, len(full_ref) + 1))
    seed_calls = call_positions(seed_reads, full_ref, full_positions, args.flank)
    cons_m = consensus(seed_calls, seed_m, full_positions, 5, 0.80)
    cons_n = consensus(seed_calls, seed_n, full_positions, 5, 0.80)

    panel = []
    for pos in sorted(set(cons_m) & set(cons_n)):
        m_base, m_n, m_fraction = cons_m[pos]
        n_base, n_n, n_fraction = cons_n[pos]
        if m_base != n_base:
            panel.append(
                {
                    "CHE_coordinate": pos,
                    "M_allele": m_base,
                    "N_allele": n_base,
                    "M_called_reads": m_n,
                    "M_modal_fraction": f"{m_fraction:.4f}",
                    "N_called_reads": n_n,
                    "N_modal_fraction": f"{n_fraction:.4f}",
                }
            )
    panel_positions = [row["CHE_coordinate"] for row in panel]
    panel_m = {row["CHE_coordinate"]: row["M_allele"] for row in panel}
    panel_n = {row["CHE_coordinate"]: row["N_allele"] for row in panel}

    candidate_reads = {read_id: reads[read_id] for read_id in informative_ids}
    candidate_calls = call_positions(candidate_reads, full_ref, panel_positions, args.flank)
    classification_rows = []
    class_by_read = {}
    for read_id in informative_ids:
        calls = candidate_calls[read_id]
        m_hits = sum(calls.get(pos) == panel_m[pos] for pos in panel_positions)
        n_hits = sum(calls.get(pos) == panel_n[pos] for pos in panel_positions)
        other = sum(
            calls.get(pos) not in {None, panel_m[pos], panel_n[pos]}
            for pos in panel_positions
        )
        discriminating = m_hits + n_hits
        if discriminating >= 5 and m_hits / discriminating >= 0.80:
            read_class = "M-like"
        elif discriminating >= 5 and n_hits / discriminating >= 0.80:
            read_class = "N-like"
        else:
            read_class = "unclassified"
        class_by_read[read_id] = read_class
        row = repeats[read_id]
        classification_rows.append(
            {
                "read_id": read_id,
                "read_length": row["read_length"],
                "sequence_background": read_class,
                "M_panel_hits": m_hits,
                "N_panel_hits": n_hits,
                "other_panel_calls": other,
                "discriminating_calls": discriminating,
                "R27": row["R27"],
                "R56": row["R56"],
                "R80": row["R80"],
                "arrays_called": row["arrays_called"],
                "total_anchor_mismatches": row["total_anchor_mismatches"],
            }
        )

    # Independently label M and N by their agreement with the assembled young NUMTs.
    # The strict common-core FASTA omits two CHE-reference bases and therefore has
    # its own coordinate axis. Call bases against that sequence directly rather
    # than reusing full-reference coordinate calls.
    common_positions = list(range(1, len(common_ref) + 1))
    common_seed_calls = call_positions(seed_reads, common_ref, common_positions, args.flank)
    common_cons_m = consensus(common_seed_calls, seed_m, common_positions, 5, 0.80)
    common_cons_n = consensus(common_seed_calls, seed_n, common_positions, 5, 0.80)
    distance_rows = []
    for background, cons in (("M", common_cons_m), ("N", common_cons_n)):
        for record_name, seq in known.items():
            compared = 0
            mismatches = 0
            for pos in common_positions:
                if pos not in cons:
                    continue
                compared += 1
                mismatches += cons[pos][0] != seq[pos - 1]
            distance_rows.append(
                {
                    "sequence_background": background,
                    "reference_sequence": record_name,
                    "compared_sites": compared,
                    "mismatches": mismatches,
                    "identity": f"{1 - mismatches / compared:.8f}" if compared else "",
                }
            )

    # Alleles at the sites shared by all three known NUMT copies versus mtDNA.
    known_shared = []
    copy_names = [name for name in known if name != "BIA_mtDNA_refprojected"]
    for pos in common_positions:
        mt_base = common_ref[pos - 1]
        copy_bases = {known[name][pos - 1] for name in copy_names}
        if len(copy_bases) == 1 and mt_base not in copy_bases:
            numt_base = next(iter(copy_bases))
            known_shared.append(
                {
                    "common_core_coordinate": pos,
                    "first_individual_mtDNA_allele": mt_base,
                    "shared_young_NUMT_allele": numt_base,
                    "M_consensus": common_cons_m.get(pos, ("", 0, 0))[0],
                    "N_consensus": common_cons_n.get(pos, ("", 0, 0))[0],
                    "M_support": common_cons_m.get(pos, ("", 0, 0))[1],
                    "N_support": common_cons_n.get(pos, ("", 0, 0))[1],
                }
            )

    write_tsv(
        Path(str(args.out_prefix) + "_diagnostic_panel.tsv"),
        panel,
        [
            "CHE_coordinate",
            "M_allele",
            "N_allele",
            "M_called_reads",
            "M_modal_fraction",
            "N_called_reads",
            "N_modal_fraction",
        ],
    )
    write_tsv(
        Path(str(args.out_prefix) + "_read_classification.tsv"),
        classification_rows,
        [
            "read_id",
            "read_length",
            "sequence_background",
            "M_panel_hits",
            "N_panel_hits",
            "other_panel_calls",
            "discriminating_calls",
            "R27",
            "R56",
            "R80",
            "arrays_called",
            "total_anchor_mismatches",
        ],
    )
    write_tsv(
        Path(str(args.out_prefix) + "_distance_to_known_sequences.tsv"),
        distance_rows,
        ["sequence_background", "reference_sequence", "compared_sites", "mismatches", "identity"],
    )
    write_tsv(
        Path(str(args.out_prefix) + "_known_NUMT_sites.tsv"),
        known_shared,
        [
            "common_core_coordinate",
            "first_individual_mtDNA_allele",
            "shared_young_NUMT_allele",
            "M_consensus",
            "N_consensus",
            "M_support",
            "N_support",
        ],
    )

    class_counts = Counter(class_by_read.values())
    repeat_combinations: dict[str, Counter] = defaultdict(Counter)
    for read_id, read_class in class_by_read.items():
        row = repeats[read_id]
        combo = (row["R27"] or "NA", row["R56"] or "NA", row["R80"] or "NA")
        repeat_combinations[read_class][combo] += 1

    summary_path = Path(str(args.out_prefix) + "_summary.txt")
    with summary_path.open("w") as out:
        out.write(f"Input reads: {len(reads)}\n")
        out.write(f"Repeat-associated candidate reads: {len(informative_ids)}\n")
        out.write(f"M seed reads (R56=4--6; >=2 arrays): {len(seed_m)}\n")
        out.write(f"N seed reads (R56=8; >=2 arrays): {len(seed_n)}\n")
        out.write(f"Diagnostic panel positions: {len(panel)}\n")
        out.write("Classified repeat-associated candidate reads:\n")
        for read_class, count in sorted(class_counts.items()):
            out.write(f"  {read_class}: {count}\n")
        out.write("Repeat-state combinations by class (R27,R56,R80):\n")
        for read_class in sorted(repeat_combinations):
            out.write(f"  {read_class}:\n")
            for combo, count in sorted(repeat_combinations[read_class].items()):
                out.write(f"    {combo}: {count}\n")
        out.write("\nDistances of group consensuses to known sequences:\n")
        for row in distance_rows:
            out.write(
                f"  {row['sequence_background']} vs {row['reference_sequence']}: "
                f"{row['mismatches']}/{row['compared_sites']} mismatches; "
                f"identity {row['identity']}\n"
            )

    print(summary_path)


if __name__ == "__main__":
    main()
