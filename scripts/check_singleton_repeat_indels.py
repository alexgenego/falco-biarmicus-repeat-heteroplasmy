#!/usr/bin/env python3
"""Verify the four singleton R56/R80 states reported in the manuscript."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
NEG = -10**9


def rc(seq: str) -> str:
    return seq.translate(COMP)[::-1]


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    name = None
    chunks: list[str] = []
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                records[name] = "".join(chunks).upper()
            name = line[1:].split()[0]
            chunks = []
        else:
            chunks.append(line)
    if name is not None:
        records[name] = "".join(chunks).upper()
    return records


def affine_global(query: str, ref: str) -> tuple[str, str, int]:
    """Global alignment; X is a gap in ref, Y is a gap in query."""
    n, m = len(query), len(ref)
    go, ge = -10, -1
    M = [[NEG] * (m + 1) for _ in range(n + 1)]
    X = [[NEG] * (m + 1) for _ in range(n + 1)]
    Y = [[NEG] * (m + 1) for _ in range(n + 1)]
    tM = [bytearray(m + 1) for _ in range(n + 1)]
    tX = [bytearray(m + 1) for _ in range(n + 1)]
    tY = [bytearray(m + 1) for _ in range(n + 1)]
    M[0][0] = 0
    for i in range(1, n + 1):
        X[i][0] = go + (i - 1) * ge
        tX[i][0] = 1
    for j in range(1, m + 1):
        Y[0][j] = go + (j - 1) * ge
        tY[0][j] = 2

    for i in range(1, n + 1):
        qi = query[i - 1]
        for j in range(1, m + 1):
            prev = (M[i - 1][j - 1], X[i - 1][j - 1], Y[i - 1][j - 1])
            state = max(range(3), key=prev.__getitem__)
            M[i][j] = prev[state] + (2 if qi == ref[j - 1] else -3)
            tM[i][j] = state

            prev_x = (M[i - 1][j] + go, X[i - 1][j] + ge, Y[i - 1][j] + go)
            state = max(range(3), key=prev_x.__getitem__)
            X[i][j] = prev_x[state]
            tX[i][j] = state

            prev_y = (M[i][j - 1] + go, X[i][j - 1] + go, Y[i][j - 1] + ge)
            state = max(range(3), key=prev_y.__getitem__)
            Y[i][j] = prev_y[state]
            tY[i][j] = state

    finals = (M[n][m], X[n][m], Y[n][m])
    state = max(range(3), key=finals.__getitem__)
    score = finals[state]
    aq: list[str] = []
    ar: list[str] = []
    i, j = n, m
    while i or j:
        if state == 0:
            aq.append(query[i - 1]); ar.append(ref[j - 1])
            state = tM[i][j]; i -= 1; j -= 1
        elif state == 1:
            aq.append(query[i - 1]); ar.append("-")
            state = tX[i][j]; i -= 1
        else:
            aq.append("-"); ar.append(ref[j - 1])
            state = tY[i][j]; j -= 1
    return "".join(reversed(aq)), "".join(reversed(ar)), score


def gap_runs(aq: str, ar: str) -> list[tuple[str, int, str, int]]:
    runs: list[tuple[str, int, str, int]] = []
    qpos = rpos = 0
    i = 0
    while i < len(aq):
        if aq[i] == "-" or ar[i] == "-":
            kind = "deletion_from_read" if aq[i] == "-" else "insertion_in_read"
            start_q, start_r = qpos, rpos
            seq: list[str] = []
            while i < len(aq) and ((aq[i] == "-") if kind == "deletion_from_read" else (ar[i] == "-")):
                if aq[i] != "-":
                    seq.append(aq[i]); qpos += 1
                if ar[i] != "-":
                    seq.append(ar[i]); rpos += 1
                i += 1
            runs.append((kind, len(seq), "".join(seq), start_r if kind == "deletion_from_read" else start_q))
            continue
        qpos += 1; rpos += 1; i += 1
    return runs


def identity(a: str, b: str) -> tuple[int, int, float]:
    matches = sum(x == y for x, y in zip(a, b) if x != "-" and y != "-")
    compared = sum(x != "-" and y != "-" for x, y in zip(a, b))
    return matches, compared, matches / compared if compared else 0.0


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--reads-341",
        type=Path,
        default=Path("data/derived/SRR34169110.candidates.union.fasta"),
    )
    parser.add_argument(
        "--calls-341",
        type=Path,
        default=Path("results/SRR34169110/repeat_calls.tsv"),
    )
    parser.add_argument(
        "--reads-257",
        type=Path,
        default=Path("data/derived/SRR25731498.candidates.fasta"),
    )
    parser.add_argument(
        "--calls-257",
        type=Path,
        default=Path("results/SRR25731498/repeat_calls.tsv"),
    )
    parser.add_argument(
        "--consensus",
        type=Path,
        default=Path("data/reference/BIA_mtDNA_consensus_full_v1.fasta"),
    )
    parser.add_argument(
        "--out-dir", type=Path, default=Path("results/indel_checks")
    )
    args = parser.parse_args()

    sample_inputs = {
        "SRR34169110": (args.reads_341, args.calls_341),
        "SRR25731498": (args.reads_257, args.calls_257),
    }
    targets = {
        "gnl|SRA|SRR34169110.3087023.1": (
            "SRR34169110", "R56", 56, 16152, 16432
        ),
        "gnl|SRA|SRR25731498.691932.1": (
            "SRR25731498", "R56", 56, 16152, 16432
        ),
        "gnl|SRA|SRR25731498.308159.1": (
            "SRR25731498", "R56", 56, 16152, 16432
        ),
        "gnl|SRA|SRR25731498.711747.1": (
            "SRR25731498", "R80", 80, 14966, 15206
        ),
    }

    reads_by_sample = {}
    calls_by_sample = {}
    for sample, (reads_path, calls_path) in sample_inputs.items():
        reads_by_sample[sample] = read_fasta(reads_path)
        with calls_path.open(newline="") as handle:
            calls_by_sample[sample] = {
                row["read_id"]: row
                for row in csv.DictReader(handle, delimiter="\t")
            }
    consensus = next(iter(read_fasta(args.consensus).values()))
    args.out_dir.mkdir(parents=True, exist_ok=True)

    report: list[str] = []
    for read_id, (sample, array, unit, left_boundary, right_boundary) in targets.items():
        row = calls_by_sample[sample][read_id]
        seq = reads_by_sample[sample][read_id]
        if row["orientation"] == "-":
            seq = rc(seq)
        left_pos = int(row[f"{array}_left_pos"])
        right_pos = int(row[f"{array}_right_pos"])
        flank = 120
        anchor_len = 40
        q0 = max(0, left_pos - flank)
        q1 = min(len(seq), right_pos + anchor_len + flank)
        r0 = left_boundary - anchor_len - flank
        r1 = right_boundary + anchor_len + flank
        query = seq[q0:q1]
        ref = consensus[r0:r1]
        aq, ar, score = affine_global(query, ref)
        runs = sorted(gap_runs(aq, ar), key=lambda x: x[1], reverse=True)
        matches, compared, ident = identity(aq, ar)
        report.extend([
            f"SAMPLE {sample}",
            f"READ {read_id}",
            f"array={array} called_copies={row[array]} orientation={row['orientation']}",
            f"query_segment={q0}:{q1} ({len(query)} bp); reference_segment={r0}:{r1} ({len(ref)} bp)",
            f"length_difference_query_minus_reference={len(query)-len(ref)} bp; alignment_score={score}",
            f"aligned_ungapped_identity={matches}/{compared} ({ident:.5f})",
            "gap_runs_len>=5:",
        ])
        for kind, length, gap_seq, pos in runs:
            if length >= 5:
                report.append(f"  {kind}\tlength={length}\tposition={pos}\tsequence={gap_seq}")
        report.append("")
        alignment_path = args.out_dir / f"{sample}_{read_id.split('.')[-2]}_{array}_alignment.txt"
        alignment_path.write_text(
            f">query\n{aq}\n>reference\n{ar}\n"
        )
    out = args.out_dir / "selected_repeat_indel_checks.txt"
    out.write_text("\n".join(report))
    print(out)
    print("\n".join(report))


if __name__ == "__main__":
    main()
