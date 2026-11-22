#!/usr/bin/env python3
"""Call R80/R56/R27 copy-number states in long Falco biarmicus reads.

The caller is deliberately anchor based.  A state is accepted only when both
non-repetitive boundaries of the corresponding array are found in one read.
For R27 the right boundary is the beginning of the linearized circular
mitogenome, so the script also handles reads rotated across that junction.
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path


COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def rc(seq: str) -> str:
    return seq.translate(COMP)[::-1]


def read_fasta(path: Path):
    name = None
    chunks = []
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks).upper()
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if name is not None:
        yield name, "".join(chunks).upper()


def edit_distance(a: str, b: str) -> int:
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i]
        for j, cb in enumerate(b, 1):
            cur.append(min(cur[-1] + 1, prev[j] + 1, prev[j - 1] + (ca != cb)))
        prev = cur
    return prev[-1]


def anchor_hits(seq: str, anchor: str, max_mm: int, min_seeds: int = 2):
    """Return (edit distance, start, end) hits from clustered exact seeds."""
    k = 10
    estimates = []
    for off in (0, 10, 20, 30):
        seed = anchor[off : off + k]
        start = 0
        while True:
            j = seq.find(seed, start)
            if j < 0:
                break
            p = j - off
            if 0 <= p < len(seq):
                estimates.append(p)
            start = j + 1
    estimates.sort()
    clusters = []
    for p in estimates:
        if not clusters or p - clusters[-1][-1] > 6:
            clusters.append([p])
        else:
            clusters[-1].append(p)
    candidates = set()
    for cluster in clusters:
        if len(cluster) < min_seeds:
            continue
        center = cluster[len(cluster) // 2]
        for shift in range(-5, 6):
            q = center + shift
            if 0 <= q < len(seq):
                candidates.add(q)
    ans = []
    for p in candidates:
        best = None
        for length in range(len(anchor) - 5, len(anchor) + 6):
            if p + length > len(seq):
                continue
            ed = edit_distance(anchor, seq[p : p + length])
            item = (ed, abs(length - len(anchor)), p + length)
            if best is None or item < best:
                best = item
        if best is not None and best[0] <= max_mm:
            ans.append((best[0], p, best[2]))
    ans.sort()
    # Collapse overlapping representations of the same local hit.
    kept = []
    for hit in ans:
        if any(abs(hit[1] - x[1]) <= 5 for x in kept):
            continue
        kept.append(hit)
    return kept


def closest_state(span: int, unit: int, lo: int, hi: int, tolerance: int):
    n = round(span / unit)
    if n < lo or n > hi:
        return None
    residual = span - n * unit
    if abs(residual) > tolerance:
        return None
    return n, residual


def best_linear_pair(left_hits, right_hits, anchor_len, unit, limits, tolerance):
    best = None
    for lmm, lp, le in left_hits:
        for rmm, rp, re in right_hits:
            span = rp - le
            state = closest_state(span, unit, *limits, tolerance)
            if state is None:
                continue
            n, residual = state
            score = (lmm + rmm, abs(residual), lmm, rmm)
            item = (score, n, residual, lp, rp, lmm, rmm)
            if best is None or item < best:
                best = item
    return best


def best_circular_pair(seq_len, left_hits, right_hits, anchor_len, unit, limits, tolerance):
    best = None
    for lmm, lp, le in left_hits:
        for rmm, rp, re in right_hits:
            span = rp - le
            wrapped = False
            if span < 0:
                span += seq_len
                wrapped = True
            state = closest_state(span, unit, *limits, tolerance)
            if state is None:
                continue
            n, residual = state
            score = (lmm + rmm, abs(residual), not wrapped, lmm, rmm)
            item = (score, n, residual, lp, rp, lmm, rmm, wrapped)
            if best is None or item < best:
                best = item
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("reads", type=Path)
    ap.add_argument(
        "--consensus",
        type=Path,
        default=Path("data/reference/BIA_mtDNA_consensus_full_v1.fasta"),
    )
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--max-mm", type=int, default=8)
    ap.add_argument("--tolerance", type=int, default=12)
    ap.add_argument(
        "--min-anchor-seeds",
        type=int,
        choices=(1, 2),
        default=2,
        help="Exact 10-mer seeds required before approximate anchor extension.",
    )
    args = ap.parse_args()

    consensus = next(read_fasta(args.consensus))[1]
    anchor_len = 40
    arrays = {
        "R80": dict(unit=80, left=consensus[14966-anchor_len:14966], right=consensus[15206:15206+anchor_len], limits=(1, 6)),
        "R56": dict(unit=56, left=consensus[16152-anchor_len:16152], right=consensus[16432:16432+anchor_len], limits=(1, 16)),
        "R27": dict(unit=27, left=consensus[17303-anchor_len:17303], right=consensus[:anchor_len], limits=(1, 80)),
    }

    rows = []
    for read_id, raw_seq in read_fasta(args.reads):
        orientation_results = []
        for orient, seq in (("+", raw_seq), ("-", rc(raw_seq))):
            result = {"orientation": orient, "seq_len": len(seq)}
            total_anchor_mm = 0
            called = 0
            for label, spec in arrays.items():
                lh = anchor_hits(
                    seq, spec["left"], args.max_mm, args.min_anchor_seeds
                )
                rh = anchor_hits(
                    seq, spec["right"], args.max_mm, args.min_anchor_seeds
                )
                # Even for R27, both boundaries must occur contiguously in the
                # raw read.  The physical ends of a raw read must never be
                # joined computationally merely because the template is
                # circular; only a read that actually traverses the selected
                # linearization junction supplies an R27 state.
                pair = best_linear_pair(lh, rh, anchor_len, spec["unit"], spec["limits"], args.tolerance)
                if pair is None:
                    result[label] = None
                    result[label + "_detail"] = None
                else:
                    result[label] = pair[1]
                    result[label + "_detail"] = pair
                    total_anchor_mm += pair[5] + pair[6]
                    called += 1
            result["called"] = called
            result["total_anchor_mm"] = total_anchor_mm
            orientation_results.append(result)
        result = min(orientation_results, key=lambda x: (-x["called"], x["total_anchor_mm"], x["orientation"] != "+"))
        row = {
            "read_id": read_id,
            "read_length": len(raw_seq),
            "orientation": result["orientation"],
            "R27": result["R27"],
            "R56": result["R56"],
            "R80": result["R80"],
            "arrays_called": result["called"],
            "total_anchor_mismatches": result["total_anchor_mm"],
        }
        for label in ("R27", "R56", "R80"):
            d = result[label + "_detail"]
            if d is not None:
                row[label + "_residual_bp"] = d[2]
                row[label + "_left_pos"] = d[3]
                row[label + "_right_pos"] = d[4]
                row[label + "_anchor_mm"] = d[5] + d[6]
            else:
                row[label + "_residual_bp"] = ""
                row[label + "_left_pos"] = ""
                row[label + "_right_pos"] = ""
                row[label + "_anchor_mm"] = ""
        rows.append(row)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"reads={len(rows)}")
    for label in ("R27", "R56", "R80"):
        c = Counter(r[label] for r in rows if r[label] is not None)
        print(label, sum(c.values()), sorted(c.items()))
    phased = Counter((r["R27"], r["R56"], r["R80"]) for r in rows if r["arrays_called"] == 3)
    print("phased", sum(phased.values()))
    for state, n in phased.most_common():
        print(state, n)


if __name__ == "__main__":
    main()
