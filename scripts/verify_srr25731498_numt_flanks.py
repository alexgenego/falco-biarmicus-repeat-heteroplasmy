#!/usr/bin/env python3
"""Verify direct linkage of SRR25731498 reads to young-NUMT nuclear flanks.

The three assembled young NUMT copies occur on JAMOAU010000287.1.  This
script downloads 1-kb flanks around each copy from NCBI, removes 21-mers that
also occur in the mitochondrial reference or in any M-like read, and counts
the remaining nuclear-specific flank 21-mers in each classified read.
"""

from __future__ import annotations

import argparse
import csv
import urllib.parse
import urllib.request
from pathlib import Path


COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
COPIES = (
    ("copy_101034_116265", 101034, 116265),
    ("copy_447310_462543", 447310, 462543),
    ("copy_465650_480881", 465650, 480881),
)


def reverse_complement(seq: str) -> str:
    return seq.translate(COMP)[::-1]


def parse_fasta_text(text: str) -> dict[str, str]:
    records: dict[str, list[str]] = {}
    name = None
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            name = line[1:].split()[0]
            records[name] = []
        elif name is not None:
            records[name].append(line.upper())
    return {name: "".join(parts) for name, parts in records.items()}


def read_fasta(path: Path) -> dict[str, str]:
    return parse_fasta_text(path.read_text())


def fetch_region(accession: str, start: int, end: int) -> str:
    params = urllib.parse.urlencode(
        {
            "db": "nuccore",
            "id": accession,
            "seq_start": start,
            "seq_stop": end,
            "rettype": "fasta",
            "retmode": "text",
        }
    )
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + params
    with urllib.request.urlopen(url, timeout=60) as response:
        records = parse_fasta_text(response.read().decode())
    if len(records) != 1:
        raise RuntimeError(f"Unexpected NCBI response for {accession}:{start}-{end}")
    return next(iter(records.values()))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--reads",
        type=Path,
        default=Path("data/derived/SRR25731498.candidates.fasta"),
    )
    parser.add_argument(
        "--classification",
        type=Path,
        default=Path("results/SRR25731498/SRR25731498_mtDNA_NUMT_read_classification.tsv"),
    )
    parser.add_argument(
        "--mitochondrial-reference",
        type=Path,
        default=Path("data/reference/BIA_mtDNA_refprojected_NC026715_v1.fasta"),
    )
    parser.add_argument("--accession", default="JAMOAU010000287.1")
    parser.add_argument("--flank-bp", type=int, default=1000)
    parser.add_argument("--k", type=int, default=21)
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("results/SRR25731498/SRR25731498_NUMT_nuclear_flank_linkage.tsv"),
    )
    args = parser.parse_args()

    reads = read_fasta(args.reads)
    mt_ref = next(iter(read_fasta(args.mitochondrial_reference).values()))
    with args.classification.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    # Exclude any sequence word observed in M-like molecules, in either strand,
    # so that a BIA-specific mitochondrial word cannot be mistaken for a
    # nuclear flank merely because it is absent from the first mtDNA reference.
    mitochondrial_words = set()
    for row in rows:
        if row["sequence_background"] != "M-like":
            continue
        seq = reads[row["read_id"]]
        for oriented in (seq, reverse_complement(seq)):
            mitochondrial_words.update(
                oriented[i : i + args.k] for i in range(len(oriented) - args.k + 1)
            )

    flank_words: dict[str, set[str]] = {}
    for copy_name, copy_start, copy_end in COPIES:
        region_start = copy_start - args.flank_bp
        region_end = copy_end + args.flank_bp
        region = fetch_region(args.accession, region_start, region_end)
        flanks = {
            "left": region[: args.flank_bp],
            "right": region[-args.flank_bp :],
        }
        for side, flank in flanks.items():
            words = {
                flank[i : i + args.k]
                for i in range(len(flank) - args.k + 1)
                if "N" not in flank[i : i + args.k]
            }
            words = {
                word
                for word in words
                if word not in mt_ref
                and reverse_complement(word) not in mt_ref
                and word not in mitochondrial_words
                and reverse_complement(word) not in mitochondrial_words
            }
            flank_words[f"{copy_name}_{side}"] = words

    output_rows = []
    for row in rows:
        read_id = row["read_id"]
        seq = reads[read_id]
        oriented_reads = (seq, reverse_complement(seq))
        counts = {}
        for flank_name, words in flank_words.items():
            counts[flank_name] = sum(
                1
                for word in words
                if word in oriented_reads[0] or word in oriented_reads[1]
            )
        max_hits = max(counts.values(), default=0)
        linked_flanks = [name for name, count in counts.items() if count >= 3]
        if linked_flanks:
            final_class = "young_NUMT_direct"
        elif row["sequence_background"] == "M-like":
            final_class = "mtDNA-like"
        else:
            final_class = "unresolved"
        output_rows.append(
            {
                "read_id": read_id,
                "sequence_background": row["sequence_background"],
                "final_classification": final_class,
                "R27": row["R27"],
                "R56": row["R56"],
                "R80": row["R80"],
                "direct_nuclear_flank_link": "yes" if linked_flanks else "no",
                "max_specific_flank_21mer_hits": max_hits,
                "linked_flanks": ",".join(linked_flanks),
                **{name + "_hits": counts[name] for name in sorted(counts)},
            }
        )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fields = list(output_rows[0])
    with args.out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)

    for background in ("M-like", "N-like", "unclassified"):
        subset = [row for row in output_rows if row["sequence_background"] == background]
        linked = [row for row in subset if row["direct_nuclear_flank_link"] == "yes"]
        print(f"{background}: {len(linked)}/{len(subset)} with direct nuclear-flank linkage")
    phased_n = [
        row
        for row in output_rows
        if row["sequence_background"] == "N-like" and row["R56"] == "8"
    ]
    phased_n_linked = [
        row for row in phased_n if row["direct_nuclear_flank_link"] == "yes"
    ]
    print(
        "N-like reads with R56=8: "
        f"{len(phased_n_linked)}/{len(phased_n)} with direct nuclear-flank linkage"
    )
    final_counts = {}
    for row in output_rows:
        final_counts[row["final_classification"]] = (
            final_counts.get(row["final_classification"], 0) + 1
        )
    print("Final classifications:", ", ".join(f"{k}={v}" for k, v in sorted(final_counts.items())))
    print(args.out)


if __name__ == "__main__":
    main()
