#!/usr/bin/env python3
"""Reassemble candidate-read FASTA files from repository-safe byte shards."""

from __future__ import annotations

import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DERIVED = ROOT / "data" / "derived"
PARTS = DERIVED / "parts"

TARGETS = {
    "SRR25731498.candidates.fasta":
        "649af028cd6b59faa7dfcb7a2466837af507adb3f16df1f26ce9180156788fd2",
    "SRR34169110.candidates.union.fasta":
        "7426a8021609080c792f212f292b4a5c1f781b2508764d09b8de692846600094",
    "SRR34169110_BLAST5000_candidates.fasta":
        "12bc0a2f21a929eceec1d30b74dad16e9933084d9f1d7731db69897cb8806cdb",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def assemble(name: str, expected: str) -> None:
    target = DERIVED / name
    if target.exists() and sha256(target) == expected:
        print(f"Verified {target.relative_to(ROOT)}")
        return

    shards = sorted(PARTS.glob(f"{name}.part*"))
    if not shards:
        raise FileNotFoundError(f"No shards found for {name}")

    temporary = target.with_suffix(target.suffix + ".tmp")
    with temporary.open("wb") as output:
        for shard in shards:
            output.write(shard.read_bytes())

    observed = sha256(temporary)
    if observed != expected:
        temporary.unlink(missing_ok=True)
        raise ValueError(
            f"SHA-256 mismatch for {name}: expected {expected}, observed {observed}"
        )
    temporary.replace(target)
    print(f"Assembled {target.relative_to(ROOT)} from {len(shards)} shards")


def main() -> None:
    for name, expected in TARGETS.items():
        assemble(name, expected)


if __name__ == "__main__":
    main()
