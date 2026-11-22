# Reproducibility notes

## Candidate discovery versus final inference

Candidate discovery was intentionally permissive. Final mitochondrial inference
requires agreement of independent evidence: distributed mtDNA alleles, absence
of the young-NUMT signature, accepted repeat-array boundaries, and, for
SRR25731498, absence of direct nuclear-flank linkage.

The 51-read SRR25731498 candidate FASTA therefore contains two N-like reads for
which the current strict caller leaves all repeat fields blank. They remain in
the panel and flank analyses but do not contribute to mitochondrial repeat
frequencies. In prose, this set should be called **repeat-associated candidate
reads**, not 51 strict repeat-informative reads.

## Circular genome and R27

The physical ends of a raw read are never joined computationally. R27 is called
only if both unique boundaries occur in the required order within one original
read. The circular-spanning reconstruction read remains evidence of molecular
circularity, but its apparent wrapped R27 span is not used as a quantitative
repeat call.

## Singleton variants

The four manuscript-linked singleton variants are checked by local affine-gap
alignment. Their rarity limits frequency estimation but not identification of
the indel sequence. The validation script requires the exact 56-bp deletion in
`SRR34169110.3087023.1` and verifies that all reported singleton checks are
present.

## Archived nuclear-flank result

The direct-flank table is distributed as a derived result because recreating it
requires an NCBI network request. `run_pipeline.py --with-ncbi-flanks` refetches
the flanks and overwrites that table; the default offline workflow validates the
archived result.
