# Falco biarmicus repeat-number heteroplasmy

Code, derived candidate reads, reference sequences, and manuscript-linked
tables for the analysis of tandem-repeat copy-number heteroplasmy in the
mitochondrial control region of the lanner falcon (*Falco biarmicus*).

The analysis separates contemporary mitochondrial molecules from young nuclear
mitochondrial DNA segments (NUMTs) before calling the R27, R56, and R80 repeat
arrays in long PacBio HiFi reads.

## Main reproducible results

- **SRR34169110:** 25 strict R27-informative molecules contain 20 states from
  29 to 58 copies; R56 contains four- and five-copy states (1/27 and 26/27);
  all 25 R80-informative molecules contain three copies.
- **SRR25731498:** 77 distributed diagnostic positions separate two sequence
  backgrounds. Twelve reads are mtDNA-like; the remaining 39 candidate reads
  are directly linked to young-NUMT nuclear flanks. The mtDNA-like reads contain
  eight R27 states (18–49 copies), R56 states with four, five, and six copies,
  and R80 states with two and three copies.
- Local affine-gap alignments confirm that the four reported singleton R56/R80
  states differ by repeat-unit sequence rather than by an unrelated local
  deletion or insertion.

## Repository layout

```text
data/reference/       reconstructed and reference-projected mtDNA; NUMT panel
data/derived/parts/   byte-exact shards of the three candidate-read FASTA files
results/reconstruction/ reconstruction evidence and the original 100-read class
results/SRR34169110/  repeat calls, strict classifications, and summary
results/SRR25731498/  repeat calls, 77-site classification, and flank linkage
results/indel_checks/ local alignments of singleton R56/R80 states
scripts/              analysis, validation, and workflow scripts
docs/                 data dictionary and reproducibility notes
```

The large primary sequencing runs remain in NCBI SRA under **SRR34169110** and
**SRR25731498**. The compact derived FASTA files are included as numbered byte
shards so that the manuscript-linked analyses can be rerun without downloading
complete WGS runs. The workflow reassembles and verifies all three FASTA files
automatically.

## Requirements

- Python 3.10 or newer;
- no third-party Python packages;
- internet access only when `--with-ncbi-flanks` is used to refetch the six
  1-kb nuclear flanks from NCBI.

## Reproduce the analysis

From the repository root:

```bash
python3 scripts/run_pipeline.py
```

This reruns repeat calling, the two mtDNA/NUMT classifications, singleton-indel
checks, and all manuscript-linked assertions. Before analysis, the three
candidate-read FASTA files are reconstructed from `data/derived/parts/` and
checked against their published SHA-256 digests. The expected final message is:

```text
All manuscript-linked validation checks passed.
```

To refetch the nuclear flanks of the three NUMT copies and recreate the linkage
table before validation:

```bash
python3 scripts/run_pipeline.py --with-ncbi-flanks
```

The flank step queries `JAMOAU010000287.1` for the NUMT copies at
101034–116265, 447310–462543, and 465650–480881.

## Two candidate-screen levels

The repository deliberately distinguishes a permissive discovery screen from
the final calls.

1. The initial screen identified 51 repeat-associated candidates in
   SRR25731498.
2. The final anchor caller assigns at least one strict repeat state to 49 of
   them. The other two are still retained for distributed-sequence and direct
   nuclear-flank classification; both belong to the young-NUMT class and do
   not enter mitochondrial repeat counts.

For SRR34169110, the original 100-read reconstruction class is extended by one
fully phased read, `SRR34169110.3087023.1`, carrying R56 = 4. A second partial
read has sufficient mtDNA diagnostic evidence and exact R56/R80 anchors but a
one-mismatch R27 anchor. It contributes only to the R56 and R80 denominators;
the inclusion decision is explicit in
`results/SRR34169110/strict_read_classification.tsv`.

## Data provenance

- Primary reads: NCBI SRA SRR34169110 and SRR25731498.
- Mitochondrial coordinate reference: *Falco cherrug* NC_026715.1.
- Additional comparison reference: *Falco rusticolus* KT989235.1.
- Young NUMTs: JAMOAU010000287.1 coordinates listed above.
- Reconstructed working mitochondrial sequence:
  `data/reference/BIA_mtDNA_consensus_full_v1.fasta` (18,088 bp).

The reconstructed sequence will be cited by its DDBJ/ENA/GenBank TPA accession
after deposition.

## Citation and permanent archive

The release citation and Zenodo DOI will be added after the public GitHub
release is archived. Preliminary citation metadata are provided in
`CITATION.cff` and `.zenodo.json`.

## Licensing

Code is released under the MIT License. Original documentation and derived
tables are released under CC BY 4.0. Candidate-read sequence subsets retain
their NCBI provenance and are not relicensed; see `DATA_LICENSE.md`.
