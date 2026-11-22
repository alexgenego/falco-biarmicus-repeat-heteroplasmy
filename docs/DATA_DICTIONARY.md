# Data dictionary

## SRR34169110 strict_read_classification.tsv

| Field | Meaning |
| --- | --- |
| `read_id` | NCBI SRA read identifier. |
| `read_length` | Full candidate-read length in base pairs. |
| `mtDNA_signature_hits` | Diagnostic positions carrying the reconstructed mtDNA allele. |
| `young_NUMT_signature_hits` | Diagnostic positions carrying the shared young-NUMT allele. |
| `other_diagnostic_calls` | Called diagnostic positions carrying neither expected allele. |
| `classification` | Conservative sequence class used by the workflow. |
| `classification_basis` | Explicit rule or source table supporting that class. |
| `R27`, `R56`, `R80` | Called repeat-copy number; blank means no accepted call. |
| `arrays_called` | Number of accepted repeat arrays in the read. |
| `total_anchor_mismatches` | Sum of mismatches across anchors of accepted arrays. |
| `<array>_residual_bp` | Difference between observed span and integer repeat count × nominal unit length. |
| `<array>_anchor_mm` | Mismatches across the two accepted anchors of that array. |
| `<array>_quantitative_subset` | Whether the call enters the manuscript denominator for that array. |
| `*_signature_positions` | Comma-separated mitochondrial-reference coordinates supporting the call. |

## SRR25731498 mtDNA/NUMT tables

| File | Contents |
| --- | --- |
| `SRR25731498_mtDNA_NUMT_diagnostic_panel.tsv` | 77 positions separating the M and N sequence backgrounds. |
| `SRR25731498_mtDNA_NUMT_read_classification.tsv` | Background class, panel support, and repeat states for all 51 candidates. |
| `SRR25731498_mtDNA_NUMT_distance_to_known_sequences.tsv` | Consensus distances to reconstructed mtDNA and the three young NUMTs. |
| `SRR25731498_mtDNA_NUMT_known_NUMT_sites.tsv` | Sites shared by all three young NUMTs relative to reconstructed mtDNA. |
| `SRR25731498_NUMT_nuclear_flank_linkage.tsv` | Direct read-level linkage to NUMT-specific nuclear-flank 21-mers. |

`sequence_background` is inferred from the 77-site panel. `final_classification`
adds direct nuclear-flank evidence. Repeat states are used for mitochondrial
heteroplasmy only when `final_classification` is `mtDNA-like`.

## Repeat-call tables

`repeat_calls.tsv` contains the raw anchor-based output for every candidate
read. A state is accepted only when both array boundaries are found in the same
physical read and the intervening span is compatible with an integer repeat
count within the specified residual tolerance.

The SRR25731498 discovery set uses one exact 10-mer seed per anchor before
approximate extension because the young NUMTs are more divergent at some
anchors. Biological interpretation still depends on the 77-site background and
direct nuclear-flank classification, not on anchor matching alone.
