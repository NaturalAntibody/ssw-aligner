# ssw-aligner

Standalone [StripedSmithWaterman (SSW)](http://www.plosone.org/article/info:doi/10.1371/journal.pone.0082138) aligner extracted from [scikit-bio](https://scikit.bio).

This micro-library provides the SSW C implementation and its Python/Cython wrapper as a drop-in replacement for `skbio.alignment.StripedSmithWaterman`.

`StripedSmithWaterman` keeps the original scikit-bio defaults for backward compatibility.
For RIOT-compatible defaults, use `NucleotideAligner` or `ProteinAligner`.

## Installation

```bash
pip install .
```

## RIOT-Compatible Wrappers

These wrappers use the same defaults as RIOT while leaving
`StripedSmithWaterman` unchanged for scikit-bio compatibility.

`NucleotideAligner` uses RIOT's `ALIGNER_PARAMS` defaults:
`match_score=1`, `mismatch_score=-1`, `gap_open_penalty=4`, `gap_extend_penalty=1`.

`ProteinAligner` uses RIOT's `AA_ALIGNER_PARAMS` defaults:
`gap_open_penalty=11`, `gap_extend_penalty=1`, `protein=True`, and `substitution_matrix=BLOSUM_62`.

```python
from ssw_aligner import NucleotideAligner, ProteinAligner

nt_aligner = NucleotideAligner(query_sequence="ACTGACTG")
nt_result = nt_aligner(target_sequence="AACTGACTGA")

protein_aligner = ProteinAligner(query_sequence="ARNDCQEGHILKMFPSTWYV")
protein_result = protein_aligner(target_sequence="ARNDAQEGHILKMFASTWYV")

print(nt_result["optimal_alignment_score"])
print(protein_result["optimal_alignment_score"])
```

Full flow with RIOT-compatible defaults and post-alignment metrics:

```python
from ssw_aligner import (
    NucleotideAligner,
    ProteinAligner,
    calculate_seq_identity,
    compute_bit_score,
    compute_evalue,
    compute_raw_score_aa,
)

database_length = 2_500_000  # Sum of lengths of all sequences in the search database.

nt_aligner = NucleotideAligner(query_sequence="ACGTACGT")
nt_result = nt_aligner(target_sequence="AACGTAGT")

nt_sequence_identity = calculate_seq_identity(
    cigar=nt_result["cigar"],
    query=nt_result["query_sequence"],
    target=nt_result["target_sequence"],
    query_start=nt_result["query_begin"],
    target_start=nt_result["target_begin"],
)
nt_bit_score = compute_bit_score(raw_score=nt_result["optimal_alignment_score"])
nt_evalue = compute_evalue(
    query_length=len(nt_result["query_sequence"]),
    db_length=database_length,
    bit_score=nt_bit_score,
)

protein_aligner = ProteinAligner(query_sequence="ARNDCQEGHILKMFPSTWYV")
protein_result = protein_aligner(target_sequence="ARNDAQEGHILKMFASTWYV")

protein_sequence_identity = calculate_seq_identity(
    cigar=protein_result["cigar"],
    query=protein_result["query_sequence"],
    target=protein_result["target_sequence"],
    query_start=protein_result["query_begin"],
    target_start=protein_result["target_begin"],
)
protein_raw_score = compute_raw_score_aa(
    query=protein_result["query_sequence"][protein_result["query_begin"] : protein_result["query_end"] + 1],
    target=protein_result["target_sequence"][protein_result["target_begin"] : protein_result["target_end_optimal"] + 1],
    cigar=protein_result["cigar"],
)
protein_bit_score = compute_bit_score(raw_score=protein_raw_score)
protein_evalue = compute_evalue(
    query_length=len(protein_result["query_sequence"]),
    db_length=database_length,
    bit_score=protein_bit_score,
)

print(nt_sequence_identity)
print(nt_bit_score)
print(nt_evalue)
print(protein_sequence_identity)
print(protein_raw_score)
print(protein_bit_score)
print(protein_evalue)
```

`calculate_seq_identity` reports the fraction of identical positions across the
alignment span encoded by the CIGAR string.

`compute_bit_score` converts a raw alignment score into a database-independent
normalized score.

`compute_evalue` estimates how often an alignment with that score would occur by
chance for a given query length and total database length. `db_length` is the
total number of residues or bases across the entire search database, not the
length of the single matched target sequence.

## Usage

```python
from ssw_aligner import BLOSUM_62, StripedSmithWaterman, calculate_seq_identity, compute_bit_score

# Nucleotide alignment
aligner = StripedSmithWaterman(query_sequence="ACTGACTG")
result = aligner(target_sequence="AACTGACTGA")
print(result["optimal_alignment_score"])
print(result["cigar"])
print(
    calculate_seq_identity(
        cigar=result["cigar"],
        query=result["query_sequence"],
        target=result["target_sequence"],
    )
)

# Protein alignment
aligner = StripedSmithWaterman(
    query_sequence="ARNDCQEGHILKMFPSTWYV",
    protein=True,
    substitution_matrix=BLOSUM_62,
    gap_open_penalty=11,
    gap_extend_penalty=1,
)
result = aligner(target_sequence="ARNDCQEGHILKMFPSTWYV")
print(compute_bit_score(raw_score=result["optimal_alignment_score"]))
```

The package also exposes pure-Python helpers for common post-alignment metrics:
`BLOSUM_62`, `compute_raw_score_aa`, `compute_bit_score`, `compute_evalue`,
and `calculate_seq_identity`.

## Metrics Examples

The helpers can also be used directly with the base `StripedSmithWaterman`
interface. For E-values, pass the total database length rather than the matched
target length.

```python
from ssw_aligner import StripedSmithWaterman, calculate_seq_identity, compute_bit_score, compute_evalue

database_length = 2_500_000

aligner = StripedSmithWaterman(query_sequence="ACGTACGT")
result = aligner(target_sequence="AACGTAGT")

sequence_identity = calculate_seq_identity(
    cigar=result["cigar"],
    query=result["query_sequence"],
    target=result["target_sequence"],
    query_start=result["query_begin"],
    target_start=result["target_begin"],
)
bit_score = compute_bit_score(raw_score=result["optimal_alignment_score"])
evalue = compute_evalue(
    query_length=len(result["query_sequence"]),
    db_length=database_length,
    bit_score=bit_score,
)

print(sequence_identity)
print(bit_score)
print(evalue)
```

```python
from ssw_aligner import (
    BLOSUM_62,
    StripedSmithWaterman,
    calculate_seq_identity,
    compute_bit_score,
    compute_evalue,
    compute_raw_score_aa,
    compute_raw_score_from_bit_score,
)

database_length = 2_500_000

aligner = StripedSmithWaterman(
    query_sequence="ARNDCQEGHILKMFPSTWYV",
    protein=True,
    substitution_matrix=BLOSUM_62,
    gap_open_penalty=11,
    gap_extend_penalty=1,
)
result = aligner(target_sequence="ARNDAQEGHILKMFASTWYV")

sequence_identity = calculate_seq_identity(
    cigar=result["cigar"],
    query=result["query_sequence"],
    target=result["target_sequence"],
    query_start=result["query_begin"],
    target_start=result["target_begin"],
)

raw_score = compute_raw_score_aa(
    query=result["query_sequence"][result["query_begin"] : result["query_end"] + 1],
    target=result["target_sequence"][result["target_begin"] : result["target_end_optimal"] + 1],
    cigar=result["cigar"],
    substitution_matrix=BLOSUM_62,
)
bit_score = compute_bit_score(raw_score=raw_score)
reconstructed_raw_score = compute_raw_score_from_bit_score(bit_score=bit_score)
evalue = compute_evalue(
    query_length=len(result["query_sequence"]),
    db_length=database_length,
    bit_score=bit_score,
)

print(sequence_identity)
print(raw_score)
print(bit_score)
print(reconstructed_raw_score)
print(evalue)
```

## Migrating from scikit-bio

Replace:
```python
from skbio.alignment import StripedSmithWaterman
```

With:
```python
from ssw_aligner import StripedSmithWaterman
```

No other code changes are needed.

## License

- SSW C library: MIT License (Copyright 2012-2015 Boston College)
- Cython wrapper: Modified BSD License (Copyright 2013-- scikit-bio development team)
- SIMDe header: MIT License
