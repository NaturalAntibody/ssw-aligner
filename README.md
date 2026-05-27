# ssw-aligner

Standalone [StripedSmithWaterman (SSW)](http://www.plosone.org/article/info:doi/10.1371/journal.pone.0082138) aligner extracted from [scikit-bio](https://scikit.bio).

This micro-library provides the SSW C implementation and its Python/Cython wrapper as a drop-in replacement for `skbio.alignment.StripedSmithWaterman`.

## Installation

```bash
pip install .
```

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

```python
from ssw_aligner import StripedSmithWaterman, calculate_seq_identity, compute_bit_score, compute_evalue

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
    db_length=len(result["target_sequence"]),
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
    compute_bit_score,
    compute_evalue,
    compute_raw_score_aa,
    compute_raw_score_from_bit_score,
)

aligner = StripedSmithWaterman(
    query_sequence="ARNDCQEGHILKMFPSTWYV",
    protein=True,
    substitution_matrix=BLOSUM_62,
    gap_open_penalty=11,
    gap_extend_penalty=1,
)
result = aligner(target_sequence="ARNDAQEGHILKMFASTWYV")

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
    db_length=len(result["target_sequence"]),
    bit_score=bit_score,
)

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
