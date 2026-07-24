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

## Usage

```python
from ssw_aligner import BLOSUM_62, StripedSmithWaterman

# Nucleotide alignment
aligner = StripedSmithWaterman(query_sequence="ACTGACTG")
result = aligner(target_sequence="AACTGACTGA")
print(result["optimal_alignment_score"])
print(result["cigar"])

# Protein alignment
aligner = StripedSmithWaterman(
    query_sequence="ARNDCQEGHILKMFPSTWYV",
    protein=True,
    substitution_matrix=BLOSUM_62,
    gap_open_penalty=11,
    gap_extend_penalty=1,
)
result = aligner(target_sequence="ARNDCQEGHILKMFPSTWYV")
print(result["optimal_alignment_score"])
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

This package is licensed as `BSD-3-Clause AND MIT`. Full text is in
[`LICENSE.txt`](LICENSE.txt):

- SSW C library (`ssw.c` / `ssw.h`): MIT License (Copyright 2012-2015 Boston College)
- Cython wrapper: Modified BSD License (Copyright 2013-- scikit-bio development team)
- SIMDe header: MIT License
