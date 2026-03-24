# ssw-aligner

Standalone [StripedSmithWaterman (SSW)](http://www.plosone.org/article/info:doi/10.1371/journal.pone.0082138) aligner extracted from [scikit-bio](https://scikit.bio).

This micro-library provides the SSW C implementation and its Python/Cython wrapper as a drop-in replacement for `skbio.alignment.StripedSmithWaterman`.

## Installation

```bash
pip install .
```

## Usage

```python
from ssw_aligner import StripedSmithWaterman

# Nucleotide alignment
aligner = StripedSmithWaterman("ACTGACTG")
result = aligner("AACTGACTGA")
print(result["optimal_alignment_score"])
print(result["cigar"])

# Protein alignment (requires a substitution matrix)
import blosum
aligner = StripedSmithWaterman(
    "ARNDCQEGHILKMFPSTWYV",
    protein=True,
    substitution_matrix=blosum.BLOSUM(62),
    gap_open_penalty=11,
    gap_extend_penalty=1,
)
result = aligner("ARNDCQEGHILKMFPSTWYV")
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
