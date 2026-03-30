# ssw-aligner

SIMD-accelerated [Smith-Waterman](https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm) local aligner for Python.

The alignment engine is extracted from [MMseqs2](https://github.com/soedinglab/MMseqs2)
and exposes a pure-Python (ctypes) API that is a drop-in replacement for
`skbio.alignment.StripedSmithWaterman`.

## Features

- **Fast** — SSE2 / AVX2 SIMD striped Smith-Waterman (byte + word modes).
- **Nucleotide & protein** — built-in match/mismatch scoring for DNA;
  arbitrary substitution matrices (e.g. BLOSUM62) for proteins.
- **Profile (PSSM) alignment** — align targets against a position-specific
  scoring matrix via `SmithWatermanProfileAligner`.
- **Gumbel statistics** — compute λ and K parameters via the
  [ALP](https://doi.org/10.1093/bioinformatics/btv575) Monte Carlo method,
  then derive bit-scores and E-values.
- **No compiled Python extension** — the shared library is plain C/C++;
  the wrapper uses `ctypes`, so there is nothing to Cython-compile or
  link against NumPy ABI.
- **Drop-in replacement** for `skbio.alignment.StripedSmithWaterman`.

## Installation

### From wheel (recommended)

```bash
pip install ssw-aligner
```

Pre-built wheels include the compiled `libssw_aligner.so` — no C++ toolchain
required.

### From source (sdist)

When installing from source, CMake ≥ 3.14 and a C++17 compiler are required.
CMake is declared as a build dependency and will be fetched automatically
by pip:

```bash
pip install ssw-aligner --no-binary ssw-aligner
```

### Development setup

```bash
git clone https://github.com/NaturalAntibody/ssw-aligner.git
cd ssw-aligner

# Poetry (recommended)
poetry install
poetry build          # produces wheel + sdist in dist/

# Or plain pip (editable)
pip install -e ".[dev]"
```

The build system uses [poetry-core](https://python-poetry.org/) as the
PEP 517 backend.  A build script (`scripts/build_ext.py`) compiles the
C++ library via CMake automatically during `poetry build` / `pip install`.

To enable AVX2 for ~2× throughput on supported CPUs, set the environment
variable before building:

```bash
CMAKE_ARGS="-DHAVE_AVX2=ON" poetry build
```

## Quick start

```python
from ssw_aligner import StripedSmithWaterman

# --- Nucleotide alignment ---
aligner = StripedSmithWaterman("ACTGACTG")
result = aligner("AACTGACTGA")
print(result.optimal_alignment_score)   # 8
print(result.cigar)                     # "8M"

# --- Protein alignment ---
import blosum
aligner = StripedSmithWaterman(
    "EVQLVESGGGLVQPGG",
    protein=True,
    substitution_matrix=blosum.BLOSUM(62),
    gap_open_penalty=11,
    gap_extend_penalty=1,
)
result = aligner("EVQLVESGGGLVQPGG")
print(result.optimal_alignment_score)   # 89
```

See [examples/basic_usage.py](examples/basic_usage.py) for a more complete
walkthrough.

## Alignment modes

### Sequence–sequence alignment

The standard mode aligns a query string against target strings using a
**global substitution matrix** — the same score table applies at every query
position.

For **nucleotides** the matrix is built automatically from `match_score` /
`mismatch_score`.  For **proteins** you supply a substitution matrix
(e.g. BLOSUM62) as a nested `dict[str, dict[str, int]]`.

```python
from ssw_aligner import StripedSmithWaterman

# Nucleotide — matrix built from match/mismatch
aligner = StripedSmithWaterman("ACTGACTG")
result = aligner("AACTGACTGA")

# Protein — explicit substitution matrix
import blosum
aligner = StripedSmithWaterman(
    "EVQLVESGGGLVQPGG",
    protein=True,
    substitution_matrix=blosum.BLOSUM(62),
    gap_open_penalty=11, gap_extend_penalty=1,
)
```

Internally the query is encoded to numeric indices (0–4 for `ACGTN`,
0–19 for `ARNDCQEGHILKMFPSTWYV`) and a SIMD-striped query profile is
built once via the C function `ssw_init`.  Each `aligner(target)` call
invokes `ssw_align` which runs the Smith-Waterman DP in byte (uint8) or
word (uint16) SIMD lanes.

### Profile (PSSM) alignment

`SmithWatermanProfileAligner` replaces the fixed substitution matrix with
a **position-specific scoring matrix** (PSSM) of shape `(20, query_length)`.
Every query position can have a different score vector, enabling:

- **HMM profile → sequence** alignment
- **MSA-derived PSSMs** (e.g. PSI-BLAST / HHblits profiles)
- **Position-specific weighting** (e.g. stronger penalties in CDR regions
  of an antibody sequence)

```python
import numpy as np
from ssw_aligner import SmithWatermanProfileAligner

query = "EVQLVESGGGLVQPGG"
pssm = np.zeros((20, len(query)), dtype=np.int8)
# … fill pssm[aa_idx, position] with per-position scores …

aligner = SmithWatermanProfileAligner(
    query, pssm,
    gap_open_penalty=11, gap_extend_penalty=1,
)
result = aligner("EVQLVESGGGLVKPGG")
```

`pssm[i, j]` is the score for amino acid `i` (in the order
`ARNDCQEGHILKMFPSTWYV`) at query position `j`.  The PSSM is passed to
the C function `ssw_init_profile`, which sets MMseqs2's sequence type to
`DBTYPE_HMM_PROFILE` and activates the `PROFILE_SEQ` SIMD template.

When the PSSM is derived from a flat substitution matrix (e.g. every column
is just the BLOSUM62 row for the query residue at that position), the profile
aligner produces **identical scores** to `StripedSmithWaterman` with the same
matrix and gap penalties.

See [examples/profile_alignment.py](examples/profile_alignment.py) for
a full walkthrough including MSA-derived profiles and position-weighted
scoring.

## API reference

### `StripedSmithWaterman`

```python
StripedSmithWaterman(
    query_sequence,
    *,
    gap_open_penalty=5,
    gap_extend_penalty=2,
    protein=False,
    match_score=2,              # nucleotide only
    mismatch_score=-3,          # nucleotide only
    substitution_matrix=None,   # required for protein
    score_only=False,
    mask_length=15,
    mask_auto=True,
    suppress_sequences=False,
    zero_index=True,
)
```

Create an aligner for `query_sequence`, then call it on target strings:

```python
result = aligner(target_sequence)   # returns AlignmentStructure
```

### `AlignmentStructure`

Returned by calling a `StripedSmithWaterman` instance.  Supports both
attribute and dict-style access.

| Property | Description |
|---|---|
| `optimal_alignment_score` | Best alignment score |
| `suboptimal_alignment_score` | Second-best score |
| `query_begin` | 0-based start in query |
| `query_end` | 0-based end in query |
| `target_begin` | 0-based start in target |
| `target_end_optimal` | 0-based end in target (best) |
| `target_end_suboptimal` | 0-based end in target (second-best) |
| `cigar` | CIGAR string |
| `aligned_query_sequence` | Query with gaps inserted |
| `aligned_target_sequence` | Target with gaps inserted |

### `SmithWatermanProfileAligner`

```python
SmithWatermanProfileAligner(
    query_sequence,
    pssm,                       # numpy.ndarray, shape (20, query_length), dtype int8
    *,
    gap_open_penalty=11,
    gap_extend_penalty=1,
    score_only=False,
    mask_length=15,
    mask_auto=True,
    zero_index=True,
)
```

Create a profile aligner, then call it on target strings:

```python
result = aligner(target_sequence)   # returns AlignmentStructure
```

The PSSM rows correspond to amino acids in the order
`ARNDCQEGHILKMFPSTWYV` (indices 0–19); columns correspond to query
positions.

### `compute_gumbel_params`

```python
from ssw_aligner import compute_gumbel_params

params = compute_gumbel_params(
    substitution_matrix,        # dict[str, dict[str, int]]
    gap_open,                   # int
    gap_extend,                 # int
    protein=True,
    use_mmseqs_aa_freqs=False,  # use Robinson-Robinson background freqs
    bg_freqs=None,              # custom background freqs dict
    max_seconds=60.0,
)

params.lambda_              # Gumbel λ
params.K                    # Gumbel K
params.bit_score(raw_score)                          # → float
params.evalue(raw_score, query_length, db_residues)  # → float
```

Estimates Gumbel extreme-value distribution parameters via ALP Monte Carlo
simulation.  Set `use_mmseqs_aa_freqs=True` to use the Robinson-Robinson
amino acid background frequencies from MMseqs2's built-in BLOSUM62.

### `MMSEQS_AA_FREQS`

```python
from ssw_aligner import MMSEQS_AA_FREQS
```

`dict[str, float]` — Robinson-Robinson background frequencies for the 20
standard amino acids, as used by MMseqs2.

## Benchmark

Protein alignment throughput measured on 1 000 amino-acid query sequences
(NGS dataset) × 327 human V-gene germline references, BLOSUM62 / gap-open 11 /
gap-extend 1 — the same parameters used by
[riot_na](https://github.com/NaturalAntibody/riot_na):

| Engine | Alignments/s | Total (327 000 aln) |
|---|--:|--:|
| **ssw-aligner** (this package) | ~27 000 | ~12 s |
| scikit-bio 0.6.2 SSW | ~39 000 | ~8 s |

Score agreement between the two engines is **93.1 %** on this workload.
The 6.9 % divergence comes from differing traceback heuristics when
byte-mode scores overflow and the word-mode fallback is engaged — the
optimal alignment scores are identical; only start-position / CIGAR
details can differ.

Reproduce with:

```bash
pytest tests/test_performance.py -v -s
```

## Migrating from scikit-bio

Replace:
```python
from skbio.alignment import StripedSmithWaterman
```
with:
```python
from ssw_aligner import StripedSmithWaterman
```

The API is compatible — no other code changes are needed.

## License

- MMseqs2 alignment engine: [GPL-3.0](https://github.com/soedinglab/MMseqs2/blob/master/LICENSE.md) (Copyright 2016 Söding Lab)
- ALP library: Public domain (Spouge/Sheetlin, NCBI)
- Python wrapper & build glue: [BSD-3-Clause](LICENSE.txt) (Copyright 2013-- scikit-bio development team)
- SIMDe header: MIT License
