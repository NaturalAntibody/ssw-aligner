# Copilot instructions for ssw-aligner

## Project overview

`ssw-aligner` is a standalone SIMD-accelerated Smith-Waterman local aligner
for Python. The alignment engine is extracted from **MMseqs2** (Söding Lab,
GPL-3.0) and wrapped in a pure-Python ctypes API that is a drop-in replacement
for `skbio.alignment.StripedSmithWaterman`.

### Why this project exists

[scikit-bio](https://github.com/scikit-bio/scikit-bio) bundles a Cython wrapper
around the original Mengyao Zhao SSW C library (Boston College, 2010). However:

- scikit-bio may deprecate its SSW module (see
  [scikit-bio#1814](https://github.com/scikit-bio/scikit-bio/issues/1814)).
- The original C library lacks features present in MMseqs2's engine:
  profile (PSSM) alignment, block-aligner fallback for large scores,
  AVX2 support, and built-in Gumbel statistics (ALP library).
- scikit-bio's SSW is compiled via Cython and tightly coupled to NumPy ABI,
  making cross-platform builds harder.

`ssw-aligner` replaces all of that with a single shared library
(`libssw_aligner.so`) built via CMake and a pure-ctypes Python wrapper
with no Cython dependency.

### Relationship with riot_na

[riot_na](https://github.com/NaturalAntibody/riot_na) (Rapid Identification
Of anTibody sequences) is an antibody numbering tool that aligns nucleotide
and amino-acid immunoglobulin sequences against germline gene databases using
Smith-Waterman, then applies numbering schemes (IMGT, Kabat, Chothia, etc.).

- **riot_na currently depends on scikit-bio==0.6.2** solely for
  `StripedSmithWaterman`. The pin is explicit: *"pinned because of the
  danger of deprecating SSW"*.
- **ssw-aligner is the intended replacement**: swapping
  `from skbio.alignment import StripedSmithWaterman` →
  `from ssw_aligner import StripedSmithWaterman` requires no other code
  changes thanks to API compatibility.
- riot_na also has a Rust component (`riot_prefiltering`) compiled via
  maturin that performs fast k-mer prefiltering to narrow the candidate
  gene list *before* handing sequences to SSW for full alignment.
- Key alignment parameters used by riot_na: `match_score=1`,
  `mismatch_score=-1`, `gap_open_penalty=4`, `gap_extend_penalty=1`
  (nucleotide mode).

### Relationship with MMseqs2

The C++ alignment engine under `src/` is extracted from
[MMseqs2](https://github.com/soedinglab/MMseqs2) (`src/alignment/` and
`src/commons/`). Key facts for developers:

- The core DP kernel lives in `StripedSmithWaterman.cpp` with templates
  `sw_sse2_byte<type>` and `sw_sse2_word<type>`, where `type` is
  `SEQ_SEQ` or `PROFILE_SEQ`.
- Profile mode is activated when `Sequence::getSequenceType()` returns
  `DBTYPE_HMM_PROFILE` (value 2).
- The engine uses SSE2 by default; AVX2 is enabled via
  `-DHAVE_AVX2=ON` in CMake, which defines the `AVX2` macro and switches
  SIMD lane widths from 128-bit to 256-bit.
- When byte-width scores overflow (>255), the engine falls back to
  16-bit word mode. When word-mode scores overflow, it falls back to
  the block-aligner (a Rust WASM-compiled SIMD aligner under
  `lib/block-aligner/c/`).
- Gumbel statistics (λ, K for E-values) are computed via the ALP library
  under `lib/alp/`.
- MMseqs2 uses `alphabetSize = 21` for amino acids (20 AAs + X sentinel).
  The C API adapter (`ssw_api.cpp`) takes `n = 20` from Python and
  internally uses `n + 1` for profile mode so the engine's X-row zeroing
  targets the padding row, not the last real amino acid.

### Relationship with scikit-bio

- The Python API mirrors `skbio.alignment.StripedSmithWaterman` and
  `skbio.alignment.AlignmentStructure` (skbio 0.6.2).
- The test suite (`tests/test_regression.py`, 41 tests) validates
  compatibility by comparing results against scikit-bio's implementation
  (scikit-bio is a dev dependency).
- Notable API differences:
  - `protein=True` parameter for protein alignment (skbio auto-detected
    from matrix keys).
  - `SmithWatermanProfileAligner` — new, not present in scikit-bio.
  - `compute_gumbel_params()` — new, not present in scikit-bio.

---

## Repository structure

```
ssw-aligner/
├── ssw_aligner/              # Python package
│   ├── __init__.py           # Public exports
│   └── _wrapper.py           # Pure ctypes wrapper (main Python API)
│
├── src/                      # C++ source (extracted from MMseqs2)
│   ├── api/
│   │   ├── ssw_api.h         # extern "C" API header
│   │   └── ssw_api.cpp       # C adapter bridging to C++ engine
│   ├── alignment/
│   │   ├── StripedSmithWaterman.cpp/.h   # Core SIMD SW engine
│   │   ├── EvalueComputation.h           # E-value computation
│   │   └── PSSMCalculator.h              # PSSM utilities
│   └── commons/
│       ├── Sequence.cpp/.h               # Sequence encoding, PROFILE_AA_SIZE=20
│       ├── SubstitutionMatrix.cpp/.h     # Matrix loading/computation
│       ├── BaseMatrix.cpp/.h             # Base class for scoring matrices
│       ├── Parameters.cpp/.h             # DBTYPE constants, config
│       └── ...                           # Debug, Util, LambdaCalculation, etc.
│
├── lib/                      # Third-party libraries
│   ├── alp/                  # ALP — Gumbel parameter estimation (public domain)
│   ├── block-aligner/c/      # Block-aligner — fallback for large scores
│   ├── simd/                 # SIMD abstraction header
│   ├── simde/                # SIMDe — portable SSE2/AVX intrinsics
│   └── fmt/                  # {fmt} formatting library
│
├── tests/
│   └── test_regression.py    # 41 regression tests (pytest)
│
├── examples/
│   ├── basic_usage.py        # Sequence-sequence alignment examples
│   └── profile_alignment.py  # Profile (PSSM) alignment examples
│
├── CMakeLists.txt            # Build system for libssw_aligner.so
├── pyproject.toml            # Python packaging
└── README.md
```

---

## Architecture

### C layer

`src/api/ssw_api.cpp` provides five `extern "C"` functions:

| Function | Purpose |
|---|---|
| `ssw_init` | Build a query profile from a flat n×n substitution matrix |
| `ssw_init_profile` | Build a query profile from an n×queryLen PSSM |
| `ssw_align` | Align a target against the stored query profile |
| `ssw_free_cigar` | Free a cigar array from `ssw_result` |
| `ssw_destroy` | Destroy a handle and free all resources |
| `compute_gumbel_params` | Estimate Gumbel λ/K via ALP Monte Carlo |

The adapter uses a `SimpleBaseMatrix` (subclass of `BaseMatrix`) to wrap
flat int8_t scoring matrices, and an `ssw_handle` struct holding
`SmithWaterman*`, `SimpleBaseMatrix*`, `Sequence*`, the flat matrix copy,
alphabet size, and query length.

**Important internal detail:** For profile mode, `ssw_init_profile` uses
`internalAlphSize = n + 1` (21 for amino acids) because the MMseqs2 engine
zeroes row `alphabetSize - 1` as the neutral 'X' state. With `n = 20`, this
would corrupt the last real amino acid (Valine). The extra row is allocated
and zeroed explicitly.

### Python layer

`ssw_aligner/_wrapper.py` is a single-file pure-ctypes wrapper:

- **`StripedSmithWaterman`** — sequence-sequence alignment (nucleotide or
  protein). Encodes the query, builds the flat substitution matrix, calls
  `ssw_init`, and provides `__call__` for alignment.
- **`SmithWatermanProfileAligner`** — profile alignment. Takes a numpy
  PSSM of shape `(20, query_length)`, calls `ssw_init_profile`.
- **`AlignmentStructure`** — result container with score, positions,
  CIGAR, and aligned sequences. Supports both attribute and dict access.
- **`GumbelParams`** — holds λ and K with `bit_score()` / `evalue()`.
- **`compute_gumbel_params()`** — calls the C function, with precomputed
  parameter lookup for common configurations (BLOSUM62/11/1,
  nucleotide/7/1, nucleotide/5/2).

### Alphabet conventions

- **Amino acids:** 20-letter alphabet `ARNDCQEGHILKMFPSTWYV` (indices 0–19).
  This is the order used in scoring matrices, PSSMs, and `_AA_TABLE`.
- **Nucleotides:** 5-letter alphabet `ACGTN` (indices 0–4).
- Unknown characters map to the last index (19 for AA, 4 for NT).

---

## Build instructions

```bash
# Prerequisites: CMake ≥ 3.14, C++17 compiler
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)

# AVX2 (faster on supported CPUs):
cmake -B build -DCMAKE_BUILD_TYPE=Release -DHAVE_AVX2=ON
cmake --build build -j$(nproc)
```

Output: `build/libssw_aligner.so` — loaded automatically by the Python
wrapper from the `build/` directory relative to the package.

---

## Testing

```bash
pip install -e ".[dev]"    # installs pytest + scikit-bio==0.6.2
pytest tests/ -q
```

The 41 regression tests compare ssw-aligner results against scikit-bio's
`StripedSmithWaterman` for nucleotide and protein alignments with various
gap penalties, masking options, and edge cases.

---

## Development guidelines

### When modifying the C++ engine

- All changes should go through `src/api/ssw_api.cpp` — avoid modifying
  the MMseqs2 engine files directly unless fixing a bug upstream.
- If you add a new C API function, declare it in `ssw_api.h` and add
  the corresponding ctypes declaration at the top of `_wrapper.py`.
- After any C++ change: `cmake --build build -j$(nproc)` then
  `pytest tests/ -q`.

### When modifying the Python wrapper

- Keep API compatibility with `skbio.alignment.StripedSmithWaterman`.
- All public names must be exported in `__init__.py`.
- Add tests to `tests/test_regression.py` for new features.

### Precomputed Gumbel parameters

`_PRECOMPUTED_PARAMS` in `_wrapper.py` stores Gumbel parameters for
common matrix/gap configurations, avoiding the slow ALP Monte Carlo.
When adding a new precomputed config:
1. Run ALP via `compute_gumbel_params(..., recalculate_gumbel_params=True)`.
2. Add the 12 doubles to `_PRECOMPUTED_PARAMS` keyed by
   `(matrix_name, gap_open, gap_extend)`.
3. If the matrix is not auto-detected, add fingerprint logic to
   `_detect_matrix_name()`.

### Key constants

- `PROFILE_AA_SIZE = 20` (in `Sequence.h`)
- `DBTYPE_HMM_PROFILE = 2` (in `Parameters.h`) — triggers profile alignment
- `DBTYPE_AMINO_ACIDS = 0` (in `Parameters.h`)
- Internal amino acid alphabet size for the engine: 21 (20 + X sentinel)
