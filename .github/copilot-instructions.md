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
  profile (PSSM) alignment, AVX2 support, and built-in Gumbel statistics
  (ALP library).
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
  16-bit word mode.
- **Block-aligner was removed** from this project. The upstream MMseqs2
  engine includes a block-aligner (Rust WASM SIMD) fallback for scores
  that overflow 16-bit word mode. In ssw-aligner this was a stub that
  always returned failure, causing a fallback to standard SW traceback
  on every word-mode alignment — adding overhead and warning noise.
  The block-aligner code, `struct s_block`, and the
  `alignStartPosBacktraceBlock` method have all been deleted. The engine
  now always uses `alignStartPosBacktrace` directly for traceback.
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
│   ├── simd/                 # SIMD abstraction header
│   ├── simde/                # SIMDe — portable SSE2/AVX intrinsics
│   └── fmt/                  # {fmt} formatting library
│
├── scripts/
│   └── build_ext.py          # Build script invoked by poetry-core (CMake)
│
├── tests/
│   ├── test_regression.py    # 41 regression tests vs scikit-bio (pytest)
│   ├── test_mmseqs_freqs.py  # 11 unit tests for MMseqs2 AA freqs & Gumbel
│   ├── test_install.py       # 14 integration tests (wheel + sdist install)
│   ├── test_performance.py   # Throughput benchmark vs scikit-bio
│   └── smoke_test.py         # Standalone smoke test (copied into clean venvs)
│
├── examples/
│   ├── basic_usage.py        # Sequence-sequence alignment examples
│   └── profile_alignment.py  # Profile (PSSM) alignment examples
│
├── CMakeLists.txt            # Build system for libssw_aligner.so
├── pyproject.toml            # Python packaging (poetry-core backend)
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

**Important:** After rebuilding, copy the `.so` to the package directory
so the editable install picks it up:

```bash
cp build/libssw_aligner.so ssw_aligner/libssw_aligner.so
```

The wrapper searches for the library in this order:
1. `ssw_aligner/libssw_aligner.so` (package directory)
2. `build/libssw_aligner.so` (build directory)

### Building with Poetry

```bash
poetry build          # produces wheel + sdist in dist/
pip install dist/ssw_aligner-*.whl
```

The build system uses [poetry-core](https://python-poetry.org/) as the
PEP 517 backend. A build script (`scripts/build_ext.py`) runs CMake
automatically during `poetry build` / `pip install`.

---

## Testing

```bash
pip install -e ".[dev]"    # installs pytest + scikit-bio==0.6.2
pytest tests/ -q
```

### Test suites

| File | Tests | Purpose |
|---|---|---|
| `test_regression.py` | 41 | Compare ssw-aligner vs scikit-bio for nucleotide/protein alignments with various gap penalties, masking, edge cases |
| `test_mmseqs_freqs.py` | 11 | Unit tests for MMseqs2 AA background frequencies, Gumbel parameter computation, bit-scores, E-values |
| `test_install.py` | 14 | Integration tests: build wheel & sdist in fresh venvs, run `smoke_test.py` |
| `test_performance.py` | 3 | Throughput benchmark vs scikit-bio (requires `riot_na` data in sibling directory) |

### Running the performance benchmark

```bash
pytest tests/test_performance.py -v -s
```

Requires the `riot_na` repository cloned as a sibling directory (for NGS
query sequences and V-gene target databases). The benchmark aligns 1 000
amino-acid queries × 327 human V-gene references using BLOSUM62/11/1
parameters.

---

## Benchmark results

Protein alignment throughput (1 000 queries × 327 V-genes, BLOSUM62/11/1):

| Engine | Alignments/s |
|---|---|
| **ssw-aligner** | ~27 000 |
| scikit-bio 0.6.2 | ~39 000 |

Score agreement: **93.1 %** (30 439 / 32 700 pairs match exactly).

The 6.9 % score divergence occurs on alignments where byte-mode scores
overflow (>255) and the word-mode fallback engages. The two engines use
different traceback implementations for word-mode, producing different
start positions and CIGARs. The optimal alignment score itself is always
correct; only traceback details differ.

### Performance gap analysis

ssw-aligner is currently ~1.4× slower than scikit-bio. The inner SIMD
DP kernel is structurally identical between the two engines — the
slowdown comes entirely from surrounding overhead.

**Estimated breakdown of the ~1.4× slowdown:**
- ~40–50 % from **ctypes / Python overhead**
- ~30–40 % from **init overhead** (surplus allocations, unused profiles)
- ~10–20 % from **per-alignment dead code** (unused backtrace string,
  coverage computation)

#### 1. ctypes vs Cython call overhead (~40–50 %)

scikit-bio's Cython wrapper makes **direct C function calls** through
`cdef extern` — zero marshalling cost. ssw-aligner uses ctypes, which
incurs per-call overhead:

| Operation | scikit-bio (Cython) | ssw-aligner (ctypes) |
|---|---|---|
| Function call dispatch | ~5 ns (direct C) | ~200–500 ns (libffi) |
| Argument marshalling | 0 ns (compile-time cast) | ~50–100 ns per arg (boxing) |
| Sequence encoding | Cython loop | Pure Python `for i, ch in enumerate(…)` loop |
| Cigar decoding | Lazy (deferred to property access) | Eager (decoded in `AlignmentStructure.__init__`) |

The `__call__` path wraps **9 arguments** through `ctypes.c_*()` boxing
on every alignment. `_encode_sequence()` is a pure Python character loop
on every target. `AlignmentStructure.__init__` eagerly decodes the full
CIGAR (scikit-bio defers this to a memoized property).

#### 2. Init overhead: ~30 allocations vs 3 (~30–40 %)

| | scikit-bio | ssw-aligner (MMseqs2) |
|---|---|---|
| Heap allocations in init | 2–3 | ~30+ |
| Profile arrays built | 2 (byte + word) | 7 (byte/word/int32 × fwd, plus linear word/int) |
| Extra structures | — | `SimpleBaseMatrix`, `Sequence`, composition bias arrays, reverse query, reverse matrix |

The MMseqs2 `SmithWaterman` constructor (`StripedSmithWaterman.cpp`
L646–L705) allocates ~26 objects: `simd_data`, 4 DP vectors, 6 profile
arrays (byte/word/int32 × fwd+rev), 2 sequence copies, 2 composition
bias arrays, 2 linear profiles, 2 matrix copies, temp buffers, and
per-column scoring arrays. Most of these are **never used** in basic
sequence-sequence alignment (int32 profiles, linear profiles, composition
bias, reverse profiles are all dead weight).

The `createQueryProfile` loop also includes an extra memory load + add
for the zeroed composition bias array on every profile element, even when
bias correction is disabled.

#### 3. Per-alignment dead code (~10–20 %)

Code that runs on **every** alignment but produces unused results:

- **`computerBacktrace()`** — iterates the entire CIGAR building a
  `std::string` by appending "M"/"I"/"D" one character at a time,
  involving heap allocation and copying. The string is populated in
  `ssw_align_private` and immediately discarded in `ssw_api.cpp`.
- **`computeCov()` × 2** — computes query and target coverage ratios
  (always called, though `covThr=0` makes the subsequent
  `hasCoverage()` check trivially pass).
- **E-value null-pointer check** — `evaluer` is always `nullptr` in our
  API, but the branch exists.

### Potential optimizations (not yet implemented)

| Optimization | Expected impact | Effort |
|---|---|---|
| AVX2 build (`-DHAVE_AVX2=ON`) | ~2× SIMD throughput | Low (cmake flag) |
| Batch alignment API | Amortize ctypes overhead across N targets | Medium |
| Strip `computerBacktrace()` call | Eliminate per-alignment string alloc | Low |
| Strip `computeCov()` + `hasCoverage()` | Remove per-alignment dead code | Low |
| Remove int32 + linear profile building | Faster init, less memory | Low |
| Remove composition bias array loads | Fewer memory accesses in profile build | Low |
| Cython wrapper (replace ctypes) | Eliminate ~40–50 % of gap | High |

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

---

## Known issues

1. **Stale `.so` after rebuild** — the Python wrapper loads
   `ssw_aligner/libssw_aligner.so` first. After running `cmake --build`,
   the new library lands in `build/`. You must copy it manually:
   ```bash
   cp build/libssw_aligner.so ssw_aligner/libssw_aligner.so
   ```
   Forgetting this step means the old binary is still loaded at runtime.

2. **~1.4× slower than scikit-bio** — see [Benchmark results](#benchmark-results)
   for root causes and planned optimizations.

3. **6.9 % score divergence on word-mode alignments** — when byte-mode
   scores overflow (>255), the engine falls back to 16-bit word mode.
   The traceback implementation differs from scikit-bio's Mengyao Zhao
   SSW library, producing different start positions and CIGARs on those
   alignments. The optimal score itself is always correct.

4. **No Windows or macOS CI** — the build and tests are only validated on
   Linux. CMake and SIMDe should support other platforms, but this is
   untested.
