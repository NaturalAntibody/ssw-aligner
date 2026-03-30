# MMseqs2 Smith-Waterman Engine Implementation

This document captures the full thinking process, technical decisions, and
code changes made while replacing the old scikit-bio SSW wrapper with an
MMseqs2-based SIMD alignment engine. It is written for someone who
understands what sequence alignment is (comparing biological sequences to
find regions of similarity) but may not be familiar with specific
algorithms, scoring parameters, or SIMD optimisation techniques.

---

## Table of contents

1. [Background: what problem are we solving?](#1-background-what-problem-are-we-solving)
2. [How Smith-Waterman alignment works](#2-how-smith-waterman-alignment-works)
3. [How SIMD acceleration works](#3-how-simd-acceleration-works)
4. [The byte / word / int cascade — why "255" matters](#4-the-byte--word--int-cascade--why-255-matters)
5. [Why we chose MMseqs2 over the original SSW library](#5-why-we-chose-mmseqs2-over-the-original-ssw-library)
6. [Architecture of ssw-aligner](#6-architecture-of-ssw-aligner)
7. [Block-aligner: what it was and why we removed it](#7-block-aligner-what-it-was-and-why-we-removed-it)
8. [Performance analysis: why we are slower than scikit-bio](#8-performance-analysis-why-we-are-slower-than-scikit-bio)
9. [Dead code removal optimisation](#9-dead-code-removal-optimisation)
10. [The 6.9 % score divergence explained](#10-the-69--score-divergence-explained)
11. [Chronological list of code changes](#11-chronological-list-of-code-changes)
12. [Remaining work and future optimisations](#12-remaining-work-and-future-optimisations)

---

## 1. Background: what problem are we solving?

### The ssw-aligner package

`ssw-aligner` is a Python library that performs **local sequence
alignment** — finding the best-matching region between two biological
sequences (DNA or protein). It is used by
[riot_na](https://github.com/NaturalAntibody/riot_na), an antibody
numbering tool that aligns immunoglobulin sequences against germline gene
databases.

### The old implementation

Previously, `ssw-aligner` was a thin wrapper around the C library
written by Mengyao Zhao (Boston College, 2010) — the same library
bundled inside [scikit-bio](https://github.com/scikit-bio/scikit-bio).
The wrapper was written in **Cython** (a language that compiles Python-like
code into C extensions).

### Why we needed to replace it

1. **scikit-bio may deprecate its SSW module**
   ([issue #1814](https://github.com/scikit-bio/scikit-bio/issues/1814)).
   riot_na pins `scikit-bio==0.6.2` explicitly because of this risk.
2. **Missing features** — the original C library lacks profile (PSSM)
   alignment, AVX2 support, and built-in statistical significance
   computation (Gumbel parameters).
3. **Build complexity** — Cython + NumPy ABI coupling makes
   cross-platform wheel builds fragile.

### The new implementation

We extracted the Smith-Waterman engine from
[MMseqs2](https://github.com/soedinglab/MMseqs2) (a state-of-the-art
sequence search tool from the Söding Lab) and wrapped it with a
pure-Python **ctypes** interface. The result is a single shared library
(`libssw_aligner.so`) built with CMake, plus a Python wrapper that is
a drop-in replacement for `skbio.alignment.StripedSmithWaterman`.

---

## 2. How Smith-Waterman alignment works

Smith-Waterman (SW) is a **dynamic programming** algorithm for finding
the best local alignment between two sequences. "Local" means it finds
the highest-scoring subsequence match — unlike global alignment
(Needleman-Wunsch), it can ignore poor-quality regions at the edges.

### The scoring model

The algorithm fills a 2D matrix where each cell `H[i][j]` represents the
best alignment score ending at position `i` in the query and position `j`
in the target. The recurrence is:

$$
H[i][j] = \max \begin{cases}
H[i-1][j-1] + S(q_i, t_j) & \text{(match/mismatch)} \\
E[i][j] & \text{(gap in query)} \\
F[i][j] & \text{(gap in target)} \\
0 & \text{(restart: local alignment)}
\end{cases}
$$

Where:
- $S(q_i, t_j)$ is the **substitution score** — how similar residue
  $q_i$ is to residue $t_j$.
- $E$ and $F$ track gaps (insertions/deletions) with **affine gap
  penalties**: opening a gap costs `gap_open + gap_extend`, extending it
  costs just `gap_extend`. This models biology — a single mutation is
  more likely to insert multiple consecutive residues than multiple
  independent single-residue insertions.
- The $0$ term means the score never goes negative — the algorithm
  "restarts" whenever continuing would produce a negative score.

### Substitution matrices

For **nucleotide** alignment (DNA), scoring is simple: matches get a
positive score (e.g. +2), mismatches get a negative score (e.g. -3).
The letter "N" (unknown base) scores 0 against everything.

For **protein** alignment, scoring uses a **substitution matrix** like
BLOSUM62. This is a 20×20 table (one row/column per amino acid) that
reflects how often each pair of amino acids is observed to substitute for
each other in evolution. For example, two chemically similar amino acids
(like Isoleucine and Leucine) get a positive score, while dissimilar
pairs get a negative score.

### Gap penalties

Gaps represent insertions or deletions (indels). The **affine gap model**
uses two parameters:
- **gap_open** — the penalty for starting a new gap. Larger values
  discourage opening many small gaps.
- **gap_extend** — the penalty for extending an existing gap by one
  position. Usually smaller than `gap_open`, so one long gap is cheaper
  than many short gaps.

For example, riot_na's nucleotide mode uses `gap_open=4`,
`gap_extend=1`. A 3-residue gap costs $4 + 1 + 1 = 6$, not $4 + 4 + 4 = 12$.

### Alignment output

The algorithm produces:
- **Optimal alignment score** — the maximum value in the DP matrix.
- **End position** — where in query and target the best alignment ends.
- **Start position** — found by running a second, reverse SW pass from
  the end position.
- **CIGAR string** — a compact encoding of the alignment path, e.g.
  `3M1I2M` means "3 matches, 1 insertion, 2 matches". Found by running a
  banded traceback between the start and end positions.

### Profile (PSSM) alignment

Instead of comparing residue-vs-residue through a substitution matrix,
**profile alignment** uses a Position-Specific Scoring Matrix (PSSM).
Each position in the query has its own row of 20 scores — one per amino
acid — reflecting the conservation pattern at that position (e.g. from a
multiple sequence alignment). This is more sensitive than sequence-sequence
alignment for detecting remote homologs.

---

## 3. How SIMD acceleration works

### The problem with naive SW

The basic SW algorithm has $O(m \times n)$ time complexity (where $m$ and
$n$ are the sequence lengths). For riot_na, aligning 1,000 antibody
sequences against 327 V-gene references means ~327,000 alignments. Each
alignment fills a matrix of ~100×300 = 30,000 cells. That is roughly 10
billion cell operations — too slow for pure scalar code.

### SIMD: Single Instruction, Multiple Data

Modern CPUs have **SIMD** instructions that operate on multiple values
simultaneously using wide registers:

| Instruction set | Register width | Year |
|---|---|---|
| SSE2 | 128 bits | 2001 |
| AVX2 | 256 bits | 2013 |

A 128-bit SSE2 register can hold **16 × 8-bit values** or **8 × 16-bit
values**. One SIMD instruction processes all of them in parallel — a
16× or 8× speedup over scalar code.

### Striped profile layout

You cannot parallelise SW trivially because each cell depends on its
neighbours. The **Farrar striped** approach (Farrar, 2007) rearranges the
query profile so that SIMD lanes compute independent segments of the
query simultaneously. The segments are "striped" — interleaved so that
data dependencies between adjacent cells are resolved in a single
correction pass per column.

For example, with 16 SIMD lanes and a query of length 48:
- Lane 0 handles positions 0, 16, 32
- Lane 1 handles positions 1, 17, 33
- Lane 2 handles positions 2, 18, 34
- etc.

This layout is pre-computed into a **query profile** — a rearranged array
of substitution scores indexed by target residue. Building this profile
is a one-time cost per query; all subsequent target alignments reuse it.

### The inner loop

The SIMD inner loop processes one column of the DP matrix per iteration.
Both scikit-bio and MMseqs2 use the same algorithm — about 12–14 SIMD
instructions per segment:

```
for each target residue j:
    for each segment i:
        load profile score for (target[j], query segment i)
        H = H + score                    // score lookup
        H = max(H, E)                    // gap in query
        H = max(H, F)                    // gap in target
        track max score
        store H
        update E (gap extension)
        update F (gap extension)
```

**The inner loop is structurally identical between scikit-bio and
ssw-aligner.** The performance difference comes entirely from overhead
*outside* this loop.

---

## 4. The byte / word / int cascade — why "255" matters

### Trading precision for speed

The number of values that fit in one SIMD register depends on the
integer width:

| Mode | Type | Bits | Values per 128-bit register | Max representable score | Relative speed |
|---|---|---|---|---|---|
| **byte** | `uint8_t` | 8 | **16** | 255 | Fastest |
| **word** | `uint16_t` | 16 | 8 | 65,535 | ~2× slower |
| **int** | `int32_t` | 32 | 4 | ~2 billion | ~4× slower |

The engine always starts in **byte mode** for maximum speed. If the
optimal score exceeds 255 (the maximum `uint8_t` value), the byte-mode
result is saturated (clipped at 255) and incorrect. The engine detects
this (`score == 255`) and **falls back to word mode** (16-bit), which is
slower but can handle scores up to 65,535.

In practice, most nucleotide alignments score well below 255. Protein
alignments with BLOSUM62 can exceed 255 for long, highly similar
sequences — this happens for ~7% of riot_na's V-gene alignments.

The int mode (32-bit) is commented out in our code because real-world
scores virtually never exceed 65,535.

### When does overflow happen?

The alignment score depends on:
- Sequence length (longer sequences accumulate more score)
- Sequence similarity (more matches = higher score)
- Substitution matrix values (BLOSUM62 has max entry +11 for W-W)
- Gap penalties (lower penalties allow longer alignments)

For riot_na's V-gene protein alignments (~100–300 amino acids,
BLOSUM62), about 7% of pairs produce scores above 255, triggering the
word-mode fallback.

---

## 5. Why we chose MMseqs2 over the original SSW library

| Feature | Mengyao Zhao SSW (scikit-bio) | MMseqs2 engine (ssw-aligner) |
|---|---|---|
| Core SIMD algorithm | Striped SW (Farrar) | Striped SW (Farrar) — identical |
| Profile (PSSM) alignment | No | Yes |
| AVX2 support | No | Yes (compile-time flag) |
| Gumbel statistics | No | Yes (ALP library) |
| Composition bias correction | No | Yes (can be disabled) |
| Build system | Compiled inline via Cython | CMake (standalone `.so`) |
| Python wrapper | Cython (`cdef extern`) | Pure ctypes |
| Codebase complexity | ~900 lines (C) | ~1,600 lines (C++) |

The MMseqs2 engine gives us profile alignment and Gumbel statistics for
free, at the cost of higher codebase complexity and some performance
overhead from unused features.

---

## 6. Architecture of ssw-aligner

The system has three layers:

```
Python application (riot_na, user code)
        │
        ▼
┌─────────────────────────────────────────┐
│  ssw_aligner/_wrapper.py  (Python)      │
│  ├─ StripedSmithWaterman                │  ← API-compatible with scikit-bio
│  ├─ SmithWatermanProfileAligner         │  ← new: PSSM alignment
│  ├─ AlignmentStructure                  │  ← result container
│  ├─ compute_gumbel_params()             │  ← new: statistical significance
│  └─ ctypes FFI calls ──────────────────────┐
└─────────────────────────────────────────┘  │
                                              ▼
┌─────────────────────────────────────────┐
│  src/api/ssw_api.cpp  (C adapter)       │
│  ├─ ssw_init()         → build profile  │
│  ├─ ssw_init_profile() → build PSSM     │
│  ├─ ssw_align()        → run alignment  │
│  ├─ ssw_free_cigar()   → free memory    │
│  ├─ ssw_destroy()      → cleanup        │
│  └─ compute_gumbel_params()             │
└────────────────┬────────────────────────┘
                 │ calls
                 ▼
┌─────────────────────────────────────────┐
│  src/alignment/StripedSmithWaterman.cpp │
│  (MMseqs2 SIMD engine, ~1,600 lines)   │
│  ├─ SmithWaterman constructor           │  ← allocates SIMD buffers
│  ├─ ssw_init()  → build query profile   │  ← rearrange scores for SIMD
│  ├─ ssw_align() → run alignment         │  ← byte → word cascade
│  ├─ sw_sse2_byte<T>()                   │  ← inner SIMD loop (8-bit)
│  ├─ sw_sse2_word<T>()                   │  ← inner SIMD loop (16-bit)
│  ├─ alignStartPosBacktrace<T>()         │  ← reverse pass + cigar
│  └─ banded_sw<T>()                      │  ← banded traceback for cigar
└─────────────────────────────────────────┘
```

### How a single alignment flows

1. **Initialisation (once per query):**
   - Python encodes the query string to numeric indices (e.g. A→0, C→1, …)
   - Python builds a flat substitution matrix (n×n array of int8 scores)
   - ctypes calls `ssw_init()` → C adapter creates a `SmithWaterman`
     object, which allocates SIMD buffers and builds the striped query
     profile

2. **Alignment (once per target):**
   - Python encodes the target string to numeric indices
   - ctypes calls `ssw_align()` with gap penalties, flag, and mask length
   - C adapter maps the flag → `alignmentMode` (0 = score only, 3 = full)
   - Engine runs forward SW in byte mode (`sw_sse2_byte`)
   - If score == 255, re-runs in word mode (`sw_sse2_word`)
   - If `alignmentMode > 0`: runs reverse SW to find start position,
     then `banded_sw` to produce CIGAR
   - Returns score, start/end positions, and CIGAR array

3. **Result decoding:**
   - Python copies the CIGAR from C memory into a Python list
   - Frees the C-allocated CIGAR array
   - Builds an `AlignmentStructure` with score, positions, CIGAR string,
     and aligned sequences

---

## 7. Block-aligner: what it was and why we removed it

### What block-aligner is

[Block-aligner](https://github.com/Daniel-Liu-c0deb0t/block-aligner) is
a Rust library by Daniel Liu that uses an adaptive banding strategy for
Smith-Waterman alignment. In upstream MMseqs2, it serves as an
alternative traceback implementation for alignments where scores overflow
word mode (>65,535) — it's designed to be more cache-friendly than the
standard banded traceback for very long sequences.

### Why it was a stub in ssw-aligner

Building the real block-aligner requires a Rust toolchain and WASM SIMD
compilation, which would significantly complicate the build. Instead,
ssw-aligner compiled a **stub** (`block_aligner_stub.c`) where
`block_aligner_align()` always returned `score = -1,000,000,000` — an
intentional failure.

### What happened at runtime

When a word-mode alignment needed traceback, the engine would:

1. Try block-aligner → stub returns failure score (-1,000,000,000)
2. Engine detects failure → prints "Block alignment failed" to stderr
3. Falls back to the standard C++ reverse SW + `banded_sw` traceback

This happened on **every word-mode alignment** (~7% of riot_na's
workload), adding:
- Unnecessary struct allocation per alignment
- A function call that always fails
- Warning messages printed to stderr
- A wasted conditional branch

### What we removed

- `lib/block-aligner/` directory (stub C files)
- `#include "block_aligner.h"` from `StripedSmithWaterman.cpp`
- `struct s_block` allocation in the `SmithWaterman` constructor
- `alignStartPosBacktraceBlock` method (~190 lines)
- The block-aligner call path in `ssw_align_private()`
- Block-aligner references in `CMakeLists.txt`

### Impact

- **+8% throughput** (25,298 → 27,295 aln/s) — from eliminating the
  stub overhead on word-mode alignments
- **Zero warnings** — no more "Block alignment failed" stderr noise
- **No functional loss** — the stub never produced a valid result
- **Score agreement unchanged** at 93.1% vs scikit-bio

---

## 8. Performance analysis: why we are slower than scikit-bio

### Benchmark setup

We benchmark by aligning 1,000 amino-acid antibody sequences (from
riot_na's NGS dataset) against 327 human V-gene germline references using
BLOSUM62 with `gap_open=11`, `gap_extend=1`. This produces 327,000
alignments and mirrors riot_na's real protein alignment workflow.

### Results

| Engine | Throughput (aln/s) |
|---|---|
| ssw-aligner (after optimisations) | ~27,750 |
| scikit-bio 0.6.2 | ~38,100 |

ssw-aligner is **~1.37× slower**. The inner SIMD DP kernel is
structurally identical — the slowdown comes entirely from surrounding
overhead, broken down into three categories.

### Category 1: ctypes vs Cython call overhead (~40–50% of the gap)

**What ctypes and Cython do differently:**

When Python code calls a C function, there is overhead for converting
Python objects to C types and dispatching the function call.

**Cython** (scikit-bio's approach) compiles to C at build time. The C
function call is a normal `function_name(arg1, arg2, …)` call — the CPU
jumps directly to the function. Cost: ~5 nanoseconds.

**ctypes** (ssw-aligner's approach) uses Python's `libffi` at runtime to
dynamically construct each function call. Every argument must be "boxed"
into a ctypes object (`ctypes.c_int32(value)`), then unpacked by libffi
into the correct register/stack position. Cost: ~200–500 nanoseconds per
call.

For `ssw_align`, we box **9 arguments** per call:

```python
result = _lib.ssw_align(
    self._handle,                                          # void*
    ref_seq.ctypes.data_as(ctypes.POINTER(ctypes.c_int8)), # pointer conversion
    ctypes.c_int32(len(target_sequence)),                   # integer boxing
    ctypes.c_uint8(self._gap_open),                         # integer boxing
    ctypes.c_uint8(self._gap_extend),                       # integer boxing
    ctypes.c_uint8(self._bit_flag),                         # integer boxing
    ctypes.c_uint16(self._score_filter),                    # integer boxing
    ctypes.c_int32(self._distance_filter),                  # integer boxing
    ctypes.c_int32(self._mask_length),                      # integer boxing
)
```

**Additional Python-side overhead:**

- `_encode_sequence()` is a **pure Python loop** that runs on every
  target sequence, converting each character to a numeric index:
  ```python
  for i, ch in enumerate(sequence):
      arr[i] = table[ord(ch)]
  ```
  In scikit-bio, this same loop is Cython-compiled (runs as fast C code).

- `AlignmentStructure.__init__` **eagerly decodes** the entire CIGAR in
  Python. scikit-bio lazily defers cigar parsing to a memoized property
  that is only computed if the user actually accesses it.

### Category 2: initialisation overhead (~30–40% of the gap)

When we create a `StripedSmithWaterman(query)`, the C adapter builds the
MMseqs2 engine, which allocates far more memory than scikit-bio's
simpler engine:

| Resource | scikit-bio allocs | ssw-aligner allocs |
|---|---|---|
| **Profile arrays** | 2 (byte + word) | 7 (byte/word/int32 × fwd+rev, plus linear) |
| **Total heap allocations** | 2–3 | ~30+ |
| **Extra structures** | — | `SimpleBaseMatrix`, `Sequence`, bias arrays, reverse query, reverse matrix, temp buffers |

Most of the extra allocations are for features we do not use:
- **Int32 profile** (fwd + rev) — for scores >65,535, which never happen
- **Linear profiles** (word + int32) — for the `scoreIdentical` method,
  which we never call
- **Composition bias arrays** (fwd + rev) — for correcting biased amino
  acid composition, which we disable (`aaBiasCorrection=false`)
- **Reverse profiles** — pre-allocated but rebuilt on demand anyway

Additionally, the `createQueryProfile` function loads and adds the
zeroed composition bias array on every profile element, even when bias
correction is disabled — adding an extra memory load and addition per
element for zero benefit.

### Category 3: per-alignment dead code (~10–20% of the gap)

Before our optimisation (commit `a43b728`), the engine ran several
operations on every alignment whose results were immediately discarded.
See [Section 9](#9-dead-code-removal-optimisation) for details and the
fix.

---

## 9. Dead code removal optimisation

### What we found

By tracing the call path from `ssw_api.cpp` through
`ssw_align_private()` and `alignStartPosBacktrace()`, we identified four
pieces of code that ran on every alignment but produced unused results.

#### 1. `computerBacktrace()` — unused string construction

After generating the CIGAR (a compact array of 32-bit integers), the
engine called `computerBacktrace()`, which iterated over the entire CIGAR
building a `std::string` by appending "M", "I", or "D" one character at a
time:

```cpp
// For a cigar like 50M2I30M, this builds the string
// "MMMMM...MMIIMMMM...M" (82 characters)
for (int32_t c = 0; c < alignment.cigarLen; ++c) {
    char letter = cigar_int_to_op(alignment.cigar[c]);
    uint32_t length = cigar_int_to_len(alignment.cigar[c]);
    for (uint32_t i = 0; i < length; ++i) {
        backtrace.append("M");  // heap allocation + copying
    }
}
```

This string was passed by reference from `ssw_api.cpp`:
```cpp
std::string backtrace;  // allocated on stack, grows on heap
s_align a = handle->sw->ssw_align(..., backtrace, ...);
// backtrace is immediately discarded here
```

Every cigar-mode alignment built and threw away a heap-allocated string.

#### 2. `computeCov()` × 2 — coverage always passes

The engine computed query and target coverage ratios:
```cpp
align.qCov = computeCov(0, align.qEndPos1, query_length);
align.tCov = computeCov(0, align.dbEndPos1, db_length);
bool hasLowerCoverage = !(Util::hasCoverage(covThr, covMode, ...));
```

Our API always passes `covThr=0.0`, which means `hasCoverage()` trivially
returns `true` and `hasLowerCoverage` is always `false`. The coverage
values are computed and never used.

#### 3. E-value null-pointer check

```cpp
align.evalue = (evaluer != NULL) ? evaluer->computeEvalue(...) : 0.0;
bool hasLowerEvalue = align.evalue > evalueThr;
```

Our API always passes `evaluer=nullptr` and `evalueThr=DBL_MAX`, so
`evalue` is always 0.0 and `hasLowerEvalue` is always `false`. Pure
overhead.

#### 4. Correlation score computation

```cpp
if (correlationScoreWeight > 0.0) {
    int correlationScore = computeCorrelationScore(scorePerCol, mStateCnt);
    r.score1 += correlationScore * correlationScoreWeight;
}
```

Our API always passes `correlationScoreWeight=0.0`, so this block is
never entered, but the surrounding work (`computerBacktrace` populating
`scorePerCol`) still runs.

### What we changed

We simplified the function signatures by removing all dead parameters:

**Before** (12 parameters):
```cpp
s_align ssw_align(const unsigned char* db_sequence, int32_t db_length,
    std::string& backtrace,
    uint8_t gap_open, uint8_t gap_extend, uint8_t alignmentMode,
    double evalueThr, EvalueComputation* evaluer,
    int covMode, float covThr, float correlationScoreWeight,
    int32_t maskLen);
```

**After** (6 parameters):
```cpp
s_align ssw_align(const unsigned char* db_sequence, int32_t db_length,
    uint8_t gap_open, uint8_t gap_extend, uint8_t alignmentMode,
    int32_t maskLen);
```

Removed from `ssw_align_private()`:
- `backtrace` string reference and the `computerBacktrace()` call
- `evalueThr`, `evaluer`, and all E-value computation
- `covMode`, `covThr`, and all coverage computation
- `correlationScoreWeight` and `computeCorrelationScore()` call

Removed from `alignStartPosBacktrace()`:
- Same parameters as above
- `computeCov()` calls after finding start position
- `hasLowerCoverage` early-exit (unreachable with `covThr=0`)

Updated `ssw_api.cpp`:
- Removed `std::string backtrace` allocation
- Removed unused `#include <limits>` and `#include <string>`
- Simplified `ssw_align()` call to pass only 6 arguments

### About alignmentMode

The original code supported 4 alignment modes (0, 1, 2, 3). Our C API
adapter (`ssw_api.cpp`) only ever maps to two:

```cpp
uint8_t alignmentMode = (flag == 0) ? 0 : 3;
```

- **Mode 0**: return score + end positions only
- **Mode 3**: full result (score + start/end positions + CIGAR)

Modes 1 and 2 were MMseqs2-specific filter modes that early-exit when
E-value or coverage thresholds are not met. Since our API always passes
`covThr=0` (always passes) and `evaluer=nullptr` (no E-value), the
mode 1/2 early-exit conditions **always evaluated to false** — so
removing them had no effect on behaviour.

### Impact

- **Before**: 27,295 aln/s
- **After**: 27,753 aln/s (+1.7%)
- No change in alignment results (52 tests pass, 93.1% score agreement)
- Cleaner, simpler code with half the parameters

The gain is modest because this category was only ~10–20% of the total
overhead. The larger contributors (ctypes calling convention and init
allocations) require more invasive changes.

---

## 10. The 6.9 % score divergence explained

### What diverges

When comparing ssw-aligner against scikit-bio on 32,700 (100 queries ×
327 targets) alignment pairs, 93.1% produce identical optimal scores.
The 6.9% that diverge are exactly the alignments where **byte-mode
scores overflow** (score > 255).

### Why it happens

Both engines get the **same optimal score** from the word-mode forward
pass — the inner SIMD kernel is identical. The divergence is in the
**traceback** (finding start positions and CIGAR), which uses different
implementations:

1. **scikit-bio** (Mengyao Zhao library):
   - Reverses the query sequence
   - Builds a reverse profile (`qP_byte` or `qP_word`)
   - Runs `sw_sse2_byte`/`sw_sse2_word` backward
   - Uses its own `banded_sw` for CIGAR generation

2. **ssw-aligner** (MMseqs2 engine):
   - Reverses the query sequence
   - Rebuilds the reverse profile via `createQueryProfile` (different
     template-based code with composition bias handling)
   - Runs `sw_sse2_byte`/`sw_sse2_word` backward
   - Uses its own `banded_sw` (different implementation with different
     edge-case handling)

The reverse profile construction and banded traceback are subtly different
between the two codebases, producing different start positions and CIGARs
for some word-mode alignments. The optimal alignment **score** itself is
always correct in both engines — only the specific alignment path
(which of potentially multiple optimal paths is chosen) differs.

### Is this a problem?

For riot_na's use case, no. The numbering algorithm uses scores for
ranking and the CIGAR for positional mapping. When scores match (93.1%
of cases), the CIGAR almost always matches too. For the 6.9% that differ,
both alignments are valid — they just represent different optimal paths
through the DP matrix. Multiple optimal alignments frequently exist for
protein alignments with substitution matrices.

---

## 11. Chronological list of code changes

### Phase 1 — Initial extraction (commit `eeb5614`)

**Goal:** Replace the Cython/C wrapper with MMseqs2 engine + ctypes.

**Files created:**
- `src/alignment/StripedSmithWaterman.cpp/.h` — SIMD engine from MMseqs2
- `src/commons/Sequence.cpp/.h`, `SubstitutionMatrix.cpp/.h`,
  `BaseMatrix.cpp/.h`, `Parameters.cpp/.h` — supporting MMseqs2 classes
- `src/api/ssw_api.h`, `ssw_api.cpp` — C adapter layer
- `ssw_aligner/_wrapper.py` — pure ctypes Python wrapper
- `lib/alp/` — Gumbel parameter estimation library
- `lib/simde/` — portable SIMD intrinsics
- `lib/block-aligner/` — stub block-aligner
- `CMakeLists.txt` — build system
- `examples/basic_usage.py`, `examples/profile_alignment.py`
- `.github/copilot-instructions.md`
- `README.md`

**Files deleted:**
- `ssw_aligner/_lib/` (old C library: ssw.c, ssw.h)
- `ssw_aligner/_ssw_wrapper.c`, `_ssw_wrapper.pyx` (old Cython wrapper)

**Key decisions:**
- Pure ctypes (no Cython) for simpler builds, at the cost of call
  overhead (~200–500 ns per call vs ~5 ns)
- 20-letter amino acid alphabet (matching MMseqs2's `PROFILE_AA_SIZE`)
  with internal `n+1 = 21` for profile mode to prevent X-row corruption
- Precomputed Gumbel parameters for common configurations to avoid
  expensive Monte Carlo computation

### Phase 2 — Build system (commits `212bbf8` → `918e3db`)

- Removed Cython/NumPy build dependencies
- Migrated to poetry-core as PEP 517 backend
- Created `scripts/build_ext.py` (auto-runs CMake during `poetry build`)
- Fixed the build script shadowing Python's `build` module
- Added integration tests (`test_install.py`, `smoke_test.py`)

### Phase 3 — Benchmarking (commit `4fea10b`)

- Created `tests/test_performance.py` — 1,000 queries × 327 V-genes
- Initial results: ssw-aligner 25,848 aln/s vs scikit-bio 38,583 aln/s
- Identified 93.1% score agreement, 6.9% divergence on word-mode
  alignments

### Phase 4 — Block-aligner removal (commit after `4fea10b`)

**Problem:** "Block alignment failed" warnings on every word-mode
alignment. Initially misattributed to scikit-bio, then correctly traced
to ssw-aligner's stale `.so` (the Python wrapper was loading an old copy
from `ssw_aligner/libssw_aligner.so` instead of the freshly built
`build/libssw_aligner.so`).

**Investigation:**
- `nm -D libssw_aligner.so | grep block` confirmed block symbols present
- `lib/block-aligner/c/block_aligner_stub.c` always returns
  `score = -1,000,000,000`
- Engine treats this as failure → falls back to C++ traceback → prints
  warning

**Changes:**
- Deleted `lib/block-aligner/` directory
- Removed from `StripedSmithWaterman.cpp`: `#include "block_aligner.h"`,
  `struct s_block`, constructor/destructor allocations,
  `alignStartPosBacktraceBlock` method (~190 lines), block-aligner call
  path in `ssw_align_private()`, template instantiations
- Removed from `StripedSmithWaterman.h`: forward decl, method decl,
  member variable
- Removed from `CMakeLists.txt`: include dir, source file

**Post-removal results:**
- 27,295 aln/s (+8% from 25,298) — pure gain from eliminating stub
  overhead
- Zero warnings
- 93.1% score agreement (unchanged)
- Library size decreased from 1,157,656 to 1,138,648 bytes

**Lesson learned:** Always copy `build/libssw_aligner.so` to
`ssw_aligner/libssw_aligner.so` after rebuilding. The wrapper loads
from the package directory first.

### Phase 5 — Performance analysis (commit `554df03`)

Detailed investigation comparing both engines' source code:
- Inner SIMD kernels are structurally identical — same operations, same
  order
- Identified three overhead categories: ctypes (40–50%), init (30–40%),
  dead code (10–20%)
- Documented in `.github/copilot-instructions.md`

### Phase 6 — Dead code removal (commit `a43b728`)

Stripped unused code from the hot path. See
[Section 9](#9-dead-code-removal-optimisation) for full details.

Results: 27,295 → 27,753 aln/s (+1.7%)

---

## 12. Remaining work and future optimisations

### Optimisations not yet implemented

| Optimisation | Expected impact | Effort | Description |
|---|---|---|---|
| **AVX2 build** | ~2× SIMD throughput | Low | Add `-DHAVE_AVX2=ON` to CMake. Doubles SIMD register width (128→256 bits), so byte mode processes 32 cells per instruction instead of 16. Requires AVX2-capable CPU (most x86 since ~2013). |
| **Batch alignment API** | Amortise ctypes overhead | Medium | New C function `ssw_align_batch(handle, targets[], N)` that aligns multiple targets in one call, reducing the ~200–500 ns ctypes overhead from per-alignment to per-batch. |
| **Remove unused profile building** | Faster init | Low | Skip int32 profiles (forward+reverse), linear profiles (word+int32), and composition bias arrays in the `SmithWaterman` constructor. These are allocated but never used in sequence-sequence alignment. |
| **Remove composition bias loads** | Fewer memory accesses | Low | In `createQueryProfile`, skip the load and add of the zeroed composition bias array when `aaBiasCorrection=false`. |
| **Cython wrapper** | Eliminate 40–50% of gap | High | Replace ctypes with Cython `cdef extern` calls. Eliminates libffi dispatch, argument boxing, and pure-Python sequence encoding. Would require Cython as a build dependency (losing one of the original motivations for switching to ctypes). |

### Known issues

1. **Stale `.so` after rebuild** — must manually copy
   `build/libssw_aligner.so` to `ssw_aligner/libssw_aligner.so`
2. **~1.37× slower than scikit-bio** — from ctypes overhead and engine
   complexity, not from the SIMD kernel itself
3. **6.9% score divergence** — different traceback implementations for
   word-mode alignments; scores are always correct
4. **No Windows/macOS CI** — only tested on Linux

### What a new developer should know

- **The inner SIMD loop is not the bottleneck.** Do not try to optimise
  `sw_sse2_byte` or `sw_sse2_word` — they are identical to scikit-bio's
  and already optimal.
- **ctypes overhead is the biggest contributor** to the performance gap.
  Any change that reduces the number of Python→C round-trips will have
  the most impact.
- **The `SmithWaterman` constructor allocates ~30 objects** but only uses
  about 10 of them for basic sequence-sequence alignment. Stripping the
  unused allocations is a safe, low-effort optimisation.
- **Always copy the `.so` after rebuilding.** The Python wrapper loads
  `ssw_aligner/libssw_aligner.so` first, not `build/libssw_aligner.so`.
  Forgetting this step means you are testing against the old binary.
- **Run `pytest tests/test_regression.py tests/test_mmseqs_freqs.py -q`
  after every C++ change** to verify correctness (52 tests).
