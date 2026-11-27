# ssw-aligner

[![PyPI](https://img.shields.io/pypi/v/ssw-aligner-na)](https://pypi.org/project/ssw-aligner-na/)
[![Push to master](https://github.com/NaturalAntibody/ssw-aligner/actions/workflows/push-master.yml/badge.svg)](https://github.com/NaturalAntibody/ssw-aligner/actions/workflows/push-master.yml)
[![License: BSD-3-Clause AND MIT](https://img.shields.io/badge/license-BSD--3--Clause%20AND%20MIT-blue)](LICENSE.txt)

Standalone [StripedSmithWaterman (SSW)](https://doi.org/10.1371/journal.pone.0082138)
aligner extracted from [scikit-bio](https://scikit.bio).

This micro-library provides the SSW C implementation and its Python/Cython
wrapper as a drop-in replacement for `skbio.alignment.StripedSmithWaterman`.

- `StripedSmithWaterman` keeps the original scikit-bio defaults for backward
  compatibility.
- `NucleotideAligner` and `ProteinAligner` use RIOT-compatible defaults.

## Installation

```bash
pip install ssw-aligner-na   # prebuilt wheels for Linux, macOS and Windows
```

> [!NOTE]
> The distribution is named `ssw-aligner-na` on PyPI (`ssw_aligner` is taken by
> an unrelated project); the import name is still `ssw_aligner`.

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

## Equal-score extension

Of several optimal alignments with the same score, SSW reports the shortest one. With
`extend_equal_score` (`"start"`, `"end"` or `"both"`, default `"none"`), the reported alignment is
extended towards the target start and/or end by the farthest path that adds exactly 0 to the score,
with at most `max_tie_insertions` (default 2) query residues beyond the target residues it covers.
Not part of scikit-bio; see [docs/equal_score_extension.md](docs/equal_score_extension.md).

```python
from ssw_aligner import ProteinAligner

aligner = ProteinAligner("GGGGSEVQLLEQSGAELVRPG", extend_equal_score="start")
result = aligner("QVQLQQSGAELVRPG")
print(result["cigar"])  # 3M1I12M from the target start (without the extension: 12M from position 4)
```

## RIOT-compatible wrappers

These wrappers use the same defaults as RIOT while leaving
`StripedSmithWaterman` unchanged for scikit-bio compatibility.

| Wrapper | RIOT parameters | Defaults |
|---|---|---|
| `NucleotideAligner` | `ALIGNER_PARAMS` | `match_score=1`, `mismatch_score=-1`, `gap_open_penalty=4`, `gap_extend_penalty=1` |
| `ProteinAligner` | `AA_ALIGNER_PARAMS` | `gap_open_penalty=11`, `gap_extend_penalty=1`, `protein=True`, `substitution_matrix=BLOSUM_62` |

```python
from ssw_aligner import NucleotideAligner, ProteinAligner

nt_aligner = NucleotideAligner(query_sequence="ACTGACTG")
nt_result = nt_aligner(target_sequence="AACTGACTGA")

protein_aligner = ProteinAligner(query_sequence="ARNDCQEGHILKMFPSTWYV")
protein_result = protein_aligner(target_sequence="ARNDAQEGHILKMFASTWYV")

print(nt_result["optimal_alignment_score"])
print(protein_result["optimal_alignment_score"])
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

No other code changes are needed.

## Development

The project is managed with [uv](https://docs.astral.sh/uv/) and built with
[scikit-build-core](https://scikit-build-core.readthedocs.io/) (CMake + Cython).
You need a C compiler; CMake, Cython and NumPy are installed into the isolated
build environment automatically.

```bash
uv sync               # create .venv, build the extension, install dev deps
uv run pytest         # run the tests (incl. regression against scikit-bio)
uv build              # build the sdist and a wheel for this machine into dist/
```

`uv sync` rebuilds the extension whenever `CMakeLists.txt` or a `.pyx`/`.c`/`.h`
file changes; the CMake build tree is kept in `build/` so rebuilds are incremental.

## Releasing

Releases are fully automated with
[python-semantic-release](https://python-semantic-release.readthedocs.io/),
following its recommended GitHub Actions workflow:

- **Pull requests** run `.github/workflows/pull-request.yml`: tests, then
  wheels and sdist for every platform.
- **Pushes to `master`** run `.github/workflows/push-master.yml`:
  1. tests;
  2. python-semantic-release bumps the version from the conventional commits
     (and refreshes `uv.lock`);
  3. cibuildwheel builds wheels for Linux (x86_64, aarch64), macOS (universal2)
     and Windows (AMD64) for CPython 3.10–3.15 and PyPy 3.11, and
     `uv build --sdist` builds the sdist;
  4. the files are attached to the GitHub release and uploaded to PyPI with
     trusted publishing.

Steps 3–4 only run when a new version was released.

## License

This package is licensed as `BSD-3-Clause AND MIT`. Full text is in
[`LICENSE.txt`](LICENSE.txt):

- SSW C library (`ssw.c` / `ssw.h`): MIT License (Copyright 2012-2015 Boston College)
- Cython wrapper: Modified BSD License (Copyright 2013-- scikit-bio development team)
- SIMDe header: MIT License
