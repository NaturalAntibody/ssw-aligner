# ssw-aligner

Standalone [StripedSmithWaterman (SSW)](http://www.plosone.org/article/info:doi/10.1371/journal.pone.0082138) aligner extracted from [scikit-bio](https://scikit.bio).

This micro-library provides the SSW C implementation and its Python/Cython wrapper as a drop-in replacement for `skbio.alignment.StripedSmithWaterman`.

`StripedSmithWaterman` keeps the original scikit-bio defaults for backward compatibility.
For RIOT-compatible defaults, use `NucleotideAligner` or `ProteinAligner`.

## Installation

```bash
pip install ssw-aligner-na   # prebuilt wheels for Linux, macOS and Windows
```

The distribution is named `ssw-aligner-na` on PyPI (`ssw_aligner` is taken by an
unrelated project); the import name is still `ssw_aligner`.

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

Pull requests run `.github/workflows/pull-request.yml` (tests, then wheels and
sdist for every platform). Every push to `master` runs
`.github/workflows/push-master.yml` (following python-semantic-release's
recommended workflow): tests, then python-semantic-release
bumps the version from the conventional commits (and refreshes `uv.lock`),
cibuildwheel builds wheels for Linux (x86_64, aarch64), macOS (universal2) and
Windows (AMD64) for CPython 3.10–3.15 and PyPy 3.11, `uv build --sdist` builds
the sdist; the files are then attached to the GitHub release and uploaded to
PyPI with trusted publishing. Build and publish only run when a new version
was released.

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
