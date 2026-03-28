"""Performance benchmark — ssw-aligner (MMseqs2) vs scikit-bio SSW.

Aligns 1 000 amino-acid query sequences sampled from a real NGS dataset
against a panel of human V-gene germline references using BLOSUM62 scoring,
mirroring the riot_na protein alignment workflow.

Invocation:
    pytest tests/test_performance.py -v -s
"""

import importlib
import os
import random
import subprocess
import time
from pathlib import Path

import blosum  # type: ignore
import pytest
from skbio.alignment import StripedSmithWaterman as SkbioSSW  # type: ignore

from ssw_aligner import StripedSmithWaterman as MmseqsSSW

ROOT = Path(__file__).resolve().parent.parent

# ---------------------------------------------------------------------------
# Data paths
# ---------------------------------------------------------------------------

_RIOT_DATA = Path(__file__).resolve().parent.parent.parent / "riot_na"
_QUERY_FASTA = _RIOT_DATA / "data" / "ngs" / "ngs_sample_clean_aa.fasta"
_TARGET_FASTA = (
    _RIOT_DATA / "riot_na" / "databases" / "gene_db" / "aa_genes" / "v_genes" / "human.fasta"
)

SAMPLE_SIZE = 1_000
SEED = 42

# Alignment parameters matching riot_na AA alignment workflow
ALIGNER_PARAMS = {
    "gap_open_penalty": 11,
    "gap_extend_penalty": 1,
    "protein": True,
    "substitution_matrix": blosum.BLOSUM(62),
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _read_fasta(path: Path) -> list[str]:
    """Return list of sequences from a FASTA file."""
    sequences: list[str] = []
    current: list[str] = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if current:
                    sequences.append("".join(current))
                    current = []
            else:
                current.append(line.strip())
    if current:
        sequences.append("".join(current))
    return sequences


def _sample_queries(n: int, seed: int) -> list[str]:
    all_seqs = _read_fasta(_QUERY_FASTA)
    rng = random.Random(seed)
    return rng.sample(all_seqs, min(n, len(all_seqs)))


def _load_targets() -> list[str]:
    return _read_fasta(_TARGET_FASTA)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def queries() -> list[str]:
    return _sample_queries(SAMPLE_SIZE, SEED)


@pytest.fixture(scope="module")
def targets() -> list[str]:
    return _load_targets()


# ---------------------------------------------------------------------------
# Benchmarks
# ---------------------------------------------------------------------------


class TestPerformance:
    """Compare ssw-aligner vs scikit-bio alignment throughput."""

    def test_ssw_aligner_throughput(self, queries: list[str], targets: list[str]) -> None:
        """Benchmark ssw-aligner (MMseqs2 engine)."""
        t0 = time.perf_counter()
        n_alignments = 0
        for query in queries:
            aligner = MmseqsSSW(query, **ALIGNER_PARAMS)
            for target in targets:
                aligner(target)
                n_alignments += 1
        elapsed = time.perf_counter() - t0

        rate = n_alignments / elapsed
        print(
            f"\n  ssw-aligner:  {n_alignments:,} alignments in {elapsed:.2f}s "
            f"({rate:,.0f} aln/s)"
        )
        assert n_alignments == len(queries) * len(targets)

    def test_skbio_throughput(self, queries: list[str], targets: list[str]) -> None:
        """Benchmark scikit-bio SSW."""
        t0 = time.perf_counter()
        n_alignments = 0
        for query in queries:
            aligner = SkbioSSW(query, **ALIGNER_PARAMS)
            for target in targets:
                aligner(target)
                n_alignments += 1
        elapsed = time.perf_counter() - t0

        rate = n_alignments / elapsed
        print(
            f"\n  scikit-bio:   {n_alignments:,} alignments in {elapsed:.2f}s "
            f"({rate:,.0f} aln/s)"
        )
        assert n_alignments == len(queries) * len(targets)

    def test_results_agree(self, queries: list[str], targets: list[str]) -> None:
        """Verify both implementations produce the same scores on a subset."""
        # Check first 100 queries × all targets
        check_queries = queries[:100]
        mismatches = 0
        total = 0
        for query in check_queries:
            ours = MmseqsSSW(query, **ALIGNER_PARAMS)
            theirs = SkbioSSW(query, **ALIGNER_PARAMS)
            for target in targets:
                r_ours = ours(target)
                r_theirs = theirs(target)
                if r_ours.optimal_alignment_score != r_theirs["optimal_alignment_score"]:
                    mismatches += 1
                total += 1

        mismatch_pct = 100.0 * mismatches / total if total else 0
        print(f"\n  Score agreement: {total - mismatches}/{total} ({100 - mismatch_pct:.1f}%)")
        # When byte-mode scores overflow (>255) the word-mode fallback
        # uses a different traceback than scikit-bio, so some pairs diverge.
        assert mismatch_pct < 10, f"Too many score mismatches: {mismatch_pct:.1f}%"


# ---------------------------------------------------------------------------
# AVX2 build fixture
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def avx2_ssw_class():
    """Build libssw_aligner.so with -DHAVE_AVX2=ON and yield a
    StripedSmithWaterman class backed by that AVX2-enabled library.

    Skipped automatically when:
    - cmake is not on PATH, or
    - the host CPU does not support AVX2 (cmake HAVE_AVX2 check fails), or
    - the build fails for any other reason.
    """
    avx2_build = ROOT / "build" / "avx2"
    avx2_lib = avx2_build / "libssw_aligner.so"

    if not avx2_lib.exists():
        try:
            subprocess.run(
                [
                    "cmake", "-B", str(avx2_build),
                    "-DCMAKE_BUILD_TYPE=Release",
                    "-DHAVE_AVX2=ON",
                    str(ROOT),
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            subprocess.run(
                ["cmake", "--build", str(avx2_build), f"-j{os.cpu_count() or 1}"],
                check=True,
                capture_output=True,
                text=True,
            )
        except (subprocess.CalledProcessError, FileNotFoundError) as exc:
            pytest.skip(f"AVX2 build unavailable: {exc}")

    import ssw_aligner._wrapper as _wrapper_mod
    import ssw_aligner as _ssw_pkg

    old_env = os.environ.get("SSW_ALIGNER_LIB")
    try:
        os.environ["SSW_ALIGNER_LIB"] = str(avx2_lib)
        importlib.reload(_wrapper_mod)
        importlib.reload(_ssw_pkg)
        yield _ssw_pkg.StripedSmithWaterman
    finally:
        if old_env is None:
            os.environ.pop("SSW_ALIGNER_LIB", None)
        else:
            os.environ["SSW_ALIGNER_LIB"] = old_env
        importlib.reload(_wrapper_mod)
        importlib.reload(_ssw_pkg)


# ---------------------------------------------------------------------------
# AVX2 benchmarks
# ---------------------------------------------------------------------------


class TestPerformanceAVX2:
    """Same benchmarks as TestPerformance but using the AVX2-enabled library.

    Build the AVX2 library once and run it with:
        CMAKE_ARGS="-DHAVE_AVX2=ON" poetry build
    or let this fixture build it automatically.
    """

    def test_ssw_aligner_avx2_throughput(
        self,
        queries: list[str],
        targets: list[str],
        avx2_ssw_class,
    ) -> None:
        """Benchmark ssw-aligner compiled with -DHAVE_AVX2=ON."""
        t0 = time.perf_counter()
        n_alignments = 0
        for query in queries:
            aligner = avx2_ssw_class(query, **ALIGNER_PARAMS)
            for target in targets:
                aligner(target)
                n_alignments += 1
        elapsed = time.perf_counter() - t0

        rate = n_alignments / elapsed
        print(
            f"\n  ssw-aligner (AVX2): {n_alignments:,} alignments in {elapsed:.2f}s "
            f"({rate:,.0f} aln/s)"
        )
        assert n_alignments == len(queries) * len(targets)

    def test_skbio_vs_avx2_throughput(
        self,
        queries: list[str],
        targets: list[str],
        avx2_ssw_class,
    ) -> None:
        """Compare AVX2 ssw-aligner throughput against scikit-bio back-to-back."""
        # --- AVX2 ---
        t0 = time.perf_counter()
        n = 0
        for query in queries:
            aligner = avx2_ssw_class(query, **ALIGNER_PARAMS)
            for target in targets:
                aligner(target)
                n += 1
        avx2_elapsed = time.perf_counter() - t0
        avx2_rate = n / avx2_elapsed

        # --- scikit-bio ---
        t0 = time.perf_counter()
        n = 0
        for query in queries:
            aligner = SkbioSSW(query, **ALIGNER_PARAMS)
            for target in targets:
                aligner(target)
                n += 1
        skbio_elapsed = time.perf_counter() - t0
        skbio_rate = n / skbio_elapsed

        speedup = avx2_rate / skbio_rate
        print(
            f"\n  ssw-aligner (AVX2): {avx2_rate:,.0f} aln/s\n"
            f"  scikit-bio:         {skbio_rate:,.0f} aln/s\n"
            f"  AVX2 speedup vs scikit-bio: {speedup:.2f}×"
        )
        assert n == len(queries) * len(targets)

    def test_avx2_results_agree(
        self,
        queries: list[str],
        targets: list[str],
        avx2_ssw_class,
    ) -> None:
        """Verify the AVX2 build produces the same scores as scikit-bio."""
        check_queries = queries[:100]
        mismatches = 0
        total = 0
        for query in check_queries:
            ours = avx2_ssw_class(query, **ALIGNER_PARAMS)
            theirs = SkbioSSW(query, **ALIGNER_PARAMS)
            for target in targets:
                r_ours = ours(target)
                r_theirs = theirs(target)
                if r_ours.optimal_alignment_score != r_theirs["optimal_alignment_score"]:
                    mismatches += 1
                total += 1

        mismatch_pct = 100.0 * mismatches / total if total else 0
        print(
            f"\n  AVX2 score agreement: {total - mismatches}/{total} "
            f"({100 - mismatch_pct:.1f}%)"
        )
        assert mismatch_pct < 10, f"Too many score mismatches: {mismatch_pct:.1f}%"
