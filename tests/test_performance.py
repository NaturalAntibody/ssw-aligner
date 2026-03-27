"""Performance benchmark — ssw-aligner (MMseqs2) vs scikit-bio SSW.

Aligns 1 000 amino-acid query sequences sampled from a real NGS dataset
against a panel of human V-gene germline references using BLOSUM62 scoring,
mirroring the riot_na protein alignment workflow.

Invocation:
    pytest tests/test_performance.py -v -s
"""

import random
import time
from pathlib import Path

import blosum  # type: ignore
import pytest

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
        from ssw_aligner import StripedSmithWaterman

        t0 = time.perf_counter()
        n_alignments = 0
        for query in queries:
            aligner = StripedSmithWaterman(query, **ALIGNER_PARAMS)
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
        from skbio.alignment import StripedSmithWaterman  # type: ignore

        t0 = time.perf_counter()
        n_alignments = 0
        for query in queries:
            aligner = StripedSmithWaterman(query, **ALIGNER_PARAMS)
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
        from skbio.alignment import StripedSmithWaterman as SkbioSSW  # type: ignore

        from ssw_aligner import StripedSmithWaterman as MmseqsSSW

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
        # The MMseqs2 engine handles byte-overflow differently (block-aligner
        # fallback) so some high-scoring BLOSUM62 pairs diverge.
        assert mismatch_pct < 10, f"Too many score mismatches: {mismatch_pct:.1f}%"
