"""Tests for use_mmseqs_aa_freqs flag and MMSEQS_AA_FREQS constant."""

import pytest
import blosum  # type: ignore

from ssw_aligner import compute_gumbel_params, MMSEQS_AA_FREQS

# Expected values computed via ALP with Robinson-Robinson background freqs
# matching the precomputed entry for (blosum62, 11, 1) in _wrapper.py
_LAMBDA_EXPECTED = 0.27359865037
_K_EXPECTED = 0.04462092066
_SIGMA_EXPECTED = 29.602445
_TAU_EXPECTED = -601.810880
_TOL = 1e-4


@pytest.fixture
def blosum62():
    return blosum.BLOSUM(62)


@pytest.fixture
def nt_matrix():
    mat = {}
    for r in "ACGTN":
        mat[r] = {}
        for c in "ACGTN":
            mat[r][c] = 1 if (r == c and r != "N") else (-1 if r != "N" and c != "N" else 0)
    return mat


class TestMmseqsAaFreqs:
    def test_lambda_matches_expected(self, blosum62):
        p = compute_gumbel_params(
            blosum62, gap_open=11, gap_extend=1, protein=True, use_mmseqs_aa_freqs=True
        )
        assert abs(p.lambda_ - _LAMBDA_EXPECTED) < _TOL

    def test_k_matches_expected(self, blosum62):
        p = compute_gumbel_params(
            blosum62, gap_open=11, gap_extend=1, protein=True, use_mmseqs_aa_freqs=True
        )
        assert abs(p.K - _K_EXPECTED) < _TOL

    def test_sigma_matches_expected(self, blosum62):
        p = compute_gumbel_params(
            blosum62, gap_open=11, gap_extend=1, protein=True, use_mmseqs_aa_freqs=True
        )
        assert abs(p.sigma - _SIGMA_EXPECTED) < _TOL

    def test_tau_matches_expected(self, blosum62):
        p = compute_gumbel_params(
            blosum62, gap_open=11, gap_extend=1, protein=True, use_mmseqs_aa_freqs=True
        )
        assert abs(p.tau - _TAU_EXPECTED) < _TOL

    def test_uniform_freqs_differ_from_mmseqs(self, blosum62):
        # Force ALP recalculation to bypass precomputed cache so both runs use
        # different background frequency assumptions
        p_uniform = compute_gumbel_params(
            blosum62, gap_open=11, gap_extend=1, protein=True,
            recalculate_gumbel_params=True
        )
        p_mmseqs = compute_gumbel_params(
            blosum62, gap_open=11, gap_extend=1, protein=True,
            use_mmseqs_aa_freqs=True, recalculate_gumbel_params=True
        )
        # Robinson-Robinson freqs should produce different lambda than uniform
        assert abs(p_uniform.lambda_ - p_mmseqs.lambda_) > 1e-6

    def test_nucleotide_raises_value_error(self, nt_matrix):
        with pytest.raises(ValueError, match="use_mmseqs_aa_freqs"):
            compute_gumbel_params(nt_matrix, 4, 1, protein=False, use_mmseqs_aa_freqs=True)


class TestMmseqsAaFreqsConstant:
    def test_has_twenty_entries(self):
        assert len(MMSEQS_AA_FREQS) == 20

    def test_all_standard_amino_acids_present(self):
        expected = set("ARNDCQEGHILKMFPSTWYV")
        assert set(MMSEQS_AA_FREQS.keys()) == expected

    def test_frequencies_sum_to_one(self):
        # Robinson-Robinson values sum to ~1.00001 due to rounding in the
        # original publication; tolerance is set accordingly
        total = sum(MMSEQS_AA_FREQS.values())
        assert abs(total - 1.0) < 1e-4

    def test_all_frequencies_positive(self):
        assert all(v > 0 for v in MMSEQS_AA_FREQS.values())

    def test_known_values(self):
        # Leucine is the most abundant AA in Robinson-Robinson
        assert MMSEQS_AA_FREQS["L"] > MMSEQS_AA_FREQS["W"]
