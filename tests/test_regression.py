"""Regression tests for ssw-aligner.

These tests verify that the extracted SSW micro-library produces identical
results to `skbio.alignment.StripedSmithWaterman` for both nucleotide and
protein alignments, covering all the usage patterns found in riot_na.
"""

import pytest
import blosum  # type: ignore
from skbio.alignment import StripedSmithWaterman as SkbioSSW  # type: ignore

from ssw_aligner import StripedSmithWaterman


# ---------------------------------------------------------------------------
# Nucleotide alignment fixtures and helpers
# ---------------------------------------------------------------------------

NT_ALIGNER_PARAMS = {
    "match_score": 1,
    "mismatch_score": -1,
    "gap_open_penalty": 4,
    "gap_extend_penalty": 1,
}

AA_ALIGNER_PARAMS = {
    "gap_open_penalty": 11,
    "gap_extend_penalty": 1,
    "protein": True,
    "substitution_matrix": blosum.BLOSUM(62),
}


# ---------------------------------------------------------------------------
# 1. Basic nucleotide alignment (matches riot_na GeneAligner / test_skbio_alignment)
# ---------------------------------------------------------------------------

class TestNucleotideAlignment:
    """Nucleotide alignment regression tests."""

    def test_basic_alignment_from_riot_na(self):
        """Exact reproduction of riot_na's test_skbio_alignment.test_align."""
        query = "CTATACTACTATGGTTCGGGGAGTTATTATAGCCTTGACTACTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAGGGAGTGCATCCGCCCCAACCTCGT"
        target = "ACTACTTTGACTACTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAG"

        aligner = StripedSmithWaterman(query, **NT_ALIGNER_PARAMS)
        res = aligner(target)

        q_start = res["query_begin"]
        q_end = res["query_end"] + 1

        assert query[q_start:q_end] == "TTGACTACTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAG"

    def test_perfect_match(self):
        """Identical query and target must produce a perfect alignment."""
        seq = "ACGTACGTACGTACGT"
        aligner = StripedSmithWaterman(seq)
        res = aligner(seq)

        assert res["query_begin"] == 0
        assert res["query_end"] == len(seq) - 1
        assert res["target_begin"] == 0
        assert res["target_end_optimal"] == len(seq) - 1
        assert res["cigar"] == f"{len(seq)}M"
        assert res["optimal_alignment_score"] > 0

    def test_short_sequence(self):
        """Short sequences should still align properly."""
        query = "ACGT"
        target = "ACGT"
        aligner = StripedSmithWaterman(query)
        res = aligner(target)

        assert res["optimal_alignment_score"] > 0
        assert res["query_begin"] == 0
        assert res["query_end"] == 3

    def test_no_match(self):
        """Completely mismatched sequences should still return a result."""
        query = "AAAAAAAAAA"
        target = "CCCCCCCCCC"
        aligner = StripedSmithWaterman(query, match_score=2, mismatch_score=-3)
        res = aligner(target)

        # Score should be 0 or very low
        assert res["optimal_alignment_score"] == 0

    def test_gap_in_alignment(self):
        """Alignment with an insertion in target produces a gap."""
        query = "ACGTACGTACGT"
        target = "ACGTAAAACGTACGT"
        aligner = StripedSmithWaterman(query)
        res = aligner(target)

        assert res["optimal_alignment_score"] > 0
        assert "M" in res["cigar"]

    def test_query_subsequence(self):
        """Target is a substring of the query."""
        query = "AAACGTACGTAAA"
        target = "CGTACGT"
        aligner = StripedSmithWaterman(query)
        res = aligner(target)

        assert res["query_begin"] == 3
        assert res["query_end"] == 9
        assert res["target_begin"] == 0
        assert res["target_end_optimal"] == 6

    def test_alignment_result_properties(self):
        """All expected properties are accessible on the result."""
        query = "ACGTACGT"
        target = "ACGTACGT"
        aligner = StripedSmithWaterman(query)
        res = aligner(target)

        # Check that all properties riot_na uses via __getitem__ are present
        assert isinstance(res["optimal_alignment_score"], int)
        assert isinstance(res["query_begin"], int)
        assert isinstance(res["query_end"], int)
        assert isinstance(res["target_begin"], int)
        assert isinstance(res["target_end_optimal"], int)
        assert isinstance(res["cigar"], str)
        assert isinstance(res["query_sequence"], str)
        assert isinstance(res["target_sequence"], str)
        assert isinstance(res["suboptimal_alignment_score"], int)
        assert isinstance(res["target_end_suboptimal"], int)

    def test_aligned_sequences(self):
        """aligned_query_sequence and aligned_target_sequence work."""
        query = "ACGTACGT"
        target = "ACGTACGT"
        aligner = StripedSmithWaterman(query)
        res = aligner(target)

        assert res.aligned_query_sequence is not None
        assert res.aligned_target_sequence is not None

    def test_repr_and_str(self):
        """__repr__ and __str__ don't crash."""
        aligner = StripedSmithWaterman("ACGT")
        res = aligner("ACGT")
        assert isinstance(repr(res), str)
        assert isinstance(str(res), str)

    def test_suppress_sequences(self):
        """suppress_sequences hides the sequence data."""
        aligner = StripedSmithWaterman("ACGT", suppress_sequences=True)
        res = aligner("ACGT")
        assert res["query_sequence"] == ""
        assert res["target_sequence"] == ""

    def test_reuse_aligner(self):
        """Aligner can be reused with multiple targets."""
        query = "ACGTACGT"
        aligner = StripedSmithWaterman(query)

        targets = ["ACGTACGT", "ACGT", "TGCATGCA", "AAACGTACGTAAA"]
        results = [aligner(t) for t in targets]

        assert all(r["optimal_alignment_score"] >= 0 for r in results)
        # First target is a perfect match
        assert results[0]["cigar"] == "8M"

    def test_default_penalty_values(self):
        """Default gap penalties (5/2) produce valid results."""
        aligner = StripedSmithWaterman("ACGTACGT")
        res = aligner("ACGTACGT")
        assert res["optimal_alignment_score"] > 0

    def test_riot_na_nt_params(self):
        """riot_na GeneAligner params produce correct results."""
        query = "ACGTACGTACGTACGT"
        target = "ACGTACGTACGTACGT"
        aligner = StripedSmithWaterman(query, **NT_ALIGNER_PARAMS)
        res = aligner(target)

        assert res["optimal_alignment_score"] == 16  # 16 matches * 1
        assert res["cigar"] == "16M"

    def test_zero_index_true(self):
        """zero_index=True yields 0-based indices."""
        aligner = StripedSmithWaterman("ACGTACGT", zero_index=True)
        res = aligner("ACGTACGT")
        assert res["query_begin"] == 0

    def test_zero_index_false(self):
        """zero_index=False yields 1-based indices."""
        aligner = StripedSmithWaterman("ACGTACGT", zero_index=False)
        res = aligner("ACGTACGT")
        assert res["query_begin"] == 1


# ---------------------------------------------------------------------------
# 2. Protein alignment (matches riot_na AA_ALIGNER_PARAMS usage)
# ---------------------------------------------------------------------------

class TestProteinAlignment:
    """Protein alignment regression tests matching riot_na aa alignment usage."""

    def test_protein_perfect_match(self):
        """Identical protein sequences produce a perfect alignment."""
        seq = "ARNDCQEGHILKMFPSTWYV"
        aligner = StripedSmithWaterman(seq, **AA_ALIGNER_PARAMS)
        res = aligner(seq)

        assert res["query_begin"] == 0
        assert res["query_end"] == len(seq) - 1
        assert res["target_begin"] == 0
        assert res["target_end_optimal"] == len(seq) - 1
        assert res["cigar"] == f"{len(seq)}M"
        assert res["optimal_alignment_score"] > 0

    def test_protein_requires_substitution_matrix(self):
        """protein=True without substitution_matrix raises."""
        with pytest.raises(Exception, match="substitution matrix"):
            StripedSmithWaterman("ARND", protein=True)

    def test_protein_v_gene_alignment(self):
        """Simulate riot_na V gene AA alignment scenario."""
        # Realistic V gene query fragment
        query = "QVQLVQSGAEVKKPGASVKVSCKASGYTFTGYYMHWVRQAPGQGLEWMGWINPNSGGTNYA"
        target = "QVQLVQSGAEVKKPGASVKVSCKASGYTFTGYYMHWVRQAPGQGLEWMGWINPNSGGTNYAQKFQG"

        aligner = StripedSmithWaterman(query, **AA_ALIGNER_PARAMS)
        res = aligner(target)

        assert res["query_begin"] == 0
        assert res["query_end"] == len(query) - 1
        assert res["target_begin"] == 0
        assert res["optimal_alignment_score"] > 0
        assert "M" in res["cigar"]

    def test_protein_j_gene_alignment(self):
        """Simulate riot_na J gene AA alignment: short protein fragment."""
        query = "WGQGTLVTVSS"
        target = "WGQGTLVTVSSASTKGPS"

        aligner = StripedSmithWaterman(query, **AA_ALIGNER_PARAMS)
        res = aligner(target)

        assert res["query_begin"] == 0
        assert res["query_end"] == len(query) - 1
        assert res["target_begin"] == 0
        assert res["optimal_alignment_score"] > 0

    def test_protein_full_heavy_chain(self):
        """Full heavy chain alignment (same sequence used in riot_na __main__)."""
        query = (
            "QVQLQQWGAGLLKPSETLSLTCAVFGGSFSGYYWSWIRQPPGKGLEWIGEINHRGNTNDNPSLKS"
            "RVTISVDTSKNQFALKLSSVTAADTAVYYCARERGYTYGNFDHWGQGTLVTVSS"
        )
        target = query  # self-alignment
        aligner = StripedSmithWaterman(query, **AA_ALIGNER_PARAMS)
        res = aligner(target)

        assert res["query_begin"] == 0
        assert res["query_end"] == len(query) - 1
        assert res["cigar"] == f"{len(query)}M"

    def test_protein_alignment_with_mismatch(self):
        """Protein alignment with mismatches produces non-zero score."""
        query = "ARNDCQEGHILKMFPSTWYV"
        # Mutate a few positions
        target = "ARNDAQEGHILKMFASTWYV"

        aligner = StripedSmithWaterman(query, **AA_ALIGNER_PARAMS)
        res = aligner(target)

        assert res["optimal_alignment_score"] > 0
        assert res["query_begin"] >= 0
        assert res["target_begin"] >= 0

    def test_protein_reuse_aligner(self):
        """AA aligner can be reused with multiple targets (riot_na pattern)."""
        query = "QVQLVQSGAEVKKPGASVKVSCKAS"
        aligner = StripedSmithWaterman(query, **AA_ALIGNER_PARAMS)

        targets = [
            "QVQLVQSGAEVKKPGASVKVSCKAS",
            "EVQLVESGGGLVQPGGSLRLSCAAS",
            "QVQLVESGGGVVQPGRSLRLSCAAS",
        ]
        results = [aligner(t) for t in targets]

        assert all(r["optimal_alignment_score"] > 0 for r in results)
        # First target is a perfect match
        assert results[0]["cigar"] == f"{len(query)}M"


# ---------------------------------------------------------------------------
# 3. Parameter validation
# ---------------------------------------------------------------------------

class TestParameterValidation:
    """Validate that parameter constraints match scikit-bio behavior."""

    def test_gap_open_penalty_zero_raises(self):
        with pytest.raises(ValueError, match="gap_open_penalty"):
            StripedSmithWaterman("ACGT", gap_open_penalty=0)

    def test_gap_extend_penalty_zero_raises(self):
        with pytest.raises(ValueError, match="gap_extend_penalty"):
            StripedSmithWaterman("ACGT", gap_extend_penalty=0)

    def test_gap_open_negative_raises(self):
        with pytest.raises(ValueError, match="gap_open_penalty"):
            StripedSmithWaterman("ACGT", gap_open_penalty=-1)

    def test_gap_extend_negative_raises(self):
        with pytest.raises(ValueError, match="gap_extend_penalty"):
            StripedSmithWaterman("ACGT", gap_extend_penalty=-1)


# ---------------------------------------------------------------------------
# 4. AlignmentStructure interface
# ---------------------------------------------------------------------------

class TestAlignmentStructure:
    """AlignmentStructure exposes the same interface as scikit-bio."""

    def test_is_zero_based_default(self):
        aligner = StripedSmithWaterman("ACGT")
        res = aligner("ACGT")
        assert res.is_zero_based()

    def test_set_zero_based(self):
        aligner = StripedSmithWaterman("ACGT")
        res = aligner("ACGT")
        res.set_zero_based(False)
        assert not res.is_zero_based()
        res.set_zero_based(True)
        assert res.is_zero_based()


# ---------------------------------------------------------------------------
# 5. Cross-validation against scikit-bio
# ---------------------------------------------------------------------------

class TestCrossValidation:
    """Compare ssw-aligner output to scikit-bio for identical inputs."""

    FIELDS = [
        "optimal_alignment_score",
        "suboptimal_alignment_score",
        "query_begin",
        "query_end",
        "target_begin",
        "target_end_optimal",
        "target_end_suboptimal",
        "cigar",
    ]

    def _compare(self, query, target, **params):
        ssw_res = StripedSmithWaterman(query, **params)(target)
        skbio_res = SkbioSSW(query, **params)(target)

        for field in self.FIELDS:
            assert ssw_res[field] == skbio_res[field], (
                f"Mismatch on '{field}': ssw_aligner={ssw_res[field]}, "
                f"skbio={skbio_res[field]}"
            )

    # NT tests
    def test_nt_perfect_match(self):
        self._compare("ACGTACGTACGTACGT", "ACGTACGTACGTACGT")

    def test_nt_partial_match(self):
        self._compare("ACGTACGTACGTACGT", "TACGTACGT")

    def test_nt_with_gap(self):
        self._compare("ACGTACGT", "ACGTAAAAACGT")

    def test_nt_riot_params(self):
        self._compare(
            "CTATACTACTATGGTTCGGGGAGTTATTATAGCCTTGACTACTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAGGGAGTGCATCCGCCCCAACCTCGT",
            "ACTACTTTGACTACTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAG",
            **NT_ALIGNER_PARAMS,
        )

    def test_nt_reverse_complement_scenario(self):
        """Cross-validate a realistic nt scenario."""
        query = "ATGGCCATGGCCCCCAGAACTGAGATCAATAAACAAA"
        target = "ATGGCCATGGCCCCCAGAACTGAG"
        self._compare(query, target, **NT_ALIGNER_PARAMS)

    # AA tests
    def test_aa_perfect_match(self):
        self._compare("ARNDCQEGHILKMFPSTWYV", "ARNDCQEGHILKMFPSTWYV", **AA_ALIGNER_PARAMS)

    def test_aa_partial_match(self):
        self._compare("ARNDCQEGHILKMFPSTWYV", "CQEGHILKMFP", **AA_ALIGNER_PARAMS)

    def test_aa_with_mutations(self):
        self._compare("ARNDCQEGHILKMFPSTWYV", "ARNDAQEGHILKMFASTWYV", **AA_ALIGNER_PARAMS)

    def test_aa_v_gene_fragment(self):
        query = "QVQLVQSGAEVKKPGASVKVSCKASGYTFTGYYMHWVRQAPGQGLEWMGWINPNSGGTNYA"
        target = "QVQLVQSGAEVKKPGASVKVSCKASGYTFTGYYMHWVRQAPGQGLEWMGWINPNSGGTNYAQKFQG"
        self._compare(query, target, **AA_ALIGNER_PARAMS)

    def test_aa_j_gene_fragment(self):
        self._compare("WGQGTLVTVSS", "WGQGTLVTVSSASTKGPS", **AA_ALIGNER_PARAMS)

    def test_aa_full_heavy_chain(self):
        query = (
            "QVQLQQWGAGLLKPSETLSLTCAVFGGSFSGYYWSWIRQPPGKGLEWIGEINHRGNTNDNPSLKS"
            "RVTISVDTSKNQFALKLSSVTAADTAVYYCARERGYTYGNFDHWGQGTLVTVSSASTKGPSVFPL"
            "APSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSL"
            "GTQTYICNVNHKPSNTKVDKKVEPKSCDKTHTCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPE"
            "VTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKV"
            "SNKALPAPIEKTISKAKGQPREPQVYTLPPSRDELTKNQVSLTCLVKGFYPSDIAVEWESNGQPE"
            "NNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPGK"
        )
        self._compare(query, query, **AA_ALIGNER_PARAMS)

    def test_multiple_targets_nt(self):
        """Reusing an aligner with multiple targets yields identical results."""
        query = "ACGTACGTACGTACGT"
        targets = ["ACGTACGT", "TACGTACG", "AAACCCGGG", "ACGTACGTACGTACGT"]

        ssw_aligner = StripedSmithWaterman(query, **NT_ALIGNER_PARAMS)
        skbio_aligner = SkbioSSW(query, **NT_ALIGNER_PARAMS)

        for target in targets:
            ssw_res = ssw_aligner(target)
            skbio_res = skbio_aligner(target)
            for field in self.FIELDS:
                assert ssw_res[field] == skbio_res[field], (
                    f"target={target!r}, field={field!r}: "
                    f"ssw_aligner={ssw_res[field]}, skbio={skbio_res[field]}"
                )

    def test_multiple_targets_aa(self):
        """Reusing an AA aligner with multiple targets yields identical results."""
        query = "QVQLVQSGAEVKKPGASVKVSCKAS"
        targets = [
            "QVQLVQSGAEVKKPGASVKVSCKAS",
            "EVQLVESGGGLVQPGGSLRLSCAAS",
            "QVQLVESGGGVVQPGRSLRLSCAAS",
        ]

        ssw_aligner = StripedSmithWaterman(query, **AA_ALIGNER_PARAMS)
        skbio_aligner = SkbioSSW(query, **AA_ALIGNER_PARAMS)

        for target in targets:
            ssw_res = ssw_aligner(target)
            skbio_res = skbio_aligner(target)
            for field in self.FIELDS:
                assert ssw_res[field] == skbio_res[field], (
                    f"target={target!r}, field={field!r}: "
                    f"ssw_aligner={ssw_res[field]}, skbio={skbio_res[field]}"
                )
