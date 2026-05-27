from ssw_aligner import NucleotideAligner, ProteinAligner, StripedSmithWaterman
from ssw_aligner.aligners import AA_ALIGNER_PARAMS, ALIGNER_PARAMS


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


def assert_alignment_fields_match(wrapper_result, base_result):
    for field in FIELDS:
        assert wrapper_result[field] == base_result[field]


def test_nucleotide_aligner_matches_riot_defaults():
    query = "CTATACTACTATGGTTCGGGGAGTTATTATAGCCTTGACTACTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAG"
    target = "ACTACTTTGACTACTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAG"

    wrapper = NucleotideAligner(query_sequence=query)
    base = StripedSmithWaterman(query_sequence=query, **ALIGNER_PARAMS)

    assert_alignment_fields_match(
        wrapper_result=wrapper(target_sequence=target),
        base_result=base(target),
    )


def test_protein_aligner_matches_riot_defaults():
    query = "QVQLVQSGAEVKKPGASVKVSCKASGYTFTGYYMHWVRQAPGQGLEWMGWINPNSGGTNYA"
    target = "QVQLVQSGAEVKKPGASVKVSCKASGYTFTGYYMHWVRQAPGQGLEWMGWINPNSGGTNYAQKFQG"

    wrapper = ProteinAligner(query_sequence=query)
    base = StripedSmithWaterman(query_sequence=query, **AA_ALIGNER_PARAMS)

    assert_alignment_fields_match(
        wrapper_result=wrapper(target_sequence=target),
        base_result=base(target),
    )


def test_nucleotide_aligner_allows_overrides():
    query = "AAAAAAAAAA"
    target = "CCCCCCCCCC"
    override_params = {"match_score": 2, "mismatch_score": -3}

    wrapper = NucleotideAligner(query_sequence=query, **override_params)
    base = StripedSmithWaterman(query_sequence=query, gap_open_penalty=4, gap_extend_penalty=1, **override_params)

    assert_alignment_fields_match(
        wrapper_result=wrapper(target_sequence=target),
        base_result=base(target),
    )


def test_protein_aligner_allows_gap_overrides():
    query = "ARNDCQEGHILKMFPSTWYV"
    target = "ARNDAQEGHILKMFASTWYV"
    override_params = {"gap_open_penalty": 9}

    wrapper = ProteinAligner(query_sequence=query, **override_params)
    base_params = dict(AA_ALIGNER_PARAMS)
    base_params.update(override_params)
    base = StripedSmithWaterman(query_sequence=query, **base_params)

    assert_alignment_fields_match(
        wrapper_result=wrapper(target_sequence=target),
        base_result=base(target),
    )