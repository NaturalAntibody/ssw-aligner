"""Equal-score extension of alignment ends (`extend_equal_score`), see docs/equal_score_extension.md."""

import random
import re

import pytest

from ssw_aligner import ProteinAligner, StripedSmithWaterman
from ssw_aligner._blosum_62 import BLOSUM_62
from ssw_aligner.aligners import AA_ALIGNER_PARAMS

# catumaxomab VH after its scFv linker, and the mouse germline IGHV it aligns to. BLOSUM62, gap 11/1: the alignment
# from germline position 4 (95M, score 485) ties with the one from germline position 1 with an insertion
# (EVQ/QVQ +11, insertion -11).
CATUMAXOMAB = (
    "GGSGGGGSGGGGSEVQLLEQSGAELVRPGTSVKISCKASGYAFTNYWLGWVKQRPGHGLEWIGDIFPGSG"
    "NIHYNEKFKGKATLTADKSSSTAYMQLSSLTFEDSAVYFCARLRNWDEPMDYWGQGTTVTVSSGGGGS"
)
GERMLINE = "QVQLQQSGAELVRPGTSVKISCKASGYAFTNYWLGWVKQRPGHGLEWIGDIYPGSGNTYYNEKFKGKATLTADKSSSTAYMQLSSLTSEDSAVYFCAR"

FIELDS = ["optimal_alignment_score", "query_begin", "query_end", "target_begin", "target_end_optimal", "cigar"]


def align(query, target, **kwargs):
    return ProteinAligner(query, **kwargs)(target)


def fields(result):
    return {field: result[field] for field in FIELDS}


def cigar_ops(cigar):
    return [(int(length), op) for length, op in re.findall(r"(\d+)([MID])", cigar)]


def rescore(result, gap_open=11, gap_extend=1):
    """Score of the reported alignment, recomputed from its CIGAR."""
    query, target = result.query_sequence, result.target_sequence
    q, t, score = result.query_begin, result.target_begin, 0
    for length, op in cigar_ops(result.cigar):
        if op == "M":
            score += sum(BLOSUM_62[query[q + k]][target[t + k]] for k in range(length))
            q += length
            t += length
        else:
            score -= gap_open + (length - 1) * gap_extend
            q += length if op == "I" else 0
            t += length if op == "D" else 0
    assert (q - 1, t - 1) == (result.query_end, result.target_end_optimal), "CIGAR disagrees with the positions"
    return score


def test_default_is_unchanged():
    assert fields(align(CATUMAXOMAB, GERMLINE)) == fields(align(CATUMAXOMAB, GERMLINE, extend_equal_score=None))


def test_start_tie_is_extended():
    plain = align(CATUMAXOMAB, GERMLINE)
    assert (plain.query_begin, plain.target_begin, plain.cigar) == (17, 3, "95M")

    extended = align(CATUMAXOMAB, GERMLINE, extend_equal_score="start")
    assert (extended.query_begin, extended.target_begin) == (13, 0)
    assert extended.cigar == "3M1I95M"
    assert extended.optimal_alignment_score == plain.optimal_alignment_score == 485
    assert extended.aligned_query_sequence.startswith("EVQLLEQS")
    assert extended.aligned_target_sequence.startswith("QVQ-LQQS")
    assert rescore(extended) == 485


def test_end_tie_is_extended():
    """The same tie, at the end of the alignment of the reversed sequences."""
    query, target = CATUMAXOMAB[::-1], GERMLINE[::-1]
    plain = align(query, target)
    extended = align(query, target, extend_equal_score="end")

    assert plain.cigar == "95M"
    assert extended.cigar == "95M1I3M"
    assert extended.target_end_optimal == len(target) - 1
    assert extended.query_end == plain.query_end + 4
    assert rescore(extended) == extended.optimal_alignment_score == plain.optimal_alignment_score


def test_start_mode_does_not_extend_the_end():
    query, target = CATUMAXOMAB[::-1], GERMLINE[::-1]
    assert fields(align(query, target, extend_equal_score="start")) == fields(align(query, target))


def test_negative_extension_is_not_taken():
    # K/Q scores 1 less than E/Q: the extension would score -1
    query = CATUMAXOMAB.replace("SEVQLL", "SKVQLL")
    assert fields(align(query, GERMLINE, extend_equal_score="both")) == fields(align(query, GERMLINE))


def test_max_tie_insertions_zero_allows_no_insertion():
    assert fields(align(CATUMAXOMAB, GERMLINE, extend_equal_score="start", max_tie_insertions=0)) == fields(
        align(CATUMAXOMAB, GERMLINE)
    )


def test_max_tie_insertions_limits_inserted_residues():
    # QIQ/QVQ +13 pays for exactly 3 inserted residues (-11 -1 -1)
    query = "GGGGSQIQAAALQQSGAELVRPGTSVKISCKASGYAFTNYWLGWVKQRPGHGLEWIG"
    plain = align(query, GERMLINE)
    assert plain.target_begin == 3

    assert fields(align(query, GERMLINE, extend_equal_score="start")) == fields(plain)  # default: 2

    extended = align(query, GERMLINE, extend_equal_score="start", max_tie_insertions=3)
    assert extended.target_begin == 0
    assert extended.cigar.startswith("3M3I")
    assert rescore(extended) == plain.optimal_alignment_score


def test_extension_stops_at_query_start():
    query = CATUMAXOMAB[CATUMAXOMAB.index("EVQLL") :]
    extended = align(query, GERMLINE, extend_equal_score="start")
    assert (extended.query_begin, extended.target_begin, extended.cigar) == (0, 0, "3M1I95M")

    query = CATUMAXOMAB[CATUMAXOMAB.index("VQLL") :]  # without the E there is no tie
    assert fields(align(query, GERMLINE, extend_equal_score="start")) == fields(align(query, GERMLINE))


def test_nothing_to_extend_when_alignment_covers_target_start():
    query = "GGGGS" + GERMLINE
    assert fields(align(query, GERMLINE, extend_equal_score="both")) == fields(align(query, GERMLINE))


@pytest.mark.parametrize("score_size", [0, 1, 2])
def test_score_sizes(score_size):
    # short enough for 8-bit scores (< 255)
    query, target = CATUMAXOMAB[:45], GERMLINE[:32]
    extended = align(query, target, extend_equal_score="start", score_size=score_size)
    assert (extended.target_begin, extended.cigar[:4]) == (0, "3M1I")
    assert rescore(extended) == extended.optimal_alignment_score


def test_score_only_has_nothing_to_extend():
    result = align(CATUMAXOMAB, GERMLINE, extend_equal_score="both", score_only=True)
    assert result.optimal_alignment_score == 485
    assert result.cigar == ""


def test_invalid_parameters():
    with pytest.raises(ValueError):
        StripedSmithWaterman(CATUMAXOMAB, extend_equal_score="n", **AA_ALIGNER_PARAMS)
    with pytest.raises(ValueError):
        StripedSmithWaterman(CATUMAXOMAB, extend_equal_score="none", **AA_ALIGNER_PARAMS)
    with pytest.raises(ValueError):
        StripedSmithWaterman(CATUMAXOMAB, max_tie_insertions=-1, **AA_ALIGNER_PARAMS)


def mutate(sequence, rng, rate=0.15):
    residues = "ACDEFGHIKLMNPQRSTVWY"
    out = []
    for residue in sequence:
        roll = rng.random()
        if roll < rate * 0.8:
            out.append(rng.choice(residues))
        elif roll < rate * 0.9:
            out.append(residue + rng.choice(residues))
        elif roll >= rate:
            out.append(residue)
    return "".join(out)


def test_random_alignments_keep_score_and_contain_the_original():
    rng = random.Random(0)
    residues = "ACDEFGHIKLMNPQRSTVWY"
    extended_count = 0
    for _ in range(500):
        flank_start = "".join(rng.choice(residues) for _ in range(rng.randint(0, 8)))
        flank_end = "".join(rng.choice(residues) for _ in range(rng.randint(0, 8)))
        core = GERMLINE[rng.randint(0, 10) : len(GERMLINE) - rng.randint(0, 10)]
        query = flank_start + mutate(core, rng) + flank_end
        plain = align(query, GERMLINE)
        extended = align(query, GERMLINE, extend_equal_score="both")

        assert extended.optimal_alignment_score == plain.optimal_alignment_score
        assert rescore(extended) == extended.optimal_alignment_score
        assert extended.query_begin <= plain.query_begin and extended.target_begin <= plain.target_begin
        assert extended.query_end >= plain.query_end and extended.target_end_optimal >= plain.target_end_optimal
        # at most max_tie_insertions (2) query residues beyond the target residues, on each side
        assert (plain.query_begin - extended.query_begin) - (plain.target_begin - extended.target_begin) <= 2
        assert (extended.query_end - plain.query_end) - (extended.target_end_optimal - plain.target_end_optimal) <= 2
        extended_count += fields(extended) != fields(plain)
    assert extended_count > 0
