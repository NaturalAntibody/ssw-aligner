#!/usr/bin/env python
"""
Profile (PSSM) alignment examples
==================================

Demonstrates :class:`SmithWatermanProfileAligner`, which aligns target
sequences against a query using a **position-specific scoring matrix**
(PSSM) rather than a single global substitution matrix.

In profile alignment every query position can have its own set of
scores for each amino acid, enabling use-cases such as:

  * HMM profile → sequence alignment
  * MSA-derived PSSMs (e.g. PSI-BLAST iterations)
  * Custom position-specific penalties (e.g. CDR vs framework regions)

Prerequisites:
    pip install numpy blosum
    cmake -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j$(nproc)
"""

import numpy as np
import blosum

from ssw_aligner import SmithWatermanProfileAligner

# Amino acid order used internally by ssw-aligner.
AA_ORDER = "ARNDCQEGHILKMFPSTWYV"


# ──────────────────────────────────────────────────────────────────────
# Helper: Convert a substitution matrix to a flat PSSM for a query
# ──────────────────────────────────────────────────────────────────────

def substitution_matrix_to_pssm(
    query: str,
    sub_matrix: dict,
    aa_order: str = AA_ORDER,
) -> np.ndarray:
    """Build a PSSM of shape (20, query_length) from a substitution matrix.

    This is equivalent to using StripedSmithWaterman with the same matrix,
    but demonstrates the PSSM interface.  In practice the PSSM would come
    from a profile HMM, MSA, or PSI-BLAST.

    Parameters
    ----------
    query : str
        Amino-acid query sequence.
    sub_matrix : dict
        A ``{row: {col: score}}`` substitution matrix (e.g. blosum.BLOSUM(62)).
    aa_order : str
        Alphabet order matching ssw-aligner internals.

    Returns
    -------
    np.ndarray
        PSSM of shape ``(len(aa_order), len(query))`` with dtype ``int8``.
    """
    n = len(aa_order)
    qlen = len(query)
    pssm = np.zeros((n, qlen), dtype=np.int8)
    for aa_idx, aa in enumerate(aa_order):
        for pos, qaa in enumerate(query):
            pssm[aa_idx, pos] = int(sub_matrix[aa][qaa])
    return pssm


# ──────────────────────────────────────────────────────────────────────
# 1. Basic profile alignment (PSSM derived from BLOSUM62)
# ──────────────────────────────────────────────────────────────────────

def basic_profile_alignment():
    """Align using a PSSM built from BLOSUM62 — equivalent to seq-seq."""
    query = "EVQLVESGGGLVQPGGSLRLSCAASGFTFS"
    target = "EVQLVESGGGLVKPGGSLRLSCAASGFTFS"  # K vs Q at pos 12

    mat = blosum.BLOSUM(62)
    pssm = substitution_matrix_to_pssm(query, mat)

    aligner = SmithWatermanProfileAligner(
        query, pssm,
        gap_open_penalty=11,
        gap_extend_penalty=1,
    )
    result = aligner(target)

    print("=== Basic profile alignment (BLOSUM62 → PSSM) ===")
    print(f"Query:  {query}")
    print(f"Target: {target}")
    print(f"Score:  {result.optimal_alignment_score}")
    print(f"CIGAR:  {result.cigar}")
    print(f"Aligned query:  {result.aligned_query_sequence}")
    print(f"Aligned target: {result.aligned_target_sequence}")
    print()


# ──────────────────────────────────────────────────────────────────────
# 2. Position-specific scoring (boosted CDR positions)
# ──────────────────────────────────────────────────────────────────────

def position_specific_scores():
    """Boost scores at specific positions (e.g. CDR residues).

    This is a key advantage of profile alignment: you can give different
    weights to different positions.  Here we double the scores at CDR-like
    positions to penalise mutations in those regions more heavily.
    """
    query = "EVQLVESGGGLVQPGGSLRLSCAASGFTFS"
    mat = blosum.BLOSUM(62)

    # Start with a standard BLOSUM62-derived PSSM
    pssm = substitution_matrix_to_pssm(query, mat)

    # "CDR" positions (0-indexed): boost scores by 2×
    cdr_positions = list(range(24, 29))  # last 5 residues as mock CDR
    pssm[:, cdr_positions] = np.clip(pssm[:, cdr_positions] * 2, -128, 127)

    aligner = SmithWatermanProfileAligner(
        query, pssm,
        gap_open_penalty=11,
        gap_extend_penalty=1,
    )

    targets = [
        ("match_fw",  "EVQLVESGGGLVKPGGSLRLSCAASGFTFS"),  # mismatch in framework
        ("match_cdr", "EVQLVESGGGLVQPGGSLRLSCAASAAAAA"),  # mismatch in CDR
    ]

    print("=== Position-specific scoring (boosted CDR) ===")
    for label, target in targets:
        result = aligner(target)
        print(f"  {label:12s}  score={result.optimal_alignment_score:4d}  cigar={result.cigar}")
    print("  (CDR mutations are penalised more heavily)")
    print()


# ──────────────────────────────────────────────────────────────────────
# 3. MSA-derived PSSM (mock frequency-based profile)
# ──────────────────────────────────────────────────────────────────────

def msa_derived_pssm():
    """Build a PSSM from a mock multiple sequence alignment.

    In practice you would derive this from a real MSA using log-odds
    scoring.  Here we demonstrate the concept with a simple frequency
    count → rounded log-odds approach.
    """
    # Mock MSA of VH germline-like sequences (just first 10 positions)
    msa = [
        "EVQLVESGGS",
        "EVQLVESGGS",
        "EVQLVQSGGS",  # Q→Q at pos 4 in 2 seqs, but also gap variation
        "EVQLVESGGS",
        "EVQLLESGGS",  # V→L at pos 4
    ]

    query = msa[0]
    qlen = len(query)
    n = len(AA_ORDER)

    # Count amino acid frequencies at each position
    freq = np.zeros((n, qlen), dtype=np.float64)
    for seq in msa:
        for pos, aa in enumerate(seq):
            idx = AA_ORDER.index(aa)
            freq[idx, pos] += 1
    freq /= len(msa)

    # Convert to log-odds scores (vs uniform background)
    bg = 1.0 / n
    with np.errstate(divide="ignore"):
        pssm_float = np.log2(freq / bg)
    pssm_float[~np.isfinite(pssm_float)] = -4  # floor for unseen AAs

    # Round and clip to int8
    pssm = np.clip(np.round(pssm_float), -128, 127).astype(np.int8)

    print("=== MSA-derived PSSM ===")
    print(f"  Query: {query}")
    print(f"  PSSM (first 6 columns):")
    for aa_idx in range(n):
        scores = " ".join(f"{pssm[aa_idx, p]:3d}" for p in range(min(6, qlen)))
        print(f"    {AA_ORDER[aa_idx]}: {scores}")

    aligner = SmithWatermanProfileAligner(
        query, pssm,
        gap_open_penalty=5,
        gap_extend_penalty=1,
    )

    test_seqs = [
        ("consensus", "EVQLVESGGS"),
        ("L_at_5",    "EVQLLESGGS"),
        ("mutant",    "DVKLVESGGS"),
    ]

    for label, target in test_seqs:
        result = aligner(target)
        print(f"  {label:12s}  score={result.optimal_alignment_score:4d}  cigar={result.cigar}")
    print()


# ──────────────────────────────────────────────────────────────────────
# 4. Reusing a profile aligner across many targets
# ──────────────────────────────────────────────────────────────────────

def multiple_targets():
    """Align one query profile against a database of targets."""
    query = "EVQLVESGGGLVQPGG"
    mat = blosum.BLOSUM(62)
    pssm = substitution_matrix_to_pssm(query, mat)

    aligner = SmithWatermanProfileAligner(
        query, pssm,
        gap_open_penalty=11,
        gap_extend_penalty=1,
    )

    database = [
        ("VH3-23",    "EVQLVESGGGLVQPGGSLRLSCAAS"),
        ("VH1-69",    "EVQLVQSGAEVKKPGASVKVSCKAS"),
        ("VH4-34",    "QVQLQESGPGLVKPSETLSLTCTVS"),
        ("unrelated", "DIQMTQSPSSLSASVGDRVTITC"),
    ]

    print("=== Profile alignment against database ===")
    for label, target in database:
        result = aligner(target)
        print(
            f"  {label:12s}  score={result.optimal_alignment_score:4d}  "
            f"cigar={result.cigar}  "
            f"target=[{result.target_begin}:{result.target_end_optimal}]"
        )
    print()


# ──────────────────────────────────────────────────────────────────────
# 5. Score-only mode (faster, no traceback)
# ──────────────────────────────────────────────────────────────────────

def score_only_profile():
    """When you only need scores, skip the traceback for speed."""
    query = "EVQLVESGGGLVQPGG"
    mat = blosum.BLOSUM(62)
    pssm = substitution_matrix_to_pssm(query, mat)

    aligner = SmithWatermanProfileAligner(
        query, pssm,
        gap_open_penalty=11,
        gap_extend_penalty=1,
        score_only=True,
    )

    result = aligner("EVQLVESGGGLVQPGG")
    print("=== Profile score-only mode ===")
    print(f"  Score: {result.optimal_alignment_score}")
    print(f"  CIGAR: '{result.cigar}'  (empty — traceback not computed)")
    print()


# ──────────────────────────────────────────────────────────────────────
# 6. Comparing sequence-level and profile-level alignment scores
# ──────────────────────────────────────────────────────────────────────

def compare_seq_vs_profile():
    """Show that a BLOSUM62-derived PSSM gives the same score as seq-seq."""
    from ssw_aligner import StripedSmithWaterman

    query = "EVQLVESGGGLVQPGGSLRLSCAASGFTFS"
    target = "EVQLVESGGGLVKPGGSLRLSCAASGFTFS"
    mat = blosum.BLOSUM(62)

    # Sequence-sequence alignment
    seq_aligner = StripedSmithWaterman(
        query,
        protein=True,
        substitution_matrix=mat,
        gap_open_penalty=11,
        gap_extend_penalty=1,
    )
    seq_result = seq_aligner(target)

    # Profile alignment with same matrix
    pssm = substitution_matrix_to_pssm(query, mat)
    prof_aligner = SmithWatermanProfileAligner(
        query, pssm,
        gap_open_penalty=11,
        gap_extend_penalty=1,
    )
    prof_result = prof_aligner(target)

    print("=== Sequence vs Profile alignment comparison ===")
    print(f"  Seq-seq  score: {seq_result.optimal_alignment_score}")
    print(f"  Profile  score: {prof_result.optimal_alignment_score}")
    print(f"  Scores match:   {seq_result.optimal_alignment_score == prof_result.optimal_alignment_score}")
    print()


# ──────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    basic_profile_alignment()
    position_specific_scores()
    msa_derived_pssm()
    multiple_targets()
    score_only_profile()
    compare_seq_vs_profile()
    print("All profile examples completed successfully.")
