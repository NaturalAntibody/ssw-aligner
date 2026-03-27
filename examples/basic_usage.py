#!/usr/bin/env python
"""
ssw-aligner usage examples
===========================

Demonstrates nucleotide alignment, protein alignment, result inspection,
and Gumbel statistical parameter estimation.

Prerequisites:
    pip install numpy blosum
    cmake -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j$(nproc)
"""

import blosum

from ssw_aligner import (
    AlignmentStructure,
    StripedSmithWaterman,
    compute_gumbel_params,
    MMSEQS_AA_FREQS,
)


# ──────────────────────────────────────────────────────────────────────
# 1. Nucleotide alignment
# ──────────────────────────────────────────────────────────────────────

def nucleotide_example():
    """Basic DNA alignment with default scoring (match=2, mismatch=-3)."""
    query = "ACTGACTGACTG"
    target = "AACTGACTGACTGA"

    aligner = StripedSmithWaterman(query)
    result: AlignmentStructure = aligner(target)

    print("=== Nucleotide alignment ===")
    print(f"Query:  {query}")
    print(f"Target: {target}")
    print(f"Score:  {result.optimal_alignment_score}")
    print(f"CIGAR:  {result.cigar}")
    print(f"Query region:  [{result.query_begin}, {result.query_end}]")
    print(f"Target region: [{result.target_begin}, {result.target_end_optimal}]")
    print(f"Aligned query:  {result.aligned_query_sequence}")
    print(f"Aligned target: {result.aligned_target_sequence}")
    print()


# ──────────────────────────────────────────────────────────────────────
# 2. Nucleotide alignment with custom scoring
# ──────────────────────────────────────────────────────────────────────

def nucleotide_custom_scoring():
    """DNA alignment with riot_na-style scoring parameters."""
    query = "CAGGTGCAGCTGGTGGAGTCTGGG"
    target = "CAGGTGCAGCTGGTGCAGTCTGGG"  # one mismatch

    aligner = StripedSmithWaterman(
        query,
        match_score=1,
        mismatch_score=-1,
        gap_open_penalty=4,
        gap_extend_penalty=1,
    )
    result = aligner(target)

    print("=== Nucleotide (custom scoring) ===")
    print(f"Score: {result.optimal_alignment_score}")
    print(f"CIGAR: {result.cigar}")
    print()


# ──────────────────────────────────────────────────────────────────────
# 3. Protein alignment
# ──────────────────────────────────────────────────────────────────────

def protein_example():
    """Protein alignment with BLOSUM62 substitution matrix."""
    # Antibody VH gene fragment
    query = "EVQLVESGGGLVQPGGSLRLSCAASGFTFS"
    target = "EVQLVESGGGLVKPGGSLRLSCAASGFTFS"  # K vs Q at position 12

    aligner = StripedSmithWaterman(
        query,
        protein=True,
        substitution_matrix=blosum.BLOSUM(62),
        gap_open_penalty=11,
        gap_extend_penalty=1,
    )
    result = aligner(target)

    print("=== Protein alignment (BLOSUM62) ===")
    print(f"Query:  {query}")
    print(f"Target: {target}")
    print(f"Score:  {result.optimal_alignment_score}")
    print(f"CIGAR:  {result.cigar}")
    print(f"Aligned query:  {result.aligned_query_sequence}")
    print(f"Aligned target: {result.aligned_target_sequence}")
    print()


# ──────────────────────────────────────────────────────────────────────
# 4. Reusing an aligner across multiple targets
# ──────────────────────────────────────────────────────────────────────

def multiple_targets():
    """Align one query against several targets efficiently."""
    query = "EVQLVESGGGLVQPGG"

    aligner = StripedSmithWaterman(
        query,
        protein=True,
        substitution_matrix=blosum.BLOSUM(62),
        gap_open_penalty=11,
        gap_extend_penalty=1,
    )

    targets = [
        ("exact",    "EVQLVESGGGLVQPGG"),
        ("mutant",   "EVQLVESGGGLVKPGG"),
        ("partial",  "ESGGGLVQPGGSLRL"),
        ("unrelated", "DIQMTQSPSSLSASVG"),
    ]

    print("=== Multiple targets ===")
    for label, target in targets:
        result = aligner(target)
        print(f"  {label:10s}  score={result.optimal_alignment_score:4d}  cigar={result.cigar}")
    print()


# ──────────────────────────────────────────────────────────────────────
# 5. Accessing results with dict-style syntax
# ──────────────────────────────────────────────────────────────────────

def dict_style_access():
    """AlignmentStructure supports both attribute and dict-style access."""
    aligner = StripedSmithWaterman("ACTGACTG")
    result = aligner("AACTGACTGA")

    # These are equivalent:
    score_attr = result.optimal_alignment_score
    score_dict = result["optimal_alignment_score"]
    assert score_attr == score_dict

    print("=== Dict-style access ===")
    print(f"  result.optimal_alignment_score      = {score_attr}")
    print(f"  result['optimal_alignment_score']    = {score_dict}")
    print(f"  result['cigar']                      = {result['cigar']}")
    print(f"  result['query_begin']                = {result['query_begin']}")
    print()


# ──────────────────────────────────────────────────────────────────────
# 6. Gumbel parameters and E-value computation
# ──────────────────────────────────────────────────────────────────────

def gumbel_params_example():
    """Compute Gumbel statistical parameters and derive E-values."""
    mat = blosum.BLOSUM(62)

    # Compute with MMseqs2-style Robinson-Robinson background frequencies
    params = compute_gumbel_params(
        mat,
        gap_open=11,
        gap_extend=1,
        protein=True,
        use_mmseqs_aa_freqs=True,
    )

    print("=== Gumbel parameters (BLOSUM62, 11/1) ===")
    print(f"  lambda = {params.lambda_:.6f}")
    print(f"  K      = {params.K:.6f}")
    print(f"  sigma  = {params.sigma:.4f}")
    print(f"  tau    = {params.tau:.4f}")
    print()

    # Run an alignment and compute statistics
    query = "EVQLVESGGGLVQPGGSLRLSCAASGFTFS"
    target = "EVQLVESGGGLVKPGGSLRLSCAASGFTFS"
    aligner = StripedSmithWaterman(
        query,
        protein=True,
        substitution_matrix=mat,
        gap_open_penalty=11,
        gap_extend_penalty=1,
    )
    result = aligner(target)

    raw_score = result.optimal_alignment_score
    bit_score = params.bit_score(raw_score)
    db_residues = 1_000_000  # total residues in the database
    e_value = params.evalue(raw_score, len(query), db_residues)

    print(f"  Raw score  = {raw_score}")
    print(f"  Bit score  = {bit_score:.2f}")
    print(f"  E-value    = {e_value:.2e}  (db_size={db_residues:,})")
    print()


# ──────────────────────────────────────────────────────────────────────
# 7. Gumbel parameters with uniform background (default)
# ──────────────────────────────────────────────────────────────────────

def gumbel_uniform_example():
    """Compare Gumbel params with uniform vs MMseqs2 background frequencies."""
    mat = blosum.BLOSUM(62)

    p_uniform = compute_gumbel_params(mat, gap_open=11, gap_extend=1, protein=True)
    p_mmseqs = compute_gumbel_params(
        mat, gap_open=11, gap_extend=1, protein=True, use_mmseqs_aa_freqs=True
    )

    print("=== Background frequency comparison ===")
    print(f"  Uniform:  lambda={p_uniform.lambda_:.6f}  K={p_uniform.K:.6f}")
    print(f"  MMseqs2:  lambda={p_mmseqs.lambda_:.6f}  K={p_mmseqs.K:.6f}")
    print()

    # You can also inspect the Robinson-Robinson frequencies directly
    print("=== MMSEQS_AA_FREQS (top 5 by frequency) ===")
    for aa, freq in sorted(MMSEQS_AA_FREQS.items(), key=lambda x: -x[1])[:5]:
        print(f"  {aa}: {freq:.5f}")
    print()


# ──────────────────────────────────────────────────────────────────────
# 8. Score-only mode (faster, no traceback)
# ──────────────────────────────────────────────────────────────────────

def score_only_example():
    """Use score_only=True when you only need the alignment score."""
    aligner = StripedSmithWaterman(
        "EVQLVESGGGLVQPGGSLRL",
        protein=True,
        substitution_matrix=blosum.BLOSUM(62),
        gap_open_penalty=11,
        gap_extend_penalty=1,
        score_only=True,
    )
    result = aligner("EVQLVESGGGLVQPGGSLRL")

    print("=== Score-only mode ===")
    print(f"  Score: {result.optimal_alignment_score}")
    print(f"  CIGAR: '{result.cigar}'  (empty — traceback not computed)")
    print()


# ──────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    nucleotide_example()
    nucleotide_custom_scoring()
    protein_example()
    multiple_targets()
    dict_style_access()
    gumbel_params_example()
    gumbel_uniform_example()
    score_only_example()
    print("All examples completed successfully.")
