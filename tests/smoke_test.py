"""Smoke test exercising the ssw_aligner public API.

This script is executed inside a clean virtual environment by the
integration tests in test_install.py.  It prints a JSON summary to
stdout so the calling test can assert on the results.
"""

import json

import numpy as np

from ssw_aligner import (
    MMSEQS_AA_FREQS,
    AlignmentStructure,
    SmithWatermanProfileAligner,
    StripedSmithWaterman,
    compute_gumbel_params,
)

results = {}

# --- 1. Nucleotide alignment ---
query = "ACGTACGTACGT"
target = "TACGTACGTACGTAA"
aligner = StripedSmithWaterman(
    query,
    match_score=2,
    mismatch_score=-1,
    gap_open_penalty=3,
    gap_extend_penalty=1,
)
aln = aligner(target)
results["nt_score"] = aln.optimal_alignment_score
results["nt_query_begin"] = aln.query_begin
results["nt_cigar"] = aln.cigar
results["nt_has_aligned"] = bool(aln.aligned_query_sequence)

# --- 2. Protein alignment (identity matrix, score 4 on diagonal) ---
aa_order = "ARNDCQEGHILKMFPSTWYV"
n = len(aa_order)
flat = [0] * (n * n)
for i in range(n):
    flat[i * n + i] = 4
mat = {}
for i, a in enumerate(aa_order):
    mat[a] = {}
    for j, b in enumerate(aa_order):
        mat[a][b] = flat[i * n + j]
prot_q = "ARNDCQEGHILKMFPSTWYV"
prot_t = "ARNDCQEGHILKMFPSTWYV"
prot_aligner = StripedSmithWaterman(
    prot_q,
    protein=True,
    substitution_matrix=mat,
    gap_open_penalty=11,
    gap_extend_penalty=1,
)
prot_aln = prot_aligner(prot_t)
results["prot_score"] = prot_aln.optimal_alignment_score
results["prot_cigar"] = prot_aln.cigar

# --- 3. Profile alignment ---
pssm = np.zeros((20, len(prot_q)), dtype=np.int8)
for col, aa in enumerate(prot_q):
    row = aa_order.index(aa)
    pssm[row, col] = 4
profile_aligner = SmithWatermanProfileAligner(
    prot_q,
    pssm,
    gap_open_penalty=11,
    gap_extend_penalty=1,
)
prof_aln = profile_aligner(prot_t)
results["profile_score"] = prof_aln.optimal_alignment_score
results["profile_scores_match"] = (
    prot_aln.optimal_alignment_score == prof_aln.optimal_alignment_score
)

# --- 4. MMSEQS_AA_FREQS ---
results["freqs_count"] = len(MMSEQS_AA_FREQS)

# --- 5. Version ---
import ssw_aligner

results["version"] = ssw_aligner.__version__

print(json.dumps(results))
