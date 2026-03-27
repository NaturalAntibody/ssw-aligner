/*
 * ssw_api.h – C-compatible wrapper around the MMseqs2 SmithWaterman engine.
 *
 * This API is deliberately modelled after the original SSW library (ssw.h)
 * so that existing Python/Cython wrappers can be ported with minimal changes.
 *
 * Typical usage:
 *   1. ssw_init()   – build a query profile
 *   2. ssw_align()  – align one or more targets against that profile
 *   3. ssw_free_cigar() – free cigar arrays when done
 *   4. ssw_destroy() – release the context
 */

#ifndef SSW_API_H
#define SSW_API_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Opaque handle returned by ssw_init. */
typedef struct ssw_handle ssw_handle;

/* Alignment result – field names mirror the original ssw.h s_align. */
typedef struct {
    uint32_t score1;        /* optimal alignment score                       */
    uint32_t score2;        /* sub-optimal alignment score                   */
    int32_t  ref_begin1;    /* 0-based target begin  (-1 if not computed)    */
    int32_t  ref_end1;      /* 0-based target end                           */
    int32_t  read_begin1;   /* 0-based query begin   (-1 if not computed)    */
    int32_t  read_end1;     /* 0-based query end                            */
    int32_t  ref_end2;      /* sub-optimal target end                       */
    uint32_t *cigar;        /* BAM-style packed cigar (NULL if not computed) */
    int32_t  cigarLen;      /* number of cigar operations                   */
} ssw_result;

/**
 * Build a query profile from a flat substitution matrix (sequence-sequence).
 *
 * @param read_num   Numerically-encoded query sequence (values 0..n-1).
 * @param readLen    Length of the query.
 * @param mat        Flat n×n substitution matrix (row-major int8_t).
 * @param n          Alphabet size (matrix width); typically 5 for NT, 20 for AA.
 * @param score_size Ignored (kept for API compat); the engine auto-selects.
 * @return           Opaque handle – must be freed with ssw_destroy().
 */
ssw_handle *ssw_init(const int8_t  *read_num,
                     int32_t        readLen,
                     const int8_t  *mat,
                     int32_t        n,
                     int8_t         score_size);

/**
 * Build a query profile from a position-specific scoring matrix (PSSM).
 *
 * The PSSM is a flat int8_t array of shape [n × readLen] (row-major),
 * where n is the alphabet size (typically 20 for amino acids) and readLen
 * is the query length.  Element pssm[aa * readLen + pos] gives the score
 * for amino acid aa at query position pos.
 *
 * @param read_num   Numerically-encoded query sequence (values 0..n-1).
 *                   Used for consensus/traceback; the PSSM carries scoring.
 * @param readLen    Length of the query.
 * @param pssm       Flat n × readLen position-specific scoring matrix.
 * @param n          Alphabet size (number of PSSM rows); typically 20.
 * @return           Opaque handle – must be freed with ssw_destroy().
 */
ssw_handle *ssw_init_profile(const int8_t  *read_num,
                             int32_t        readLen,
                             const int8_t  *pssm,
                             int32_t        n);

/**
 * Align a target against the stored query profile.
 *
 * @param handle     Handle returned by ssw_init().
 * @param ref_num    Numerically-encoded target sequence.
 * @param refLen     Length of the target.
 * @param gap_open   Gap-open penalty (absolute value, > 0).
 * @param gap_extend Gap-extend penalty (absolute value, > 0).
 * @param flag       Bit-flag controlling output (old ssw.h convention):
 *                     0 = score + end positions only;
 *                     non-zero = compute start positions + cigar.
 * @param filters    Score filter threshold (used when flag & 0x02).
 * @param filterd    Distance filter threshold (used when flag & 0x04).
 * @param maskLen    Mask length for sub-optimal score (>= 15 recommended).
 * @return           ssw_result by value. Caller owns cigar – free with
 *                   ssw_free_cigar().
 */
ssw_result ssw_align(ssw_handle    *handle,
                     const int8_t  *ref_num,
                     int32_t        refLen,
                     uint8_t        gap_open,
                     uint8_t        gap_extend,
                     uint8_t        flag,
                     uint16_t       filters,
                     int32_t        filterd,
                     int32_t        maskLen);

/**
 * Free a cigar array returned in ssw_result.
 * Safe to call with NULL.
 */
void ssw_free_cigar(uint32_t *cigar);

/**
 * Destroy a handle created by ssw_init().
 * Safe to call with NULL.
 */
void ssw_destroy(ssw_handle *handle);

/* ------------------------------------------------------------------ */
/*  Gumbel statistical parameters (ALP library)                       */
/* ------------------------------------------------------------------ */

/**
 * Gumbel parameters for E-value computation. These are the 12 parameters
 * from the ALP (Ascending Ladder Points) Monte Carlo simulation that
 * describe the score distribution of gapped local alignments.
 *
 * E-value ≈ K * exp(-lambda * score) * area(query_len, db_size)
 * bit_score = (lambda * score - ln(K)) / ln(2)
 */
typedef struct {
    double lambda;       /* Gumbel distribution scale parameter           */
    double K;            /* Gumbel distribution prefactor                 */
    double a_I;          /* Length adjustment slope for sequence I         */
    double b_I;          /* Length adjustment intercept for sequence I     */
    double a_J;          /* Length adjustment slope for sequence J         */
    double b_J;          /* Length adjustment intercept for sequence J     */
    double alpha_I;      /* Variance slope for sequence I                 */
    double beta_I;       /* Variance intercept for sequence I             */
    double alpha_J;      /* Variance slope for sequence J                 */
    double beta_J;       /* Variance intercept for sequence J             */
    double sigma;        /* Aggregate variance slope                      */
    double tau;          /* Aggregate variance intercept                  */
    int    valid;        /* 1 if computation succeeded, 0 otherwise       */
} gumbel_params;

/**
 * Compute Gumbel statistical parameters for a given substitution matrix
 * and gap penalties using the ALP library.
 *
 * @param mat            Flat n×n substitution matrix (row-major int8_t).
 * @param n              Alphabet size (matrix width).
 * @param bg_freqs       Background frequencies for each symbol (length n).
 *                       If NULL, uniform frequencies (1/n) are used.
 * @param gap_open       Gap opening penalty (>= 0).
 * @param gap_extend     Gap extension penalty (>= 0).
 * @param max_seconds    Maximum time for Monte Carlo simulation (60.0 typical).
 * @return               gumbel_params struct. Check .valid == 1 before use.
 */
gumbel_params compute_gumbel_params(const int8_t  *mat,
                                    int32_t        n,
                                    const double  *bg_freqs,
                                    int32_t        gap_open,
                                    int32_t        gap_extend,
                                    double         max_seconds);

#ifdef __cplusplus
}  /* extern "C" */
#endif

#endif /* SSW_API_H */
