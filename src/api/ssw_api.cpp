/*
 * ssw_api.cpp – C adapter that bridges the old ssw.h-style API
 *               to the MMseqs2 SmithWaterman C++ engine.
 */

#include "ssw_api.h"

#include <algorithm>
#include <cstring>
#include <cmath>

#include "BaseMatrix.h"
#include "Parameters.h"
#include "Sequence.h"
#include "StripedSmithWaterman.h"
#include "SubstitutionMatrix.h"
#include "sls_alignment_evaluer.hpp"

/* ------------------------------------------------------------------ */
/*  Minimal BaseMatrix initialised from a flat int8_t scoring matrix  */
/* ------------------------------------------------------------------ */
class SimpleBaseMatrix : public BaseMatrix {
public:
    SimpleBaseMatrix(const int8_t *flatMatrix, int alphSize) {
        alphabetSize = alphSize;
        allocatedAlphabetSize = alphSize;
        initMatrixMemory(alphSize);

        /* Fill subMatrix (short**) from flat int8_t */
        for (int i = 0; i < alphSize; i++) {
            for (int j = 0; j < alphSize; j++) {
                subMatrix[i][j] = flatMatrix[i * alphSize + j];
            }
        }

        /* Uniform background – sufficient for basic SSW without bias correction. */
        for (int i = 0; i < alphSize; i++) {
            pBack[i] = 1.0 / alphSize;
        }

        /* Identity numeric mapping (caller already encoded to 0..n-1). */
        for (int i = 0; i < 256; i++) {
            aa2num[i] = static_cast<unsigned char>(alphSize - 1);
        }
        for (int i = 0; i < alphSize && i < 256; i++) {
            aa2num[i] = static_cast<unsigned char>(i);
            num2aa[i] = '?';
        }

        matrixName = "custom";

        /* Compute lambda from the scoring matrix */
        lambda = 0.0;
        double **dmat = new double *[alphSize];
        for (int i = 0; i < alphSize; i++) {
            dmat[i] = new double[alphSize];
            for (int j = 0; j < alphSize; j++) {
                dmat[i][j] = subMatrix[i][j];
            }
        }
        /* Very basic computation – just use pBack to estimate. */
        double sum = 0.0;
        for (int i = 0; i < alphSize; i++)
            for (int j = 0; j < alphSize; j++)
                sum += pBack[i] * pBack[j] * std::exp(dmat[i][j]);
        if (sum > 0.0 && sum != 1.0)
            lambda = std::log(sum);

        for (int i = 0; i < alphSize; i++)
            delete[] dmat[i];
        delete[] dmat;
    }

    ~SimpleBaseMatrix() override = default;
};

/* ------------------------------------------------------------------ */
/*  Opaque handle                                                     */
/* ------------------------------------------------------------------ */
struct ssw_handle {
    SmithWaterman  *sw;
    SimpleBaseMatrix *baseMat;
    Sequence       *seq;
    int8_t         *mat;          /* copy of flat matrix */
    int             alphabetSize;
    int             queryLen;
};

/* ------------------------------------------------------------------ */
/*  C API                                                             */
/* ------------------------------------------------------------------ */
extern "C" {

ssw_handle *ssw_init(const int8_t *read_num,
                     int32_t       readLen,
                     const int8_t *mat,
                     int32_t       n,
                     int8_t /* score_size – ignored */) {
    auto *h = new ssw_handle();
    h->alphabetSize = n;
    h->queryLen     = readLen;

    /* Copy flat substitution matrix */
    h->mat = new int8_t[n * n];
    std::memcpy(h->mat, mat, n * n * sizeof(int8_t));

    /* Build a BaseMatrix from the flat array */
    h->baseMat = new SimpleBaseMatrix(mat, n);

    /* Allocate SmithWaterman –
     * maxSequenceLength must cover both query and any future target.
     * We pick a generous default here; ssw_align will check. */
    size_t maxLen = std::max(static_cast<size_t>(readLen), static_cast<size_t>(4096));
    /* aaBiasCorrection = false  (the old C ssw never did it) */
    h->sw = new SmithWaterman(maxLen, n,
                              /*aaBiasCorrection=*/false,
                              /*aaBiasCorrectionScale=*/1.0f,
                              reinterpret_cast<SubstitutionMatrix *>(h->baseMat));

    /* Build a Sequence and populate with the query's numeric encoding */
    h->seq = new Sequence(static_cast<size_t>(readLen) + 1,
                          Parameters::DBTYPE_AMINO_ACIDS,  /* works for NT too */
                          h->baseMat,
                          /*kmerSize=*/0,
                          /*spaced=*/false,
                          /*aaBiasCorrection=*/false,
                          /*shouldAddPC=*/false);

    /* Directly fill numSequence (already numerically encoded by caller) */
    for (int32_t i = 0; i < readLen; i++) {
        h->seq->numSequence[i] = static_cast<unsigned char>(read_num[i]);
    }
    h->seq->L = readLen;

    /* Initialise the SIMD query profile */
    h->sw->ssw_init(h->seq, h->mat, h->baseMat);

    return h;
}

ssw_handle *ssw_init_profile(const int8_t *read_num,
                             int32_t       readLen,
                             const int8_t *pssm,
                             int32_t       n) {
    auto *h = new ssw_handle();
    h->queryLen     = readLen;

    /*
     * The engine's ssw_init() for profile mode zeroes out the last row
     * (index alphabetSize-1) in the PSSM as the neutral 'X' state.
     * MMseqs2 uses alphabetSize = 21 (20 AAs + X), so the zeroing
     * targets row 20 — a padding row that is not part of the real PSSM.
     *
     * We must mirror this: use n+1 as the internal alphabet size so
     * that the engine zeroes row n (the extra X row) instead of row
     * n-1 (the last real amino acid, Valine).
     */
    const int32_t internalAlphSize = n + 1;
    h->alphabetSize = internalAlphSize;

    /* Allocate PSSM with one extra row for the X state (zeroed). */
    int32_t pssmSize = internalAlphSize * readLen;
    h->mat = new int8_t[pssmSize];
    std::memcpy(h->mat, pssm, n * readLen * sizeof(int8_t));
    /* Zero the extra X row — engine will also zero it, but be explicit. */
    std::memset(h->mat + n * readLen, 0, readLen * sizeof(int8_t));

    /* We still need a BaseMatrix for internal engine bookkeeping.
     * Build a minimal identity-like matrix – it is not used for scoring
     * in profile mode, but the engine needs alphabetSize etc. */
    int8_t *identMat = new int8_t[internalAlphSize * internalAlphSize];
    std::memset(identMat, 0, internalAlphSize * internalAlphSize * sizeof(int8_t));
    for (int i = 0; i < internalAlphSize; i++)
        identMat[i * internalAlphSize + i] = 1;
    h->baseMat = new SimpleBaseMatrix(identMat, internalAlphSize);
    delete[] identMat;

    size_t maxLen = std::max(static_cast<size_t>(readLen), static_cast<size_t>(4096));
    h->sw = new SmithWaterman(maxLen, internalAlphSize,
                              /*aaBiasCorrection=*/false,
                              /*aaBiasCorrectionScale=*/1.0f,
                              reinterpret_cast<SubstitutionMatrix *>(h->baseMat));

    /* Build a Sequence with HMM_PROFILE type so the engine uses PROFILE_SEQ. */
    h->seq = new Sequence(static_cast<size_t>(readLen) + 1,
                          Parameters::DBTYPE_HMM_PROFILE,
                          h->baseMat,
                          /*kmerSize=*/0,
                          /*spaced=*/false,
                          /*aaBiasCorrection=*/false,
                          /*shouldAddPC=*/false);

    for (int32_t i = 0; i < readLen; i++) {
        h->seq->numSequence[i] = static_cast<unsigned char>(read_num[i]);
    }
    h->seq->L = readLen;

    /* Initialise the SIMD query profile with the PSSM */
    h->sw->ssw_init(h->seq, h->mat, h->baseMat);

    return h;
}

ssw_result ssw_align(ssw_handle   *handle,
                     const int8_t *ref_num,
                     int32_t       refLen,
                     uint8_t       gap_open,
                     uint8_t       gap_extend,
                     uint8_t       flag,
                     uint16_t   /* filters – not used by new engine */,
                     int32_t    /* filterd – not used by new engine */,
                     int32_t       maskLen) {

    /* Map old bit-flag to MMseqs2 alignmentMode:
     *   0 → 0  (score + end positions only)
     *   anything else → 3  (always compute start positions + cigar)
     */
    uint8_t alignmentMode = (flag == 0) ? 0 : 3;

    s_align a = handle->sw->ssw_align(
        reinterpret_cast<const unsigned char *>(ref_num),
        refLen,
        gap_open,
        gap_extend,
        alignmentMode,
        maskLen);

    /* Translate to ssw_result */
    ssw_result r;
    r.score1     = a.score1;
    r.score2     = a.score2;
    r.ref_begin1 = a.dbStartPos1;
    r.ref_end1   = a.dbEndPos1;
    r.read_begin1= a.qStartPos1;
    r.read_end1  = a.qEndPos1;
    r.ref_end2   = a.ref_end2;
    r.cigar      = a.cigar;      /* ownership transferred to caller */
    r.cigarLen   = a.cigarLen;
    return r;
}

void ssw_free_cigar(uint32_t *cigar) {
    if (cigar) {
        delete[] cigar;
    }
}

void ssw_destroy(ssw_handle *handle) {
    if (!handle) return;
    delete handle->sw;
    delete handle->seq;
    delete handle->baseMat;
    delete[] handle->mat;
    delete handle;
}

gumbel_params compute_gumbel_params(const int8_t  *mat,
                                    int32_t        n,
                                    const double  *bg_freqs,
                                    int32_t        gap_open,
                                    int32_t        gap_extend,
                                    double         max_seconds) {
    gumbel_params result;
    std::memset(&result, 0, sizeof(result));
    result.valid = 0;

    /* Build temporary long** matrix from flat int8_t (ALP requires long**) */
    long **tmpMat = new long *[n];
    long *tmpMatData = new long[n * n];
    for (int i = 0; i < n; i++) {
        tmpMat[i] = &tmpMatData[i * n];
        for (int j = 0; j < n; j++) {
            tmpMat[i][j] = static_cast<long>(mat[i * n + j]);
        }
    }

    /* Background frequencies – use supplied or uniform */
    double *freqs = nullptr;
    bool ownFreqs = false;
    if (bg_freqs != nullptr) {
        freqs = const_cast<double *>(bg_freqs);
    } else {
        freqs = new double[n];
        ownFreqs = true;
        for (int i = 0; i < n; i++) {
            freqs[i] = 1.0 / n;
        }
    }

    const double lambdaTolerance = 0.01;
    const double kTolerance = 0.05;
    const double maxMegabytes = 500;
    const long randomSeed = 42;

    Sls::AlignmentEvaluer evaluer;

    bool isGapped = (gap_open > 0 || gap_extend > 0);
    /* Use full alphabet size – caller is responsible for passing only
       the symbols that should participate in the simulation. */
    int alpSize = n;

    try {
        if (isGapped) {
            evaluer.initGapped(
                alpSize,
                (const long *const *)tmpMat,
                freqs, freqs,
                gap_open, gap_extend,
                gap_open, gap_extend,
                false,
                lambdaTolerance, kTolerance,
                max_seconds, maxMegabytes, randomSeed);
        } else {
            evaluer.initGapless(
                alpSize,
                (const long *const *)tmpMat,
                freqs, freqs,
                max_seconds);
        }

        if (evaluer.isGood()) {
            const Sls::ALP_set_of_parameters &p = evaluer.parameters();
            result.lambda   = p.lambda;
            result.K        = p.K;
            result.a_I      = p.a_I;
            result.b_I      = p.b_I;
            result.a_J      = p.a_J;
            result.b_J      = p.b_J;
            result.alpha_I  = p.alpha_I;
            result.beta_I   = p.beta_I;
            result.alpha_J  = p.alpha_J;
            result.beta_J   = p.beta_J;
            result.sigma    = p.sigma;
            result.tau      = p.tau;
            result.valid    = 1;
        }
    } catch (...) {
        /* ALP can throw on convergence failure – return invalid */
        result.valid = 0;
    }

    delete[] tmpMatData;
    delete[] tmpMat;
    if (ownFreqs) delete[] freqs;

    return result;
}

}  /* extern "C" */
