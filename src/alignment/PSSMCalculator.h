// Stub PSSMCalculator.h - minimal declaration to satisfy include from Sequence.cpp
// The actual PSSMCalculator in MMseqs2 includes Matcher.h which creates a circular
// dependency back to StripedSmithWaterman.h. This stub breaks that cycle.
#ifndef MMSEQS_PSSM_H
#define MMSEQS_PSSM_H

#include <cstddef>
#include <cstdint>
#include <string>

class BaseMatrix;
class Sequence;

class PSSMCalculator {
public:
    struct Profile {
        char *pssm;
        float *prob;
        const float *neffM;
#ifdef GAP_POS_SCORING
        const uint8_t *gDel;
        const uint8_t *gIns;
#endif
        unsigned char *consensus;

#ifdef GAP_POS_SCORING
        Profile(char *pssm, float *prob, float *neffM, const uint8_t *gDel, const uint8_t *gIns, unsigned char *consensus)
            : pssm(pssm), prob(prob), neffM(neffM), gDel(gDel), gIns(gIns), consensus(consensus) {}
#else
        Profile(char *pssm, float *prob, float *neffM, unsigned char *consensus)
            : pssm(pssm), prob(prob), neffM(neffM), consensus(consensus) {}
#endif
    };
};

#endif // MMSEQS_PSSM_H
