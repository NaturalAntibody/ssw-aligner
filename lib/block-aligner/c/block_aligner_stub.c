/* Stub implementation for block_aligner.
   All functions are no-ops or return failure values, causing StripedSmithWaterman
   to fall back to the standard banded SW traceback path. */

#include "block_aligner.h"
#include <string.h>

/* Internal stub structures - minimal to satisfy alloc/free */
struct AAMatrix { int8_t dummy; };
struct AAProfile {
    int8_t* pos_aa;
    int16_t* aa_pos;
    size_t curr_len;
    size_t str_len;
    int8_t gap_extend;
};
struct Cigar { size_t len; };
struct PaddedBytes { int8_t dummy; };
struct PosBias { int8_t dummy; };

/* Single static block handle - block alignment always "fails" */
static int _block_dummy = 0;

/* --- AAMatrix stubs --- */
struct AAMatrix *block_new_simple_aamatrix(int8_t match_score, int8_t mismatch_score) {
    (void)match_score; (void)mismatch_score;
    return (struct AAMatrix*)calloc(1, sizeof(struct AAMatrix));
}
void block_set_aamatrix(struct AAMatrix *m, uint8_t a, uint8_t b, int8_t score) { (void)m;(void)a;(void)b;(void)score; }
void block_set_aamatrix_num(struct AAMatrix *m, int8_t a, int8_t b, int8_t score) { (void)m;(void)a;(void)b;(void)score; }
void block_free_aamatrix(struct AAMatrix *m) { free(m); }

/* --- AAProfile stubs --- */
struct AAProfile *block_new_aaprofile(uintptr_t str_len, uintptr_t block_size, int8_t gap_extend) {
    (void)block_size;
    struct AAProfile* p = (struct AAProfile*)calloc(1, sizeof(struct AAProfile));
    if (p) {
        /* +2 for padding on both sides, *32 for alphabet stride */
        p->str_len = str_len;
        p->curr_len = str_len + 2;
        p->gap_extend = gap_extend;
        p->pos_aa = (int8_t*)calloc((str_len + 2) * 32, sizeof(int8_t));
        p->aa_pos = (int16_t*)calloc(32 * (str_len + 2), sizeof(int16_t));
    }
    return p;
}
uintptr_t block_len_aaprofile(const struct AAProfile *profile) { return profile ? profile->str_len : 0; }
void block_clear_aaprofile(struct AAProfile *profile, uintptr_t str_len) { (void)profile; (void)str_len; }
void block_set_aaprofile(struct AAProfile *profile, uintptr_t i, uint8_t b, int8_t score) { (void)profile;(void)i;(void)b;(void)score; }
void block_set_all_aaprofile(struct AAProfile *p, const uint8_t *order, uintptr_t order_len,
                             const int8_t *scores, uintptr_t scores_len,
                             uintptr_t left_shift, uintptr_t right_shift) {
    (void)p;(void)order;(void)order_len;(void)scores;(void)scores_len;(void)left_shift;(void)right_shift;
}
void block_set_all_rev_aaprofile(struct AAProfile *p, const uint8_t *order, uintptr_t order_len,
                                 const int8_t *scores, uintptr_t scores_len,
                                 uintptr_t left_shift, uintptr_t right_shift) {
    (void)p;(void)order;(void)order_len;(void)scores;(void)scores_len;(void)left_shift;(void)right_shift;
}
int8_t* aaprofile_pos_aa(struct AAProfile *profile) { return profile ? profile->pos_aa : NULL; }
int16_t* aaprofile_aa_pos(struct AAProfile *profile) { return profile ? profile->aa_pos : NULL; }
void block_set_gap_open_C_aaprofile(struct AAProfile *p, uintptr_t i, int8_t gap) { (void)p;(void)i;(void)gap; }
void block_set_gap_close_C_aaprofile(struct AAProfile *p, uintptr_t i, int8_t gap) { (void)p;(void)i;(void)gap; }
void block_set_gap_open_R_aaprofile(struct AAProfile *p, uintptr_t i, int8_t gap) { (void)p;(void)i;(void)gap; }
void block_set_all_gap_open_C_aaprofile(struct AAProfile *p, int8_t gap) { (void)p;(void)gap; }
void block_set_all_gap_close_C_aaprofile(struct AAProfile *p, int8_t gap) { (void)p;(void)gap; }
void block_set_all_gap_open_R_aaprofile(struct AAProfile *p, int8_t gap) { (void)p;(void)gap; }
int8_t block_get_aaprofile(const struct AAProfile *p, uintptr_t i, uint8_t b) { (void)p;(void)i;(void)b; return 0; }
int8_t block_get_gap_extend_aaprofile(const struct AAProfile *p) { return p ? p->gap_extend : 0; }
size_t block_get_curr_len_aaprofile(const struct AAProfile *p) { return p ? p->curr_len : 0; }
void block_free_aaprofile(struct AAProfile *p) {
    if (p) { free(p->pos_aa); free(p->aa_pos); free(p); }
}

/* --- Cigar stubs --- */
struct Cigar *block_new_cigar(uintptr_t query_len, uintptr_t reference_len) {
    (void)query_len; (void)reference_len;
    return (struct Cigar*)calloc(1, sizeof(struct Cigar));
}
struct OpLen block_get_cigar(const struct Cigar *cigar, uintptr_t i) {
    (void)cigar; (void)i;
    struct OpLen op = { Sentinel, 0 };
    return op;
}
uintptr_t block_len_cigar(const struct Cigar *cigar) { return cigar ? cigar->len : 0; }
void block_free_cigar(struct Cigar *cigar) { free(cigar); }

/* --- PaddedBytes stubs --- */
struct PaddedBytes *block_new_padded_aa(uintptr_t len, uintptr_t max_size) {
    (void)len; (void)max_size;
    return (struct PaddedBytes*)calloc(1, sizeof(struct PaddedBytes));
}
void block_set_bytes_padded_aa(struct PaddedBytes *p, const uint8_t *s, uintptr_t len, uintptr_t max_size) {
    (void)p;(void)s;(void)len;(void)max_size;
}
void block_set_bytes_padded_aa_numsequence(struct PaddedBytes *p, const uint8_t *s, uintptr_t len, uintptr_t max_size) {
    (void)p;(void)s;(void)len;(void)max_size;
}
void block_free_padded_aa(struct PaddedBytes *p) { free(p); }

/* --- PosBias stubs --- */
struct PosBias *block_new_pos_bias(uintptr_t len, uintptr_t max_size) {
    (void)len; (void)max_size;
    return (struct PosBias*)calloc(1, sizeof(struct PosBias));
}
void block_set_pos_bias(struct PosBias *bias, const int16_t *b, uintptr_t len) { (void)bias;(void)b;(void)len; }
void block_free_pos_bias(struct PosBias *bias) { free(bias); }

/* --- Block alignment stubs --- */
/* All alignment functions are no-ops; block_res always returns a score that won't match,
   causing SSW to fall back to the standard SW traceback. */
BlockHandle block_new_aa_trace_xdrop(uintptr_t query_len, uintptr_t reference_len, uintptr_t max_size) {
    (void)query_len; (void)reference_len; (void)max_size;
    return (BlockHandle)&_block_dummy;
}
void block_align_aa_trace_xdrop(BlockHandle b, const struct PaddedBytes *q, const struct PaddedBytes *r,
                                const struct AAMatrix *m, struct Gaps g, struct SizeRange s, int32_t x) {
    (void)b;(void)q;(void)r;(void)m;(void)g;(void)s;(void)x;
}
void block_align_profile_aa_trace_xdrop(BlockHandle b, const struct PaddedBytes *q, const struct AAProfile *r,
                                        struct SizeRange s, int32_t x) {
    (void)b;(void)q;(void)r;(void)s;(void)x;
}
void block_align_aa_trace_xdrop_posbias(BlockHandle b, const struct PaddedBytes *q, const struct PosBias *q_bias,
                                        const struct PaddedBytes *r, const struct PosBias *r_bias,
                                        const struct AAMatrix *m, struct Gaps g, struct SizeRange s, int32_t x) {
    (void)b;(void)q;(void)q_bias;(void)r;(void)r_bias;(void)m;(void)g;(void)s;(void)x;
}
struct AlignResult block_res_aa_trace_xdrop(BlockHandle b) {
    (void)b;
    /* Return a very negative score so SSW never accepts block aligner result */
    struct AlignResult res = { -1000000000, 0, 0 };
    return res;
}
void block_cigar_aa_trace_xdrop(BlockHandle b, uintptr_t query_idx, uintptr_t reference_idx, struct Cigar *cigar) {
    (void)b;(void)query_idx;(void)reference_idx;(void)cigar;
}
void block_free_aa_trace_xdrop(BlockHandle b) { (void)b; }
