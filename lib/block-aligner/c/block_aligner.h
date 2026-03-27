#ifndef block_aligner_h
#define block_aligner_h

/* Stub header for block_aligner - provides types and no-op function declarations.
   The actual Rust block-aligner library is not linked; StripedSmithWaterman will
   fall back to the standard SW traceback path when block alignment returns a
   mismatched score. */

#include <stdarg.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#define ALIGNED(n) __attribute__ ((aligned(n)))

enum Operation
#ifdef __cplusplus
  : uint8_t
#endif
{
  Sentinel = 0,
  M = 1,
  Eq = 2,
  X = 3,
  I = 4,
  D = 5,
};
#ifndef __cplusplus
typedef uint8_t Operation;
#endif

typedef struct AAMatrix AAMatrix;
typedef struct AAProfile AAProfile;
typedef struct Cigar Cigar;
typedef struct NucMatrix NucMatrix;
typedef struct PaddedBytes PaddedBytes;
typedef struct PosBias PosBias;

typedef struct OpLen {
  Operation op;
  uintptr_t len;
} OpLen;

typedef void *BlockHandle;

typedef struct Gaps {
  int8_t open;
  int8_t extend;
} Gaps;

typedef struct SizeRange {
  uintptr_t min;
  uintptr_t max;
} SizeRange;

typedef struct AlignResult {
  int32_t score;
  uintptr_t query_idx;
  uintptr_t reference_idx;
} AlignResult;

typedef struct ByteMatrix {
  int8_t match_score;
  int8_t mismatch_score;
} ByteMatrix;

#ifdef __cplusplus
extern "C" {
#endif

/* --- AAMatrix --- */
struct AAMatrix *block_new_simple_aamatrix(int8_t match_score, int8_t mismatch_score);
void block_set_aamatrix(struct AAMatrix *matrix, uint8_t a, uint8_t b, int8_t score);
void block_set_aamatrix_num(struct AAMatrix *matrix, int8_t a, int8_t b, int8_t score);
void block_free_aamatrix(struct AAMatrix *matrix);

/* --- AAProfile --- */
struct AAProfile *block_new_aaprofile(uintptr_t str_len, uintptr_t block_size, int8_t gap_extend);
uintptr_t block_len_aaprofile(const struct AAProfile *profile);
void block_clear_aaprofile(struct AAProfile *profile, uintptr_t str_len);
void block_set_aaprofile(struct AAProfile *profile, uintptr_t i, uint8_t b, int8_t score);
void block_set_all_aaprofile(struct AAProfile *profile, const uint8_t *order,
                             uintptr_t order_len, const int8_t *scores,
                             uintptr_t scores_len, uintptr_t left_shift, uintptr_t right_shift);
void block_set_all_rev_aaprofile(struct AAProfile *profile, const uint8_t *order,
                                 uintptr_t order_len, const int8_t *scores,
                                 uintptr_t scores_len, uintptr_t left_shift, uintptr_t right_shift);
int8_t* aaprofile_pos_aa(struct AAProfile *profile);
int16_t* aaprofile_aa_pos(struct AAProfile *profile);
void block_set_gap_open_C_aaprofile(struct AAProfile *profile, uintptr_t i, int8_t gap);
void block_set_gap_close_C_aaprofile(struct AAProfile *profile, uintptr_t i, int8_t gap);
void block_set_gap_open_R_aaprofile(struct AAProfile *profile, uintptr_t i, int8_t gap);
void block_set_all_gap_open_C_aaprofile(struct AAProfile *profile, int8_t gap);
void block_set_all_gap_close_C_aaprofile(struct AAProfile *profile, int8_t gap);
void block_set_all_gap_open_R_aaprofile(struct AAProfile *profile, int8_t gap);
int8_t block_get_aaprofile(const struct AAProfile *profile, uintptr_t i, uint8_t b);
int8_t block_get_gap_extend_aaprofile(const struct AAProfile *profile);
size_t block_get_curr_len_aaprofile(const struct AAProfile *profile);
void block_free_aaprofile(struct AAProfile *profile);

/* --- Cigar --- */
struct Cigar *block_new_cigar(uintptr_t query_len, uintptr_t reference_len);
struct OpLen block_get_cigar(const struct Cigar *cigar, uintptr_t i);
uintptr_t block_len_cigar(const struct Cigar *cigar);
void block_free_cigar(struct Cigar *cigar);

/* --- PaddedBytes --- */
struct PaddedBytes *block_new_padded_aa(uintptr_t len, uintptr_t max_size);
void block_set_bytes_padded_aa(struct PaddedBytes *padded, const uint8_t *s, uintptr_t len, uintptr_t max_size);
void block_set_bytes_padded_aa_numsequence(struct PaddedBytes *padded, const uint8_t *s, uintptr_t len, uintptr_t max_size);
void block_free_padded_aa(struct PaddedBytes *padded);

/* --- PosBias --- */
struct PosBias *block_new_pos_bias(uintptr_t len, uintptr_t max_size);
void block_set_pos_bias(struct PosBias *bias, const int16_t *b, uintptr_t len);
void block_free_pos_bias(struct PosBias *bias);

/* --- Block alignment (trace + xdrop) --- */
BlockHandle block_new_aa_trace_xdrop(uintptr_t query_len, uintptr_t reference_len, uintptr_t max_size);
void block_align_aa_trace_xdrop(BlockHandle b, const struct PaddedBytes *q, const struct PaddedBytes *r,
                                const struct AAMatrix *m, struct Gaps g, struct SizeRange s, int32_t x);
void block_align_profile_aa_trace_xdrop(BlockHandle b, const struct PaddedBytes *q, const struct AAProfile *r,
                                        struct SizeRange s, int32_t x);
void block_align_aa_trace_xdrop_posbias(BlockHandle b, const struct PaddedBytes *q, const struct PosBias *q_bias,
                                        const struct PaddedBytes *r, const struct PosBias *r_bias,
                                        const struct AAMatrix *m, struct Gaps g, struct SizeRange s, int32_t x);
struct AlignResult block_res_aa_trace_xdrop(BlockHandle b);
void block_cigar_aa_trace_xdrop(BlockHandle b, uintptr_t query_idx, uintptr_t reference_idx, struct Cigar *cigar);
void block_free_aa_trace_xdrop(BlockHandle b);

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* block_aligner_h */
