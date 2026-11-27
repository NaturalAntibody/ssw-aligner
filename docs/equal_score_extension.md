# Equal-score extension of alignment ends

An optional post-processing step: after SSW returns the optimal local alignment, extend its start and/or end over
neighbouring residues that add exactly 0 to the score. The score stays the same, the alignment stays optimal, it just
covers more of the target (germline).

Motivation: RIOT (antibody numbering) aligns a query to V/J germline genes and then *forces* the uncovered germline
ends as ungapped matches (`force_n_terminus_matches`, `force_c_terminus_matches`). When the optimal alignment has an
equal-score variant that reaches further, SSW reports the shortest one, and the forcing then numbers the wrong
residues.

## Example: catumaxomab VH

BLOSUM62, gap open 11, gap extend 1. Two optimal alignments with the same score (485):

```
A (SSW reports this one)            B (same score, covers the germline start)
query     L L E Q S …               query     E V Q L L E Q S …
germline  L V E S G …               germline  Q V Q - L V E S G …
          germline position 4                 germline position 1, 3M1I…

B - A = E/Q +2, V/V +4, Q/Q +5, insertion -11 = 0
```

RIOT then forces `VQL` onto germline `QVQ` (ungapped) and leaves `E` out of the domain. With B the numbering would
be E1 V2 Q3 L3.1 L4.

## Where SSW breaks ties (`ssw_aligner/_lib/ssw.c`)

The query is the `read`, the germline is the `ref`. SSW consistently picks the **shortest** equal-score alignment:

| End | Code | Picks |
|---|---|---|
| End (C-terminal) | `sw_sse2_byte`/`_word`: `if (temp > max)`, then the smallest `end_read` equal to max | the first ref column reaching the max, the smallest read position in it |
| Start (N-terminal) | `ssw_align`: reverse pass terminating at `maxColumn[i] == terminate` | the first column (scanning back from the end) reaching `score1`, i.e. the latest start |
| Gap placement inside | `banded_sw`: `if (temp1 <= temp2)` prefers the diagonal | out of scope, doesn't change coverage |

Why not change these lines:

1. The SIMD forward pass enters the update block only when a lane max grows (`cmp != 0xffff`), so an equal score in a
   later column isn't even detected.
2. The reverse pass is a free local alignment in `[0..ref_end] × [0..read_end]`, not anchored at the end. Scanning
   further for the farthest start could pick a different equal-score alignment (e.g. another domain copy in a
   multidomain query).
3. `banded_sw` floors scores at 0, so with a zero-score extension the traceback could stop before the corner and
   the CIGAR would disagree with the begin/end positions.

## Design

A new C function, run after `ssw_align`, leaving the SSW kernels untouched:

```c
void ssw_extend_ties(const s_profile* prof, const int8_t* ref, int32_t refLen,
                     uint8_t weight_gapO, uint8_t weight_gapE, s_align* r,
                     uint8_t sides /* 1: start, 2: end, 3: both */, int32_t max_insertions);
```

For each enabled side:

- **Anchored DP** (scalar, int32, affine gaps / Gotoh, **no zero floor**) from the alignment edge outward: from
  `(read_begin1, ref_begin1)` backwards for the start, from `(read_end1, ref_end1)` forwards for the end. Every path
  starts at the edge, so the extension stays attached to the reported alignment.
- **Window:** `g` = ref residues outside the alignment on that side (germline positions before the start / after the
  end), and at most `min(read residues on that side, g + max_insertions)` read residues. The extension reaches at most
  the germline start/end, as RIOT's forcing does (`min(query, germline)`), plus up to `max_insertions` inserted query
  residues. `max_insertions = 0` allows deletions only. Matrix size is at most `(g + max_insertions + 1) × (g + 1)`.
- **Band:** the path stays in `i - j <= max_insertions` (`i` read, `j` ref residues from the edge): at no point more
  than `max_insertions` read residues beyond the ref residues covered so far. Only cells in the band are computed,
  about half of the matrix.
- **Candidates:** cells (in the band) with `H == 0` whose last step is a match. The best extension score is 0: a positive one would mean SSW's alignment wasn't optimal. Positive
  cells are ignored.
- **Choice:** most ref (germline) residues covered, then fewest read residues (fewer insertions, less linker or other
  non-germline sequence pulled in).
- **CIGAR:** traceback of the extension path, prepended (start) or appended (end), adjacent operations merged. Neither
  the alignment nor the extension can end with a gap (removing it would raise the score), so they always join M to M.
- **Update** `read_begin1`/`ref_begin1` or `read_end1`/`ref_end1` and the CIGAR. `score1` is unchanged.

Not covered: equal-score alignments that are not extensions of the reported one (a different path, not a longer
version of the same).

### Python API

`StripedSmithWaterman(..., extend_equal_score=None | "start" | "end" | "both", max_tie_insertions=2)`.
Default `None`: outputs identical to the previous version (`tests/test_regression.py` unchanged).

### Noise

The step only adds zero-score extensions; anything scoring > 0 SSW already includes. Gapped extensions are the ones
to watch: matched residues have to pay the whole gap, e.g. `EVQ` (+14) followed by 4 inserted residues (-14) ties.
`max_tie_insertions` bounds the number of such residues (2 by default).

## Performance

Requirement: negligible overhead in RIOT.

- No DP when there is nothing to extend (`g == 0` or no read residues on that side), the common case for complete
  sequences.
- Otherwise the band of a `(g + max_insertions + 1) × (g + 1)` matrix, scalar; the ref residues' substitution matrix
  rows are looked up once per call. Stack buffer up to ~32 × 32 cells, `malloc` above.

Measured (x86-64, Python 3.11, minimum of 7–9 runs):

| Benchmark | Calls | DP runs | Changed | Overhead |
|---|---|---|---|---|
| therapeutic queries vs their top 10 V germlines, `"start"` | 3,000 | 377 | 122 | none measurable (-0.2%) |
| therapeutic queries vs their top 5 J germlines, `"end"` | 1,500 | 41 | 23 | none measurable (-0.2%) |

Worst case, every call runs the DP: V germline truncated by `g` residues, with `g + 10` random residues before it
(~45 µs per call without the extension):

| `g` | added per call |
|---|---|
| 5 | +0.7 µs (+1.5%) |
| 15 | +0.4 µs (+0.8%) |
| 30 | +1.2 µs (+2.7%) |
| 60 | +4.2 µs (+10.5%) |

(Before the band and the ref row lookup: +8.8 µs at `g` 60.) The numbers around 1 µs are at the level of noise.

RIOT, multidomain AA numbering of 564 therapeutics (every 4th), V `"start"`, J `"end"`: 3.35 s off vs 3.32 s on, no
measurable difference.

## Results in RIOT

All 2,256 therapeutics (2,438 domains), extension off vs on (V `"start"`, J `"end"`):

- 63 V or J alignments change: 47 V starts (germline start 2→1, 3→1, 4→3, 4→2, 4→1, 10→9), 16 J ends (one more
  germline residue). Scores unchanged. (65 domains in total: the 2 VL domains preceding catumaxomab's and solitomab's
  VH end one residue earlier.)
- Numbering (regions, scheme mapping, gene calls, domain boundaries) changes only for catumaxomab and solitomab, the
  same VH: `3M1I95M` from germline position 1, numbered E1 V2 Q3 L3.1 L4 (before: V1 Q2 L3 L4 E5…, with `E` left in
  the preceding VL domain). The other 63 are ungapped extensions RIOT already forced (`force_n/c_terminus_matches`).
- Only those two used an insertion (1); no extension used 2. No extension took linker residues.
- RIOT tests pass.

All 63 changed alignments, tie off vs on: [equal_score_extension_therapeutics.md](equal_score_extension_therapeutics.md).

### NGS sample

Single domain AA numbering (`create_riot_aa`) of 10,000 sequences sampled (seed 0) from RIOT's
`data/ngs/ngs_sample_clean_aa.fasta`; 9,637 got a V alignment.

- 80 sequences have a different alignment: 58 V starts, 22 J ends. Most V changes are reads starting inside FR1 with
  one or two residues before the alignment start (germline start 16→15 12×, 14→13 9×, 17→16 6×; 2→1 10×).
- No extension used an insertion: all are ungapped residue pairs scoring 0 in total.
- Numbering is unchanged in 79 (residues RIOT already forced). One sequence changes, see below.
- Runtime, 2,500 sequences: 7.27 s off vs 7.29 s on (+0.3%, noise).

### Side effect: alignment length filters

The extension happens before RIOT filters the alignments, so a filter on the alignment's germline length
(`t_end - t_start`) sees the extended alignment. One sequence of the NGS sample, a broken read (stop codon at
position 23, no J motif in frame), gets a J only with the extension:

```
SRR8365476.30165943  …EDFAVYYYQQYGTSVRTLLARGPGWR, IGKJ3*01 = FTFGPGTKVDIK

Tie off: 3M, germline 4–6, score 19
query     QQYGTSVRTLLARGPGWR
germline  -------------GPG--
match                  |||
region    ·············JJJ··
                         ^ SW J end (91, germline position 6 of 12, 3M)

Tie on: 5M, germline 4–8, score 19
query     QQYGTSVRTLLARGPGWR
germline  -------------GPGTK
match                  |||..
region    ·············JJJJJ
                           ^ SW J end (93, germline position 8 of 12, 5M)

G/G +6, P/P +7, G/G +6 = 19, then W/T -2, R/K +2 = 19 (the end of the read)
```

| | Tie off | Tie on |
|---|---|---|
| J germline length (RIOT minimum: 5) | 3, filtered out | 5, kept |
| E-value (RIOT maximum: 0.3) | 0.252 | 0.252 |
| `j_call` | none | `IGKJ3*01` |
| CDR3 / FR4 | `QQYGTS` / none | `QQYGTSVRTLLA` / `RGPGWR` |

Had the 3-residue alignment passed, RIOT's forced C-terminal extension would have numbered the same residues
(`WR` onto `TK`): the extension doesn't change the numbering of an accepted alignment, only whether this one is accepted.

Decision: not a problem, kept as is. It is a borderline J of a broken read (the e-value threshold of J, 0.3, is very
lenient), 1 sequence in 10,000. Filtering on the alignment before the extension would need the pre-extension positions
from ssw-aligner, not worth it for now.

## Plan

Status: 1–4 implemented. ssw-aligner: branch `feat/extend_equal_score`. RIOT: branch `multidomain_alignement`
(uncommitted); it needs the released version for amino acid numbering with `extend_alignment` (the old one rejects
the new keyword argument).

1. `ssw.c`/`ssw.h`: `ssw_extend_ties`.
2. Wrapper (`_ssw_wrapper.pyx`, `.pyi`): `extend_equal_score`, `max_tie_insertions`. Tests: catumaxomab-like start
   case (`3M1I…`), the same at the end, a -1 extension not taken, no tie → identical result, `max_tie_insertions=0`
   blocks the catumaxomab case, 3 insertions blocked at the default 2, stops at sequence edges, both `score_size`
   paths (8 and 16 bit). Release as `feat:` (minor bump).
3. RIOT: amino acid numbering (`RiotNumberingAA`, `MultidomainRiotNumberingAA`) with `extend_alignment` only:
   `"start"` for V, `"end"` for J, C not extended. Without `extend_alignment` alignments are plain Smith-Waterman, as
   before. Not used for nucleotide numbering for now (small integer DNA scores make zero-sum extensions much more
   frequent). `force_n/c_terminus_matches` stay, they just force less.
4. Evaluation on RIOT's therapeutics: every alignment whose begin/end/CIGAR changes (per side), added insertions taking
   `G`/`S` linker residues or crossing domain boundaries, catumaxomab → E1 V2 Q3 L3.1 L4, how many extensions use 1 or
   2 insertions, runtime (target < 5% total).
