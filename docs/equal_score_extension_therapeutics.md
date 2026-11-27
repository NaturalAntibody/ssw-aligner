# Equal-score extension: RIOT therapeutics, Tie break off vs on

Every V start and J end in `riot_na/data/therapeutics.csv` whose alignment changes when the equal-score extension
([equal_score_extension.md](equal_score_extension.md)) is switched on: Tie break off (`extend_equal_score="none"`) vs Tie break on
(V `"start"`, J `"end"`, `max_tie_insertions=2`).

State: ssw-aligner branch `feat/extend_equal_score` (uncommitted), RIOT branch `multidomain_alignement` (commit
`b95a439` plus the params change), multidomain AA numbering (`create_multidomain_riot_aa`, IMGT), data as of
2026-11-26: 2,256 chains, 2,438 domains.

## Overview

- 63 alignments change: 47 V starts, 16 J ends. Scores are unchanged in all of them.
- Numbering changes in 2 (catumaxomab, solitomab; the same VH). In the other 61 the extension covers residues that
  RIOT already forced as matches (`force_n_terminus_matches`, `force_c_terminus_matches`), so the regions and the
  scheme mapping stay the same; only the AIRR alignment fields (`v_sequence_start_aa`, `v_germline_start_aa`,
  `v_cigar_aa`, `j_sequence_end_aa`, `j_germline_end_aa`, `j_cigar_aa`, identities) change.
- Every extension adds mismatches, scoring 0 together (BLOSUM62):
  - 56 add a single residue pair that scores 0: V `D/Q` 24, `E/S` 4, `Q/D` 4, `D/S` 3, `S/Q` 3, `V/A` 1, `H/Q` 1;
    J `K/N` 8, `K/S` 4, `N/K` 2, `S/K` 1, `T/A` 1 (query/germline).
  - 5 add two pairs scoring +2 and -2: danvilostomig, umizortamig (`LV`/`IQ`), regdanvimab (`EL`/`QS`), tregalizumab
    (`EE`/`QV`), vatelizumab (`DF`/`QI`).
  - 2 add `EVQ` + an insertion (`L`): catumaxomab, solitomab (+11 -11).

  So, with the extension, an alignment can start or end with a mismatch, which a Smith-Waterman alignment never does
  (an edge pair scoring <= 0 is dropped, as SSW does on ties).
- No extension used 2 insertions, none took linker residues (`G`/`S` repeats).

## How to read the views

Each view shows the alignment edge with the Tie break off, then on, over the same query positions:

```
query     the query residues (amino acids)
germline  the aligned germline residues, '-' where there is no alignment (or a gap)
match     '|' identical, '.' different, ' ' not aligned or gap
region    V / J : residue in the Smith-Waterman V / J alignment
          e     : residue covered only by RIOT's forced extension (force_n_terminus_matches /
                  force_c_terminus_matches), matched to the germline residue shown
          ·     : not aligned
          ‖     : domain boundary
          ^     : Smith-Waterman V alignment start / J alignment end (1-based query position, germline position)
```

Positions are 1-based positions in the chain (`molecule_sequence`). Each example is named by the therapeutic, the chain
(`molecule_name` in therapeutics.csv) and the domain number in the chain. Gene names like `IGKV-ND3J` are germlines of
the mouse database (IDs instead of IMGT names), `vhh3-1` and `ighJ-5` alpaca.

## Summary

| Therapeutic | Chain | Domain | Gene | Germline | Change | Numbering |
|---|---|---|---|---|---|---|
| catumaxomab | chain | 2 | V | IGHV-5O75 | germline start 4 → 1 | segment_start, fwr1_aa, scheme_residue_mapping, preceding domain end 129 → 128 |
| solitomab | chain | 2 | V | IGHV-5O75 | germline start 4 → 1 | segment_start, fwr1_aa, scheme_residue_mapping, preceding domain end 129 → 128 |
| amatuximab | anti-MSLN light chain 1 | 1 | V | IGKV-ND3J | germline start 2 → 1 | unchanged |
| blinatumomab | chain | 3 | V | IGHV-22XF | germline start 2 → 1 | unchanged |
| blinatumomab | chain | 4 | V | IGKV-O6WO | germline start 2 → 1 | unchanged |
| brivekimig | Heavy chain | 1 | V | vhh3-1 | germline start 2 → 1 | unchanged |
| catumaxomab | chain | 3 | V | IGHV-22XF | germline start 2 → 1 | unchanged |
| clivatuzumab tetraxetan | anti-MUC1 light chain 1 | 1 | V | IGKV-3GON | germline start 2 → 1 | unchanged |
| cosfroviximab | anti-Zaire Ebolavirus GP heavy chain 1 | 1 | V | IGHV-KRVC | germline start 2 → 1 | unchanged |
| danvilostomig | Heavy chain | 3 | V | IGKV1-27*01 | germline start 4 → 2 | unchanged |
| elipovimab | anti-HIV-1 gp120 light chain 1 | 1 | V | IGLV3-21*02 | germline start 10 → 9 | unchanged |
| emfizatamab | anti-CD3E heavy chain 1 | 6 | V | IGHV3-66*01 | germline start 4 → 3 | unchanged |
| gontivimab | scVH-VH-VH chain | 1 | V | vhh3-1 | germline start 2 → 1 | unchanged |
| gumokimab | heavy chain 1 | 1 | V | IGHV4-30-4*01 | germline start 2 → 1 | unchanged |
| isecarosmab | scVH-VH' chain | 1 | V | vhh3-1 | germline start 2 → 1 | unchanged |
| ivicentamab | light chain 2 | 1 | V | IGKV1-5*03 | germline start 4 → 3 | unchanged |
| izuralimab | Chain IG scFv-h-CH2-CH3 (anti-PDCD1) | 1 | V | IGKV-EAXP | germline start 2 → 1 | unchanged |
| lacnotuzumab | anti-CSF1 heavy chain 1 | 1 | V | IGHV-O4YA | germline start 2 → 1 | unchanged |
| lancastotug | heavy chain 1 | 1 | V | IGHV4-30-4*01 | germline start 2 → 1 | unchanged |
| lucorafusp alfa | heavy chain 1 | 1 | V | IGHV4-30-4*01 | germline start 2 → 1 | unchanged |
| lunsekimig | Heavy chain | 1 | V | vhh3-1 | germline start 2 → 1 | unchanged |
| mipletamig | chain | 4 | V | IGKV-HQL3 | germline start 2 → 1 | unchanged |
| misitatug | light chain 1 | 1 | V | IGKV-GLVF | germline start 2 → 1 | unchanged |
| naxitamab | anti-Ganglioside GD2 light chain 1 | 1 | V | IGKV-2UKG | germline start 2 → 1 | unchanged |
| nivatrotamab | light chain | 1 | V | IGKV-2UKG | germline start 2 → 1 | unchanged |
| regdanvimab | anti-SARS-CoV-2 Spike RBD light chain 1 | 1 | V | IGLV1-51*01 | germline start 3 → 1 | unchanged |
| resugosbart | light chain 1 | 1 | V | IGKV-2UKG | germline start 2 → 1 | unchanged |
| rimteravimab | heavy chain 1 | 1 | V | vhh3-1 | germline start 2 → 1 | unchanged |
| solitomab | chain | 3 | V | IGHV1-46*04 | germline start 2 → 1 | unchanged |
| solitomab | chain | 4 | V | IGKV-O6WO | germline start 2 → 1 | unchanged |
| sonelokimab | chain | 1 | V | vhh3-1 | germline start 2 → 1 | unchanged |
| tilogotamab | anti-TNFRSF10B light chain 1 | 1 | V | IGKV-K25A | germline start 2 → 1 | unchanged |
| timolumab | anti-AOC3 light chain 1 | 1 | V | IGKV1D-13*02 | germline start 2 → 1 | unchanged |
| tislelizumab | anti-PDCD1 light chain 1 | 1 | V | IGKV-2UKG | germline start 2 → 1 | unchanged |
| tregalizumab | anti-CD4 heavy chain 1 | 1 | V | IGHV-MG5V | germline start 3 → 1 | unchanged |
| trontinemab | Heavy chain anti-APP A beta 1-40/42 and anti-TFRC | 2 | V | IGHV4-4*08 | germline start 4 → 3 | unchanged |
| umizortamig | heavy chain 1 | 6 | V | IGHV3-66*01 | germline start 4 → 3 | unchanged |
| umizortamig | heavy chain 1 | 7 | V | IGKV1-5*03 | germline start 4 → 2 | unchanged |
| vatelizumab | anti-ITGA2 light chain 1 | 1 | V | IGKV-O6WO | germline start 3 → 1 | unchanged |
| veltuzumab | anti-MS4A1 light chain 1 | 1 | V | IGKV-27KA | germline start 2 → 1 | unchanged |
| vudalimab | Chain scFv-h-CH2-CH3 | 1 | V | IGKV-EAXP | germline start 2 → 1 | unchanged |
| yttrium y-90 clivatuzumab tetraxetan | anti-MUC1 light chain 1 | 1 | V | IGKV-3GON | germline start 2 → 1 | unchanged |
| zeleciment | heavy chain 1 | 1 | V | IGHV-FMN7 | germline start 2 → 1 | unchanged |
| zeleciment basivarsen | heavy chain 1 | 1 | V | IGHV-FMN7 | germline start 2 → 1 | unchanged |
| zeleciment rostudirsen | heavy chain 1 | 1 | V | IGHV-FMN7 | germline start 2 → 1 | unchanged |
| zifibancimig | heavy chain 1 | 1 | V | IGHV3-21*07 | germline start 4 → 3 | unchanged |
| zuberitamab | anti-MS4A1 light chain 1 | 1 | V | IGKV-27KA | germline start 2 → 1 | unchanged |
| docaravimab | anti-Rabies Virus Strain ERA GP Ectodomain Epitope G-III light chain 1 | 1 | J | IGKJ-WRAP | germline end 11 → 12 (of 12) | unchanged |
| enristomig | Heavy chain | 1 | J | ighJ-5 | germline end 14 → 15 (of 15) | unchanged |
| enristomig | Heavy chain | 2 | J | ighJ-5 | germline end 14 → 15 (of 15) | unchanged |
| futuximab | anti-EGFR light chain 1 | 1 | J | IGKJ-WRAP | germline end 11 → 12 (of 12) | unchanged |
| lecanemab | anti-APP light chain 1 | 1 | J | IGKJ-WRAP | germline end 11 → 12 (of 12) | unchanged |
| lemalesomab | light chain 1 | 1 | J | IGKJ-WRAP | germline end 11 → 12 (of 12) | unchanged |
| muromonab-cd3 | anti-CD3E light chain 1 | 1 | J | IGKJ-Z5J4 | germline end 11 → 12 (of 12) | unchanged |
| naptumomab estafenatox | Light chain | 1 | J | IGKJ-WRAP | germline end 11 → 12 (of 12) | unchanged |
| ordesekimab | anti-IL15 light chain 1 | 1 | J | IGKJ2*01 | germline end 11 → 12 (of 12) | unchanged |
| ozekibart | heavy chain 1 | 1 | J | ighJ-4 | germline end 13 → 14 (of 15) | unchanged |
| ozekibart | heavy chain 1 | 2 | J | ighJ-4 | germline end 13 → 14 (of 15) | unchanged |
| ruplizumab | anti-CD40LG light chain 1 | 1 | J | IGKJ-WRAP | germline end 11 → 12 (of 12) | unchanged |
| tegoprubart | light chain 1 | 1 | J | IGKJ-WRAP | germline end 11 → 12 (of 12) | unchanged |
| tomuzotuximab | anti-EGFR heavy chain 1 | 1 | J | IGHJ-KZ3G | germline end 14 → 15 (of 15) | unchanged |
| veligrotug | light chain 1 | 1 | J | IGKJ-WRAP | germline end 11 → 12 (of 12) | unchanged |
| zolimomab | anti-CD3E light chain 1 | 1 | J | IGKJ-Z5J4 | germline end 11 → 12 (of 12) | unchanged |

## Numbering changes (2)

The VH after the scFv linker `(G4S)3`: the extension moves `E` from the linker into the VH and numbers it.
### catumaxomab (chain), domain 2: V start

- V gene: IGHV-5O75; domain 129–253 (off: 130–253)
- germline start 4 → 1, `95M` → `3M1I95M`, extension 3 → 0; score 176.5 bits (unchanged)
- Numbering: `segment_start` `130` → `129`; `fwr1_aa` `VQLLEQSGAELVRPGTSVKISCKAS` → `EVQLLEQSGAELVRPGTSVKISCKAS`; preceding domain end 129 → 128; IMGT 1–5 1:V 2:Q 3:L 4:L 5:E → 1:E 2:V 3:Q 3.1:L 4:L

Tie break off:

```
query     GGGGSGGGGSE‖VQLLEQSGAELVRPGT
germline  -----------‖QVQLQQSGAELVRPGT
match                ‖...|.|||||||||||
region    ···········‖eeeVVVVVVVVVVVVV
                         ^ SW V start (133, germline position 4, 95M)
```

Tie break on:

```
query     GGGGSGGGGS‖EVQLLEQSGAELVRPGT
germline  ----------‖QVQ-LQQSGAELVRPGT
match               ‖.|| |.|||||||||||
region    ··········‖VVVVVVVVVVVVVVVVV
                     ^ SW V start (129, germline position 1, 3M1I95M)
```

### solitomab (chain), domain 2: V start

- V gene: IGHV-5O75; domain 129–253 (off: 130–253)
- germline start 4 → 1, `95M` → `3M1I95M`, extension 3 → 0; score 176.5 bits (unchanged)
- Numbering: `segment_start` `130` → `129`; `fwr1_aa` `VQLLEQSGAELVRPGTSVKISCKAS` → `EVQLLEQSGAELVRPGTSVKISCKAS`; preceding domain end 129 → 128; IMGT 1–5 1:V 2:Q 3:L 4:L 5:E → 1:E 2:V 3:Q 3.1:L 4:L

Tie break off:

```
query     GGGGSGGGGSE‖VQLLEQSGAELVRPGT
germline  -----------‖QVQLQQSGAELVRPGT
match                ‖...|.|||||||||||
region    ···········‖eeeVVVVVVVVVVVVV
                         ^ SW V start (133, germline position 4, 95M)
```

Tie break on:

```
query     GGGGSGGGGS‖EVQLLEQSGAELVRPGT
germline  ----------‖QVQ-LQQSGAELVRPGT
match               ‖.|| |.|||||||||||
region    ··········‖VVVVVVVVVVVVVVVVV
                     ^ SW V start (129, germline position 1, 3M1I95M)
```

## V starts, numbering unchanged (45)

### amatuximab (anti-MSLN light chain 1), domain 1: V start

- V gene: IGKV-ND3J; domain 1–213
- germline start 2 → 1, `94M` → `95M`, extension 1 → 0; score 166.6 bits (unchanged)
- Numbering: unchanged

Tie break off:

```
query     DIELTQSPAIMSAS
germline  QIVLTQSPAIMSAS
match     .|.|||||||||||
region    eVVVVVVVVVVVVV
           ^ SW V start (2, germline position 2, 94M)
```

Tie break on:

```
query     DIELTQSPAIMSAS
germline  QIVLTQSPAIMSAS
match     .|.|||||||||||
region    VVVVVVVVVVVVVV
          ^ SW V start (1, germline position 1, 95M)
```

### blinatumomab (chain), domain 3: V start

- V gene: IGHV-22XF; domain 256–392
- germline start 2 → 1, `97M` → `98M`, extension 1 → 0; score 177.2 bits (unchanged)
- Numbering: unchanged

Tie break off:

```
query     VTVSSGGGGS‖DIKLQQSGAELARP
germline  ----------‖QVQLQQSGAELARP
match               ‖...|||||||||||
region    ··········‖eVVVVVVVVVVVVV
                      ^ SW V start (257, germline position 2, 97M)
```

Tie break on:

```
query     VTVSSGGGGS‖DIKLQQSGAELARP
germline  ----------‖QVQLQQSGAELARP
match               ‖...|||||||||||
region    ··········‖VVVVVVVVVVVVVV
                     ^ SW V start (256, germline position 1, 98M)
```

### blinatumomab (chain), domain 4: V start

- V gene: IGKV-O6WO; domain 393–504
- germline start 2 → 1, `94M` → `95M`, extension 1 → 0; score 168.7 bits (unchanged)
- Numbering: unchanged

Tie break off:

```
query     GGSGGSGGVD‖DIQLTQSPAIMSAS
germline  ----------‖QIVLTQSPAIMSAS
match               ‖.|.|||||||||||
region    ··········‖eVVVVVVVVVVVVV
                      ^ SW V start (394, germline position 2, 94M)
```

Tie break on:

```
query     GGSGGSGGVD‖DIQLTQSPAIMSAS
germline  ----------‖QIVLTQSPAIMSAS
match               ‖.|.|||||||||||
region    ··········‖VVVVVVVVVVVVVV
                     ^ SW V start (393, germline position 1, 95M)
```

### brivekimig (Heavy chain), domain 1: V start

- V gene: vhh3-1; domain 1–133
- germline start 2 → 1, `30M1I67M` → `31M1I67M`, extension 1 → 0; score 159.9 bits (unchanged)
- Numbering: unchanged

Tie break off:

```
query     DVQLVESGGGVVQP
germline  QVQLVESGGGLVQA
match     .|||||||||.||.
region    eVVVVVVVVVVVVV
           ^ SW V start (2, germline position 2, 30M1I67M)
```

Tie break on:

```
query     DVQLVESGGGVVQP
germline  QVQLVESGGGLVQA
match     .|||||||||.||.
region    VVVVVVVVVVVVVV
          ^ SW V start (1, germline position 1, 31M1I67M)
```

### catumaxomab (chain), domain 3: V start

- V gene: IGHV-22XF; domain 254–378
- germline start 2 → 1, `97M` → `98M`, extension 1 → 0; score 177.2 bits (unchanged)
- Numbering: unchanged

Tie break off:

```
query     VTVSSGGGGS‖DIKLQQSGAELARP
germline  ----------‖QVQLQQSGAELARP
match               ‖...|||||||||||
region    ··········‖eVVVVVVVVVVVVV
                      ^ SW V start (255, germline position 2, 97M)
```

Tie break on:

```
query     VTVSSGGGGS‖DIKLQQSGAELARP
germline  ----------‖QVQLQQSGAELARP
match               ‖...|||||||||||
region    ··········‖VVVVVVVVVVVVVV
                     ^ SW V start (254, germline position 1, 98M)
```

### clivatuzumab tetraxetan (anti-MUC1 light chain 1), domain 1: V start

- V gene: IGKV-3GON; domain 1–215
- germline start 2 → 1, `95M` → `96M`, extension 1 → 0; score 152.8 bits (unchanged)
- Numbering: unchanged

Tie break off:

```
query     DIQLTQSPSSLSAS
germline  QIVLTQSPAIMSAS
match     .|.|||||...|||
region    eVVVVVVVVVVVVV
           ^ SW V start (2, germline position 2, 95M)
```

Tie break on:

```
query     DIQLTQSPSSLSAS
germline  QIVLTQSPAIMSAS
match     .|.|||||...|||
region    VVVVVVVVVVVVVV
          ^ SW V start (1, germline position 1, 96M)
```

### cosfroviximab (anti-Zaire Ebolavirus GP heavy chain 1), domain 1: V start

- V gene: IGHV-KRVC; domain 1–452
- germline start 2 → 1, `99M` → `100M`, extension 1 → 0; score 153.5 bits (unchanged)
- Numbering: unchanged

Tie break off:

```
query     DVKLLESGGGLVQP
germline  QVTLKESGPGILKP
match     .|.|.|||.|...|
region    eVVVVVVVVVVVVV
           ^ SW V start (2, germline position 2, 99M)
```

Tie break on:

```
query     DVKLLESGGGLVQP
germline  QVTLKESGPGILKP
match     .|.|.|||.|...|
region    VVVVVVVVVVVVVV
          ^ SW V start (1, germline position 1, 100M)
```

### danvilostomig (Heavy chain), domain 3: V start

- V gene: IGKV1-27*01; domain 593–701
- germline start 4 → 2, `86M` → `88M`, extension 3 → 1; score 144.0 bits (unchanged)
- Numbering: unchanged

Tie break off:

```
query     GGGGSGGGGS‖ALVMTQSPSSLSASVG
germline  ----------‖DIQMTQSPSSLSASVG
match               ‖...|||||||||||||
region    ··········‖eeeVVVVVVVVVVVVV
                        ^ SW V start (596, germline position 4, 86M)
```

Tie break on:

```
query     GGGGSGGGGS‖ALVMTQSPSSLSASVG
germline  ----------‖DIQMTQSPSSLSASVG
match               ‖...|||||||||||||
region    ··········‖eVVVVVVVVVVVVVVV
                      ^ SW V start (594, germline position 2, 88M)
```

### elipovimab (anti-HIV-1 gp120 light chain 1), domain 1: V start

- V gene: IGLV3-21*02; domain 1–211
- germline start 10 → 9, `54M3I29M` → `55M3I29M`, extension 2 → 1; score 118.2 bits (unchanged)
- Numbering: unchanged

Tie break off:

```
query     SDISVAPGETARISC
germline  PSVSVAPGQTARITC
match     ...|||||.||||.|
region    eeVVVVVVVVVVVVV
            ^ SW V start (3, germline position 10, 54M3I29M)
```

Tie break on:

```
query     SDISVAPGETARISC
germline  PSVSVAPGQTARITC
match     ...|||||.||||.|
region    eVVVVVVVVVVVVVV
           ^ SW V start (2, germline position 9, 55M3I29M)
```

### emfizatamab (anti-CD3E heavy chain 1), domain 6: V start

- V gene: IGHV3-66*01; domain 984–1119
- germline start 4 → 3, `94M` → `95M`, extension 3 → 2; score 144.3 bits (unchanged)
- Numbering: unchanged

Tie break off:

```
query     KGGGGSGGGG‖SRSLVESGGGLVQPGG
germline  ----------‖EVQLVESGGGLVQPGG
match               ‖...|||||||||||||
region    ··········‖eeeVVVVVVVVVVVVV
                        ^ SW V start (987, germline position 4, 94M)
```

Tie break on:

```
query     KGGGGSGGGG‖SRSLVESGGGLVQPGG
germline  ----------‖EVQLVESGGGLVQPGG
match               ‖...|||||||||||||
region    ··········‖eeVVVVVVVVVVVVVV
                       ^ SW V start (986, germline position 3, 95M)
```

### gontivimab (scVH-VH-VH chain), domain 1: V start

- V gene: vhh3-1; domain 1–141
- germline start 2 → 1, `97M` → `98M`, extension 1 → 0; score 140.1 bits (unchanged)
- Numbering: unchanged

Tie break off:

```
query     DVQLVESGGGLVQA
germline  QVQLVESGGGLVQA
match     .|||||||||||||
region    eVVVVVVVVVVVVV
           ^ SW V start (2, germline position 2, 97M)
```

Tie break on:

```
query     DVQLVESGGGLVQA
germline  QVQLVESGGGLVQA
match     .|||||||||||||
region    VVVVVVVVVVVVVV
          ^ SW V start (1, germline position 1, 98M)
```

### gumokimab (heavy chain 1), domain 1: V start

- V gene: IGHV4-30-4*01; domain 1–450
- germline start 2 → 1, `30M1D67M` → `31M1D67M`, extension 1 → 0; score 162.7 bits (unchanged)
- Numbering: unchanged

Tie break off:

```
query     DVQLQESGPGLVKP
germline  QVQLQESGPGLVKP
match     .|||||||||||||
region    eVVVVVVVVVVVVV
           ^ SW V start (2, germline position 2, 30M1D67M)
```

Tie break on:

```
query     DVQLQESGPGLVKP
germline  QVQLQESGPGLVKP
match     .|||||||||||||
region    VVVVVVVVVVVVVV
          ^ SW V start (1, germline position 1, 31M1D67M)
```

### isecarosmab (scVH-VH' chain), domain 1: V start

- V gene: vhh3-1; domain 1–159
- germline start 2 → 1, `97M` → `98M`, extension 1 → 0; score 158.8 bits (unchanged)
- Numbering: unchanged

Tie break off:

```
query     DVQLVESGGGVVQP
germline  QVQLVESGGGLVQA
match     .|||||||||.||.
region    eVVVVVVVVVVVVV
           ^ SW V start (2, germline position 2, 97M)
```

Tie break on:

```
query     DVQLVESGGGVVQP
germline  QVQLVESGGGLVQA
match     .|||||||||.||.
region    VVVVVVVVVVVVVV
          ^ SW V start (1, germline position 1, 98M)
```

### ivicentamab (light chain 2), domain 1: V start

- V gene: IGKV1-5*03; domain 1–217
- germline start 4 → 3, `92M` → `93M`, extension 3 → 2; score 148.2 bits (unchanged)
- Numbering: unchanged

Tie break off:

```
query     AYDMTQSPSTLSASVG
germline  DIQMTQSPSTLSASVG
match     ...|||||||||||||
region    eeeVVVVVVVVVVVVV
             ^ SW V start (4, germline position 4, 92M)
```

Tie break on:

```
query     AYDMTQSPSTLSASVG
germline  DIQMTQSPSTLSASVG
match     ...|||||||||||||
region    eeVVVVVVVVVVVVVV
            ^ SW V start (3, germline position 3, 93M)
```

### izuralimab (Chain IG scFv-h-CH2-CH3 (anti-PDCD1)), domain 1: V start

- V gene: IGKV-EAXP; domain 1–127
- germline start 2 → 1, `94M` → `95M`, extension 1 → 0; score 139.1 bits (unchanged)
- Numbering: unchanged

Tie break off:

```
query     EIVLTQSPATLSAS
germline  SIVMTQTPKFLPVS
match     .||.||.|..|..|
region    eVVVVVVVVVVVVV
           ^ SW V start (2, germline position 2, 94M)
```

Tie break on:

```
query     EIVLTQSPATLSAS
germline  SIVMTQTPKFLPVS
match     .||.||.|..|..|
region    VVVVVVVVVVVVVV
          ^ SW V start (1, germline position 1, 95M)
```

### lacnotuzumab (anti-CSF1 heavy chain 1), domain 1: V start

- V gene: IGHV-O4YA; domain 1–448
- germline start 2 → 1, `96M` → `97M`, extension 1 → 0; score 164.8 bits (unchanged)
- Numbering: unchanged

Tie break off:

```
query     QVQLQESGPGLVKP
germline  DVQLQESGPGLVKP
match     .|||||||||||||
region    eVVVVVVVVVVVVV
           ^ SW V start (2, germline position 2, 96M)
```

Tie break on:

```
query     QVQLQESGPGLVKP
germline  DVQLQESGPGLVKP
match     .|||||||||||||
region    VVVVVVVVVVVVVV
          ^ SW V start (1, germline position 1, 97M)
```

### lancastotug (heavy chain 1), domain 1: V start

- V gene: IGHV4-30-4*01; domain 1–451
- germline start 2 → 1, `30M1D67M` → `31M1D67M`, extension 1 → 0; score 162.4 bits (unchanged)
- Numbering: unchanged

Tie break off:

```
query     DVQLQESGPGLVKP
germline  QVQLQESGPGLVKP
match     .|||||||||||||
region    eVVVVVVVVVVVVV
           ^ SW V start (2, germline position 2, 30M1D67M)
```

Tie break on:

```
query     DVQLQESGPGLVKP
germline  QVQLQESGPGLVKP
match     .|||||||||||||
region    VVVVVVVVVVVVVV
          ^ SW V start (1, germline position 1, 31M1D67M)
```

### lucorafusp alfa (heavy chain 1), domain 1: V start

- V gene: IGHV4-30-4*01; domain 1–604
- germline start 2 → 1, `30M1D67M` → `31M1D67M`, extension 1 → 0; score 162.4 bits (unchanged)
- Numbering: unchanged

Tie break off:

```
query     DVQLQESGPGLVKP
germline  QVQLQESGPGLVKP
match     .|||||||||||||
region    eVVVVVVVVVVVVV
           ^ SW V start (2, germline position 2, 30M1D67M)
```

Tie break on:

```
query     DVQLQESGPGLVKP
germline  QVQLQESGPGLVKP
match     .|||||||||||||
region    VVVVVVVVVVVVVV
          ^ SW V start (1, germline position 1, 31M1D67M)
```

### lunsekimig (Heavy chain), domain 1: V start

- V gene: vhh3-1; domain 1–158
- germline start 2 → 1, `97M` → `98M`, extension 1 → 0; score 160.2 bits (unchanged)
- Numbering: unchanged

Tie break off:

```
query     DVQLVESGGGVVQP
germline  QVQLVESGGGLVQA
match     .|||||||||.||.
region    eVVVVVVVVVVVVV
           ^ SW V start (2, germline position 2, 97M)
```

Tie break on:

```
query     DVQLVESGGGVVQP
germline  QVQLVESGGGLVQA
match     .|||||||||.||.
region    VVVVVVVVVVVVVV
          ^ SW V start (1, germline position 1, 98M)
```

### mipletamig (chain), domain 4: V start

- V gene: IGKV-HQL3; domain 650–757
- germline start 2 → 1, `94M` → `95M`, extension 1 → 0; score 146.1 bits (unchanged)
- Numbering: unchanged

Tie break off:

```
query     GGGGSGGGGS‖DIQMTQSPSTLSAS
germline  ----------‖QIVLTQSPAIMSAS
match               ‖.|..||||...|||
region    ··········‖eVVVVVVVVVVVVV
                      ^ SW V start (651, germline position 2, 94M)
```

Tie break on:

```
query     GGGGSGGGGS‖DIQMTQSPSTLSAS
germline  ----------‖QIVLTQSPAIMSAS
match               ‖.|..||||...|||
region    ··········‖VVVVVVVVVVVVVV
                     ^ SW V start (650, germline position 1, 95M)
```

### misitatug (light chain 1), domain 1: V start

- V gene: IGKV-GLVF; domain 1–213
- germline start 2 → 1, `94M` → `95M`, extension 1 → 0; score 158.8 bits (unchanged)
- Numbering: unchanged

Tie break off:

```
query     DVVMTQSPAFLSVT
germline  QIVLTQSPAIMSAS
match     ..|.|||||..|..
region    eVVVVVVVVVVVVV
           ^ SW V start (2, germline position 2, 94M)
```

Tie break on:

```
query     DVVMTQSPAFLSVT
germline  QIVLTQSPAIMSAS
match     ..|.|||||..|..
region    VVVVVVVVVVVVVV
          ^ SW V start (1, germline position 1, 95M)
```

### naxitamab (anti-Ganglioside GD2 light chain 1), domain 1: V start

- V gene: IGKV-2UKG; domain 1–211
- germline start 2 → 1, `93M` → `94M`, extension 1 → 0; score 151.1 bits (unchanged)
- Numbering: unchanged

Tie break off:

```
query     EIVMTQTPATLSVS
germline  SIVMTQTPKFLLVS
match     .|||||||..|.||
region    eVVVVVVVVVVVVV
           ^ SW V start (2, germline position 2, 93M)
```

Tie break on:

```
query     EIVMTQTPATLSVS
germline  SIVMTQTPKFLLVS
match     .|||||||..|.||
region    VVVVVVVVVVVVVV
          ^ SW V start (1, germline position 1, 94M)
```

### nivatrotamab (light chain), domain 1: V start

- V gene: IGKV-2UKG; domain 1–228
- germline start 2 → 1, `93M` → `94M`, extension 1 → 0; score 151.1 bits (unchanged)
- Numbering: unchanged

Tie break off:

```
query     EIVMTQTPATLSVS
germline  SIVMTQTPKFLLVS
match     .|||||||..|.||
region    eVVVVVVVVVVVVV
           ^ SW V start (2, germline position 2, 93M)
```

Tie break on:

```
query     EIVMTQTPATLSVS
germline  SIVMTQTPKFLLVS
match     .|||||||..|.||
region    VVVVVVVVVVVVVV
          ^ SW V start (1, germline position 1, 94M)
```

### regdanvimab (anti-SARS-CoV-2 Spike RBD light chain 1), domain 1: V start

- V gene: IGLV1-51*01; domain 1–216
- germline start 3 → 1, `96M` → `98M`, extension 2 → 0; score 184.6 bits (unchanged)
- Numbering: unchanged

Tie break off:

```
query     ELVLTQPPSVSAAPG
germline  QSVLTQPPSVSAAPG
match     ..|||||||||||||
region    eeVVVVVVVVVVVVV
            ^ SW V start (3, germline position 3, 96M)
```

Tie break on:

```
query     ELVLTQPPSVSAAPG
germline  QSVLTQPPSVSAAPG
match     ..|||||||||||||
region    VVVVVVVVVVVVVVV
          ^ SW V start (1, germline position 1, 98M)
```

### resugosbart (light chain 1), domain 1: V start

- V gene: IGKV-2UKG; domain 1–214
- germline start 2 → 1, `94M` → `95M`, extension 1 → 0; score 150.7 bits (unchanged)
- Numbering: unchanged

Tie break off:

```
query     DIVMTQSPSSLSAS
germline  SIVMTQTPKFLLVS
match     .|||||.|..|..|
region    eVVVVVVVVVVVVV
           ^ SW V start (2, germline position 2, 94M)
```

Tie break on:

```
query     DIVMTQSPSSLSAS
germline  SIVMTQTPKFLLVS
match     .|||||.|..|..|
region    VVVVVVVVVVVVVV
          ^ SW V start (1, germline position 1, 95M)
```

### rimteravimab (heavy chain 1), domain 1: V start

- V gene: vhh3-1; domain 1–361
- germline start 2 → 1, `97M` → `98M`, extension 1 → 0; score 174.7 bits (unchanged)
- Numbering: unchanged

Tie break off:

```
query     DVQLVESGGGLVQP
germline  QVQLVESGGGLVQA
match     .||||||||||||.
region    eVVVVVVVVVVVVV
           ^ SW V start (2, germline position 2, 97M)
```

Tie break on:

```
query     DVQLVESGGGLVQP
germline  QVQLVESGGGLVQA
match     .||||||||||||.
region    VVVVVVVVVVVVVV
          ^ SW V start (1, germline position 1, 98M)
```

### solitomab (chain), domain 3: V start

- V gene: IGHV1-46*04; domain 254–390
- germline start 2 → 1, `97M` → `98M`, extension 1 → 0; score 153.2 bits (unchanged)
- Numbering: unchanged

Tie break off:

```
query     VTVSSGGGGS‖DVQLVQSGAEVKKP
germline  ----------‖QVQLVQSGAEVKKP
match               ‖.|||||||||||||
region    ··········‖eVVVVVVVVVVVVV
                      ^ SW V start (255, germline position 2, 97M)
```

Tie break on:

```
query     VTVSSGGGGS‖DVQLVQSGAEVKKP
germline  ----------‖QVQLVQSGAEVKKP
match               ‖.|||||||||||||
region    ··········‖VVVVVVVVVVVVVV
                     ^ SW V start (254, germline position 1, 98M)
```

### solitomab (chain), domain 4: V start

- V gene: IGKV-O6WO; domain 391–502
- germline start 2 → 1, `94M` → `95M`, extension 1 → 0; score 153.5 bits (unchanged)
- Numbering: unchanged

Tie break off:

```
query     GGSGGSGGAD‖DIVLTQSPATLSLS
germline  ----------‖QIVLTQSPAIMSAS
match               ‖.||||||||..|.|
region    ··········‖eVVVVVVVVVVVVV
                      ^ SW V start (392, germline position 2, 94M)
```

Tie break on:

```
query     GGSGGSGGAD‖DIVLTQSPATLSLS
germline  ----------‖QIVLTQSPAIMSAS
match               ‖.||||||||..|.|
region    ··········‖VVVVVVVVVVVVVV
                     ^ SW V start (391, germline position 1, 95M)
```

### sonelokimab (chain), domain 1: V start

- V gene: vhh3-1; domain 1–132
- germline start 2 → 1, `97M` → `98M`, extension 1 → 0; score 157.1 bits (unchanged)
- Numbering: unchanged

Tie break off:

```
query     DVQLVESGGGLVQP
germline  QVQLVESGGGLVQA
match     .||||||||||||.
region    eVVVVVVVVVVVVV
           ^ SW V start (2, germline position 2, 97M)
```

Tie break on:

```
query     DVQLVESGGGLVQP
germline  QVQLVESGGGLVQA
match     .||||||||||||.
region    VVVVVVVVVVVVVV
          ^ SW V start (1, germline position 1, 98M)
```

### tilogotamab (anti-TNFRSF10B light chain 1), domain 1: V start

- V gene: IGKV-K25A; domain 1–213
- germline start 2 → 1, `94M` → `95M`, extension 1 → 0; score 154.6 bits (unchanged)
- Numbering: unchanged

Tie break off:

```
query     DIQLTQSPSSLSAS
germline  QIVLTQSPAILSAS
match     .|.|||||..||||
region    eVVVVVVVVVVVVV
           ^ SW V start (2, germline position 2, 94M)
```

Tie break on:

```
query     DIQLTQSPSSLSAS
germline  QIVLTQSPAILSAS
match     .|.|||||..||||
region    VVVVVVVVVVVVVV
          ^ SW V start (1, germline position 1, 95M)
```

### timolumab (anti-AOC3 light chain 1), domain 1: V start

- V gene: IGKV1D-13*02; domain 1–214
- germline start 2 → 1, `94M` → `95M`, extension 1 → 0; score 171.9 bits (unchanged)
- Numbering: unchanged

Tie break off:

```
query     VIQLTQSPSSLSAS
germline  AIQLTQSPSSLSAS
match     .|||||||||||||
region    eVVVVVVVVVVVVV
           ^ SW V start (2, germline position 2, 94M)
```

Tie break on:

```
query     VIQLTQSPSSLSAS
germline  AIQLTQSPSSLSAS
match     .|||||||||||||
region    VVVVVVVVVVVVVV
          ^ SW V start (1, germline position 1, 95M)
```

### tislelizumab (anti-PDCD1 light chain 1), domain 1: V start

- V gene: IGKV-2UKG; domain 1–214
- germline start 2 → 1, `94M` → `95M`, extension 1 → 0; score 140.1 bits (unchanged)
- Numbering: unchanged

Tie break off:

```
query     DIVMTQSPDSLAVS
germline  SIVMTQTPKFLLVS
match     .|||||.|..|.||
region    eVVVVVVVVVVVVV
           ^ SW V start (2, germline position 2, 94M)
```

Tie break on:

```
query     DIVMTQSPDSLAVS
germline  SIVMTQTPKFLLVS
match     .|||||.|..|.||
region    VVVVVVVVVVVVVV
          ^ SW V start (1, germline position 1, 95M)
```

### tregalizumab (anti-CD4 heavy chain 1), domain 1: V start

- V gene: IGHV-MG5V; domain 1–454
- germline start 3 → 1, `97M` → `99M`, extension 2 → 0; score 149.3 bits (unchanged)
- Numbering: unchanged

Tie break off:

```
query     EEQLVESGGGLVKPG
germline  QVQLVETGGGLVRPG
match     ..||||.|||||.||
region    eeVVVVVVVVVVVVV
            ^ SW V start (3, germline position 3, 97M)
```

Tie break on:

```
query     EEQLVESGGGLVKPG
germline  QVQLVETGGGLVRPG
match     ..||||.|||||.||
region    VVVVVVVVVVVVVVV
          ^ SW V start (1, germline position 1, 99M)
```

### trontinemab (Heavy chain anti-APP A beta 1-40/42 and anti-TFRC), domain 2: V start

- V gene: IGHV4-4*08; domain 473–703
- germline start 4 → 3, `8M1I59M2D25M` → `9M1I59M2D25M`, extension 3 → 2; score 123.9 bits (unchanged)
- Numbering: unchanged

Tie break off:

```
query     SGGGGSGGGG‖SQSMQESGPGLVVKPS
germline  ----------‖QVQLQESGPGL-VKPS
match               ‖....||||||| ||||
region    ··········‖eeeVVVVVVVVVVVVV
                        ^ SW V start (476, germline position 4, 8M1I59M2D25M)
```

Tie break on:

```
query     SGGGGSGGGG‖SQSMQESGPGLVVKPS
germline  ----------‖QVQLQESGPGL-VKPS
match               ‖....||||||| ||||
region    ··········‖eeVVVVVVVVVVVVVV
                       ^ SW V start (475, germline position 3, 9M1I59M2D25M)
```

### umizortamig (heavy chain 1), domain 6: V start

- V gene: IGHV3-66*01; domain 979–1119
- germline start 4 → 3, `30M1I19M2I45M` → `31M1I19M2I45M`, extension 3 → 2; score 153.2 bits (unchanged)
- Numbering: unchanged

Tie break off:

```
query     KGGGGSGGGG‖SQSLVESGGGLVQPGG
germline  ----------‖EVQLVESGGGLVQPGG
match               ‖...|||||||||||||
region    ··········‖eeeVVVVVVVVVVVVV
                        ^ SW V start (982, germline position 4, 30M1I19M2I45M)
```

Tie break on:

```
query     KGGGGSGGGG‖SQSLVESGGGLVQPGG
germline  ----------‖EVQLVESGGGLVQPGG
match               ‖...|||||||||||||
region    ··········‖eeVVVVVVVVVVVVVV
                       ^ SW V start (981, germline position 3, 31M1I19M2I45M)
```

### umizortamig (heavy chain 1), domain 7: V start

- V gene: IGKV1-5*03; domain 1120–1231
- germline start 4 → 2, `86M` → `88M`, extension 3 → 1; score 142.9 bits (unchanged)
- Numbering: unchanged

Tie break off:

```
query     GGGGSGGGGS‖ALVMTQSPSTLSASVG
germline  ----------‖DIQMTQSPSTLSASVG
match               ‖...|||||||||||||
region    ··········‖eeeVVVVVVVVVVVVV
                        ^ SW V start (1123, germline position 4, 86M)
```

Tie break on:

```
query     GGGGSGGGGS‖ALVMTQSPSTLSASVG
germline  ----------‖DIQMTQSPSTLSASVG
match               ‖...|||||||||||||
region    ··········‖eVVVVVVVVVVVVVVV
                      ^ SW V start (1121, germline position 2, 88M)
```

### vatelizumab (anti-ITGA2 light chain 1), domain 1: V start

- V gene: IGKV-O6WO; domain 1–213
- germline start 3 → 1, `93M` → `95M`, extension 2 → 0; score 144.7 bits (unchanged)
- Numbering: unchanged

Tie break off:

```
query     DFVMTQSPAFLSVTP
germline  QIVLTQSPAIMSASP
match     ..|.|||||..|..|
region    eeVVVVVVVVVVVVV
            ^ SW V start (3, germline position 3, 93M)
```

Tie break on:

```
query     DFVMTQSPAFLSVTP
germline  QIVLTQSPAIMSASP
match     ..|.|||||..|..|
region    VVVVVVVVVVVVVVV
          ^ SW V start (1, germline position 1, 95M)
```

### veltuzumab (anti-MS4A1 light chain 1), domain 1: V start

- V gene: IGKV-27KA; domain 1–213
- germline start 2 → 1, `94M` → `95M`, extension 1 → 0; score 148.9 bits (unchanged)
- Numbering: unchanged

Tie break off:

```
query     DIQLTQSPSSLSAS
germline  QIVLSQSPAILSAS
match     .|.|.|||..||||
region    eVVVVVVVVVVVVV
           ^ SW V start (2, germline position 2, 94M)
```

Tie break on:

```
query     DIQLTQSPSSLSAS
germline  QIVLSQSPAILSAS
match     .|.|.|||..||||
region    VVVVVVVVVVVVVV
          ^ SW V start (1, germline position 1, 95M)
```

### vudalimab (Chain scFv-h-CH2-CH3), domain 1: V start

- V gene: IGKV-EAXP; domain 1–127
- germline start 2 → 1, `94M` → `95M`, extension 1 → 0; score 139.1 bits (unchanged)
- Numbering: unchanged

Tie break off:

```
query     EIVLTQSPATLSAS
germline  SIVMTQTPKFLPVS
match     .||.||.|..|..|
region    eVVVVVVVVVVVVV
           ^ SW V start (2, germline position 2, 94M)
```

Tie break on:

```
query     EIVLTQSPATLSAS
germline  SIVMTQTPKFLPVS
match     .||.||.|..|..|
region    VVVVVVVVVVVVVV
          ^ SW V start (1, germline position 1, 95M)
```

### yttrium y-90 clivatuzumab tetraxetan (anti-MUC1 light chain 1), domain 1: V start

- V gene: IGKV-3GON; domain 1–215
- germline start 2 → 1, `95M` → `96M`, extension 1 → 0; score 152.8 bits (unchanged)
- Numbering: unchanged

Tie break off:

```
query     DIQLTQSPSSLSAS
germline  QIVLTQSPAIMSAS
match     .|.|||||...|||
region    eVVVVVVVVVVVVV
           ^ SW V start (2, germline position 2, 95M)
```

Tie break on:

```
query     DIQLTQSPSSLSAS
germline  QIVLTQSPAIMSAS
match     .|.|||||...|||
region    VVVVVVVVVVVVVV
          ^ SW V start (1, germline position 1, 96M)
```

### zeleciment (heavy chain 1), domain 1: V start

- V gene: IGHV-FMN7; domain 1–227
- germline start 2 → 1, `97M` → `98M`, extension 1 → 0; score 165.5 bits (unchanged)
- Numbering: unchanged

Tie break off:

```
query     QVQLQESGPGLVKP
germline  DVQLQESGPGLVKP
match     .|||||||||||||
region    eVVVVVVVVVVVVV
           ^ SW V start (2, germline position 2, 97M)
```

Tie break on:

```
query     QVQLQESGPGLVKP
germline  DVQLQESGPGLVKP
match     .|||||||||||||
region    VVVVVVVVVVVVVV
          ^ SW V start (1, germline position 1, 98M)
```

### zeleciment basivarsen (heavy chain 1), domain 1: V start

- V gene: IGHV-FMN7; domain 1–227
- germline start 2 → 1, `97M` → `98M`, extension 1 → 0; score 165.5 bits (unchanged)
- Numbering: unchanged

Tie break off:

```
query     QVQLQESGPGLVKP
germline  DVQLQESGPGLVKP
match     .|||||||||||||
region    eVVVVVVVVVVVVV
           ^ SW V start (2, germline position 2, 97M)
```

Tie break on:

```
query     QVQLQESGPGLVKP
germline  DVQLQESGPGLVKP
match     .|||||||||||||
region    VVVVVVVVVVVVVV
          ^ SW V start (1, germline position 1, 98M)
```

### zeleciment rostudirsen (heavy chain 1), domain 1: V start

- V gene: IGHV-FMN7; domain 1–227
- germline start 2 → 1, `97M` → `98M`, extension 1 → 0; score 165.5 bits (unchanged)
- Numbering: unchanged

Tie break off:

```
query     QVQLQESGPGLVKP
germline  DVQLQESGPGLVKP
match     .|||||||||||||
region    eVVVVVVVVVVVVV
           ^ SW V start (2, germline position 2, 97M)
```

Tie break on:

```
query     QVQLQESGPGLVKP
germline  DVQLQESGPGLVKP
match     .|||||||||||||
region    VVVVVVVVVVVVVV
          ^ SW V start (1, germline position 1, 98M)
```

### zifibancimig (heavy chain 1), domain 1: V start

- V gene: IGHV3-21*07; domain 1–223
- germline start 4 → 3, `24M1D70M` → `25M1D70M`, extension 3 → 2; score 128.8 bits (unchanged)
- Numbering: unchanged

Tie break off:

```
query     SEHLVESGGGLVKPGG
germline  EVQLVESGGGLVKPGG
match     ...|||||||||||||
region    eeeVVVVVVVVVVVVV
             ^ SW V start (4, germline position 4, 24M1D70M)
```

Tie break on:

```
query     SEHLVESGGGLVKPGG
germline  EVQLVESGGGLVKPGG
match     ...|||||||||||||
region    eeVVVVVVVVVVVVVV
            ^ SW V start (3, germline position 3, 25M1D70M)
```

### zuberitamab (anti-MS4A1 light chain 1), domain 1: V start

- V gene: IGKV-27KA; domain 1–213
- germline start 2 → 1, `94M` → `95M`, extension 1 → 0; score 173.6 bits (unchanged)
- Numbering: unchanged

Tie break off:

```
query     DIELSQSPAILSAS
germline  QIVLSQSPAILSAS
match     .|.|||||||||||
region    eVVVVVVVVVVVVV
           ^ SW V start (2, germline position 2, 94M)
```

Tie break on:

```
query     DIELSQSPAILSAS
germline  QIVLSQSPAILSAS
match     .|.|||||||||||
region    VVVVVVVVVVVVVV
          ^ SW V start (1, germline position 1, 95M)
```

## J ends, numbering unchanged (16)

### docaravimab (anti-Rabies Virus Strain ERA GP Ectodomain Epitope G-III light chain 1), domain 1: J end

- J gene: IGKJ-WRAP; domain 1–213
- germline end 11 → 12 (of 12), `11M` → `12M`, extension 1 → 0; score 26.1 bits (unchanged)
- Numbering: unchanged

Tie break off:

```
query     LPPTFGGGTKLEIKRADAAPTVS
germline  --PTFGGGTKLEIN---------
match       |||||||||||.         
region    ··JJJJJJJJJJJe·········
                      ^ SW J end (106, germline position 11 of 12, …11M)
```

Tie break on:

```
query     LPPTFGGGTKLEIKRADAAPTVS
germline  --PTFGGGTKLEIN---------
match       |||||||||||.         
region    ··JJJJJJJJJJJJ·········
                       ^ SW J end (107, germline position 12 of 12, …12M)
```

### enristomig (Heavy chain), domain 1: J end

- J gene: ighJ-5; domain 1–117
- germline end 14 → 15 (of 15), `12M` → `13M`, extension 1 → 0; score 29.3 bits (unchanged)
- Numbering: unchanged

Tie break off:

```
query     VFEYWGQGTLVTVKPGGSGGS‖EV
germline  -FEYWGQGTLVTVS-------‖--
match      ||||||||||||.       ‖  
region    ·JJJJJJJJJJJJe·······‖··
                      ^ SW J end (109, germline position 14 of 15, …12M)
```

Tie break on:

```
query     VFEYWGQGTLVTVKPGGSGGS‖EV
germline  -FEYWGQGTLVTVS-------‖--
match      ||||||||||||.       ‖  
region    ·JJJJJJJJJJJJJ·······‖··
                       ^ SW J end (110, germline position 15 of 15, …13M)
```

### enristomig (Heavy chain), domain 2: J end

- J gene: ighJ-5; domain 118–466
- germline end 14 → 15 (of 15), `14M` → `15M`, extension 1 → 0; score 27.1 bits (unchanged)
- Numbering: unchanged

Tie break off:

```
query     SVAYWGQGTLVTVKPGGGGDKTH
germline  QFEYWGQGTLVTVS---------
match     ...||||||||||.         
region    JJJJJJJJJJJJJe·········
                      ^ SW J end (236, germline position 14 of 15, …14M)
```

Tie break on:

```
query     SVAYWGQGTLVTVKPGGGGDKTH
germline  QFEYWGQGTLVTVS---------
match     ...||||||||||.         
region    JJJJJJJJJJJJJJ·········
                       ^ SW J end (237, germline position 15 of 15, …15M)
```

### futuximab (anti-EGFR light chain 1), domain 1: J end

- J gene: IGKJ-WRAP; domain 1–214
- germline end 11 → 12 (of 12), `11M` → `12M`, extension 1 → 0; score 26.1 bits (unchanged)
- Numbering: unchanged

Tie break off:

```
query     VPPTFGGGTKLEIKRTVAAPSVF
germline  --PTFGGGTKLEIN---------
match       |||||||||||.         
region    ··JJJJJJJJJJJe·········
                      ^ SW J end (106, germline position 11 of 12, …11M)
```

Tie break on:

```
query     VPPTFGGGTKLEIKRTVAAPSVF
germline  --PTFGGGTKLEIN---------
match       |||||||||||.         
region    ··JJJJJJJJJJJJ·········
                       ^ SW J end (107, germline position 12 of 12, …12M)
```

### lecanemab (anti-APP light chain 1), domain 1: J end

- J gene: IGKJ-WRAP; domain 1–219
- germline end 11 → 12 (of 12), `11M` → `12M`, extension 1 → 0; score 23.3 bits (unchanged)
- Numbering: unchanged

Tie break off:

```
query     VPPTFGPGTKLEIKRTVAAPSVF
germline  --PTFGGGTKLEIN---------
match       ||||.||||||.         
region    ··JJJJJJJJJJJe·········
                      ^ SW J end (111, germline position 11 of 12, …11M)
```

Tie break on:

```
query     VPPTFGPGTKLEIKRTVAAPSVF
germline  --PTFGGGTKLEIN---------
match       ||||.||||||.         
region    ··JJJJJJJJJJJJ·········
                       ^ SW J end (112, germline position 12 of 12, …12M)
```

### lemalesomab (light chain 1), domain 1: J end

- J gene: IGKJ-WRAP; domain 1–218
- germline end 11 → 12 (of 12), `11M` → `12M`, extension 1 → 0; score 26.1 bits (unchanged)
- Numbering: unchanged

Tie break off:

```
query     VPPTFGGGTKLEIKRADAAPTVS
germline  --PTFGGGTKLEIN---------
match       |||||||||||.         
region    ··JJJJJJJJJJJe·········
                      ^ SW J end (111, germline position 11 of 12, …11M)
```

Tie break on:

```
query     VPPTFGGGTKLEIKRADAAPTVS
germline  --PTFGGGTKLEIN---------
match       |||||||||||.         
region    ··JJJJJJJJJJJJ·········
                       ^ SW J end (112, germline position 12 of 12, …12M)
```

### muromonab-cd3 (anti-CD3E light chain 1), domain 1: J end

- J gene: IGKJ-Z5J4; domain 1–106
- germline end 11 → 12 (of 12), `11M` → `12M`, extension 1 → 0; score 25.0 bits (unchanged)
- Numbering: unchanged

Tie break off:

```
query     NPFTFGSGTKLEIN
germline  --FTFGSGTKLEIK
match       |||||||||||.
region    ··JJJJJJJJJJJe
                      ^ SW J end (105, germline position 11 of 12, …11M)
```

Tie break on:

```
query     NPFTFGSGTKLEIN
germline  --FTFGSGTKLEIK
match       |||||||||||.
region    ··JJJJJJJJJJJJ
                       ^ SW J end (106, germline position 12 of 12, …12M)
```

### naptumomab estafenatox (Light chain), domain 1: J end

- J gene: IGKJ-WRAP; domain 1–214
- germline end 11 → 12 (of 12), `11M` → `12M`, extension 1 → 0; score 26.1 bits (unchanged)
- Numbering: unchanged

Tie break off:

```
query     SPPTFGGGTKLEIKRADAAPTVS
germline  --PTFGGGTKLEIN---------
match       |||||||||||.         
region    ··JJJJJJJJJJJe·········
                      ^ SW J end (106, germline position 11 of 12, …11M)
```

Tie break on:

```
query     SPPTFGGGTKLEIKRADAAPTVS
germline  --PTFGGGTKLEIN---------
match       |||||||||||.         
region    ··JJJJJJJJJJJJ·········
                       ^ SW J end (107, germline position 12 of 12, …12M)
```

### ordesekimab (anti-IL15 light chain 1), domain 1: J end

- J gene: IGKJ2*01; domain 1–214
- germline end 11 → 12 (of 12), `11M` → `12M`, extension 1 → 0; score 24.0 bits (unchanged)
- Numbering: unchanged

Tie break off:

```
query     SSHTFGQGTKLEISRTVAAPSVF
germline  --YTFGQGTKLEIK---------
match       .||||||||||.         
region    ··JJJJJJJJJJJe·········
                      ^ SW J end (106, germline position 11 of 12, …11M)
```

Tie break on:

```
query     SSHTFGQGTKLEISRTVAAPSVF
germline  --YTFGQGTKLEIK---------
match       .||||||||||.         
region    ··JJJJJJJJJJJJ·········
                       ^ SW J end (107, germline position 12 of 12, …12M)
```

### ozekibart (heavy chain 1), domain 1: J end

- J gene: ighJ-4; domain 1–130
- germline end 13 → 14 (of 15), `13M` → `14M`, extension 2 → 1; score 28.2 bits (unchanged)
- Numbering: unchanged

Tie break off:

```
query     KYDYWGQGTLVTVKPGGSGGSE‖VQ
germline  EYDYWGQGTQVTVSS-------‖--
match     .||||||||.|||..       ‖  
region    JJJJJJJJJJJJJee·······‖··
                      ^ SW J end (121, germline position 13 of 15, …13M)
```

Tie break on:

```
query     KYDYWGQGTLVTVKPGGSGGSE‖VQ
germline  EYDYWGQGTQVTVSS-------‖--
match     .||||||||.|||..       ‖  
region    JJJJJJJJJJJJJJe·······‖··
                       ^ SW J end (122, germline position 14 of 15, …14M)
```

### ozekibart (heavy chain 1), domain 2: J end

- J gene: ighJ-4; domain 131–481
- germline end 13 → 14 (of 15), `13M` → `14M`, extension 2 → 1; score 28.2 bits (unchanged)
- Numbering: unchanged

Tie break off:

```
query     KYDYWGQGTLVTVKPGGGGDKTHT
germline  EYDYWGQGTQVTVSS---------
match     .||||||||.|||..         
region    JJJJJJJJJJJJJee·········
                      ^ SW J end (251, germline position 13 of 15, …13M)
```

Tie break on:

```
query     KYDYWGQGTLVTVKPGGGGDKTHT
germline  EYDYWGQGTQVTVSS---------
match     .||||||||.|||..         
region    JJJJJJJJJJJJJJe·········
                       ^ SW J end (252, germline position 14 of 15, …14M)
```

### ruplizumab (anti-CD40LG light chain 1), domain 1: J end

- J gene: IGKJ-WRAP; domain 1–218
- germline end 11 → 12 (of 12), `11M` → `12M`, extension 1 → 0; score 26.1 bits (unchanged)
- Numbering: unchanged

Tie break off:

```
query     IPPTFGGGTKLEIKRTVAAPSVF
germline  --PTFGGGTKLEIN---------
match       |||||||||||.         
region    ··JJJJJJJJJJJe·········
                      ^ SW J end (110, germline position 11 of 12, …11M)
```

Tie break on:

```
query     IPPTFGGGTKLEIKRTVAAPSVF
germline  --PTFGGGTKLEIN---------
match       |||||||||||.         
region    ··JJJJJJJJJJJJ·········
                       ^ SW J end (111, germline position 12 of 12, …12M)
```

### tegoprubart (light chain 1), domain 1: J end

- J gene: IGKJ-WRAP; domain 1–218
- germline end 11 → 12 (of 12), `11M` → `12M`, extension 1 → 0; score 26.1 bits (unchanged)
- Numbering: unchanged

Tie break off:

```
query     IPPTFGGGTKLEIKRTVAAPSVF
germline  --PTFGGGTKLEIN---------
match       |||||||||||.         
region    ··JJJJJJJJJJJe·········
                      ^ SW J end (110, germline position 11 of 12, …11M)
```

Tie break on:

```
query     IPPTFGGGTKLEIKRTVAAPSVF
germline  --PTFGGGTKLEIN---------
match       |||||||||||.         
region    ··JJJJJJJJJJJJ·········
                       ^ SW J end (111, germline position 12 of 12, …12M)
```

### tomuzotuximab (anti-EGFR heavy chain 1), domain 1: J end

- J gene: IGHJ-KZ3G; domain 1–449
- germline end 14 → 15 (of 15), `13M` → `14M`, extension 1 → 0; score 30.3 bits (unchanged)
- Numbering: unchanged

Tie break off:

```
query     FAYWGQGTLVTVSTASTKGPSVF
germline  FAYWGQGTLVTVSA---------
match     |||||||||||||.         
region    JJJJJJJJJJJJJe·········
                      ^ SW J end (118, germline position 14 of 15, …13M)
```

Tie break on:

```
query     FAYWGQGTLVTVSTASTKGPSVF
germline  FAYWGQGTLVTVSA---------
match     |||||||||||||.         
region    JJJJJJJJJJJJJJ·········
                       ^ SW J end (119, germline position 15 of 15, …14M)
```

### veligrotug (light chain 1), domain 1: J end

- J gene: IGKJ-WRAP; domain 1–219
- germline end 11 → 12 (of 12), `11M` → `12M`, extension 1 → 0; score 26.1 bits (unchanged)
- Numbering: unchanged

Tie break off:

```
query     VPPTFGGGTKLEIKRTVAAPSVF
germline  --PTFGGGTKLEIN---------
match       |||||||||||.         
region    ··JJJJJJJJJJJe·········
                      ^ SW J end (111, germline position 11 of 12, …11M)
```

Tie break on:

```
query     VPPTFGGGTKLEIKRTVAAPSVF
germline  --PTFGGGTKLEIN---------
match       |||||||||||.         
region    ··JJJJJJJJJJJJ·········
                       ^ SW J end (112, germline position 12 of 12, …12M)
```

### zolimomab (anti-CD3E light chain 1), domain 1: J end

- J gene: IGKJ-Z5J4; domain 1–106
- germline end 11 → 12 (of 12), `11M` → `12M`, extension 1 → 0; score 25.0 bits (unchanged)
- Numbering: unchanged

Tie break off:

```
query     NPFTFGSGTKLEIN
germline  --FTFGSGTKLEIK
match       |||||||||||.
region    ··JJJJJJJJJJJe
                      ^ SW J end (105, germline position 11 of 12, …11M)
```

Tie break on:

```
query     NPFTFGSGTKLEIN
germline  --FTFGSGTKLEIK
match       |||||||||||.
region    ··JJJJJJJJJJJJ
                       ^ SW J end (106, germline position 12 of 12, …12M)
```
