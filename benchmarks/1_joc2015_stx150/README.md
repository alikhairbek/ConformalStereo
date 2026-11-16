# ConformalStereo: cross-check of published SI data, version 3

14 November 2026. Version 3 replaces versions 1 and 2. New in it: DP4+ and MM-DP4+ as coded in DP4+App 2.2.3,
with its parameter bank; the probabilities printed in the SI of the 2023 and 2025 papers recomputed from their
printed shieldings; two inputs of DP4+App that no SI prints (pairs of exchangeable nuclei, class of each nucleus)
recovered where the printed probabilities need them; a test that allows for the rounding of the printed shieldings;
and, with all that, a new status for every STX150 record.

Four pieces of work:

1. The STX150 validation set (110 cases, J. Chem. Inf. Model. 2026, 66, 11062) tested against the SI of the two
   earlier papers its cases come from, and against itself.
2. The probabilities of the two earlier papers recomputed from their SI (36 levels of the 2023 paper, three levels
   of the 2025 paper).
3. A status for every STX150 (case, level) record: verified, close, corrected, unverified or excluded.
4. An audit of the probabilities published with the DP4-AI paper (unchanged from version 1).

## Run it

    pdftotext -layout np3c00566_si_001.pdf np3c00566_si_001.txt      # SI of J. Nat. Prod. 2023, 86, 2360
    pdftotext -layout d4ob01959h1_suppl.pdf d4ob01959h1_suppl.txt    # SI part I of Org. Biomol. Chem. 2025, 23, 2480
    sh run_all.sh np3c00566_si_001.txt d4ob01959h1_suppl.txt out [python]

Python 3.9 or later. Steps 11 and 21 to 27 need numpy and scipy; the rest uses the standard library. Steps 21 to 24
draw random numbers with fixed seeds; the whole run takes about four minutes. The two SI files are not in this
archive. `inputs/` holds what was read earlier from the STX150 SI (shieldings of every printed
candidate at six levels, the candidate sets, the published table, the probabilities recomputed from the printed
shieldings), the DP4+ parameters of the 2015 and 2026 papers, and `app_parameters.csv`: the parameter bank of
DP4+App 2.2.3 (24 QM, 36 MM and 3 xTB levels), read from its spreadsheets with `scripts/app_parameters.py`
(DP4+App: copyright 2022 M. M. Zanardi and A. M. Sarotti, MIT licence, github.com/RosarioCCLab/DP4plus-App).
`lib/` holds the DP4+ code. `tables/` and `logs/` are the
output of `run_all.sh`; a run from a clean copy gives the same files byte for byte.

## Method

Cases are matched between SI files through their experimental shifts, never through their names. Four tests need
no parameters:

- **Identity.** If two SI files print the same calculation, the numbers are identical.
- **Difference test.** The difference between the shieldings of two candidates is a property of the two
  structures. At two levels it must be nearly the same vector (correlation near +1). Near -1 means the two
  candidates are exchanged at one of the levels. Weak positive values also occur for flexible molecules and are
  not counted as defects.
- **Copied columns.** The same numbers under two levels, or under two candidates.
- **Fit of a candidate.** After linear scaling a stereoisomer fits the experimental shifts within a few ppm
  (median RMS error 2.1 ppm for 13C and 0.22 ppm for 1H over all printed alternatives). A candidate with a 13C
  error above 6 ppm or a 1H error above 0.6 ppm, and more than three times the error of the correct isomer, does
  not have the atom order of the table; the case tests nothing at that level.

One test needs the DP4+ parameters: recompute the probability from the printed data and compare with the printed
probability.

| Log | Script | What it does |
|---|---|---|
| 01, 02 | `jnp_tables.py`, `jnp_probs.py` | 2023 SI: 122 shielding tables (sections 7, 9), probabilities of 84 cases (section 10) |
| 03 | `stx_vs_jnp.py` | STX150 table columns MMFF_mPW, MMFF_wB9, B3L_mPW against the 2023 values |
| 04, 05 | `stx_jnp_shieldings.py`, `stx_jnp_identify.py` | which isomers of the 2023 SI the printed STX150 candidates are |
| 06 | `stx_jnp_89.py` | case 89a at MMFF_PBE |
| 07 | `jnp_table_order.py` | 2023 SI, section 9: which printed candidate fits the experimental column |
| 08 | `si_internal_consistency.py` | difference test inside the 2023 SI |
| 09 | `obc_tables.py` | 2025 SI: all 39 compounds (five table layouts), misprinted nucleus labels and cells |
| 10 | `obc_probs.py` | 2025 probabilities at Levels A, B, C against the STX150 table; numbering; the repeated column |
| 11 | `obc_dp4plus.py` | DP4+ at Level B recomputed from the 2025 SI against its nine printed probabilities per compound |
| 12 | `obc_internal.py` | difference test inside the 2025 SI |
| 13 | `stx_obc_identify.py` | STX150 cases above 121 against the 2025 shieldings (identity with Level D; order of the candidates) |
| 14, 15 | `stx_copied_columns.py`, `stx_column_checks.py` | columns printed twice in the STX150 SI |
| 16, 17 | `stx_internal_consistency.py`, `stx_exchanged_candidates.py` | difference test inside the STX150 SI (105 candidate pairs, 15 level pairs each) |
| 18 | `stx_candidate_fit.py` | fit of every printed STX150 candidate at every level |
| 19, 20 | `jnp_tables_s7.py`, `jnp_section8.py` | 2023 SI: the 102 tables of section 7 (34 molecules, 36 levels) and the probabilities of section 8 |
| 21 | `jnp_verify_s7.py` | 2023 SI: the 36 printed probabilities of each case recomputed |
| 22 | `jnp_verify_s9.py` | 2023 SI: compounds 101 to 121, two formulations at three levels |
| 23 | `obc_verify.py` | 2025 SI: nine printed probabilities per compound at Levels A, B and C |
| 24 | `stx_verified.py` | status of every STX150 record; writes `stx_records.csv` and `stx_candidates.csv` |
| 25 | `stx_differences.py` | recomputed against published, by the defect that bears on the record |
| 26 | `stx_preview.py` | DP4+ confidence and split-conformal sets on the records |
| 27 | `stx_sp2_sensitivity.py` | does the class (sp2 or sp3) of single carbons change the unexplained entries? |

## How DP4+App computes, and what the SI files leave out

`lib/mmdp4.py` follows `dp4_module.py` of DP4+App 2.2.3: unscaled shift = reference - shielding; scaled shift from
the straight line of the unscaled shifts on the experimental ones, 13C and 1H apart; every error goes through a
Student t distribution with the mean, width and degrees of freedom of its class (sp2 or sp3, scaled or unscaled) at
that level; the products are normalised over the candidates. DP4+, MM-DP4+ and xTB-DP4+ differ in the parameter bank
only.

Two columns of the DP4+App input sheet are printed in no SI: `sp2` (the class of each nucleus) and `exchange`
(pairs of nuclei marked with a letter; for each candidate separately the larger experimental shift goes to the
smaller shielding). `lib/infer.py` recovers them where the printed probabilities need them: one change at a time
(the class of a carbon between 85 and 125 ppm with its protons, or of a proton between 4.3 and 7.5 ppm; one more
exchangeable pair among geminal protons, protons within 0.6 ppm or carbons within 3 ppm), kept if the printed
probabilities are reproduced better.

The SI files print shieldings with two or three decimals. For probabilities between the extremes that rounding
alone moves the result by 5 to 13 points. `lib/verify.py` draws the shieldings inside their rounding interval
(200 draws) and reports a printed probability as reproduced (within one point), consistent (inside the range of the
draws) or not reproduced.

## Results: the STX150 table and SI

- **Earlier results are reused.** For the 38 cases above 121 the table's MMFF_wB9 and B3L_mPW columns are Levels A
  and B of the 2025 paper in 38 of 38 (STX number = 2025 compound + 121 up to 143a, + 120 from 144a; compound 23 is
  not in STX150). The xTB_mPW shieldings printed in the STX150 SI are the Level D numbers of the 2025 SI for the
  twelve cases that have both (383 of 411 values). For cases 73 to 121 the MMFF_mPW column equals the 2023 values
  (formulation ScTMS+Uns) in 66 of 69, MMFF_wB9 in 58 of 69.
- **One table column is repeated.** MMFF_mPW equals MMFF_wB9 in all 38 rows above 121. The 2025 paper reports that
  level (its Level C) separately: 16 cases differ, 8 fall on the other side of 50% (122a, 124a, 126a, 134a, 139a,
  146a, 148a, 155a). The number of correct assignments is 30 of 38 either way.
- **Isomers exchanged between levels.** 127a, 137a, 154a: the two candidates are exchanged between the MMFF/xTB
  columns and the B3LYP-geometry columns (r from -0.80 to -1.00 for 13C and 1H). For 127a and 137a the 2025 SI
  fixes which block is the correct isomer; the published 0% at B3L_PBE and B3L_wB9 (four entries) follow from the
  exchange. 154a is undecided: the 2025 SI has the same exchange for that compound (its Levels A and C).
- **Copied columns.** Twelve: MMFF_PBE repeats B3L_PBE for 110a, 111a, 112a, 112b, 113a, 114i, 117a, 118a, 119c,
  119d, 120a and repeats xTB_PBE for 153a. Two candidate pairs share one column (90b and its alternative at
  xTB_wB9; 150a and 151a at xTB_PBE).
- **Candidates without a credible fit.** The alternative of 106a at B3L_PBE and B3L_wB9 (13C error 35 ppm), of
  115a at all six levels (1H error 0.9 ppm against 0.15 for the correct isomer) and of 142a at xTB_mPW and
  xTB_wB9 (13C 20 ppm, inherited from Level D of the 2025 SI); 117b at MMFF_PBE (55 ppm).
- **Fewer candidates printed than the source paper used.** 77a, 77b, 77d (3 of 4); 106a, 114i, 116a, 119c, 119d,
  141a, 146a, 149a (2 of 4); 121a (2 of 8). Level B DP4+ from the 2025 SI reproduces the nine printed values of
  compounds 3, 18 and 20 within 0.11 points with all four isomers; for compound 26 (146a) the correct isomer gets
  94.9% with four isomers (published 94.1%) and 98.2% with the two printed in the STX150 SI.
- **All 654 records with both values:** 487 within 1 point, 69 between 1 and 10, 51 between 10 and 50, 47 above 50.
  Of the 47, 7 are exchanges, 6 copied columns, 1 a candidate without a credible fit, 8 fewer candidates, 12 an
  unprinted input that was recovered, 13 without a cause found.
- **Flat contradictions.** Table 0-1% against recomputed 95% or more: 13 entries, 12 with a cause found, 1 without
  (111a xTB_PBE). The entry 89a at MMFF_PBE is explained: with H5 and H6 (2.05 and 1.74 ppm; H8a and H8b in the
  2023 SI) marked as exchangeable the six table values come out as 0, 95, 2, 0, 100, 100 against the printed
  0, 93, 1, 0, 100, 98; the same pair is recovered independently from the 36 values of the 2023 SI.
  Table 99-100% against recomputed 5% or less: 9 entries, 7 with a cause found.
- 29 of the 109 cases with experimental data are touched by a defect of the printed data.

## Results: records that can be used (`tables/stx_records.csv`, `tables/stx_candidates.csv`)

A record is one case at one level: experimental shifts, the candidates printed in the SI and their shieldings.

| Status | Rule | Records |
|---|---|---|
| verified | no defect touches the record, and the published probability is reproduced within 1 point (481) or lies inside the range that the rounding of the printed shieldings allows (24) | 505 |
| close | within 10 points of the published value, outside the rounding range | 42 |
| corrected | 127a and 137a, relabelled from the 2025 SI | 12 |
| unverified | more than 10 points from the published value: 13 with candidates missing from the SI, 50 without explanation (37 at 10 to 50 points, 13 above 50) | 63 |
| excluded | copied or shared column, candidate without a credible fit, 154a, no experimental data printed for 153a | 38 |

By level, verified / close: MMFF_PBE 67 / 8, xTB_PBE 82 / 9, xTB_mPW 88 / 10, xTB_wB9 87 / 7, B3L_PBE 92 / 2,
B3L_wB9 89 / 6. For 12 cases (66 records) the inputs were recovered from the six table values with at most three
changes: 77d, 82a, 86c, 87c, 89a, 103a, 106a, 117a, 130a, 143a, 145a, 158a. 47 cases are verified at all six
levels; 58 are verified, close or corrected at all six. Verified means consistent with the published table and free
of the defects that the tests can see. Where both probabilities are 100% the agreement proves little; the tests on
the printed data carry the weight there.

On the 559 verified, close and corrected records (probabilities recomputed for the printed candidates): the correct
isomer is first in 517 (92.5%); a call of 99% or more is right in 470 of 490 (95.9%); the correct isomer is at 1% or
less in 20 (3.6%). Split-conformal sets with the probability as score, calibrated on half of the candidate sets:
coverage 91.0 to 91.6% at a 90% target and 80.7 to 81.0% at an 80% target, at every level. On the 517 verified and
corrected records alone: first in 486 (94.0%), calls right in 460 of 478 (96.2%). With the 63 unverified records
added (622): first in 555 (89.2%), calls right in 481 of 505 (95.2%), coverage 90.7 to 91.3% at the 90% target.

## Results: the probabilities of the earlier papers recomputed

- **2023 SI, section 7 against section 8 (36 MMFF-geometry levels).** 34 molecules; 28 cases from pairs of isomers
  (81a, 82a, 83a, 84a, 91a, 93a have no printed alternative). 21 cases verified, 3 partly (74a, 89a, 90b),
  2 unverified (74b, 90a), 2 excluded (95a, 95b: the SI prints the same shieldings and experimental values for both
  isomers). Of 936 records, 802 are reproduced within 1 point (85.7%), 108 are consistent with the rounding of the
  printed shieldings, 26 are not reproduced (2.8%). Inputs recovered for five cases: 74a C6/C9, 74b H4/H5,
  89a H8a/H8b, 90a H4a/H7b and C3/C1, 100b the two aromatic carbons.
- **2023 SI, compounds 101 to 121 (section 9 against section 10; two formulations, three levels).** 25 cases:
  11 verified (102a, 103a, 108a, 113a, 114a, 115a, 116a, 117a, 117b, 119a, 119b), 5 partly (111a, 112a, 112b, 121a,
  121b), 9 excluded. For 101a, 118a and 120a neither printed candidate fits the experimental column. For 104a, 105a,
  106a, 107a, 109a and 110a the alternative misses the experimental column by 22 to 45 ppm: the published 100% is
  reproduced, and says nothing. Of the 96 values of the cases that are not excluded, 82 are reproduced, 6 consistent,
  8 not.
- **2025 SI, Levels A, B and C (nine printed probabilities per compound and level).** 39 compounds, 117 records:
  65 verified (Level A 19, B 21, C 25), 12 partly, 29 unverified, 11 with misprinted cells. Of 1,053 printed
  values 836 are reproduced. 17 compounds are verified at all three levels, 4 at two. The correct isomer is not the
  first printed block for compound 3 (Level C), 5 (A), 9 (A and B), 12 (C) and 31 (A). Inputs recovered for 12
  compounds: the class of one carbon set to sp3 for 7 (C4), 8 (C6), 25 (C10), 26 (C13), 35 and 36 (C5); exchangeable
  pairs for 4, 5, 32, 33, 38 and 39.

## Results: the 2025 SI by itself

- All 39 compounds are read: 29 standard tables, 5 continuations (four isomers: compounds 3, 18, 20, 26, 29), one
  side-by-side table (4 and 5) and four pair tables (30/31, 32/33, 35/36, 37/38).
- Misprints: the nucleus label of 15 rows in 10 compounds (protons printed as C or the reverse); 123 cells in 32
  of the 292 shielding columns that cannot be the shielding of their nucleus (compounds 3, 21, 22, 23, 30, 31,
  32, 33; in several the Level C column or the second isomer is shifted against the rows); the experimental
  columns of compounds 32 and 33 repeat those of 30 and 31.
- Difference test: isomers exchanged between levels for compounds 9 (Levels A, B against C, D) and 34.
- DP4+ at Level B with the inputs taken by rule (step 11): all nine printed probabilities within one point for 19
  compounds, within 10 for 6 more, not for 14. With recovered inputs (step 23): 21 verified at Level B. For
  compound 9 the first printed block gets 0.01% and the second 99.99% (published 99.97%), which says the two blocks
  are exchanged at that level.

## Limits

- The SI files print neither the class of each nucleus nor the exchangeable pairs. By rule, carbons at or above
  100 ppm are sp2 and protons follow the carbon with the same label; no pair is exchangeable. Recovered inputs are
  the smallest change that reproduces the printed probabilities, not the authors' input sheets.
- The xTB sheets of the DP4+App parameter bank hold copies of QM and MM sets, not the xTB sets printed in the
  STX150 paper; the STX150 levels are recomputed with the sets printed in that paper.
- The 2023 SI gives shieldings for 52 of the 72 cases numbered up to 121. In its section 9 the second candidate of
  several tables (104, 105, 106, 107, 109, 110) has errors of 22 to 45 ppm against the experimental column, so
  those tables cannot tell which STX150 candidate is which; the tables of 101, 118 and 120 do not fit either
  candidate. Its table for isomers 121e-h has shifted row labels.
- MMFF_mPW, MMFF_wB9, B3L_mPW and xTB_PBE_STX150 have no shieldings in the STX150 SI, so those four columns of the
  table can be compared with the earlier papers but not recomputed.
- The relabelling of 127a and 137a rests on the statement of the 2025 SI that its first isomer is the correct one.

## Published DP4-AI probabilities

`dp4ai_published/` holds the output of three scripts that are not part of `run_all.sh`:

    python scripts/dp4ai_published.py d0sc00442a1_suppl.xlsx <out>            # needs openpyxl
    python scripts/dp4ai_published_truth.py <out> <res/core of Kaggle run 1> <scripts dir of the DP4-AI archive>
    python scripts/dp4ai_published_audit.py <out> <tables dir of the DP4-AI archive>

The spreadsheet lists every candidate of 47 molecules at three levels and by two routes (automatic; hand-made peak
lists paired in order of shift). Calls of 99% or more that are right: automatic 28 of 50 (56.0%), hand-made 41 of
49 (83.7%). Correct isomer at 1% or less: 46 of 141 and 26 of 141. Against the recomputation of the DP4-AI archive
the hand-made route has the same top isomer in 41 of 44 (B3LYP energies), 37 of 44 (M06-2X energies) and 34 of 43
(MMFF geometries). In the sheet for B3LYP geometries the column marking the correct isomer sits one row low.
