# ConformalStereo: one calibration for everything, or one per kind of calculation?

14 November 2026. All five benchmarks in one table, and conformal calibration across them.

| Benchmark | What | Method | Records | Families |
|---|---|---|---|---|
| JOC2015 | 48 cases x 24 NMR levels, recomputed from the printed shieldings | DP4+ on B3LYP geometries | 1,152 | 18 |
| STX150 | the usable records of the 110 validation cases at six levels | DP4+, MM-DP4+, xTB-DP4+ | 559 | 80 |
| DP4-AI | 46 molecules x 4 protocols, unassigned shifts | DP4 (single Gaussian) | 182 | 40 |
| Cam2017 | 42 compounds x 12 NMR levels on MMFF geometries | DP4 | 504 | 30 |
| Cam2019 | 25 of the same experiments x 30 workflows on six kinds of geometry | DP4 | 703 | 17 |

3,100 records, 14,437 candidates, 152 families. A family is a set of molecules that share candidates, or that occur
in two benchmarks (JOC2015 and STX150 share 16 families; Cam2017 and Cam2019 share all 17 of the latter). Families,
not records, are the exchangeable units: a calibration score is weighted by one over the number of records of its
family in the calibration pool, and coverage and set size are averaged per test family (hierarchical conformal
prediction). The score is 1 minus the probability of the candidate.

## Results (`logs/04_pooled_calibration.txt`)

**Every method is overconfident on every benchmark.** The temperature is the factor by which the log-likelihoods
would have to be divided for honest probabilities.

| Benchmark | Group | Actual structure first | Calls of 99% or more: right / all | Temperature |
|---|---|---|---|---|
| JOC2015 | DP4+, B3LYP geometries | 857 of 1,152 (74.4%) | 613 / 698 (87.8%) | 5.5 |
| STX150 | DP4+, B3LYP geometries | 187 of 193 (96.9%) | 171 / 174 (98.3%) | 2.7 |
| STX150 | MM-DP4+, MMFF geometries | 66 of 77 (85.7%) | 61 / 65 (93.8%) | 5.0 |
| STX150 | xTB-DP4+, xTB geometries | 264 of 289 (91.3%) | 238 / 251 (94.8%) | 4.2 |
| DP4-AI | DP4, B3LYP geometries | 74 of 138 (53.6%) | 50 / 62 (80.6%) | 7.1 |
| DP4-AI | DP4, MMFF geometries | 22 of 44 (50.0%) | 15 / 19 (78.9%) | 6.0 |
| Cam2017 | DP4, MMFF geometries | 326 of 504 (64.7%) | 215 / 258 (83.3%) | 8.4 |
| Cam2019 | DP4, HF or DFT geometries | 246 of 373 (66.0%) | 122 / 147 (83.0%) | 5.1 |
| Cam2019 | DP4, MMFF geometries | 62 of 110 (56.4%) | 21 / 31 (67.7%) | 9.6 |
| Cam2019 | DP4, semi-empirical geometries | 111 of 220 (50.5%) | 60 / 86 (69.8%) | 19.6 |

A stated probability of 99.9% or more for the top candidate is right in 92.6% of the DP4+ records on B3LYP
geometries (651 of 703), 94.6% for MM-DP4+, 94.8% for xTB-DP4+, 87.4% for DP4 on MMFF geometries, 84.1% on HF or
DFT geometries and 65.3% on semi-empirical ones.

**Coverage at a 90% target, by calibration design** (300 random halves of the families; mean set size in brackets):

| Benchmark | Own level | By benchmark | By group | One for all |
|---|---|---|---|---|
| JOC2015 | 98.3% (2.48) | 98.9% (2.65) | 83.0% (1.04) | 93.9% (1.52) |
| STX150 | 90.5% (1.01) | 91.8% (1.02) | 93.3% (1.07) | 96.7% (1.25) |
| DP4-AI | 92.7% (3.64) | 93.0% (3.35) | 93.1% (3.28) | 86.3% (2.25) |
| Cam2017 | 92.2% (5.14) | 93.1% (5.24) | 86.5% (3.02) | 80.1% (2.05) |
| Cam2019 | 97.4% (5.88) | 98.6% (6.83) | 93.2% (3.68) | 83.4% (2.10) |
| All | 92.2% (2.73) | 93.1% (2.75) | 90.9% (2.09) | 90.5% (1.70) |

- *One for all* reaches 90.5% on average and misses where it matters: 80 to 86% for the three DP4 benchmarks, 94 to
  97% for the two DP4+ ones. At an 80% target the DP4 benchmarks get 66 to 68%.
- *By group* (method and kind of geometry) is not enough either: JOC2015 gets 83.0%. Its 24 NMR levels include weak
  ones, the STX150 levels of the same group are strong, and one threshold cannot serve both.
- *Own level* and *by benchmark* hold everywhere. With few families they are conservative (JOC2015 has 18, Cam2019
  17): a threshold from nine families can only be their largest score.

**Transfer** (threshold from the other benchmarks, families of the test benchmark left out; 90% target):

| Test on | From every other record | From the same groups only |
|---|---|---|
| JOC2015 | 94.4% | 72.7% |
| STX150 | 98.7% | 100.0% (DP4+ records only, sets of 2.1) |
| DP4-AI | 85.0% | 96.9% |
| Cam2017 | 76.2% | 85.0% |
| Cam2019 | 81.1% | 89.6% |

A calibration made on one benchmark does not carry over to another, in either direction, even within one method.

**Energy-averaged score under the same family-weighted calibration** (`logs/05`, `logs/06`; the 1,108 DP4 records
that have both scores): at an 80% target the sets hold 2.14 structures instead of 2.37 with a threshold per level,
and 1.82 instead of 2.04 with a threshold per group; at 90%, 4.37 instead of 4.51 and 3.19 instead of 3.32.
Coverage is unchanged.

## What follows for the method

The calibration set has to be made with the workflow that is used on the new molecule: same method, same kind of
geometry, same NMR level. A calibration shipped with the paper is valid for its workflow only, and the paper has to
say so. Probabilities stated by DP4 or DP4+ should not be read as frequencies at any level examined here.

## Limits

- The benchmarks differ in their molecules as well as in their methods, so a failed transfer cannot be assigned to
  the method alone.
- A set can be empty (no candidate passes the threshold); mean sizes below one reflect that. In use the top candidate
  would always be kept, which cannot lower the coverage.
- Earlier previews in this project averaged over records without weighting by family; the numbers here weight by
  family and are not comparable with them one to one.
- JOC2015: recomputed from printed shieldings; every candidate is within 0.5 point of the published probability in
  1,040 of 1,152 records. STX150: the 559 records with status verified, close or corrected.

## Files

- `tables/pooled_candidates.csv`: every candidate: benchmark, method, geometry, level, case, family, whether it is
  the actual structure, number of nuclei, log-likelihood relative to the best candidate (and its energy-averaged
  version where there is one).
- `tables/pooled_overview*.csv`, `pooled_designs*.csv`: the tables behind the logs.
- `tables/joc_loglik.csv`, `stx_loglik.csv`: the exports of the two DP4+ benchmarks.
- `inputs/`, `lib/`: the 2015 data and the DP4+ routines of the first archive.

## To rerun

    sh run_all.sh <ConformalStereo_SI_crosscheck (v3)> <ConformalStereo_DP4AI_energy_test/tables> <ConformalStereo_Cambridge_benchmarks> <out>

Python 3.8 or later with numpy and scipy; about two minutes.
