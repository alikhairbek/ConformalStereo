# ConformalStereo: the energy-averaged likelihood on independent molecules (DP4-AI deposit)

14 November 2026. Two additions to the DP4-AI benchmark (`ConformalStereo_DP4AI_benchmark.zip`; data of Howarth,
Ermanis and Goodman, Chem. Sci. 2020, 11, 4351, doi 10.17863/CAM.47721):

1. a pairing rule for experimental descriptions that list fewer or more values than there are atoms (7 of the 46
   molecules);
2. the test of the energy-averaged likelihood on these 46 molecules. None of them is in the two Cambridge benchmarks
   on which that score was developed, and the noise level, 4 kJ/mol, was fixed before this test.

## Method

- **Records.** 46 molecules at 4 protocols (MMFF geometries; B3LYP geometries with B3LYP, mPW1PW91 or M06-2X
  energies): 182 records, 2 to 32 candidates each. The reported structure is isomer 1 in every case. 141 records are
  clean: the geometry check of the DP4-AI archive is passed and the description has one value per atom.
- **Pairing (`scripts/pairing.py`).** The shift lists are unassigned. With one value per atom the sorted lists are
  paired in order, as in the 2020 code. Otherwise the 2020 code pairs from the top and leaves out the atoms (or the
  signals) with the lowest shifts. Here the order is kept and the items left out are those whose omission gives the
  highest DP4 likelihood (dynamic programme, repeated with the internal scaling, three starting points).
- **Energy averaging (`scripts/dp4ai_energy.py`).** Each conformer energy gets independent Gaussian noise; in each of
  1,000 draws the shieldings are averaged with the new Boltzmann weights, sorted, paired and scored as in DP4; the
  likelihood of the candidate is the mean over the draws.

## Results

**Pairing.** For descriptions with one value per atom nothing changes (log-likelihoods equal to 2.5e-9). For the 28
records of the 7 other molecules the reported structure is first in 13 under the 2020 rule and in 15 under the new
rule. The changes go both ways: JB7 moves from rank 7 to 14 to rank 1 at three protocols, TS3B from rank 1 to rank 2
at two. The dynamic programme equals exhaustive search in 300 of 300 random problems (`logs/00`).

**Energy averaging**, all 182 records (`logs/02`, `logs/03`):

| Noise on the energies | Reported structure first | Calls of 99% or more: right / all | Reported structure at 1% or less | Temperature | Set size at 90% / 80% |
|---|---|---|---|---|---|
| none | 96 (52.7%) | 65 / 81 (80.2%) | 35 | 6.8 | 3.56 / 2.13 |
| 2 kJ/mol | 103 (56.6%) | 63 / 77 (81.8%) | 33 | 5.8 | 3.50 / 2.08 |
| 4 kJ/mol | 109 (59.9%) | 62 / 73 (84.9%) | 30 | 5.0 | 3.41 / 1.96 |
| 8 kJ/mol | 110 (60.4%) | 67 / 76 (88.2%) | 30 | 4.7 | 3.45 / 1.96 |
| 16 kJ/mol | 106 (58.2%) | 65 / 72 (90.3%) | 27 | 4.6 | 3.53 / 2.10 |

Temperature: the factor by which the log-likelihoods would have to be divided for honest probabilities. Set size:
split conformal prediction with the probability as score, calibration and test separated by candidate-set family,
one threshold per protocol; coverage is 91.4 to 91.6% at the 90% target throughout.

On the 141 clean records the reported structure is first in 75 without noise and in 86 at 4 kJ/mol; calls of 99% or
more are right in 51 of 63 and in 48 of 57; the sets shrink from 4.57 to 4.17 structures at the 90% target and from
2.23 to 2.12 at 80%.

Record by record at 4 kJ/mol: 20 records are right only with the noise, 7 only without it. By molecule: 13 have
more protocols right, 6 fewer (sign test, two-sided p = 0.17).

## Limits

- The direction is the one seen on the two Cambridge benchmarks, but with 46 molecules the gain is not significant
  at the level of molecules. The four protocols of a molecule are not independent.
- The pairing rule is a rule, not a fact: for incomplete lists nobody knows which signals are missing. Methyl protons
  are not averaged and exchangeable protons are not removed, as in the 2020 code. The missing values do not match the
  number of exchangeable protons in any of the five molecules concerned.
- The numbers move a little with the random draws: with 200 draws the count at 4 kJ/mol was 108 to 112 over three
  seeds; 1,000 draws are used here.

## Files

- `tables/dp4ai_records.csv`: one line per record: probability and rank of the reported structure under the 2020
  rule, under the new pairing, and at each noise level.
- `tables/dp4ai_energy_marginal.csv` (all records), `dp4ai_energy_marginal_clean.csv`: log-likelihood of every
  candidate under the 2020 rule and at each noise level. `tables/dp4ai_truth.csv`: position of the reported structure.
- `inputs/`: two tables of the DP4-AI archive (reported structures, geometry check).
- `logs/`, `scripts/`. To rerun: `sh run_all.sh <res/core of Kaggle run 1> <output dir>`; Python 3.8 or later with
  numpy and scipy; under a minute.
