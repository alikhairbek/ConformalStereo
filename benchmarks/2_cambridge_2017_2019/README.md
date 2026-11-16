# ConformalStereo: two benchmarks from the Cambridge deposits, and a first test of conformal scores

14 November 2026. This archive replaces `ConformalStereo_Cambridge2017_benchmark.zip`: the tables of the 2017
deposit are the same byte for byte, the scripts now serve both deposits, and the 2019 deposit and the score tests
are new.

| | 2017 deposit | 2019 deposit |
|---|---|---|
| doi | 10.17863/CAM.13222 | 10.17863/CAM.32463 |
| Paper | Org. Biomol. Chem. 2017, 15, 8998 | the group's follow-up on geometries and energies |
| Experiments with a known actual structure | 42 | 25 (all among the 42) |
| Levels | 12 NMR levels on MMFF geometries | 3 NMR levels on 6 kinds of geometry, 2 choices of energies: 30 workflows |
| Records | 504 | 703 |
| Candidates per experiment | 2 to 32 | 2 to 32 |
| Tables | `tables17/` | `tables19/` |

## Method (the same for both)

1. Every PyDP4 output of the deposit is read. For every structure the shieldings are Boltzmann-averaged over its
   conformers (298.15 K, SCF energies of the NMR job).
2. Each isomer of each output is matched to a structure: `shift = (reference - shielding) / (1 - reference * 1e-6)`,
   so `(shift + shielding) / (1 + shift * 1e-6)` must be one constant per element. This also shows which energies
   stand behind an output.
3. DP4 is recomputed with one rule: atoms and pairing as in PyDP4 (2015-2017), single Gaussian error model
   (2.269 ppm for 13C, 0.187 ppm for 1H) on internally scaled shifts.
4. Actual structures. 2017: the candidate that shows the six model-free numbers printed per compound in section S-4
   of the paper's SI. 2019: taken over from 2017 after checking that the candidates are the same.

## 2017 deposit: results

- Outputs of the deposit: 3,346 of 3,480 isomer blocks reproduced with the energies of their own level, 113 only
  with another level's energies, 21 not. 455 outputs are comparable in full; the recomputed probabilities are within
  0.1 point of every printed value in all 455.
- All four tables of S-4 are reproduced: the six numbers are within the printed precision for 37, 39, 36 and 39 of
  the 42 compounds, and within twice that for 41, 42, 42 and 40.
- Actual structure first in 326 of 504 records (64.7%); 24 to 31 of 42 by level. Calls of 99% or more: 258, of which
  215 right (83.3%). Sets built from the probabilities to sum to 90% hold the actual structure in 74.8%.
- Split conformal sets (by molecule archive, 2,000 splits): 91.5% at a 90% target with 3.73 structures on average.

## 2019 deposit: results

- Levels: MMFF, PM6, PM7, HF/6-31G(d,p), B3LYP/6-31G(d,p) and M06-2X/6-31G(d,p) geometries, each with NMR at
  B3LYP/6-31G(d,p), mPW1PW91/6-311G(d) and M06-2X/6-31G(d,p). 22 experiments have all 18; artemisinin and the two
  aspergillides have the HF and DFT geometries only.
- Outputs of the deposit: of 2,412 isomer blocks, 1,802 reproduced with the energies of their own level, 344 with
  another level's (mostly mPW1PW91 shifts weighted with M06-2X energies). Five output files in mPW1PW91 folders
  hold the B3LYP result (64 blocks). 202 blocks on HF and DFT geometries are not reproduced, 108 of them by less
  than 0.5 ppm for carbon; the cause was not found. 207 outputs are comparable in full; the recomputed probabilities
  are within 0.1 point of every printed value in all 207.
- Actual structures: same structure names in the same order for all 25 experiments; averaged shieldings at
  B3LYP/6-31G(d,p) on MMFF geometries identical in both deposits for 20 of the 22 experiments that have that level
  (ezetimibe and rosuvastatin use other conformer sets).
- By geometry, with the energies of the NMR level (actual structure first; calls of 99% or more that are right):
  MMFF 35 of 66 (14 of 20), PM6 28 of 66 (13 of 24), PM7 37 of 66 (24 of 29), HF 49 of 75 (28 of 31),
  B3LYP 51 of 74 (24 of 28), M06-2X 45 of 75 (20 of 26).
- The factor by which the log-likelihoods would have to be divided for honest probabilities: HF 3.9, B3LYP 4.7,
  M06-2X 6.8, MMFF 9.6, PM7 13.9, PM6 25.8 (2017 benchmark: 8.4).
- Raw files: 391 of 1,312 PM6 optimisations did not converge (60 of 82 structures); the HF and DFT optimisations
  converged in all but two jobs. One record is left out (artemisinin, mPW1PW91 on B3LYP geometries: outputs missing).

## Conformal scores: first results (`logs/23` to `29`)

Every score orders the candidates; the set keeps those below a threshold. Compared at matched coverage and under
split conformal prediction by molecule archive.

- Rescaling the probability (temperature, log-likelihood ratio, per nucleus): at a 90% target on the 2017 benchmark
  the sets hold 2.72 structures with the probability and 2.55 with the log ratio per square root of the number of
  nuclei (3.75 and 3.39 under split conformal). Tempering alone changes nothing.
- Heavier-tailed error models (Student t, Laplace, capped Gaussian): actual structure first in 62 to 65% of the
  records for every model, set sizes 2.6 to 2.9. No gain.
- A consensus over NMR levels: no gain (29 or 30 of 42 first, against 31 for the best single level).
- Averaging the likelihood over errors in the conformer energies (independent Gaussian noise on each conformer
  energy, 200 draws): a gain of about a tenth in set size at the same coverage on both benchmarks, levelling off
  between 4 and 8 kJ/mol.

## Files

- `tables17/`, `tables19/`: `*_benchmark.csv` (one line per record), `*_candidates.csv`, `*_shifts.csv` (every
  assigned atom), `*_loglik.csv` (log-likelihood of every candidate, 13C and 1H), `*_truth.csv`, `*_records.csv`,
  `*_qc.csv`, `*_reference_constants.csv`, `*_experiments.json`, `*_scores.csv`, `*_energy_marginal.csv`;
  `dep17_si_check.csv`, `dep17_error_models.csv`, `dep19_workflows.csv`.
- `inputs/obc2017_si_s4_*.csv`: the four tables of section S-4 of the 2017 SI, typed in from the PDF.
- `logs/`: the printed summary of each of the 29 steps. `scripts/`: the code.

## To rerun

`sh run_all.sh <res/core of Kaggle run 2> <res/core of Kaggle run 1> <output dir>`. Python 3.8 or later; steps 25 to
29 need numpy and scipy. About two minutes.

## Limits

- The actual structures are those the 2017 paper treats as such; they were not checked against the original
  structure papers.
- One error model (2015 parameters) at all levels. The models the papers fit per level are not recomputed.
- The two benchmarks share compounds, so nothing found on one is confirmed independently by the other.
- With 22 to 25 experiments per workflow the differences between geometries in the 2019 benchmark are descriptive.
