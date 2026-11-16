# ConformalStereo: coupling constants of the 2017 Cambridge deposit

Third Kaggle run (notebook 1.2, 14 November 2026) and what its data show.

The deposit "Doubling the power of DP4" (Ermanis, Parkes, Agback and Goodman, 2017; doi 10.17863/CAM.13222, CC BY 4.0)
holds 3,797 Gaussian jobs with proton-proton coupling constants, at B3LYP/6-31G(d,p) on MMFF geometries, in 25 folders
whose names end in `_J`. The run read all of them: the matrix of total couplings and its Fermi-contact part for every
conformer. This package checks the matrices, reproduces what the depositors' PyDP4 outputs print about couplings, and
asks whether couplings make conformal prediction sets smaller.

## Answers in short

1. **The data are sound.** 3,797 outputs, 11,391 matrices, none with an unexpected layout; the text kept for 150
   matrices, read with the reader of PyDP4, gives the stored numbers.
2. **The depositors' outputs are reproduced from the raw data**: every one of the 212 candidate blocks of the 35
   outputs with experimental couplings, for the shifts (own shieldings and energies of the coupling jobs), for the
   assignment of the protons, and for the calculated couplings (1,462 of 1,468 protons as printed; the other 6 differ
   only in the order of two equally good matchings). The outputs used the total coupling, not its Fermi-contact part.
3. **Three properties of the depositors' coupling term make it depend on the numbering or order of the atoms.**
   - In 9 outputs the candidates are numbered differently from each other. PyDP4 renumbered the shieldings but used
     one list of hydrogen numbers, that of the last candidate, for the couplings of all. In at least 19 candidate
     blocks the printed couplings therefore belong to other protons.
   - Calculated couplings were matched to experimental ones by the smallest summed absolute error. That sum is the
     same for every matching whenever all calculated values lie on one side of the experimental ones; 118 of 1,468
     protons have several equally good matchings, and the choice among them (the first in the order of the atoms)
     changes a Gaussian coupling term by more than a factor of ten for 58 of them.
   - Signals with equal shifts received the couplings of the first such signal (3 signals in 2 descriptions).
4. **Two statistical models are behind the printed coupling terms.** Six outputs (Rosuvastatin twice, Swerilactone A
   and B, Testosterone twice) follow the Gaussian of the PyDP4 source exactly. The other 29 follow one common model
   that is not in the source history: corrected coupling = (calculated + 0.651) / 1.055, and a density of the
   absolute error with a peak at zero and a shoulder out to 4 Hz. It was recovered here from the printed terms
   themselves; with it all 182 usable blocks are matched to the rounding of the printed numbers, and 161 of 181 when
   the output in question is left out of the fit.
5. **Couplings do not buy smaller sets with any confidence.** On the 33 compounds (22 families), with the couplings on
   the right atoms and a matching that does not depend on atom order:

| Data and model | Actual first | Calls at 99% or more: right | Fitted temperature | Sets at 90%: size (random halves) | Size (one family out) |
| --- | --- | --- | --- | --- | --- |
| Shifts (C and H) | 21 of 33 | 16 of 17 | 6.2 | 3.10 | 2.72 |
| Shifts + J, Gaussian of the source | 22 | 21 of 29 | 23.4 | 6.67 | 7.00 |
| Shifts + J, Fermi-contact part | 22 | 18 of 24 | 18.2 | 6.08 | 6.80 |
| Shifts + J, recovered density | 22 | 17 of 20 | 6.6 | 3.26 | 2.34 |
| Shifts, energies averaged (4 kJ/mol) | 23 | 14 of 15 | 4.0 | 2.49 | 2.49 |
| Shifts + J recovered, energies averaged | 22 | 16 of 16 | 4.4 | 2.33 | 2.27 |

   The Gaussian of the source (standard deviation 0.79 Hz) is far too narrow: of the 427 couplings of the actual
   structures, 23% miss by more than two standard deviations and 12% by more than 3 Hz. With it the couplings
   overrule the shifts, eight calls at 99% or more are wrong, and the conformal sets double in size. With the
   heavy-tailed recovered density the couplings neither help nor harm beyond noise: the two schemes disagree in sign
   at 90% (+0.16 and -0.38 structures) and at 80% (-0.01 and +0.22). The gain that holds in both schemes is the one
   found before, averaging over conformer energies: 2.49 against 3.10 or 2.72 at 90%, and a little more with couplings
   added (2.33, 2.27).
6. **The depositors' own numbers look better than that**: actual structure first for 21 compounds with shifts and 25
   with couplings. Recomputed with each output's own model, two of those four gains (Callipeltin fragment SRR,
   Tetrahydrofuran 3) are not there.

## Contents

| Path | Content |
| --- | --- |
| `run_all.sh` | the whole chain, seven steps |
| `scripts/jcheck.py` | checks of the matrices of run 3 against run 2 and against the reader of PyDP4 |
| `scripts/jsource.py` | which structure, shieldings and energies stand behind every candidate block |
| `scripts/jrepro.py` | couplings and coupling term of the deposit outputs, reproduced step by step |
| `scripts/jmodel.py` | scaling and density behind the coupling term of 29 outputs, with a leave-one-output-out test |
| `scripts/jpreview.py` | the depositors' printed numbers, with and without couplings (the earlier "first look") |
| `scripts/jbench.py` | shift and coupling likelihoods of all 208 candidates from the raw data, three coupling models, energy averaging |
| `scripts/jeval.py` | ranking, confidence, conformal sets, printed against recomputed |
| `lib/` | readers of the deposit (from the Cambridge benchmark package), `jlib.py` (matrices; PyDP4's coupling steps), `jref.py` (PyDP4's matrix reader, used as an independent check) |
| `inputs/` | experiments, actual structures and level-1 records of the 2017 benchmark |
| `tables/`, `logs/` | results of the run that made this package |
| `LICENSE_PyDP4.txt` | licence of PyDP4, whose routines `jlib.py` and `jref.py` follow |

Main tables: `jbench_candidates.csv` (one row per compound and candidate: log-likelihoods for C, H and the three
coupling models, plain and energy-averaged), `jbench_couplings.csv` (every matched coupling: proton, partner,
experimental, calculated, corrected), `jeval_overview.csv`, `jeval_conformal.csv`, `jeval_compounds.csv`,
`jeval_printed_vs_recomputed.csv`, `jrepro_blocks.csv`, `jmodel_recovered.json`, `jmodel_density.csv`.

## How to run

    sh run_all.sh <res/core of Kaggle run 2> <res/core of Kaggle run 3> <out folder> [python]

Run 2 is the unpacked `cs11_core_NN.zip` set, run 3 the unpacked `cs12_core_01.zip`. Python 3.9 or later with numpy
and scipy; about four minutes. Random numbers are seeded; a rerun gives the same tables.

## What the third run found in the files

- With `NMR=ReadAtoms` and `atoms=H`, Gaussian 09 prints each matrix for all atoms and leaves zeros in the rows and
  columns of the others. In the first job step the diagonal holds one constant (385 GHz) on the hydrogen atoms.
- Off the diagonal the Fermi-contact matrix of the second step repeats that of the first.
- 1,444,853 proton pairs in 3,797 conformers; 11.6% are above 1 Hz in size. For those the total coupling differs from
  its Fermi-contact part by 0.38 Hz on average (total = 0.983 x FC - 0.075 Hz).
- In the outputs of Swerilactone A and B the third candidate block is a copy of the fourth: seven distinct candidates
  were compared there, not eight. The recomputation uses all eight.

## The coupling term used here (`jbench.py`)

Kept from PyDP4: couplings below 0.2 Hz are not partners; equivalent protons are averaged and their mutual couplings
set to zero; calculated couplings are corrected linearly and the absolute corrected error is scored.

Changed, so that the term does not depend on the order or numbering of the atoms:

- candidates are brought to the numbering of the first candidate for the couplings as well as for the shieldings;
- every experimental signal keeps its own couplings;
- the matching of calculated to experimental couplings is the one with the largest likelihood under the scoring model;
- an experimental coupling without a calculated partner is compared with zero (PyDP4 skipped a proton without
  calculated couplings, and set all calculated values of a proton to zero when they were fewer than the experimental).

Models: `gauss`, slope 1.160, intercept -0.366 Hz, Gaussian with mean -0.131 and standard deviation 0.793 Hz (source);
`fc`, the Fermi-contact matrix with slope 1.170, intercept -0.438, mean -0.154, standard deviation 0.921 Hz (source);
`rec`, slope 1.055, intercept -0.651 Hz and the recovered density (beyond 10 Hz it is continued at the rate of decay
found between 8 and 10 Hz). Energy averaging: independent Gaussian noise of 2, 4 or 8 kJ/mol on the conformer
energies, 500 draws, the same draws for shifts and couplings, pairings fixed at the unperturbed calculation.

Conformal sets: score 1 - probability, families (molecule archives) as exchangeable units with weights as in the
pooled calibration; random halves of the families (500 splits) and each family against all the others.

## Limits

- 33 compounds in 22 families at one level of theory. Differences of a few tenths of a structure in set size are
  within noise; the two calibration schemes disagree in sign for the couplings.
- The recovered density was built by the depositors, very likely from these same compounds; its good behaviour here is
  not an out-of-sample result.
- The actual structures are those of the 2017 benchmark, identified from the statistics printed in the SI of the paper.
- Where the recomputed coupling term differs from the printed one, the cause is one of the four changes listed above;
  the printed terms themselves are reproduced with PyDP4's own rules in `jrepro.py`.
- Nine of the 42 compounds have no output with experimental couplings and are not covered.
