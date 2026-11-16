# ConformalStereo: the DP4-AI deposit as a third benchmark

DP4 recomputed from the raw Gaussian outputs of the DP4-AI deposit (Howarth, Ermanis and Goodman,
Chem. Sci. 2020, 11, 4351; data: https://doi.org/10.17863/CAM.47721, CC BY 4.0), with a quality check of the
deposit, an audit of the stated confidence and a first conformal test.

Written 14 November 2026. Everything here is rebuilt by one command from the result files of the Kaggle
notebook (run 1, the `cs_core_NN.zip` files).

## How to rebuild

```
sh scripts/run_all.sh <folder with the unpacked cs_core zips> <work folder> <output folder> [python with RDKit]
```

Python 3.9 or later. RDKit is needed for two steps (`truth.py`, `swap_test.py`); the other steps use the standard
library only. The run takes about two minutes. A fresh run reproduces the tables in `tables/` byte for byte.

## What the deposit holds

| | |
| --- | --- |
| Molecules (archives) | 46 (the paper describes a test set of 47) |
| Constitutions | 40; six pairs of molecules share one candidate set (JB2/JB3, JB13A/B, NL1A/B, NL2A/B, NP3A/B, TS3A/B) |
| Candidate diastereomers | 344, from 2 to 32 per molecule |
| Conformer NMR jobs | 8,562 on MMFF geometries, 8,655 and 8,670 in the two DFT-geometry folders |
| Gaussian job steps read | 53,114 |
| Reported structure | isomer 1 in all 46 cases (45 from the input structure file, NP5 from the number of 13C signals) |
| Experimental data | one description file per molecule (`<case>NMR`), unassigned shift lists; 4 molecules have 1H only |

## Protocols

NMR is mPW1PW91 throughout (6-311G(d), def2 sets for a few molecules), with PCM.

| Name | Geometries | Conformer energies | Folder |
| --- | --- | --- | --- |
| `MM` | MMFF | energy of the NMR job | `MM` |
| `DFT_EB3LYP` | B3LYP/6-31G(d,p) | B3LYP energy of the optimisation | `Opt` |
| `DFT_EmPW` | same | energy of the NMR job | `Opt` |
| `DFT_EM062X` | same | M06-2X/def2-TZVP single point | `Opt_E` |

Statistics and pairing follow the January 2020 code of PyDP4 (single Gaussian, sigma 2.269 ppm for 13C and
0.187 ppm for 1H, internal scaling; the i-th largest experimental shift goes to the atom with the i-th largest
calculated shift). The published evaluation used a three-Gaussian model fitted to this test set; its parameters
are not in the code repository, so the numbers below are those of the default model, the one the stored outputs use.

## Checks

1. **The original outputs are reproduced.** The deposit keeps 131 DP4 output files. For 109 of them the calculated
   shifts are reproduced to 0.01 ppm (25,492 values) once the right shieldings and conformer energies are found,
   and for 107 of those the pairing and the probabilities as well (largest difference 0.05 percentage points,
   24,856 atom-peak pairs). The other 22 belong to earlier runs (other isomer sets or calculations no longer in
   the deposit).
2. **Which energies the stored outputs used.** Of the 28 reproduced outputs computed on DFT geometries, 26 are
   reproduced only with the energy of the NMR job, one with the B3LYP energy and one with the M06-2X energy.
   Another 18 outputs stored in the `Opt` and `Opt_E` folders hold MM-level results. The stored outputs are
   therefore not a record of the three published levels; the levels have to be recomputed, as done here.
3. **Geometry check.** For 214 conformers in `Opt` (215 in `Opt_E`) the NMR job did not use the final optimised
   geometry: NP3A 66 of 83, NP3B 83 of 83, NP4 41 of 47, NP5 3 of 3, OD1 21 of 82. For 193 of them the NMR job
   ran on the starting geometry, finished before the optimisation did, and has the shieldings of the `MM` folder.
   In `Opt_E`, 297 single points are on another geometry than the NMR job (the same 215 conformers, 81 of 675
   for TP2 and one for KE2). A case is called *clean* at a protocol when no isomer has more than 1% of its
   population on such conformers.
4. **Exchanged labels.** For the six pairs that share a candidate set, no analysis prefers the partner's
   structure in both directions (`logs/12_swap_test.txt`).

## First results (description file `<case>NMR`, combined 13C and 1H where both exist)

"Clean" = geometry check passed and the description lists exactly one value per carbon and per proton
(39 of 46 molecules have such a description).

| Protocol | Cases (all / clean) | Reported structure ranked first | Calls at 99% or more: right / all | Reported structure at 1% or less |
| --- | --- | --- | --- | --- |
| `MM` | 44 / 37 | 21 / 19 | 15 of 19 / 13 of 16 | 11 / 8 |
| `DFT_EB3LYP` | 46 / 35 | 23 / 18 | 17 of 21 / 12 of 14 | 10 / 5 |
| `DFT_EmPW` | 46 / 35 | 25 / 20 | 17 of 20 / 12 of 14 | 10 / 6 |
| `DFT_EM062X` | 46 / 34 | 25 / 18 | 19 of 25 / 14 of 19 | 9 / 7 |

Pooled over the four protocols on clean cases (141 case-protocol pairs, not independent): calls at 99% or more
are right in 51 of 63 (81.0%); the reported structure gets 1% or less in 26 pairs (18.4%); nominal 90% sets
cover it in 104 pairs (73.8%).

The authors' own stored outputs show the same picture without any recomputation: among 66 distinct reproduced
outputs with real experimental data, isomer 1 is ranked first in 39 (59%), and calls at 99% or more are right
in 26 of 34.

**Sensitivity.** Each molecule has up to eight versions of the same analysis (four protocols, two description
files, `<case>NMR` and the older `<case>NMRold`): 337 analyses of 46 molecules. The reported structure is ranked
first in every version for 13 molecules, in none for 8, in some for 25. In 7 molecules two versions each give 99%
or more to a different isomer.

**Conformal sets** (score = 1 - P(reported structure), 2,000 random half/half splits by candidate-set family).
At a 90% target the coverage holds (91 to 99%), and the sets keep 64 to 97% of the candidates: with these inputs
DP4 can rule out little at that level. At an 80% target the coverage is 82 to 84%, the sets hold 2.1 to 3.2
candidates on average, and one candidate is left and right in about a quarter of the cases
(`logs/11_audit_and_conformal.txt`).

## What is not settled

- The pairing rule when a description lists fewer or more values than atoms (7 molecules), exchangeable
  protons and methyl averaging: the 2020 code has no special handling, and neither has this run.
- The three-Gaussian model of the paper is not reproduced.
- MMFF energies for a classical DP4 weighting need the MacroModel files (run 2 of the notebook for NP1).

## Files

| File | Content |
| --- | --- |
| `tables/dp4ai_cases.csv` | one row per molecule, protocol and data type (all, C, H): probability and rank of the reported structure, top isomer, flags |
| `tables/dp4ai_probabilities.csv` | the probability of every candidate |
| `tables/*_NMRold.csv` | the same with the older description files |
| `tables/dp4ai_conformers.csv` | every conformer: energy, population, job flags (level definitions of the original code) |
| `tables/dp4ai_truth.json` | how the reported structure was identified, InChI, candidate-set family |
| `tables/dp4files_explained.json`, `dp4files_validated.json` | every stored DP4 output: which shieldings and energies reproduce it, and how well |
| `tables/dp4ai_descriptions.json` | every version of the experimental description files |
| `logs/` | the text output of each step, numbered in run order |
