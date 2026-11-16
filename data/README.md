# Data

## `pooled_candidates.csv`

One row for every candidate structure of every assignment: 14,437 candidates in 3,100 assignments (records) of
152 families of related molecules, computed with 76 workflows in five benchmarks.

| Column | Meaning |
|---|---|
| `benchmark` | JOC2015, STX150, DP4-AI, Cam2017 or Cam2019 |
| `method` | DP4, DP4+, MM-DP4+ or xTB-DP4+ |
| `geometry` | method of the geometries (MMFF, xTB, PM6, PM7, HF, B3LYP, M06-2X) |
| `level` | NMR level, labelled as in the source; with the three columns before it, this names the workflow |
| `case` | the assignment: one experimental data set with its candidates |
| `family` | molecules that share a skeleton or come from one series; the unit of calibration |
| `clean` | 0 for the 83 records with a known data problem, which the paper reports separately; 1 otherwise |
| `candidate` | name of the candidate structure |
| `is_true` | 1 for the actual structure, 0 otherwise |
| `n_nuclei` | number of shifts N in the likelihood |
| `log_ratio` | natural logarithm of the likelihood of the candidate, minus that of the best candidate of the record (zero for the best) |
| `log_ratio_energy_4` | the same for the likelihood averaged over a noise of 4 kJ/mol on the conformer energies; empty where the source holds no conformer data (filled for 1,108 records) |

The stated probability of a candidate is `exp(log_ratio)` divided by the sum of `exp(log_ratio)` over its record.
The rows of one workflow are a valid calibration table for the companion program (page *Calibrate*).

The log-likelihoods were computed again for this work: DP4+, MM-DP4+ and xTB-DP4+ from the shieldings printed in the
Supporting Information of the two source papers, with the published parameters of each level; DP4 from the raw
Gaussian outputs of the three Cambridge deposits, with the single Gaussian of PyDP4 (2.269 ppm for carbon, 0.187 ppm
for proton, internally scaled shifts). Table S1 of the paper gives the rules and the agreement with the published
numbers; the code is in `benchmarks/`.

## `coupling_candidates.csv`

The 33 compounds of the Cam2017 deposit with calculated proton-proton coupling constants (B3LYP/6-31G(d,p) on MMFF
geometries), one row per candidate: log-likelihoods of the carbon shifts, of the proton shifts and of the couplings
under three error models (`logL_J_gauss`, the Gaussian of the PyDP4 source; `logL_J_fc`, the same for the
Fermi-contact part; `logL_J_rec`, the density recovered from the deposited outputs), and the likelihoods averaged
over noise of 2, 4 and 8 kJ/mol on the conformer energies (`E2_...`, `E4_...`, `E8_...`). Used for Table S12.

## Sources and terms

The three Cambridge deposits are published under CC BY 4.0 (identifiers in the main README). The shieldings of
JOC2015 and STX150 were read from the Supporting Information of the two papers; they are not copied here. What is
stored here are quantities computed from them.
