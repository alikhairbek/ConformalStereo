# From the public data to the pooled table

`data/pooled_candidates.csv` was built in the steps below. Each folder holds the scripts, the logs of the run and
the description that was written with it. The folders do not hold their inputs and intermediate tables:

- the raw Gaussian outputs of the three Cambridge deposits are several gigabytes; the notebooks of step 0 fetch and
  parse them (they were run on Kaggle);
- the shieldings of JOC2015 and STX150 are printed in the Supporting Information of the two papers and are not
  copied here; the scripts of step 1 read text extracted from those files.

The paper itself needs only the pooled table: `analysis/run_all.py` starts from it.

| Step | Folder | What it does |
|---|---|---|
| 0 | `0_fetch_deposits_kaggle/` | downloads the deposits of DP4-AI, Cam2017 and Cam2019 from the repository of the University of Cambridge and parses the Gaussian outputs (three runs; the third adds the coupling constants) |
| 1 | `1_joc2015_stx150/` | DP4+, MM-DP4+ and xTB-DP4+ for JOC2015 and STX150 from the printed shieldings, checked against the published probabilities |
| 2 | `2_cambridge_2017_2019/` | DP4 for Cam2017 and Cam2019 from the parsed outputs, with the single Gaussian of PyDP4; checked against the deposited outputs |
| 3 | `3_dp4ai/` | DP4 for DP4-AI: assignment of the unassigned shift lists, four protocols |
| 4 | `4_dp4ai_conformer_energies/` | the likelihood averaged over errors in the conformer energies, tested on DP4-AI |
| 5 | `5_coupling_constants/` | proton-proton coupling constants of 33 compounds of Cam2017 (`data/coupling_candidates.csv`) |
| 6 | `6_pooled_table/` | puts the five benchmarks into one table of candidates |

Sources of the data and their identifiers are listed in the main README.

Third-party material: the parameters of the DP4+ levels are those of DP4+App (Franco, Luciano, Sarotti and Zanardi,
*J. Nat. Prod.* 2023, 86, 2360; package `dp4plus-app`, MIT licence). The treatment of the coupling constants follows
PyDP4 (Ermanis, Parkes, Agback and Goodman; MIT licence, notice in `5_coupling_constants/LICENSE_PyDP4.txt`).
