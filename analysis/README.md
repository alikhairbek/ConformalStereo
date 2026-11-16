# Analysis

Every table and figure of the paper is produced from one input, `data/pooled_candidates.csv`.

    pip install numpy scipy matplotlib
    python analysis/run_all.py              # tables and logs, about two minutes
    python analysis/run_all.py --figures    # also draws Figures 1 to 7 and S1 to S4 into analysis/figures

Python 3.9 or newer. The seeds are fixed, so a rerun gives the same files. The tables carry eight decimals; every
number in the paper is rounded once, from these values.

## Scripts

| Script | What it does | Writes |
|---|---|---|
| `cslib.py` | the scores, the family-weighted cut-off, the sets, the resampling of families | (library) |
| `final_method.py` | the method on the five benchmarks: reliability, sets, origin of the cut-off, stress tests, classes, kinds of calculation, conformer energies | `tables/final_*.csv`, `logs/01_final_method.txt` |
| `temperature_size.py` | temperature against the number of shifts | `tables/final_temperature_size.csv`, `logs/02_temperature_size.txt` |
| `check_final.py` | checks of the procedure: labels drawn at random, invariances, consistency of the tables | `logs/03_checks.txt` |
| `explore_candidates.py` | the score divided also by the number of candidates (not adopted) | `logs/04_explore_candidates.txt` |
| `exhibits.py` | the worked example, temperatures by thirds, kinds of sets, coverage of single workflows, misses by family, coupling constants, sources | `tables/exhibit_*.csv`, `logs/05_exhibits.txt` |
| `paper_tables.py` | Tables 1 to 5 and S1 to S13 as printed | `tables/paper/` |
| `make_cutoffs.py` | the cut-offs of the 76 workflows for the companion program, computed with the program's own functions and compared with `final_workflows.csv` | `ConformalStereo_App/app/cutoffs.csv` |
| `make_examples.py` | the example files of the companion program | `ConformalStereo_App/app/examples/` |
| `figures.py`, `si_figures.py` | the figures | `figures/` (not stored in the repository) |

## Where each table of the paper comes from

`tables/paper/` holds every table as plain CSV (`Table_1.csv` ... `Table_S13.csv`), and `tables/paper/index.csv`
repeats this list with titles, captions and notes.

| Table | Content | Built from |
|---|---|---|
| Table 1 | the five benchmarks | `final_benchmarks.csv` |
| Table 2 | calibrated sets and their price | `final_scores.csv`, `exhibit_set_kinds.csv` |
| Table 3 | where the cut-off may come from | `final_designs.csv`; number of workflows below the level from `logs/01_final_method.txt` |
| Table 4 | stress tests | `final_stress.csv` |
| Table 5 | kinds of calculation after calibration | `final_kinds.csv`, `final_scores.csv` |
| Table S1 | sources and reproduction | `exhibit_sources.csv` |
| Table S2 | families or records as units | `final_units.csv` |
| Table S3 | checks of the procedure | `logs/03_checks.txt`, `final_random_halves.csv` |
| Table S4 | reliability of the stated probabilities | `final_reliability.csv` |
| Table S5 | temperature and size of the molecule | `final_temperature_size.csv`, `exhibit_temperature_bins.csv` |
| Table S6 | scores and baselines | `final_scores.csv` |
| Table S7 | cut-offs of the 76 workflows | `final_workflows.csv` |
| Table S8 | coverage by class | `final_conditional.csv` |
| Table S9 | misses by family | `final_misses.csv` |
| Table S10 | kinds of calculation compared | `final_kinds_paired.csv` |
| Table S11 | uncertain conformer energies | `final_energy.csv` |
| Table S12 | coupling constants | `exhibit_couplings.csv` |
| Table S13 | the number of candidates in the score | `logs/04_explore_candidates.txt` |

## Where each figure comes from

| Figure | Built from |
|---|---|
| Figure 1 | `exhibit_example.csv` |
| Figure 2 | `final_reliability.csv` |
| Figure 3 | `exhibit_temperature_bins.csv`, `final_temperature_size.csv` |
| Figure 4 | `exhibit_set_kinds.csv` |
| Figure 5 | `exhibit_workflow_coverage.csv` |
| Figure 6 | `final_stress.csv`, `exhibit_family_misses.csv` |
| Figure 7 | `final_workflows.csv` |
| Figure S1 | `final_records.csv` |
| Figure S2 | `final_workflows.csv` |
| Figure S3 | `data/pooled_candidates.csv` (one curve per workflow) |
| Figure S4 | `data/pooled_candidates.csv` (columns `log_ratio` and `log_ratio_energy_4`) |

Other tables in `tables/`: `final_records.csv` lists every one of the 3,100 assignments with its sets;
`final_confident_calls.csv` and `final_scores_paired.csv` hold numbers quoted in the text.
