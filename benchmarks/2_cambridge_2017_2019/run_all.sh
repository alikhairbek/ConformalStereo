#!/bin/sh
# Two Cambridge deposits of the Goodman group, recomputed from the raw Gaussian outputs:
#   2017  doi 10.17863/CAM.13222  42 compounds, twelve NMR levels on MMFF geometries      (Org. Biomol. Chem. 2017, 15, 8998)
#   2019  doi 10.17863/CAM.32463  25 of the same experiments, three NMR levels on six kinds of geometry
# and a first comparison of conformal scores on both.
# usage: sh run_all.sh <res/core of Kaggle run 2> <res/core of Kaggle run 1> <output dir> [python]
#   run 2 (notebook ConformalStereo_Cambridge_fetch v1.1) holds gaussian_dp4x2_*, small_dp4x2_001.zip, small_index_dp4x2.tsv.gz
#   run 1 (notebook v1.0) holds gaussian_conf2018_*, small_conf2018_001.zip, small_index_conf2018.tsv.gz
#   python 3.8+; steps 25 to 29 need numpy and scipy, the rest uses the standard library; about two minutes in all
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
C17=$1; C19=$2; OUT=$3; PY=${4:-python3}
S=$HERE/scripts; IN=$HERE/inputs; T=$OUT/tables17; U=$OUT/tables19; L=$OUT/logs; W=$OUT/work17; V=$OUT/work19
mkdir -p "$T" "$U" "$L" "$W" "$V"
# --- 2017 deposit ---
export CS_DEPOSIT=dp4x2
$PY $S/runs17.py "$C17" "$W/runs.json"                                                              > "$L/01_2017_outputs_of_the_deposit.txt"
$PY $S/build17.py "$C17" "$W/structures.jsonl.gz"                                                   > "$L/02_2017_structures.txt"
$PY $S/validate17.py "$C17" "$W/structures.jsonl.gz" "$W/runs.json" "$W/validate.json"              > "$L/03_2017_shifts_of_all_outputs.txt"
$PY $S/experiments17.py "$W/runs.json" "$W/validate.json" "$T/dep17_experiments.json"               > "$L/04_2017_experiments.txt"
$PY $S/energy_source17.py "$C17" "$T/dep17_experiments.json" "$W/energy_source.json"                > "$L/05_2017_energies_behind_the_outputs.txt"
$PY $S/dp4_17.py "$C17" "$T/dep17_experiments.json" "$W/energy_source.json" "$W/structures.jsonl.gz" "$T" > "$L/06_2017_dp4_recomputed.txt"
$PY $S/truth17.py "$IN/obc2017_si_s4_b3lyp.csv" "$T/dep17_shifts.csv" "$T/dep17_candidates.csv" BF1 "$T/dep17_truth.csv" > "$L/07_2017_actual_structures.txt"
$PY $S/si_check17.py "$C17" "$T/dep17_experiments.json" "$W/energy_source.json" "$W/structures.jsonl.gz" "$T/dep17_truth.csv" "$IN" "$T/dep17_si_check.csv" > "$L/08_2017_si_tables.txt"
$PY $S/benchmark17.py "$T/dep17_candidates.csv" "$T/dep17_truth.csv" "$T" 2000                      > "$L/09_2017_benchmark.txt"
$PY $S/qc17.py "$C17" "$T/dep17_experiments.json" "$T/dep17_truth.csv" "$T/dep17_qc.csv"            > "$L/10_2017_quality_checks.txt"
$PY $S/loglik.py "$C17" "$T/dep17_experiments.json" "$W/energy_source.json" "$W/structures.jsonl.gz" "$T/dep17_loglik.csv" > "$L/11_2017_log_likelihoods.txt"
# --- 2019 deposit ---
export CS_DEPOSIT=conf2018
$PY $S/runs17.py "$C19" "$V/runs.json"                                                              > "$L/12_2019_outputs_of_the_deposit.txt"
$PY $S/build17.py "$C19" "$V/structures.jsonl.gz"                                                   > "$L/13_2019_structures.txt"
$PY $S/validate17.py "$C19" "$V/structures.jsonl.gz" "$V/runs.json" "$V/validate.json"              > "$L/14_2019_shifts_of_all_outputs.txt"
$PY $S/experiments17.py "$V/runs.json" "$V/validate.json" "$U/dep19_experiments.json"               > "$L/15_2019_experiments.txt"
$PY $S/energy_source17.py "$C19" "$U/dep19_experiments.json" "$V/energy_source.json"                > "$L/16_2019_energies_behind_the_outputs.txt"
$PY $S/dp4_17.py "$C19" "$U/dep19_experiments.json" "$V/energy_source.json" "$V/structures.jsonl.gz" "$U" > "$L/17_2019_dp4_recomputed.txt"
$PY $S/truth19.py "$T" "$U" "$U/dep19_truth.csv"                                                    > "$L/18_2019_actual_structures.txt"
$PY $S/workflows19.py "$C19" "$U/dep19_experiments.json" "$V/energy_source.json" "$V/structures.jsonl.gz" "$U/dep19_workflows.csv" > "$L/19_2019_workflows.txt"
$PY $S/loglik.py "$C19" "$U/dep19_experiments.json" "$V/energy_source.json" "$V/structures.jsonl.gz" "$U/dep19_loglik.csv" > "$L/20_2019_log_likelihoods.txt"
$PY $S/benchmark19.py "$U/dep19_workflows.csv" "$U/dep19_truth.csv" "$U" 2000                       > "$L/21_2019_benchmark.txt"
$PY $S/qc19.py "$C19" "$U/dep19_experiments.json" "$U/dep19_qc.csv"                                 > "$L/22_2019_quality_checks.txt"
# --- conformal scores ---
$PY $S/scores.py "$T/dep17_loglik.csv" "$T/dep17_truth.csv" "$T/dep17_scores.csv" 1000              > "$L/23_2017_scores.txt"
$PY $S/scores.py "$U/dep19_loglik.csv" "$U/dep19_truth.csv" "$U/dep19_scores.csv" 500               > "$L/24_2019_scores.txt"
$PY $S/error_models.py "$T/dep17_shifts.csv" "$T/dep17_truth.csv" "$T/dep17_error_models.csv" 300   > "$L/25_2017_error_models.txt"
CS_DEPOSIT=dp4x2 $PY $S/energy_marginal.py "$C17" "$T/dep17_experiments.json" "$W/energy_source.json" "$W/structures.jsonl.gz" "$T/dep17_truth.csv" "$T/dep17_energy_marginal.csv" 200 > "$L/26_2017_energy_noise.txt"
$PY $S/energy_marginal_eval.py "$T/dep17_energy_marginal.csv" "$T/dep17_truth.csv" 300             > "$L/27_2017_energy_noise_sets.txt"
$PY $S/energy_marginal.py "$C19" "$U/dep19_experiments.json" "$V/energy_source.json" "$V/structures.jsonl.gz" "$U/dep19_truth.csv" "$U/dep19_energy_marginal.csv" 200 > "$L/28_2019_energy_noise.txt"
$PY $S/energy_marginal_eval.py "$U/dep19_energy_marginal.csv" "$U/dep19_truth.csv" 300             > "$L/29_2019_energy_noise_sets.txt"
echo "done: tables in $T and $U, logs in $L"
