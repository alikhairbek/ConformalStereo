#!/bin/sh
# DP4-AI deposit (doi 10.17863/CAM.47721): a pairing rule for incomplete shift lists, and the energy-averaged
# likelihood tested on molecules that the two Cambridge benchmarks do not contain.
# usage: sh run_all.sh <res/core of Kaggle run 1> <output dir> [python]
#   <res/core> holds gaussian_dp4ai_*.jsonl.gz and small_index_dp4ai.tsv.gz; inputs/ holds two tables of the archive
#   ConformalStereo_DP4AI_benchmark (reported structures, geometry check). python 3.8+ with numpy and scipy; under a minute
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
CORE=$1; OUT=$2; PY=${3:-python3}
S=$HERE/scripts; T=$OUT/tables; L=$OUT/logs
mkdir -p "$T" "$L"
$PY $S/pairing_test.py                                                                      > "$L/00_pairing_test.txt"
$PY $S/dp4ai_energy.py "$CORE" "$HERE/inputs" "$T" 1000                                     > "$L/01_dp4ai_energy.txt"
echo "== all records"                                                                       > "$L/02_energy_noise_sets.txt"
$PY $S/energy_marginal_eval.py "$T/dp4ai_energy_marginal.csv" "$T/dp4ai_truth.csv" 500     >> "$L/02_energy_noise_sets.txt"
echo "== clean records with one value per atom"                                            >> "$L/02_energy_noise_sets.txt"
$PY $S/energy_marginal_eval.py "$T/dp4ai_energy_marginal_clean.csv" "$T/dp4ai_truth.csv" 500 >> "$L/02_energy_noise_sets.txt"
$PY $S/dp4ai_summary.py "$T/dp4ai_energy_marginal.csv" "$T/dp4ai_truth.csv" "$T/dp4ai_records.csv" > "$L/03_dp4ai_summary.txt"
echo "done: tables in $T, logs in $L"
