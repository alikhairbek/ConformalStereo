#!/bin/sh
# Rebuilds every table and log of the DP4-AI benchmark from the result files of the Kaggle notebook.
# Usage: sh run_all.sh <folder with the unpacked cs_core zips> <work folder> <output folder> [python with RDKit]
# The data files are read as data only; Python is started with -I so that nothing in those folders is imported.
set -e
CORE=$1; WORK=$2; OUT=$3; PY=${4:-python3}
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$WORK" "$OUT" "$OUT/logs"
$PY -I "$HERE/census_dp4ai.py"      "$CORE"                > "$OUT/logs/01_census.txt"
$PY -I "$HERE/qc_geometry.py"       "$CORE" "$WORK"        > "$OUT/logs/02_qc_geometry.txt"
$PY -I "$HERE/truth.py"             "$CORE" "$WORK"        > "$OUT/logs/03_reported_structure.txt"
$PY -I "$HERE/desc_inventory.py"    "$CORE" "$WORK"        > "$OUT/logs/04_nmr_descriptions.txt"
$PY -I "$HERE/check_model.py"       "$CORE"                > "$OUT/logs/05_error_model_check.txt"
$PY -I "$HERE/explain_dp4files.py"  "$CORE" "$WORK"        > "$OUT/logs/06_original_outputs_shifts.txt"
$PY -I "$HERE/validate_dp4.py"      "$CORE" "$WORK"        > "$OUT/logs/07_original_outputs_dp4.txt"
$PY -I "$HERE/build_dp4ai.py"       "$CORE" "$WORK"        > "$OUT/logs/08_isomer_table.txt"
$PY -I "$HERE/dp4ai_run.py"         "$CORE" "$WORK" "$OUT"        > "$OUT/logs/09_dp4_NMR.txt"
$PY -I "$HERE/dp4ai_run.py"         "$CORE" "$WORK" "$OUT" NMRold > "$OUT/logs/10_dp4_NMRold.txt"
$PY -I "$HERE/dp4ai_audit.py"       "$OUT" 2000            > "$OUT/logs/11_audit_and_conformal.txt"
$PY -I "$HERE/swap_test.py"         "$CORE" "$OUT"         > "$OUT/logs/12_swap_test.txt"
echo "done: tables in $OUT, logs in $OUT/logs"
