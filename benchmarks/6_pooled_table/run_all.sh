#!/bin/sh
# All five benchmarks in one table, and conformal calibration across them: per level, per benchmark, per kind of
# calculation, one for all, and transfer from one benchmark to another.
# usage: sh run_all.sh <folder of ConformalStereo_SI_crosscheck (v3)> <tables folder of ConformalStereo_DP4AI_energy_test>
#                      <folder of ConformalStereo_Cambridge_benchmarks> <output dir> [python]
#   python 3.8+ with numpy and scipy (steps 01 and 02); about two minutes
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
X=$1; AI=$2; CAM=$3; OUT=$4; PY=${5:-python3}
S=$HERE/scripts; T=$OUT/tables; L=$OUT/logs; W=$OUT/work
mkdir -p "$T" "$L" "$W"
$PY $S/export_joc.py "$HERE/lib" "$HERE/inputs" "$T/joc_loglik.csv"                                    > "$L/01_export_2015.txt"
$PY $S/stx_loglik.py "$X/lib" "$X/lib" "$X/inputs" "$X/tables/stx_internal_consistency.csv" "$X/tables/stx_candidate_fit.csv" "$W" 200 > "$L/02_export_stx150.txt"
cp "$W/stx_candidates.csv" "$T/stx_loglik.csv"
$PY $S/pool.py "$T/joc_loglik.csv" "$T/stx_loglik.csv" "$AI" "$CAM/tables17" "$CAM/tables19" "$T"     > "$L/03_pool.txt"
$PY $S/pooled_calibration.py "$T/pooled_candidates.csv" "$T" 300                                       > "$L/04_pooled_calibration.txt"
$PY $S/pooled_calibration.py "$T/pooled_candidates.csv" "$T" 300 20261114 log_ratio log_ratio_energy_4 _dp4_plain > "$L/05_dp4_records_plain_score.txt"
$PY $S/pooled_calibration.py "$T/pooled_candidates.csv" "$T" 300 20261114 log_ratio_energy_4 log_ratio_energy_4 _dp4_energy > "$L/06_dp4_records_energy_averaged_score.txt"
echo "done: tables in $T, logs in $L"
