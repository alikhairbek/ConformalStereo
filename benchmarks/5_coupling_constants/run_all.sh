#!/bin/sh
# Usage: run_all.sh <res/core folder of Kaggle run 2> <res/core folder of Kaggle run 3> <out folder> [python]
# Python 3.9 or later with numpy and scipy. About five minutes.
set -e
C2="$1"; C3="$2"; OUT="$3"; PY="${4:-python3}"
HERE="$(cd "$(dirname "$0")" && pwd)"
L="$HERE/lib"; S="$HERE/scripts"; I="$HERE/inputs"; T="$OUT/tables"; G="$OUT/logs"
mkdir -p "$T" "$G"
step() { n="$1"; shift; echo "== $n"; "$PY" -I -B "$@" > "$G/$n.txt" 2>&1 || { tail -5 "$G/$n.txt"; exit 1; }; tail -2 "$G/$n.txt"; }
step 1_jcheck   "$S/jcheck.py"   "$C3" "$C2" "$L" "$T"
step 2_jsource  "$S/jsource.py"  "$C2" "$L" "$I/dep17_experiments.json" "$T/jsource.json"
step 3_jrepro   "$S/jrepro.py"   "$C2" "$C3" "$L" "$I/dep17_experiments.json" "$T/jsource.json" "$T"
step 4_jmodel   "$S/jmodel.py"   "$T/jrepro_couplings.json" "$T"
step 5_jpreview "$S/jpreview.py" "$C2" "$L" "$I/dep17_experiments.json" "$I/dep17_truth.csv" "$T" "$I/dep17_records_level1.csv"
step 6_jbench   "$S/jbench.py"   "$C2" "$C3" "$L" "$T/jsource.json" "$I/dep17_truth.csv" "$T/jmodel_recovered.json" "$T" 500 20261114
step 7_jeval    "$S/jeval.py"    "$T/jbench_candidates.csv" "$T/jbench_couplings.csv" "$T/jpreview_outputs.csv" "$T/jmodel_recovered.json" "$T" 500 20261114
echo "done: tables in $T, logs in $G"
