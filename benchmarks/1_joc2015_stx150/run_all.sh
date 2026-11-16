#!/bin/sh
# STX150 validation set: cross-check against the SI of the two earlier papers and against itself, DP4+ of the earlier
# papers recomputed from their printed shieldings, and the table of STX150 records that can be used without the
# authors' files.
# usage: sh run_all.sh <np3c00566_si_001.txt> <d4ob01959h1_suppl.txt> <output dir> [python]
#   the two text files are made with:  pdftotext -layout <file>.pdf <file>.txt
#   python needs numpy and scipy for steps 11 and 21-27; steps 21-24 draw random numbers with fixed seeds and take a few minutes
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
JNP=$1; OBC=$2; OUT=$3; PY=${4:-python3}
S=$HERE/scripts; IN=$HERE/inputs; T=$OUT/tables; L=$OUT/logs
mkdir -p "$T" "$L"
# --- SI of the 2023 paper (J. Nat. Prod. 2023, 86, 2360) ---
$PY $S/jnp_tables.py "$JNP" "$T"                                                            > "$L/01_tables_2023.txt"
$PY $S/jnp_probs.py "$JNP" "$T"                                                             > "$L/02_probabilities_2023.txt"
$PY $S/stx_vs_jnp.py "$T" "$IN"                                                             > "$L/03_stx_table_vs_2023.txt"
$PY $S/stx_jnp_shieldings.py "$T/jnp_tables.json" "$IN" "$T/stx_vs_jnp_shieldings.csv"      > "$L/04_stx_shieldings_vs_2023.txt"
for LEV in MMFF_PBE xTB_PBE B3L_PBE; do
$PY $S/stx_jnp_identify.py "$T/jnp_tables.json" "$IN" "$T/stx_jnp_identify_$LEV.csv" $LEV MPW_B_PCM > "$L/05_stx_candidates_in_2023_$LEV.txt"
done
$PY $S/stx_jnp_89.py "$T/jnp_tables.json" "$IN"                                             > "$L/06_case_89a.txt"
$PY $S/jnp_table_order.py "$T/jnp_tables.json" "$T/jnp_probabilities.csv" "$IN" "$T/stx_jnp_identify_MMFF_PBE.csv" > "$L/07_columns_of_2023_tables.txt"
$PY $S/si_internal_consistency.py jnp "$T/jnp_tables.json"                                  > "$L/08_internal_test_2023.txt"
# --- SI of the 2025 paper (Org. Biomol. Chem. 2025, 23, 2480) ---
$PY $S/obc_tables.py "$OBC" "$T"                                                            > "$L/09_tables_2025.txt"
$PY $S/obc_probs.py "$OBC" "$IN/stx110_pcorrect.csv" "$T/stx_vs_obc_probabilities.csv"      > "$L/10_stx_table_vs_2025.txt"
$PY $S/obc_dp4plus.py "$HERE/lib" "$IN/joc_dp4plus_parameters.csv" "$T/obc_compounds.json" "$OBC" "$T" > "$L/11_dp4plus_levelB_2025.txt"
$PY $S/obc_internal.py "$T/obc_compounds.json"                                              > "$L/12_internal_test_2025.txt"
$PY $S/stx_obc_identify.py "$T/obc_compounds.json" "$IN" "$T/stx_obc_identify.csv"          > "$L/13_stx_shieldings_vs_2025.txt"
# --- the STX150 SI against itself ---
$PY $S/stx_copied_columns.py "$IN"                                                          > "$L/14_stx_copied_columns.txt"
$PY $S/stx_column_checks.py "$IN"                                                           > "$L/15_stx_column_distances.txt"
$PY $S/stx_internal_consistency.py "$IN" "$T/stx_internal_consistency.csv"                  > "$L/16_stx_internal_test.txt"
$PY $S/stx_exchanged_candidates.py "$T/stx_internal_consistency.csv"                        > "$L/17_stx_exchanged_candidates.txt"
$PY $S/stx_candidate_fit.py "$IN" "$T/stx_candidate_fit.csv"                                > "$L/18_stx_candidate_fit.txt"
# --- DP4+ of the earlier papers recomputed from their printed shieldings (MM-DP4+ / DP4+ as coded in DP4+App 2.2.3) ---
mkdir -p "$T/mm"
$PY $S/jnp_tables_s7.py "$JNP" "$T/mm"                                                      > "$L/19_tables_2023_sections_7_and_9.txt"
$PY $S/jnp_section8.py "$JNP" "$T/jnp_section8.csv"                                         > "$L/20_published_2023_section_8.txt"
$PY $S/jnp_verify_s7.py "$HERE/lib" "$HERE/lib" "$IN/app_parameters.csv" "$T/mm/jnp_tables.json" "$T/jnp_section8.csv" "$T" > "$L/21_verify_2023_36_levels.txt"
$PY $S/jnp_verify_s9.py "$HERE/lib" "$HERE/lib" "$IN/app_parameters.csv" "$T/mm/jnp_tables.json" "$T/jnp_probabilities.csv" "$T" > "$L/22_verify_2023_compounds_101_121.txt"
$PY $S/obc_verify.py "$HERE/lib" "$HERE/lib" "$IN/app_parameters.csv" "$T/obc_compounds.json" "$OBC" "$T" > "$L/23_verify_2025_levels_A_B_C.txt"
# --- STX150 records that can be used ---
$PY $S/stx_verified.py "$HERE/lib" "$HERE/lib" "$IN" "$T/stx_internal_consistency.csv" "$T/stx_candidate_fit.csv" "$T" > "$L/24_stx_records.txt"
$PY $S/stx_differences.py "$T"                                                              > "$L/25_stx_differences.txt"
$PY $S/stx_preview.py "$T" 2000                                                             > "$L/26_stx_preview.txt"
: > "$L/27_sp2_sensitivity.txt"
for C in 89a 111a 115a 106a; do $PY $S/stx_sp2_sensitivity.py "$HERE/lib" "$IN" $C >> "$L/27_sp2_sensitivity.txt"; done
echo "done: tables in $T, logs in $L"
