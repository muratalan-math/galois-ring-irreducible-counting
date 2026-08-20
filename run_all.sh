#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
BUILD="$ROOT/_build"
OUT="$ROOT/_run_outputs"
CXX=${CXX:-g++}
PYTHON=${PYTHON:-python3}

rm -rf "$BUILD" "$OUT"
mkdir -p "$BUILD" "$OUT"

compile() {
  "$CXX" -O3 -std=c++17 "$1" -o "$2"
}

compile "$ROOT/ground_truth/local_direct_enumeration.cpp" "$BUILD/local_direct_enumeration"
compile "$ROOT/ground_truth/global_z8_direct.cpp" "$BUILD/global_z8_direct"
compile "$ROOT/ground_truth/z8_local_all_degrees.cpp" "$BUILD/z8_local_all_degrees"
compile "$ROOT/finite_maps/mixed_precision_fibers.cpp" "$BUILD/mixed_precision_fibers"
compile "$ROOT/finite_maps/nonzero_curve_state.cpp" "$BUILD/nonzero_curve_state"
compile "$ROOT/finite_maps/weighted_origin.cpp" "$BUILD/weighted_origin"

{
  for args in '2 3 2 2' '2 3 2 3' '2 4 2 2' '2 4 2 3' \
              '5 2 1 4' '5 3 1 4' '7 2 1 4' '7 3 1 4' '5 2 2 4' '5 4 1 4'; do
    echo "### $args"
    "$BUILD/local_direct_enumeration" $args
  done
} > "$OUT/local_direct_enumeration.txt"

"$BUILD/global_z8_direct" > "$OUT/global_z8_direct.txt"
"$BUILD/z8_local_all_degrees" > "$OUT/z8_local_all_degrees.txt"

{
  for s in 4 5 6; do
    "$BUILD/local_direct_enumeration" 2 "$s" 1 3
  done
  for s in 5 6; do
    "$BUILD/local_direct_enumeration" 3 "$s" 1 3
  done
} > "$OUT/cubic_direct.txt"

"$BUILD/mixed_precision_fibers" > "$OUT/mixed_precision_fibers.txt"
"$BUILD/nonzero_curve_state" > "$OUT/nonzero_curve_state.txt"
"$BUILD/weighted_origin" > "$OUT/weighted_origin.txt"

"$PYTHON" "$ROOT/finite_maps/gf25_residual_classification.py" > "$OUT/gf25_residual_classification.txt"
"$PYTHON" "$ROOT/proof_boundaries/proof_boundary_checks.py" > "$OUT/proof_boundary_checks.txt"
"$PYTHON" "$ROOT/symbolic/denominator_check.py" > "$OUT/denominator_check.txt"
"$PYTHON" "$ROOT/symbolic/factor_start_index_check.py" > "$OUT/factor_start_index_check.txt"
"$PYTHON" "$ROOT/symbolic/closed_formula_checks.py" > "$OUT/closed_formula_checks.txt"
"$PYTHON" "$ROOT/symbolic/quartic_automaton.py" > "$OUT/quartic_automaton.txt"

for f in \
  local_direct_enumeration \
  global_z8_direct \
  z8_local_all_degrees \
  cubic_direct \
  mixed_precision_fibers \
  nonzero_curve_state \
  weighted_origin \
  gf25_residual_classification \
  proof_boundary_checks \
  denominator_check \
  factor_start_index_check \
  closed_formula_checks \
  quartic_automaton; do
  diff -u "$ROOT/expected/$f.txt" "$OUT/$f.txt"
done

echo "ALL REPRODUCIBILITY CHECKS PASSED"
