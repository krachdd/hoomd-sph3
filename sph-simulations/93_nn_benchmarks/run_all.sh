#!/bin/bash
# Run all non-Newtonian benchmarks (num_length=20, default step counts).
# Each benchmark creates its own geometry, runs, and reports L2 errors.

set -e
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
export PYTHONPATH="$ROOT/hoomd-blue/build:$PYTHONPATH"

NL=${1:-20}   # particles per reference length (resolution)
echo "════════════════════════════════════════════════════════════"
echo " Non-Newtonian SPH Benchmarks  (num_length=$NL)"
echo "════════════════════════════════════════════════════════════"

echo ""
echo "── 01 Power-Law Plane Poiseuille ──────────────────────────"
(cd "$(dirname "$0")/01_powerlaw_poiseuille" && python3 run.py $NL)

echo ""
echo "── 02 Bingham Plane Poiseuille ─────────────────────────────"
(cd "$(dirname "$0")/02_bingham_poiseuille"  && python3 run.py $NL)

echo ""
echo "── 03 Two-Phase Power-Law Couette ──────────────────────────"
(cd "$(dirname "$0")/03_twophase_powerlaw_couette" && python3 run.py $NL)

echo ""
echo "── 04 Bingham Slump Test ───────────────────────────────────"
(cd "$(dirname "$0")/04_bingham_slump_test" && python3 run.py $NL)

echo ""
echo "════════════════════════════════════════════════════════════"
echo " All benchmarks complete."
echo "════════════════════════════════════════════════════════════"
