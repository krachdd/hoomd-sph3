#!/usr/bin/env bash
# ==============================================================================
# run_all_re.sh  –  BCC sphere packing: Darcy–Forchheimer permeability sweep
#
# Creates a 100×100×100 BCC geometry once, then runs one SPH simulation per
# body force, sweeping Re_grain from ≈ 0.01 to ≈ 10 000 to capture the
# Darcy–Forchheimer transition (onset of inertia).
#
# PHYSICAL SETUP
# --------------
#   Geometry : BCC sphere packing, 4×4×4 unit cells, a = 25 vox, R = 10 vox
#   Domain   : 100 × 100 × 100 mm  (vsize = 1 mm)
#   Grains   : d = 20 mm  (sphere diameter)
#   Fluid    : ρ = 1000 kg/m³,  μ = 0.001 Pa·s  (water at 20 °C)
#   Porosity : φ ≈ 0.464
#   Boundaries: all periodic  (body force drives flow in x-direction)
#
# DIMENSIONLESS NUMBERS
# ---------------------
#   Re_grain = ρ φ <u> d / μ    (grain-scale Reynolds number)
#   k        = μ φ <u> / (ρ fx) (apparent permeability [m²])
#   k_norm   = k / k_KC          (1 in Darcy regime, < 1 in Forchheimer regime)
#
#   Kozeny–Carman estimate: k_KC ≈ 7.7 × 10⁻⁷ m²  (≈ 0.78 Darcy)
#   Ergun β coefficient   : β   ≈ 469 m⁻¹
#
# DARCY–REGIME PREDICTION (Re_grain, linear)
# -------------------------------------------
#   Re_grain ≈ 15 440 × fx    (from k_KC and Darcy's law)
#
#   fx [m/s²]   Re_grain (Darcy)   Regime
#   ---------   ----------------   ------
#   1.0e-6      0.015              Darcy
#   2.0e-6      0.031              Darcy
#   5.0e-6      0.077              Darcy
#   1.0e-5      0.154              Darcy
#   2.0e-5      0.309              Darcy
#   5.0e-5      0.772              Darcy
#   1.0e-4      1.54               onset of inertia zone ↑
#   2.0e-4      3.09               onset of inertia
#   5.0e-4      7.72               onset of inertia  ← expect k_norm < 1.0
#   1.0e-3      15.4               early Forchheimer (actual Re ~ 10–12)
#   2.0e-3      30.9               Forchheimer       (actual Re ~ 15–20)
#   5.0e-3      77.2               Forchheimer       (actual Re ~ 25–40)
#   1.0e-2      154                Forchheimer       (actual Re ~ 35–60)
#   2.0e-2      309                Forchheimer       (actual Re ~ 50–100)
#   5.0e-2      772                Forchheimer       (actual Re ~ 70–150)
#   1.0e-1      1 540              Forchheimer       (actual Re ~ 100–200)
#   3.0e-1      4 630              Forchheimer       (actual Re ~ 150–400)
#   1.0e+0      15 400             inertia-dominant  (actual Re ~ 200–700)
#   3.0e+0      46 200             inertia-dominant  (actual Re ~ 300–1 500)
#   1.0e+1      154 000            turbulent-like    (actual Re ~ 700–3 000)
#   3.0e+1      462 000            turbulent-like    (actual Re ~ 1 200–6 000)
#   1.0e+2      1 540 000          turbulent-like    (actual Re ~ 2 000–10 000)
#
# NOTE on high-Re cases:
#   For fx > 0.1 m/s² the Forchheimer resistance is dominant.  SPH models
#   steady inertial effects correctly up to Re ~ 100–200; beyond that the
#   flow is turbulent and WCSPH gives only time-averaged Forchheimer-regime
#   behaviour.  All cases are included to map the full curve.
#
# MPI: 8 ranks → 2×2×2 decomposition of the 100³ cube (no z constraint).
#
# Usage:  bash run_all_re.sh [steps] [damp]
#   steps : simulation steps per run   (default: 50001)
#   damp  : body-force ramp time [steps] (default: 5000)
#           The body force is smoothly ramped over the first `damp` steps to
#           avoid acoustic shock at startup; does NOT affect steady-state k.
# ==============================================================================
set -euo pipefail

STEPS=${1:-100001}
DAMP=${2:-5000}
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

banner() {
    echo
    echo "══════════════════════════════════════════════════════════════"
    printf "  %-58s\n" "$*"
    echo "══════════════════════════════════════════════════════════════"
}

banner "BCC Sphere Packing — Darcy–Forchheimer sweep  |  steps=${STEPS}"

# ─── Step 1: create geometry (only if init file does not yet exist) ───────────
banner "Create BCC geometry (100³ voxels, a=25, R=10)"
pushd "${SCRIPT_DIR}" > /dev/null
if [[ ! -f bcc100_init.gsd ]]; then
    python3 create_bcc_geometry.py
else
    echo "  bcc100_init.gsd already exists – skipping geometry creation."
fi
INIT="bcc100_init.gsd"
popd > /dev/null

run_case() {
    local FX="$1"
    local LABEL="$2"
    echo
    echo "── fx = ${FX} m/s²   (${LABEL}) ──"
    mpirun -np 8 ./run_bcc_permeability.py \
        "${INIT}" "${FX}" "${STEPS}" "${DAMP}"
}

# ─── Step 2: Darcy regime  (Re_grain ≈ 0.01 – 0.8) ──────────────────────────
banner "Darcy regime  (fx = 1e-6 … 5e-5,  Re ≈ 0.015 – 0.77)"
run_case 1.0e-6  "Re_Darcy ≈ 0.015"
run_case 2.0e-6  "Re_Darcy ≈ 0.031"
run_case 5.0e-6  "Re_Darcy ≈ 0.077"
run_case 1.0e-5  "Re_Darcy ≈ 0.15"
run_case 2.0e-5  "Re_Darcy ≈ 0.31"
run_case 5.0e-5  "Re_Darcy ≈ 0.77"

# ─── Step 3: onset of inertia  (Re_grain ≈ 1 – 10) ──────────────────────────
banner "Onset of inertia  (fx = 1e-4 … 5e-4,  Re ≈ 1.5 – 8)"
run_case 1.0e-4  "Re_Darcy ≈ 1.5  — onset zone"
run_case 2.0e-4  "Re_Darcy ≈ 3.1  — onset zone"
run_case 5.0e-4  "Re_Darcy ≈ 7.7  — expect k_norm < 1.0"

# ─── Step 4: early Forchheimer  (Re_grain ≈ 10 – 100) ───────────────────────
banner "Early Forchheimer  (fx = 1e-3 … 5e-2,  Re ≈ 10 – 150)"
run_case 1.0e-3  "Re ≈ 10–12"
run_case 2.0e-3  "Re ≈ 15–20"
run_case 5.0e-3  "Re ≈ 25–40"
run_case 1.0e-2  "Re ≈ 35–60"
run_case 2.0e-2  "Re ≈ 50–100"
run_case 5.0e-2  "Re ≈ 70–150"

# ─── Step 5: deep Forchheimer / high inertia  (Re_grain ≈ 100 – 10 000) ─────
banner "Deep Forchheimer  (fx = 1e-1 … 1e+2,  Re ≈ 100 – 10 000)"
run_case 1.0e-1  "Re ≈ 100–200"
run_case 3.0e-1  "Re ≈ 150–400"
run_case 1.0e+0  "Re ≈ 200–700"
run_case 3.0e+0  "Re ≈ 300–1500"
run_case 1.0e+1  "Re ≈ 700–3000"
run_case 3.0e+1  "Re ≈ 1200–6000"
run_case 1.0e+2  "Re ≈ 2000–10000"

# ─── Step 6: collect and sort the summary file ───────────────────────────────
banner "Collecting results"
SUMMARY="${SCRIPT_DIR}/bcc100_permeability_summary.dat"
if [[ -f "${SUMMARY}" ]]; then
    echo "  Summary file: ${SUMMARY}"
    # Sort by body force (first column), skip comment lines
    TMPFILE="$(mktemp)"
    grep -v '^#' "${SUMMARY}" | sort -k1 -g > "${TMPFILE}"
    # Rebuild with header
    {
        echo "# fx[m/s2]  Re_grain  k[m2]  k_norm  k_KC[m2]  phi"
        cat "${TMPFILE}"
    } > "${SUMMARY}.sorted"
    mv "${SUMMARY}.sorted" "${SUMMARY}"
    rm -f "${TMPFILE}"

    echo "  Body force     Re_grain      k [m²]        k/k_KC"
    echo "  ----------     --------      ------         ------"
    awk 'NF==6 && !/^#/ {
        printf "  %10.3e     %10.3e    %10.3e     %8.4f\n", $1, $2, $3, $4
    }' "${SUMMARY}"
else
    echo "  No summary file found — individual run logs are in ${SCRIPT_DIR}/"
fi

banner "All BCC permeability runs complete."
echo "  GSD trajectories : bcc100_fx*_run.gsd"
echo "  Per-run logs     : bcc100_fx*_run.log"
echo "  Summary table    : bcc100_permeability_summary.dat"
echo
echo "  To identify onset of inertia: plot k_norm vs Re_grain."
echo "  Onset ≈ Re_grain where k_norm first drops below 0.95."
