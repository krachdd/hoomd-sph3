#!/bin/bash
# ----------------------------------------------------------
# maintainer: dkrach, david.krach@mib.uni-stuttgart.de
#
# Rayleigh-number sweep for the disordered-pack fingering benchmark.
#
#   Ra = drho * g * K * H / (phi * mu * D)
#
# Workflow:
#   1. build geometry (pack + periodic permeability cell, same seed)
#   2. measure K on the periodic cell
#   3. paste the printed K below (K_M2) and run the D sweep
#
# The three D values target Ra ~ 200 / 1000 / 5000 for the full pack
# (H = 12 mm, phi ~ 0.7, K ~ 1e-9..1e-8 measured). Recompute them after
# the permeability run:  D = drho*g*K*H/(phi*mu*Ra_target).
#
# RESOLUTION: Pe_grid = Ra*phi*dx/H regardless of mu and D, so
# Ra_max ~ 4H/(phi*dx): at 10 dx/diameter Ra_max ~ 1600 — run the Ra=5000
# member at 20 dx/diameter (or accept Pe ~ 12 knowingly).
#
# Usage:  ./run_all_ra.sh <num_diameter> [seed]
# ----------------------------------------------------------
set -e

ND=${1:-10}
SEED=${2:-42}
NP=${NP:-4}                      # MPI ranks for the fingering runs
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-8}

K_M2=${K_M2:?"Set K_M2 (measured permeability, m^2) in the environment first:
  python3 create_input_geometry.py $ND $SEED --perm
  python3 run_permeability.py porous_ra_perm_${ND}_seed${SEED}_init.gsd
  K_M2=<printed K> ./run_all_ra.sh $ND $SEED"}

INIT=porous_ra_pack_${ND}_seed${SEED}_init.gsd
[ -f "$INIT" ] || python3 create_input_geometry.py "$ND" "$SEED"

# D values for Ra ~ 200 / 1000 / 5000  (recomputed from K_M2; H=12mm phi=0.7 mu=1e-3)
for RA in 200 1000 5000; do
    D=$(python3 -c "print(f'{50*9.81*${K_M2}*0.012/(0.7*1e-3*${RA}):.3e}')")
    echo "── Ra = $RA  →  D = $D m²/s ──"
    mpirun -np "$NP" --bind-to none python3 run_fingering_ra.py "$INIT" "$K_M2" "$D"
done

# stable reference (light on top) at the middle Ra
D=$(python3 -c "print(f'{50*9.81*${K_M2}*0.012/(0.7*1e-3*1000):.3e}')")
mpirun -np "$NP" --bind-to none python3 run_fingering_ra.py "$INIT" "$K_M2" "$D" --stable
