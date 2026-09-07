#!/bin/bash
# WP2 body-force-driven invasion sweep on the wolfgang A100 nodes: one sbatch
# job per (Y, g/g_crit) case (grid as used for the CPU sweep, res=12, sigma=0.01,
# vcap=1.0 and damp=1000 defaults).  Usage:  ./sweep_bodyforce_gpu.sh [--dry-run]
#   case      tau_y  g       steps   m_gd
cases=(
 "Y12_g06   150    112.5    50000  1"
 "Y12_g09   150    168.75   50000  1"
 "Y12_g14   150    253.1    50000  1"
 "Y12_g20   150    375.0    50000  1"
 "Y48_g06   600    450     111000  1"
 "Y48_g09   600    675     111000  1"
 "Y48_g14   600    1012.5  111000  1"
 "Y48_g20   600    1500    111000  1"
)
RES=12; SIGMA=0.01
cd "$(dirname "${BASH_SOURCE[0]}")"
if [[ ! -f bodyforce_domain_res${RES}_init.gsd && "${1:-}" != "--dry-run" ]]; then
    mpirun -np 1 python3 create_piston_domain.py --res ${RES} --no-piston --out bodyforce_domain_res${RES}_init.gsd
fi
for c in "${cases[@]}"; do
    read -r name tau_y g steps m_gd <<< "$c"
    cmd="sbatch --job-name=inv_${name} invasion_bodyforce_gpu.sh ${tau_y} ${g} ${SIGMA} ${RES} ${steps} --m_gd ${m_gd}"
    echo "$cmd"
    [[ "${1:-}" == "--dry-run" ]] || $cmd
done
