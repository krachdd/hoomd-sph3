#!/bin/bash
# WP2 piston-driven invasion sweep on the wolfgang A100 nodes: one sbatch job
# per (Ca, Y) case, each on its own GPU (grid as used for the CPU sweep, res=12,
# resA_mm=12, ramp=2000 defaults).  Usage:  ./sweep_piston_gpu.sh [--dry-run]
#   case        tau_y  U_p   sigma  steps   m_gd
cases=(
 "Ca02_Y0      0      0.02  0.01   130000  1"
 "Ca02_Y3      37.5   0.02  0.01   431000  1"
 "Ca1_Y0       0      0.04  0.004   92000  1"
 "Ca1_Y3       15     0.04  0.004  112000  1"
 "Ca1_Y12      60     0.04  0.004  181000  1"
 "Ca1_Y48      240    0.04  0.004  590000  1"
 "Ca4_Y0       0      0.04  0.001   92000  1"
 "Ca4_Y3       3.75   0.04  0.001   98000  1"
 "Ca4_Y12      15     0.04  0.001  112000  1"
 "Ca4_Y12_m3   15     0.04  0.001  148000  3"
 "Ca4_Y48      60     0.04  0.001  181000  1"
)
RES=12
cd "$(dirname "${BASH_SOURCE[0]}")"
# shared init domain (U_p only affects the printed travel estimate, not the geometry)
if [[ ! -f piston_domain_res${RES}_init.gsd && "${1:-}" != "--dry-run" ]]; then
    mpirun -np 1 python3 create_piston_domain.py --U_p 0.04 --res ${RES} --out piston_domain_res${RES}_init.gsd
fi
for c in "${cases[@]}"; do
    read -r name tau_y U_p sigma steps m_gd <<< "$c"
    cmd="sbatch --job-name=inv_${name} invasion_piston_gpu.sh ${tau_y} ${U_p} ${sigma} ${RES} ${steps} --m_gd ${m_gd}"
    echo "$cmd"
    [[ "${1:-}" == "--dry-run" ]] || $cmd
done
