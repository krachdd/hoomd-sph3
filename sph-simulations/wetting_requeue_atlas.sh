#!/bin/bash
# Wetting-benchmark requeue for atlas-101 (2x10-core E5-2680 v2, 2 NUMA nodes,
# no scheduler). Reruns the wetting-affected WP1 cases against the FIXED build:
#   c0..c6    capillary rise, NL=10, theta = 30..150 deg
#   c2_hires  capillary rise, NL=15, theta = 60
#   b1 / b2   sessile droplet shear / snapback
#
# WETTING FIX VALIDATED 2026-08-23 (curvature-form CSF): cleared to launch.
# POSTPROCESSING NOTE: contact-line mobility under no-slip makes the rise
# slow — fit h(t) (Washburn/exponential) for h_inf instead of reading the
# endpoint, and compare against Jurin at the EFFECTIVE angle (sessile
# harness map). Safe to re-run: .done markers skip finished cases.
#
# Layout: FIVE concurrent jobs x 4 OpenMP threads (20 physical cores).
# 4-thread OpenMP runs at ~80% parallel efficiency per core vs ~65% at 8
# threads, so wider-but-slower packing maximizes total throughput; jobs are
# ordered longest-first to minimize the makespan (~9 days total).
#
# TOTAL ESTIMATED WALL TIME ON ATLAS: ~9 days.
#
S30=545000   # theta 30/150 : 3*tau = 0.43 s / dt 7.9e-7
S45=450000   # theta 45/135 : 0.35 s
S60=320000   # theta 60/120 : 0.25 s
S90=100000   # theta 90     : reference, no rise expected
SHIRES=150000  # NL=15 theta 60: 1*tau only -- needs exponential fit in
               # postprocessing; 3*tau (475k) would be ~12 days.

set -u
DIR="$(cd "$(dirname "$0")" && pwd)"       # sph-simulations/
REPO="$(dirname "$DIR")"                   # hoomd-sph3/
RUNDIR="${WETTING_RUNDIR:-$HOME/wetting_requeue}"
export PYTHONPATH="$REPO/hoomd-blue/build:$REPO/dependencies/gsd-sph/gsd/build:$REPO/helper_modules"
export OMP_NUM_THREADS=4

CR="$DIR/01_twophaseflow_benchmarks/03_capillary_rise"
B6="$DIR/01_twophaseflow_benchmarks/06_sessile_droplet_shear"
B7="$DIR/01_twophaseflow_benchmarks/07_sessile_droplet_snapback"

mkdir -p "$RUNDIR"/{caprise,b1,b2}
stamp() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$RUNDIR/status.log"; }

run_case() {   # run_case <name> <workdir> <cmd...>
    local name="$1" wd="$2"; shift 2
    if [ -f "$wd/.done_$name" ]; then stamp "SKIP  $name"; return 0; fi
    stamp "START $name"
    ( cd "$wd" && "$@" > "$wd/${name}.out" 2>&1 )
    local rc=$?
    if [ $rc -eq 0 ]; then touch "$wd/.done_$name"; stamp "DONE  $name"
    else stamp "FAIL  $name (rc=$rc)"; fi
}

DC="$RUNDIR/caprise"

# geometries (minutes)
run_case create10 "$DC" python3 "$CR/create_capillary_geometry.py" 10   # ~2 min
run_case create15 "$DC" python3 "$CR/create_capillary_geometry.py" 15   # ~5 min
INIT10=$(ls -t "$DC"/caprise_40_*_init.gsd | head -1)
INIT15=$(ls -t "$DC"/caprise_60_*_init.gsd | head -1)

# Wave A: the five longest jobs, 4 threads each (20 cores busy)
run_case c0       "$DC" python3 "$CR/run_capillary_rise.py" 10 "$INIT10" 0 $S30   &  # ~8900 min (6.2 d)
run_case c6       "$DC" python3 "$CR/run_capillary_rise.py" 10 "$INIT10" 6 $S30   &  # ~8900 min (6.2 d)
run_case c1       "$DC" python3 "$CR/run_capillary_rise.py" 10 "$INIT10" 1 $S45   &  # ~7350 min (5.1 d)
run_case c5       "$DC" python3 "$CR/run_capillary_rise.py" 10 "$INIT10" 5 $S45   &  # ~7350 min (5.1 d)
run_case c2_hires "$DC" python3 "$CR/run_capillary_rise.py" 15 "$INIT15" 2 $SHIRES &  # ~9400 min (6.5 d; 1*tau -> exp. fit)

# Wave B: remaining jobs start as Wave-A slots free up (simple: wait for the
# two 5.1 d jobs, then launch; kernel scheduler balances the overlap)
wait  # all of wave A (c2/c4 then run on a lightly loaded machine, faster than est.)

run_case c2 "$DC" python3 "$CR/run_capillary_rise.py" 10 "$INIT10" 2 $S60 &  # ~5200 min (3.6 d, with c4/c3/b1/b2)
run_case c4 "$DC" python3 "$CR/run_capillary_rise.py" 10 "$INIT10" 4 $S60 &  # ~5200 min (3.6 d)
run_case c3 "$DC" python3 "$CR/run_capillary_rise.py" 10 "$INIT10" 3 $S90 &  # ~1600 min (1.1 d)
(
  run_case create_b1 "$RUNDIR/b1" python3 "$B6/create_input_geometry.py" 20
  INIT6=$(ls -t "$RUNDIR"/b1/sessile_shear_*_init.gsd | head -1)
  run_case run_b1 "$RUNDIR/b1" python3 "$B6/run_sessile_droplet_shear.py" 20 "$INIT6"   # ~400 min
) &
(
  run_case create_b2 "$RUNDIR/b2" python3 "$B7/create_input_geometry.py" 20
  INIT7=$(ls -t "$RUNDIR"/b2/*_init.gsd | head -1)
  run_case run_b2 "$RUNDIR/b2" python3 "$B7/run_sessile_droplet_snapback.py" 20 "$INIT7"   # ~400 min
) &
wait

stamp "════ wetting requeue finished ════"
