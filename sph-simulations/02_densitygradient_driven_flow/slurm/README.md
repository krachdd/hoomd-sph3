# SLURM scripts — dd_pm production campaign (wolfgang)

Production runs for the pore-scale density-driven fingering benchmarks
(cases `04_porous_fingering_single_mode` and `05_porous_fingering_rayleigh`).
All scripts follow the sph3-conda-env launch convention: `module purge`,
`PYTHONPATH` to the repo builds, and the **env's own `mpirun`** (a mismatched
launcher fractures the MPI world into singletons — the "GSD: Not a GSD file"
race).

Before first use, adjust in every script if your layout differs:
`REPO=/data/work/ac126015/imb_test_conda/hoomd-sph3` and the
`MPIRUN`/`PY` paths into the sph3 conda env.

## Scripts & dependency chain

| # | Script | What | Size / steps | Layout | Wall time / partition |
|---|--------|------|--------------|--------|-----------------------|
| 01 | `01_stageA_perm.sh` | Stage A permeability cells, n_d = 10 **and** 20 (creates geometry if missing) | 27 k / 212 k particles, 20 k steps each | 4×8 | 4 h, `cpu` |
| 02 | `02_stageA_finger_nd10.sh` | Stage A unstable (165 k steps, σ·t ≈ 3) **+** stable control (50 k) | 354 k particles | 8×8 | 24 h, `cpu` |
| 03 | `03_stageA_finger_nd20.sh` | Stage A convergence member, n_d = 20 | 2.83 M particles, 330 k steps | **2 nodes**, 32×8 | 2 d, `cpu` |
| 04 | `04_stageB_perm.sh` | Stage B permeability cells (seed 42), n_d = 10 and 20 | 170 k / 1.36 M particles, 20 k steps | 4×8 | 6 h, `cpu` |
| 05 | `05_stageB_ra_sweep.sh` | **Array job** 0–2: Ra = 200, Ra = 1000, stable control | 1.05 M particles, 100 k / 60 k steps | 16×8 each | 2 d, `cpu` |
| 06 | `06_stageB_ra5000_nd20.sh` | Ra = 5000 at n_d = 20 (required by Pe_grid = Ra·φ·dx/H) | ~8.4 M particles, 200 k steps | **3 nodes**, 48×8 | 2 d, `cpu` |

**Partition limits (wolfgang):** `cpu` caps at 2-00:00:00 (12 nodes);
`cpu-long` (5 d, wolfgang-cpu[09-12] only) is currently **drained**, so the
long runs compress into the 2-day cap by scaling out instead: job 03 uses
two nodes (~1.5 d estimated), job 06 three nodes (~1.3 d, ~175 k
particles/rank). Both headers carry the fallbacks (fewer nodes / shorter
σ·t horizon / cpu-long variants for when it returns). **Shake down the
multi-node MPI first** — run job 06 once with `STEPS=2000` before
committing the full horizon; per-rank subdomains stay far above the kernel
cutoff in every configuration, so the decomposition itself is safe. Job 05
requests the full 2-day cap — trim `--time` once real TPS numbers are in so
the scheduler can backfill it sooner.

The fingering jobs read the **measured permeability automatically** from the
`*_permeability_summary.dat` written by the corresponding perm job (last
line, column 2; porosity from column 4) and compute the diffusivity for the
target Ra from it — no manual K hand-off needed.

## Submission

```bash
cd <this directory>
j1=$(sbatch --parsable 01_stageA_perm.sh)
j4=$(sbatch --parsable 04_stageB_perm.sh)
sbatch --dependency=afterok:$j1 02_stageA_finger_nd10.sh
sbatch --dependency=afterok:$j1 03_stageA_finger_nd20.sh
sbatch --dependency=afterok:$j4 05_stageB_ra_sweep.sh
sbatch --dependency=afterok:$j4 06_stageB_ra5000_nd20.sh
```

## Notes

- Wall times are sized from the measured local throughput
  (~0.8 M particle·steps/s on 8 cores) scaled to 64–128 cores with margin;
  trim them once the first cluster TPS numbers are in.
- For the mixing-width scaling W(t) (paper TODO), extend the Ra-sweep
  members to 237 000 steps (~0.3 t_c) and raise `--time`.
- Outputs land next to the geometries in the case directories: `*_run.gsd`
  (dump), `*_run.log` (time series), `*.npz` (Stage B post-processing),
  inline analysis in the job log. Expect ~2 GB (Stage A) to ~6 GB (Stage B)
  per run, ~45 GB for run 06.
- Job 06 has no restart mechanism yet — if it exceeds the partition limit,
  split into legs by re-initialising from the last dump frame.
