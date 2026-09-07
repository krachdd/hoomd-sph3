#!/usr/bin/env python3
"""
Copyright (c) 2025-2026 David Krach, Daniel Rostan.
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its contributors
   may be used to endorse or promote products derived from this software without
   specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.

maintainer: dkrach, david.krach@mib.uni-stuttgart.de

Viscoplastic invasion benchmark (05) — post-processing.

Compares a Newtonian-defender and a Bingham-defender run:
  Panel A : invasion front position x_f(t) (99th percentile of invader x)
  Panels B, C : final mid-plane slice (|z| < 1.5 dx) phase maps; defender
                particles slower than 5% of the invader velocity scale are
                marked as stagnant.

Usage:
    python3 postprocess.py [base] [tauy_newt] [tauy_bing]
      base       : GSD basename prefix (default "packing")
      tauy_newt  : tau_y label of the baseline run   (default 0)
      tauy_bing  : tau_y label of the Bingham run    (default 15)

Expects <base>_tauy<X>_run.gsd and <base>_tauy<X>_dt.txt written by
run_invasion.py. Output: invasion_comparison.png.
"""

import sys
import numpy as np
import gsd.hoomd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

base = sys.argv[1] if len(sys.argv) > 1 else "packing"
ty0  = sys.argv[2] if len(sys.argv) > 2 else "0"
ty1  = sys.argv[3] if len(sys.argv) > 3 else "15"

DX = 2.0e-4
U_SCALE = 0.08
V_STAG = 0.05 * U_SCALE

INK = "#182126"
INK_SOFT = "#46555c"
GRID = "#d8e0e2"
C_NEWT = "#0b7e9d"     # teal   — Newtonian baseline / invader
C_BING = "#c2571d"     # orange — Bingham defender case
C_SOLID = "#c7ced1"    # neutral gray — solid grains (recessive background)
C_DEF = "#d99a5b"      # warm  — mobile defender
C_STAG = "#8a4315"     # dark warm — stagnant (unyielded) defender

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 9,
    "text.color": INK,
    "axes.edgecolor": GRID,
    "axes.labelcolor": INK_SOFT,
    "xtick.color": INK_SOFT,
    "ytick.color": INK_SOFT,
})


def front_series(fname):
    ts, fronts = [], []
    with gsd.hoomd.open(fname, "r") as traj:
        for s in traj:
            A = s.particles.typeid == 0
            ts.append(s.configuration.step)
            fronts.append(np.percentile(s.particles.position[A, 0], 99))
    return np.array(ts, dtype=float), np.array(fronts)


def final_slice(fname):
    with gsd.hoomd.open(fname, "r") as traj:
        s = traj[-1]
    sl = np.abs(s.particles.position[:, 2]) < 1.5 * DX
    return (s.particles.position[sl], s.particles.typeid[sl],
            s.particles.velocity[sl])


f0name = f"{base}_tauy{ty0}_run.gsd"
f1name = f"{base}_tauy{ty1}_run.gsd"
t0_steps, f0 = front_series(f0name)
t1_steps, f1 = front_series(f1name)
DT0 = float(np.loadtxt(f"{base}_tauy{ty0}_dt.txt"))
DT1 = float(np.loadtxt(f"{base}_tauy{ty1}_dt.txt"))
t0 = t0_steps * DT0 * 1e3   # ms
t1 = t1_steps * DT1 * 1e3

with gsd.hoomd.open(f0name, "r") as traj:
    L = traj[0].configuration.box[0]

p0, id0, v0 = final_slice(f0name)
p1, id1, v1 = final_slice(f1name)

# ─── figure ──────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(9.2, 3.4), dpi=170, facecolor="white")
gs = fig.add_gridspec(1, 3, width_ratios=[1.35, 1, 1], wspace=0.28,
                      left=0.07, right=0.985, top=0.86, bottom=0.16)

ax = fig.add_subplot(gs[0])
ax.set_facecolor("white")
for spine in ("top", "right"):
    ax.spines[spine].set_visible(False)
ax.grid(axis="y", color=GRID, linewidth=0.7)
ax.set_axisbelow(True)
ax.plot(t0, (f0 + L/2) * 1e3, color=C_NEWT, lw=2, solid_capstyle="round")
ax.plot(t1, (f1 + L/2) * 1e3, color=C_BING, lw=2, solid_capstyle="round")
ax.scatter([t0[-1]], [(f0[-1] + L/2) * 1e3], s=18, color=C_NEWT, zorder=3)
ax.scatter([t1[-1]], [(f1[-1] + L/2) * 1e3], s=18, color=C_BING, zorder=3)
ax.annotate("Newtonian defender", xy=(t0[-1], (f0[-1] + L/2) * 1e3),
            xytext=(-10, -4), textcoords="offset points", ha="right", va="top",
            color=INK, fontsize=8.5)
ax.annotate(f"Bingham defender\n(τ$_y$ = {ty1} Pa)",
            xy=(t1[-1], (f1[-1] + L/2) * 1e3),
            xytext=(-8, -22), textcoords="offset points", ha="right",
            color=INK, fontsize=8.5)
ax.set_xlabel("time  [ms]")
ax.set_ylabel("invasion front  [mm]")
ax.set_title("A   Invasion front position", loc="left", fontsize=10,
             color=INK, fontweight="bold")

for k, (p, tid, v, ttl) in enumerate((
        (p0, id0, v0, "B   Newtonian defender — final"),
        (p1, id1, v1, "C   Bingham defender — final"))):
    axs = fig.add_subplot(gs[k + 1])
    axs.set_facecolor("white")
    axs.set_aspect("equal")
    for spine in axs.spines.values():
        spine.set_visible(False)
    axs.set_xticks([]); axs.set_yticks([])
    speed = np.linalg.norm(v, axis=1)
    mm = 1e3
    sel = {"sol": tid == 2,
           "dfm": (tid == 1) & (speed >= V_STAG),
           "dfs": (tid == 1) & (speed < V_STAG),
           "inv": tid == 0}
    for key, col in (("sol", C_SOLID), ("dfm", C_DEF), ("dfs", C_STAG), ("inv", C_NEWT)):
        m = sel[key]
        axs.scatter(p[m, 0]*mm, p[m, 1]*mm, s=4.5, c=col, marker="s", lw=0)
    axs.set_title(ttl, loc="left", fontsize=10, color=INK, fontweight="bold")
    axs.set_xlim(-L/2*mm, L/2*mm); axs.set_ylim(-L/2*mm, L/2*mm)

handles = [
    plt.Line2D([], [], marker="s", ls="", ms=7, color=C_NEWT, label="invader (Newtonian)"),
    plt.Line2D([], [], marker="s", ls="", ms=7, color=C_DEF, label="defender, flowing"),
    plt.Line2D([], [], marker="s", ls="", ms=7, color=C_STAG, label="defender, stagnant (|v| < 0.05 U)"),
    plt.Line2D([], [], marker="s", ls="", ms=7, color=C_SOLID, label="solid grains"),
]
fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False,
           fontsize=8, bbox_to_anchor=(0.62, -0.015))

fig.savefig("invasion_comparison.png", dpi=170, facecolor="white",
            bbox_inches="tight", pad_inches=0.15)
print("saved invasion_comparison.png")

for name, tid_, v_, f_ in (("Newtonian", id0, v0, f0), ("Bingham", id1, v1, f1)):
    speed = np.linalg.norm(v_, axis=1)
    dfs = (tid_ == 1) & (speed < V_STAG)
    frac = dfs.sum() / max((tid_ == 1).sum(), 1)
    print(f"{name}: front advance = {(f_[-1]-f_[0])*1e3:.2f} mm, "
          f"stagnant defender fraction (slice) = {frac*100:.0f}%")
