#!/usr/bin/env python3
"""Supplementary Figure S5 - per-column conservation for the four families
that have a profile (IspG, IspH_uniq, NadA, ThiC).

One stacked panel per family: conservation along the alignment as a filled
profile, gap fraction as a grey band underneath, and the most conserved
columns marked. Ligand-binding cysteines are highlighted when
figure_data/panelE_coordination_cys.csv names them.

Run from /media/budheswar-dehury/DATADRIVE0/Fe-S/ :
    python3 build_figS5.py
"""
import glob, os, re
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PHY = "step13_phylo"
INK, MUTE, GRID = "#1a1a1a", "#5c5c5c", "#ebebeb"
BLUE, GREY, CYS = "#2F5C8A", "#C9D2D8", "#B5533C"
plt.rcParams.update({"font.family": "sans-serif",
                     "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
                     "pdf.fonttype": 42, "ps.fonttype": 42})


def read_any(p):
    return pd.read_csv(p, sep=None, engine="python")


files = sorted(glob.glob(f"{PHY}/*_conservation.tsv"))
if not files:
    raise SystemExit(f"no *_conservation.tsv under {PHY}/")
print("profiles found:", [os.path.basename(f) for f in files])

fig, axes = plt.subplots(len(files), 1, figsize=(7.2, 1.62 * len(files)),
                         dpi=200, sharex=False)
axes = np.atleast_1d(axes)
summary = []

for ax, f in zip(axes, files):
    fam = re.sub(r"_conservation\.tsv$", "", os.path.basename(f))
    d = read_any(f)
    d.columns = [c.strip().lower() for c in d.columns]
    col = d["col"].values
    cons = pd.to_numeric(d["conservation"], errors="coerce").values
    gap = 1 - pd.to_numeric(d["pct_gapless"], errors="coerce").values / 100
    ok = ~np.isnan(cons)

    # gap fraction as a grey band along the bottom, on the same 0-1 scale
    ax.fill_between(col, 0, gap, color=GREY, alpha=.85, lw=0, zorder=1)
    ax.fill_between(col[ok], 0, cons[ok], color=BLUE, alpha=.42, lw=0, zorder=2)
    ax.plot(col[ok], cons[ok], color=BLUE, lw=.65, zorder=3)

    # conservation is near-saturated in these alignments, so the informative
    # columns are the VARIABLE ones, not the conserved majority
    invariant = ok & (cons >= 0.999)
    thr_var = np.nanpercentile(cons, 5)
    var = ok & (cons <= thr_var)
    ax.scatter(col[var], cons[var], s=5.0, color="#B5533C", zorder=5, lw=0)

    # conserved ligand cysteines, the residues the Fe-S cluster hangs off
    ncys = 0
    if "top_aa" in d.columns:
        iscys = (ok & (cons >= .95) &
                 (d["top_aa"].astype(str).str.upper().str[0] == "C").values)
        ncys = int(iscys.sum())
        if ncys:
            ax.scatter(col[iscys], np.full(ncys, 1.08), s=16, marker="v",
                       color="#D9983E", edgecolors="white", linewidths=.4,
                       zorder=6)

    ax.set_xlim(col.min(), col.max())
    ax.set_ylim(0, 1.17)
    ax.set_yticks([0, .5, 1])
    ax.set_ylabel("conservation", fontsize=6.8)
    ax.text(.004, 1.02, fam.replace("_", " "), transform=ax.transAxes,
            fontsize=8.0, fontweight="bold", color=INK, va="bottom")
    ax.text(.998, 1.02,
            f"{int(ok.sum())} columns · mean {np.nanmean(cons):.2f} · "
            f"{100*invariant.sum()/max(ok.sum(),1):.0f} % invariant · "
            f"{int(var.sum())} variable" + (f" · {ncys} Cys" if ncys else ""),
            transform=ax.transAxes, ha="right", va="bottom",
            fontsize=5.6, color=MUTE)
    summary.append(dict(family=fam, columns=int(ok.sum()),
                        mean_conservation=round(float(np.nanmean(cons)), 3),
                        pct_invariant=round(100*float(invariant.sum())/max(int(ok.sum()),1), 1),
                        median_gap=round(float(np.nanmedian(gap)), 3),
                        n_variable=int(var.sum()), n_conserved_cys=ncys))

axes[-1].set_xlabel("Alignment column", fontsize=8.5)
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
axes[-1].legend(handles=[Patch(facecolor=BLUE, alpha=.42, label="conservation"),
                        Patch(facecolor=GREY, label="gap fraction"),
                        Line2D([], [], marker="o", ls="", mfc="#B5533C",
                               mec="#B5533C", ms=3.5,
                               label="most variable 5 % of columns"),
                        Line2D([], [], marker="v", ls="", mfc="#D9983E",
                               mec="white", ms=5, label="conserved Cys")],
                fontsize=6.0, frameon=False, ncol=4, loc="upper center",
                bbox_to_anchor=(.5, -.42), handlelength=1.2,
                handletextpad=.5, columnspacing=1.8)

for a in axes:
    for s in ("top", "right"):
        a.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        a.spines[s].set_linewidth(.7)
        a.spines[s].set_color("#4d4d4d")
    a.tick_params(labelsize=6.8, colors=MUTE, width=.7)
    a.yaxis.grid(True, lw=.5, color=GRID, zorder=0)
    a.set_axisbelow(True)

fig.tight_layout(h_pad=1.25)
fig.savefig("FigureS5_conservation.png", dpi=400, bbox_inches="tight")
fig.savefig("FigureS5_conservation.pdf", bbox_inches="tight")
fig.savefig("FigureS5_conservation.tiff", dpi=600, bbox_inches="tight",
            pil_kwargs={"compression": "tiff_lzw"})
s = pd.DataFrame(summary)
s.to_csv("figS5_conservation_summary.csv", index=False)
print(s.to_string(index=False))
print("\nwrote FigureS5_conservation.png/.pdf/.tiff, "
      "figS5_conservation_summary.csv")
