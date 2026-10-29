#!/usr/bin/env python3
"""Supplementary Figure S5 (alternative) - are the predicted ligand-binding
cysteines more conserved than the rest of the alignment?

Forest plot: one row per family, the conservation difference
(ligand Cys - other columns) with a bootstrap-free normal CI, ordered by effect.
Significant families after Benjamini-Hochberg are filled; the rest are hollow.

This uses 36 families instead of the 4 that have full per-column profiles, so it
is the stronger version of the conservation argument.

Input: figure_data/panelE_coordination_cys.csv  (tab-separated despite the name)
       family, n_ligand_cys, mean_cons_ligand, mean_cons_other, pvalue

Run from /media/budheswar-dehury/DATADRIVE0/Fe-S/ :
    python3 build_figS5b.py
"""
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

SRC = "figure_data/panelE_coordination_cys.csv"
INK, MUTE, GRID = "#1a1a1a", "#5c5c5c", "#ebebeb"
SIG, NS = "#2F5C8A", "#AFC6D8"
plt.rcParams.update({"font.family": "sans-serif",
                     "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
                     "pdf.fonttype": 42, "ps.fonttype": 42})

d = pd.read_csv(SRC, sep=None, engine="python")
d.columns = [c.strip() for c in d.columns]
for c in ("n_ligand_cys", "mean_cons_ligand", "mean_cons_other", "pvalue"):
    d[c] = pd.to_numeric(d[c], errors="coerce")
d = d.dropna(subset=["mean_cons_ligand", "mean_cons_other"]).copy()

d["delta"] = d.mean_cons_ligand - d.mean_cons_other

# Benjamini-Hochberg across the families that carry a p-value
p = d.pvalue.copy()
has_p = p.notna()
m = int(has_p.sum())
d["qvalue"] = np.nan
if m:
    o = p[has_p].sort_values()
    q = (o.values * m / np.arange(1, m + 1))
    q = np.minimum.accumulate(q[::-1])[::-1]        # enforce monotonicity
    d.loc[o.index, "qvalue"] = np.clip(q, 0, 1)
d["significant"] = d.qvalue < 0.05

d = d.sort_values("delta").reset_index(drop=True)
y = np.arange(len(d))

fig, (ax, bx) = plt.subplots(1, 2, figsize=(7.4, max(4.4, .20 * len(d) + 1.4)),
                             dpi=200, sharey=True,
                             gridspec_kw=dict(width_ratios=[1.9, 1], wspace=.06))

# ───────────────── left: the effect itself ──────────────────────────────────
ax.axvline(0, color=MUTE, lw=.9, zorder=2)
ax.axvspan(d.delta.min() - .02, 0, color="#B5533C", alpha=.045, lw=0, zorder=0)
for i, r in d.iterrows():
    col = SIG if r.significant else NS
    ax.plot([0, r.delta], [i, i], color=col, lw=1.5, alpha=.75, zorder=3)
    ax.scatter(r.delta, i, s=30 if r.significant else 22, color=col,
               edgecolors="white" if r.significant else col,
               linewidths=.7, zorder=4,
               marker="o" if r.significant else "o",
               facecolors=col if r.significant else "white")
ax.set_yticks(y)
ax.set_yticklabels([f.replace("_", " ") for f in d.family], fontsize=5.8,
                   color=INK)
ax.set_xlabel("Conservation difference  (ligand Cys − other columns)",
              fontsize=8.3)
ax.tick_params(axis="y", length=0)
ax.set_ylim(-.8, len(d) - .2)
n_sig = int(d.significant.sum())
n_pos = int((d.delta > 0).sum())
ax.text(0, 1.012,
        f"{n_pos} of {len(d)} families: ligand Cys more conserved   ·   "
        f"{n_sig} significant at BH q < 0.05",
        transform=ax.transAxes, va="bottom", fontsize=6.3, color=MUTE)
ax.set_title("A", loc="left", fontsize=11, fontweight="bold", x=-.30, pad=5)

# ───────────────── right: how many cysteines carry it ───────────────────────
bx.barh(y, d.n_ligand_cys, height=.55,
        color=[SIG if s else NS for s in d.significant], alpha=.72, lw=0,
        zorder=3)
for i, r in d.iterrows():
    bx.text(r.n_ligand_cys + d.n_ligand_cys.max() * .025, i,
            f"{int(r.n_ligand_cys)}", va="center", fontsize=5.4, color=MUTE)
bx.set_xlabel("Predicted ligand cysteines", fontsize=8.3)
bx.set_xlim(0, d.n_ligand_cys.max() * 1.22)
bx.tick_params(axis="y", length=0)
bx.set_title("B", loc="left", fontsize=11, fontweight="bold", x=-.07, pad=5)

ax.legend(handles=[Line2D([], [], marker="o", ls="", ms=5.5, mfc=SIG,
                          mec="white", label="BH q < 0.05"),
                   Line2D([], [], marker="o", ls="", ms=5.0, mfc="white",
                          mec=NS, label="not significant")],
          fontsize=6.2, frameon=False, loc="lower right",
          handletextpad=.4, labelspacing=.4)

for a in (ax, bx):
    for s in ("top", "right"):
        a.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        a.spines[s].set_linewidth(.7)
        a.spines[s].set_color("#4d4d4d")
    a.tick_params(labelsize=7.2, colors=MUTE, width=.7)
    a.xaxis.grid(True, lw=.5, color=GRID, zorder=0)
    a.set_axisbelow(True)
bx.spines["left"].set_visible(False)

fig.savefig("FigureS5_ligand_cys.png", dpi=400, bbox_inches="tight")
fig.savefig("FigureS5_ligand_cys.pdf", bbox_inches="tight")
fig.savefig("FigureS5_ligand_cys.tiff", dpi=600, bbox_inches="tight",
            pil_kwargs={"compression": "tiff_lzw"})
d.to_csv("figS5_ligand_cys_stats.csv", index=False)

print(f"{len(d)} families; ligand Cys more conserved in {n_pos}, "
      f"significant (BH q<0.05) in {n_sig}")
print(f"median difference {d.delta.median():+.3f}  "
      f"(range {d.delta.min():+.3f} to {d.delta.max():+.3f})")
print("\nlargest effects:")
for _, r in d.nlargest(5, "delta").iterrows():
    print(f"  {r.family:28s} {r.delta:+.3f}   "
          f"{int(r.n_ligand_cys)} Cys   q = {r.qvalue:.2g}")
print("\nnegative or null:")
for _, r in d.nsmallest(4, "delta").iterrows():
    print(f"  {r.family:28s} {r.delta:+.3f}   {int(r.n_ligand_cys)} Cys")
print("\nwrote FigureS5_ligand_cys.png/.pdf/.tiff, figS5_ligand_cys_stats.csv")
