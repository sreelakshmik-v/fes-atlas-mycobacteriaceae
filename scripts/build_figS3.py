#!/usr/bin/env python3
"""Supplementary Figure S3 - is the Fe-S complement just a function of genome size?

A  Fe-S count vs proteome size, OLS fit with a 95 % prediction band, outliers named
B  Studentised residuals against fitted value, |t| > 3 flagged - the formal outlier test

Needs only fes_vs_proteome.csv (genome, fes_count, proteome_size, fes_pct,
resid, student, fitted). Species labels are optional: supply genome_labels.csv
with columns accession,species and they are used for the outlier annotation.
"""
import os
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

INK, MUTE, GRID = "#1a1a1a", "#5c5c5c", "#ebebeb"
BLUE, RED, BAND = "#2F5C8A", "#B5533C", "#AFC6D8"
plt.rcParams.update({"font.family": "sans-serif",
                     "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
                     "pdf.fonttype": 42, "ps.fonttype": 42})

d = pd.read_csv("fes_vs_proteome.csv")
if "genome_labels.csv" in os.listdir("."):
    lab = pd.read_csv("genome_labels.csv")
    d = d.merge(lab, left_on="genome", right_on=lab.columns[0], how="left")
    d["label"] = d[lab.columns[1]].fillna(d.genome)
else:
    d["label"] = d.genome

x, y = d.proteome_size.values, d.fes_count.values
sl, ic, r, p, se = stats.linregress(x, y)
xs = np.linspace(x.min() * .97, x.max() * 1.03, 200)
yhat = ic + sl * xs
n = len(x)
resid_sd = np.sqrt(((y - (ic + sl * x)) ** 2).sum() / (n - 2))
se_pred = resid_sd * np.sqrt(1 + 1 / n + (xs - x.mean()) ** 2 /
                             ((x - x.mean()) ** 2).sum())
tcrit = stats.t.ppf(.975, n - 2)

out = d[d.student.abs() > 3].copy()

fig, (ax, bx) = plt.subplots(1, 2, figsize=(7.4, 3.5), dpi=200,
                             gridspec_kw=dict(width_ratios=[1.12, 1], wspace=.30))

# ───────────────────────────── A ────────────────────────────────────────────
ax.fill_between(xs, yhat - tcrit * se_pred, yhat + tcrit * se_pred,
                color=BAND, alpha=.38, lw=0, zorder=1,
                label="95 % prediction band")
ax.plot(xs, yhat, color=BLUE, lw=1.5, zorder=3)
keep = d.student.abs() <= 3
ax.scatter(d.loc[keep, "proteome_size"], d.loc[keep, "fes_count"], s=11,
           color=BLUE, alpha=.42, edgecolors="none", zorder=2)
ax.scatter(out.proteome_size, out.fes_count, s=26, color=RED,
           edgecolors="white", linewidths=.6, zorder=4)
OFF = [(10, -10, "left"), (10, 9, "left"), (-10, 9, "right"),
       (-10, -11, "right"), (12, 0, "left"), (-12, 0, "right")]
for k, (_, rw) in enumerate(out.sort_values("proteome_size").iterrows()):
    dx_, dy_, ha_ = OFF[k % len(OFF)]
    ax.annotate(str(rw.label).replace("_", " "),
                (rw.proteome_size, rw.fes_count),
                textcoords="offset points", xytext=(dx_, dy_), ha=ha_,
                va="center", fontsize=5.4, color=RED,
                arrowprops=dict(arrowstyle="-", color="#d8b5ac", lw=.4,
                                shrinkA=0, shrinkB=3))
ax.set_xlabel("Proteins in the genome", fontsize=8.5)
ax.set_ylabel("Fe–S proteins", fontsize=8.5)
ax.text(.035, .965,
        f"$R^2$ = {r**2:.2f}   slope = {sl*1000:.1f} per 1,000 proteins\n"
        f"n = {n} genomes,  p = {p:.1e}",
        transform=ax.transAxes, va="top", fontsize=6.2, color=MUTE,
        linespacing=1.4)
ax.legend(fontsize=6.0, frameon=False, loc="lower right", handlelength=1.4)
ax.set_title("A", loc="left", fontsize=11, fontweight="bold", x=-.135, pad=5)

# ───────────────────────────── B ────────────────────────────────────────────
bx.axhspan(-3, 3, color="#9AA7AD", alpha=.07, lw=0, zorder=0)
for v in (-3, 3):
    bx.axhline(v, color=RED, lw=.8, ls=(0, (4, 3)), zorder=2)
bx.axhline(0, color=MUTE, lw=.7, zorder=2)
bx.scatter(d.loc[keep, "fitted"], d.loc[keep, "student"], s=11, color=BLUE,
           alpha=.42, edgecolors="none", zorder=3)
bx.scatter(out.fitted, out.student, s=26, color=RED, edgecolors="white",
           linewidths=.6, zorder=4)
for k, (_, rw) in enumerate(out.sort_values("fitted").iterrows()):
    dx_, dy_, ha_ = OFF[k % len(OFF)]
    bx.annotate(str(rw.label).replace("_", " "), (rw.fitted, rw.student),
                textcoords="offset points", xytext=(dx_, dy_), ha=ha_,
                va="center", fontsize=5.4, color=RED,
                arrowprops=dict(arrowstyle="-", color="#d8b5ac", lw=.4,
                                shrinkA=0, shrinkB=3))
bx.text(.985, .035, "|t| > 3", transform=bx.transAxes, ha="right",
        fontsize=5.8, color=RED, style="italic")
bx.set_xlabel("Fe–S proteins expected from genome size", fontsize=8.5)
bx.set_ylabel("Studentised residual", fontsize=8.5)
bx.set_title("B", loc="left", fontsize=11, fontweight="bold", x=-.155, pad=5)

for a in (ax, bx):
    for s in ("top", "right"):
        a.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        a.spines[s].set_linewidth(.7)
        a.spines[s].set_color("#4d4d4d")
    a.tick_params(labelsize=7.5, colors=MUTE, width=.7)
    a.grid(True, lw=.5, color=GRID, zorder=0)
    a.set_axisbelow(True)

fig.savefig("FigureS3_genome_size.png", dpi=400, bbox_inches="tight")
fig.savefig("FigureS3_genome_size.pdf", bbox_inches="tight")
fig.savefig("FigureS3_genome_size.tiff", dpi=600, bbox_inches="tight",
            pil_kwargs={"compression": "tiff_lzw"})
out.sort_values("student").to_csv("figS3_outliers.csv", index=False)

print(f"R2 = {r**2:.3f}, slope = {sl*1000:.2f} Fe-S per 1,000 proteins, p = {p:.3e}")
print(f"{len(out)} genomes with |studentised residual| > 3:")
for _, rw in out.sort_values("student").iterrows():
    print(f"  {rw.label:22s} {int(rw.fes_count):4d} observed vs "
          f"{rw.fitted:6.1f} expected   t = {rw.student:+.2f}")
print("\nwrote FigureS3_genome_size.png/.pdf/.tiff, figS3_outliers.csv")
