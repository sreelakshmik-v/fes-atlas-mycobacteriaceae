#!/usr/bin/env python3
"""Figure 4 - structural diversity of the Fe-S proteome.

A  TM-score matrix, complete-linkage clustered, dendrogram on top, fold groups boxed
B  pLDDT vs pTM with marginal histograms - model reliability for all 44 families
C  Slope plot: within-fold-group vs between-group mean TM-score, Wilcoxon

Fold groups are defined by cutting the complete-linkage dendrogram at
TM-score 0.5, the conventional same-fold threshold (Xu & Zhang 2010), not at
an arbitrary k. Complete linkage means every pair inside a group clears 0.5.

Run from /media/budheswar-dehury/DATADRIVE0/Fe-S/ :
    python3 build_figure4.py                  # cut at TM 0.5
    python3 build_figure4.py --tm 0.45        # try a looser cut
"""
import argparse, os
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform
from scipy.stats import wilcoxon

S = "step14_structural"
TM, PLD, PTM = (f"{S}/tmalign_similarity_matrix.csv",
                f"{S}/model_plddt.tsv", f"{S}/model_confidence.tsv")

INK, MUTE, GRID = "#1a1a1a", "#5c5c5c", "#ebebeb"
UNIQ = "#9AA7AD"                       # families that share a fold with nobody
FG = ["#2F5C8A", "#B5533C", "#3E8E6F", "#D9983E", "#6B4E9B", "#386B7A",
      "#8C4F5E", "#B07C3F", "#527B4E", "#9B5B8A"]
plt.rcParams.update({"font.family": "sans-serif",
                     "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
                     "pdf.fonttype": 42, "ps.fonttype": 42})


def nice(f):
    return f.replace("_", " ")


def load():
    d = pd.read_csv(TM, index_col=0)
    d.index, d.columns = [str(i) for i in d.index], [str(c) for c in d.columns]
    fams = [f for f in d.index if f in set(d.columns)]
    M = d.loc[fams, fams].astype(float).values
    M = np.fmax(M, M.T)
    np.fill_diagonal(M, 1.0)
    T = pd.DataFrame(np.nan_to_num(M), index=fams, columns=fams)
    p = pd.read_csv(PLD, sep=None, engine="python")
    q = pd.read_csv(PTM, sep=None, engine="python")
    m = p.merge(q, on="family", how="outer")
    m["family"] = m.family.astype(str)
    return T, m


def build(tm_cut):
    T, qual = load()
    fams = list(T.index)
    n = len(fams)
    D = np.clip(1 - T.values, 0, None)
    np.fill_diagonal(D, 0)
    Z = linkage(squareform(D, checks=False), method="complete")
    raw = fcluster(Z, 1 - tm_cut, criterion="distance")

    # multi-member clusters keep an identity; singletons are pooled as "unique"
    sizes = pd.Series(raw).value_counts()
    multi = [c for c in sizes.index if sizes[c] > 1]
    # order the real fold groups largest first so FG1 is the biggest
    multi = sorted(multi, key=lambda c: -sizes[c])
    gid = {c: i + 1 for i, c in enumerate(multi)}
    grp_of = np.array([gid.get(c, 0) for c in raw])          # 0 == unique fold
    colour = [FG[(g - 1) % len(FG)] if g else UNIQ for g in grp_of]

    order = dendrogram(Z, no_plot=True)["leaves"]
    ofams = [fams[i] for i in order]
    ogrp = grp_of[order]

    info = pd.DataFrame({"family": fams, "fold_group": grp_of}).merge(
        qual, on="family", how="left")
    info["fold_group_name"] = info.fold_group.map(
        lambda g: f"FG{g}" if g else "unique fold")
    info.sort_values(["fold_group", "family"]).to_csv(
        "fig4_fold_groups.csv", index=False)

    fig = plt.figure(figsize=(7.4, 8.9), dpi=200)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.42, 1], width_ratios=[1.1, 1],
                          hspace=.26, wspace=.46,
                          left=.20, right=.95, top=.955, bottom=.065)

    # ───────────────── A: clustered TM matrix (full top row) ────────────────
    axA = fig.add_subplot(gs[0, :])
    b = axA.get_position()
    cbw = .016
    mw = b.width - cbw - .022                      # leave room for the colour bar
    axA.set_position([b.x0, b.y0, mw, b.height * .775])
    axD = fig.add_axes([b.x0, b.y0 + b.height * .815, mw, b.height * .150])
    cax = fig.add_axes([b.x0 + mw + .016, b.y0 + b.height * .14,
                        cbw, b.height * .50])
    dendrogram(Z, ax=axD, no_labels=True,
               link_color_func=lambda _: "#9f9f9f")
    axD.axhline(1 - tm_cut, color="#B5533C", lw=.8, ls=(0, (4, 3)))
    axD.text(axD.get_xlim()[1], 1 - tm_cut, f" TM {tm_cut:g}", va="center",
             ha="left", fontsize=5.4, color="#B5533C")
    axD.set_axis_off()

    R = T.values[np.ix_(order, order)]
    im = axA.imshow(R, cmap="magma_r", vmin=0, vmax=1, aspect="auto")
    start = 0
    for g in ogrp:
        pass
    i = 0
    while i < n:
        g = ogrp[i]
        j = i
        while j + 1 < n and ogrp[j + 1] == g:
            j += 1
        if g:                                  # box real fold groups only
            w = j - i + 1
            axA.add_patch(Rectangle((i - .5, i - .5), w, w, fill=False,
                                    edgecolor=FG[(g - 1) % len(FG)], lw=1.6,
                                    zorder=5))
            axA.text(i + w / 2 - .5, -0.95, f"FG{g}", ha="center", va="bottom",
                     fontsize=6.2, fontweight="bold",
                     color=FG[(g - 1) % len(FG)])
        i = j + 1

    axA.set_xticks(range(n))
    axA.set_xticklabels([nice(f) for f in ofams], rotation=90, fontsize=4.3)
    axA.set_yticks(range(n))
    axA.set_yticklabels([nice(f) for f in ofams], fontsize=4.3)
    for t, f in zip(axA.get_xticklabels(), ofams):
        t.set_color(colour[fams.index(f)] if grp_of[fams.index(f)] else INK)
    for t, f in zip(axA.get_yticklabels(), ofams):
        t.set_color(colour[fams.index(f)] if grp_of[fams.index(f)] else INK)
    axA.tick_params(length=0)
    for s in axA.spines.values():
        s.set_visible(False)
    cb = fig.colorbar(im, cax=cax)
    cb.set_label("TM-score", fontsize=6.6, color=INK)
    cb.ax.tick_params(labelsize=5.8, colors=MUTE, width=.6, length=2)
    cb.outline.set_visible(False)
    axD.set_title("A", loc="left", fontsize=11, fontweight="bold",
                  x=-.185, pad=1)

    # ───────────────── B: pLDDT vs pTM, with marginals ──────────────────────
    sub = fig.add_subplot(gs[1, 0]).get_position()
    fig.axes[-1].remove()
    pad, mh = .012, .080
    axB = fig.add_axes([sub.x0, sub.y0, sub.width - mh - pad,
                        sub.height - mh - pad])
    axBt = fig.add_axes([sub.x0, sub.y0 + sub.height - mh,
                         sub.width - mh - pad, mh], sharex=axB)
    axBr = fig.add_axes([sub.x0 + sub.width - mh, sub.y0, mh,
                         sub.height - mh - pad], sharey=axB)

    q = info.dropna(subset=["mean_plddt", "ptm"])
    axB.scatter(q.mean_plddt, q.ptm, s=26,
                c=[FG[(g - 1) % len(FG)] if g else UNIQ for g in q.fold_group],
                alpha=.82, edgecolors="white", linewidths=.6, zorder=3)
    axB.axvline(90, color=MUTE, lw=.7, ls=(0, (3, 3)), zorder=2)
    axB.axhline(.80, color=MUTE, lw=.7, ls=(0, (3, 3)), zorder=2)
    axB.text(axB.get_xlim()[1], .80, " pTM 0.80", va="bottom", ha="right",
             fontsize=5.2, color=MUTE)
    # name the weakest models - a reviewer will ask about exactly these
    weak = q.nsmallest(3, "ptm").sort_values("mean_plddt")
    for k, (_, r) in enumerate(weak.iterrows()):
        _o = [(11, 8, "left"), (-9, 10, "right"), (11, 7, "left")][k % 3]
        dxy, ha_ = _o[:2], _o[2]
        axB.annotate(nice(r.family), (r.mean_plddt, r.ptm),
                     textcoords="offset points", xytext=dxy,
                     ha=ha_, va="center", fontsize=5.2, color=MUTE,
                     arrowprops=dict(arrowstyle="-", color="#c8c8c8", lw=.4,
                                     shrinkA=0, shrinkB=3))
    axB.set_xlabel("Mean per-model pLDDT", fontsize=8.3)
    axB.set_ylabel("pTM", fontsize=8.3)
    axBt.hist(q.mean_plddt, bins=14, color="#2F5C8A", alpha=.55, lw=0)
    axBr.hist(q.ptm, bins=14, orientation="horizontal", color="#2F5C8A",
              alpha=.55, lw=0)
    for a in (axBt, axBr):
        a.set_axis_off()
    axBt.set_title("B", loc="left", fontsize=11, fontweight="bold",
                   x=-.165, pad=3)
    axB.text(.03, .96, f"all {len(q)} families pLDDT $\\geq$ "
             f"{q.mean_plddt.min():.0f}", transform=axB.transAxes,
             fontsize=5.6, color=MUTE, va="top")

    # ───────────────── C: within vs between fold group ──────────────────────
    axC = fig.add_subplot(gs[1, 1])
    rows = []
    for i, f in enumerate(fams):
        if grp_of[i] == 0:
            continue                               # no within-group partner
        same = (grp_of == grp_of[i]) & (np.arange(n) != i)
        rows.append(dict(family=f, fold_group=int(grp_of[i]),
                         tm_within=float(T.values[i][same].mean()),
                         tm_between=float(T.values[i][grp_of != grp_of[i]].mean())))
    sl = pd.DataFrame(rows)
    sl.to_csv("fig4_within_between_tm.csv", index=False)

    if len(sl):
        for _, r in sl.iterrows():
            col = FG[(int(r.fold_group) - 1) % len(FG)]
            axC.plot([0, 1], [r.tm_between, r.tm_within], color=col, lw=.9,
                     alpha=.6, zorder=2, solid_capstyle="round")
            axC.scatter([0, 1], [r.tm_between, r.tm_within], s=15, color=col,
                        edgecolors="white", linewidths=.5, zorder=3)
        axC.plot([0, 1], [sl.tm_between.mean(), sl.tm_within.mean()], color=INK,
                 lw=2.1, zorder=5)
        axC.scatter([0, 1], [sl.tm_between.mean(), sl.tm_within.mean()], s=38,
                    color=INK, edgecolors="white", linewidths=.9, zorder=6)
        try:
            ptxt = ("p < 0.001" if wilcoxon(sl.tm_within, sl.tm_between)[1] < 1e-3
                    else f"p = {wilcoxon(sl.tm_within, sl.tm_between)[1]:.3g}")
        except Exception:
            ptxt = "n too small"
        axC.text(.5, sl[["tm_within", "tm_between"]].values.max() + .045,
                 f"Wilcoxon signed-rank, {ptxt}\n"
                 f"{len(sl)} families in {sl.fold_group.nunique()} fold groups",
                 ha="center", va="bottom", fontsize=5.8, color=MUTE,
                 linespacing=1.35)
        axC.set_ylim(0, sl[["tm_within", "tm_between"]].values.max() + .19)
    axC.set_xticks([0, 1])
    axC.set_xticklabels(["between\nfold groups", "within\nfold group"],
                        fontsize=7.0, color=INK)
    axC.set_xlim(-.32, 1.32)
    axC.set_ylabel("Mean TM-score to other families", fontsize=8.3)
    axC.set_title("C", loc="left", fontsize=11, fontweight="bold",
                  x=-.28, pad=5)
    axC.legend(handles=[Line2D([], [], marker="o", ls="", ms=5,
                               mfc=FG[(g - 1) % len(FG)], mec="white",
                               label=f"FG{g} (n={int((grp_of == g).sum())})")
                        for g in sorted(set(grp_of) - {0})],
               fontsize=5.8, frameon=False, loc="lower right",
               handletextpad=.3, labelspacing=.35,
               bbox_to_anchor=(1.03, -.02))

    for a in (axB, axC):
        for s in ("top", "right"):
            a.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            a.spines[s].set_linewidth(.7)
            a.spines[s].set_color("#4d4d4d")
        a.tick_params(labelsize=7.4, colors=MUTE, width=.7)
        a.grid(True, lw=.5, color=GRID, zorder=0)
        a.set_axisbelow(True)
    axC.xaxis.grid(False)
    axC.tick_params(axis="x", length=0, labelcolor=INK)

    fig.savefig("Figure4.png", dpi=400, bbox_inches="tight")
    fig.savefig("Figure4.pdf", bbox_inches="tight")

    print(f"cut at TM {tm_cut:g}, complete linkage")
    print(f"{len(multi)} multi-family fold groups, "
          f"{int((grp_of == 0).sum())} structurally unique families\n")
    for g in sorted(set(grp_of) - {0}):
        m = info[info.fold_group == g].sort_values("family")
        print(f"  FG{g}  n={len(m)}  pLDDT {m.mean_plddt.min():.1f}-"
              f"{m.mean_plddt.max():.1f}   {', '.join(m.family)}")
    if len(sl):
        print(f"\nwithin {sl.tm_within.mean():.3f} vs between "
              f"{sl.tm_between.mean():.3f}   {ptxt}")
    print(f"\npLDDT {info.mean_plddt.min():.1f}-{info.mean_plddt.max():.1f}, "
          f"median {info.mean_plddt.median():.1f};  "
          f"pTM median {info.ptm.median():.2f}, "
          f"{int((info.ptm >= .8).sum())}/{len(info)} at pTM >= 0.80")
    print("\nwrote Figure4.png/.pdf, fig4_fold_groups.csv, "
          "fig4_within_between_tm.csv")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--tm", type=float, default=0.50)
    build(ap.parse_args().tm)
