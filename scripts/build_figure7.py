#!/usr/bin/env python3
"""Figure 7 - the Fe-S vulnerability of Mycobacterium tuberculosis.

A  Waffle, one square per protein, grouped into labelled blocks by
   essentiality call so each category is identifiable without a legend.
B  Treemap of the families those proteins belong to, area by protein count.
C  Mosaic of the 2x2 behind the essentiality enrichment, area by count.
D  Essentiality against identity to the closest human protein; proteins with
   no human homologue sit on the zero axis in the shaded selective band.

Inputs, in order of preference:
  selectivity.csv                      locus,gene,family,Rv,DeJesus_call,
                                       pident,has_human,ESGD
  or rebuilt from
  step11_function/mtb_fes_dejesus.tsv  + mtb_vs_human.tsv

    python3 build_figure7.py
    python3 build_figure7.py --dejesus ess.tsv --human hits.tsv
"""
import argparse, os, sys
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
from scipy.stats import fisher_exact

INK, MUTE, GRID = "#1a1a1a", "#5c5c5c", "#ebebeb"

# one palette shared by panels A, C and D so the colours mean the same thing
CALL_COL = {"ES": "#C0392B", "ESD": "#E07B39", "GD": "#F0A202",
            "NE": "#2E86AB", "GA": "#7B4FA3", "Uncertain": "#8E9AA3",
            "unmapped": "#CFD6DB"}
CALL_LAB = {"ES": "essential", "ESD": "essential domain",
            "GD": "growth defect", "NE": "non-essential",
            "GA": "growth advantage", "Uncertain": "uncertain",
            "unmapped": "not mapped"}
PAL = ["#2E86AB", "#2A9D8F", "#C0392B", "#F0A202", "#7B4FA3", "#1B6B7A",
       "#B5495B", "#D08C34", "#4C9A5E", "#A8577E", "#4C7FA8", "#E07B39",
       "#5B8C5A", "#9B5DE5", "#E85D75", "#3A8EBA", "#C77DBB", "#6AA84F"]
SEL = "#2A9D8F"

plt.rcParams.update({"font.family": "sans-serif",
                     "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
                     "pdf.fonttype": 42, "ps.fonttype": 42})


# ─────────────────────────────── data ───────────────────────────────────────
def load_data(a):
    if os.path.exists(a.selectivity):
        d = pd.read_csv(a.selectivity)
    else:
        if not os.path.exists(a.dejesus):
            sys.exit(f"need {a.selectivity} or {a.dejesus}; see --help")
        d = pd.read_csv(a.dejesus, sep=None, engine="python")
        d.columns = [c.strip() for c in d.columns]
        for c in ("locus", "gene", "family", "Rv", "DeJesus_call"):
            if c not in d.columns:
                sys.exit(f"{a.dejesus} is missing column '{c}'")
        d["pident"], d["has_human"] = 0.0, False
        if os.path.exists(a.human):
            cols = ["qseqid", "sseqid", "pident", "length", "mismatch",
                    "gapopen", "qstart", "qend", "sstart", "send", "evalue",
                    "bitscore", "qcovhsp", "stitle"]
            h = pd.read_csv(a.human, sep="\t", header=None, names=cols,
                            engine="python", on_bad_lines="skip")
            for c in ("pident", "evalue", "qcovhsp", "bitscore"):
                h[c] = pd.to_numeric(h[c], errors="coerce")
            keep = h[(h.evalue <= a.evalue) & (h.pident >= a.hom_id) &
                     (h.qcovhsp >= a.hom_cov)]
            best = (keep.sort_values("bitscore", ascending=False)
                        .drop_duplicates("qseqid").set_index("qseqid"))
            d["pident"] = d.locus.map(best.pident).fillna(0.0)
            d["has_human"] = d.locus.isin(best.index)
            print(f"rebuilt from {a.dejesus} + {a.human}: "
                  f"{int(d.has_human.sum())}/{len(d)} with a human homologue "
                  f"at E<={a.evalue:g}")
        else:
            print(f"warning: {a.human} not found - no homologues called")
        d["ESGD"] = d.DeJesus_call.isin(["ES", "ESD", "GD"])
        d.to_csv(a.selectivity, index=False)
        print(f"wrote {a.selectivity}")
    d["has_human"] = d.has_human.astype(str).str.lower().isin(
        ["true", "yes", "1"])
    d["pident"] = pd.to_numeric(d.pident, errors="coerce").fillna(0.0)
    return d


def squarify(values, x, y, w, h):
    """Minimal squarified treemap; one rectangle per value."""
    vals, rects = list(values), []
    cx, cy, cw, ch = x, y, w, h
    i = 0
    while i < len(vals):
        remaining = sum(vals[i:])
        if remaining <= 0:
            break
        horizontal = cw >= ch
        side = ch if horizontal else cw
        row, best, j = [], None, i
        while j < len(vals):
            trial = row + [vals[j]]
            s = sum(trial)
            thick = (cw if horizontal else ch) * s / remaining
            worst = (max(max(thick / (side * v / s), (side * v / s) / thick)
                         for v in trial) if thick > 0 else 1e9)
            if best is None or worst <= best:
                best, row = worst, trial
                j += 1
            else:
                break
        s = sum(row)
        thick = (cw if horizontal else ch) * s / remaining
        off = 0.0
        for v in row:
            ln = side * v / s
            rects.append((cx, cy + off, thick, ln) if horizontal
                         else (cx + off, cy, ln, thick))
            off += ln
        if horizontal:
            cx += thick; cw -= thick
        else:
            cy += thick; ch -= thick
        i += len(row)
    return rects


# ─────────────────────────────── figure ─────────────────────────────────────
def main(a):
    d = load_data(a)
    n_mapped = len(d)
    n_unmapped = max(a.n_total - n_mapped, 0)

    fig = plt.figure(figsize=(7.4, 7.6), dpi=200)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.05, 1], width_ratios=[1, 1.08],
                          hspace=.30, wspace=.26,
                          left=.055, right=.965, top=.945, bottom=.075)

    # ───────────── A: essentiality breakdown, directly labelled ─────────────
    ax = fig.add_subplot(gs[0, 0])
    order = ["ES", "ESD", "GD", "NE", "GA", "Uncertain"]
    counts = [(c, int((d.DeJesus_call == c).sum())) for c in order]
    counts = [(c, n) for c, n in counts if n]
    if n_unmapped:
        counts.append(("unmapped", n_unmapped))

    labs = [CALL_LAB[c] for c, _ in counts]
    vals = [n for _, n in counts]
    cols = [CALL_COL[c] for c, _ in counts]
    y = np.arange(len(vals))[::-1]

    ax.barh(y, vals, height=.62, color=cols, edgecolor="white", linewidth=.8,
            zorder=3)
    for yy, v, c in zip(y, vals, cols):
        ax.text(v + max(vals) * .022, yy, f"{v}", va="center", ha="left",
                fontsize=8.0, color=c, fontweight="bold")
        ax.text(v + max(vals) * .105, yy, f"{100*v/a.n_total:.1f} %",
                va="center", ha="left", fontsize=6.3, color=MUTE)
    ax.set_yticks(y)
    ax.set_yticklabels(labs, fontsize=7.6)
    for t, c in zip(ax.get_yticklabels(), cols):
        t.set_color(c); t.set_fontweight("bold")
    ax.set_xlim(0, max(vals) * 1.32)
    ax.set_ylim(-.65, len(vals) - .35)
    ax.set_xlabel("Proteins in M. tuberculosis H37Rv", fontsize=8.2)
    ax.text(0, len(vals) - .18, f"{a.n_total} Fe–S proteins in H37Rv",
            fontsize=8.4, color=INK, fontweight="bold", va="bottom")
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)
    for s_ in ("left", "bottom"):
        ax.spines[s_].set_linewidth(.7); ax.spines[s_].set_color("#4d4d4d")
    ax.tick_params(labelsize=7.2, colors=MUTE, width=.7)
    ax.tick_params(axis="y", length=0)
    ax.xaxis.grid(True, lw=.5, color=GRID, zorder=0); ax.set_axisbelow(True)
    ax.set_title("A", loc="left", fontsize=11, fontweight="bold", x=-.30, pad=6)

    # ───────────── B: treemap of families ───────────────────────────────────
    bx = fig.add_subplot(gs[0, 1])
    fam = d.family.value_counts()
    for k, ((x0, y0, w, h), (name, n)) in enumerate(
            zip(squarify(fam.values, 0, 0, 1, 1), fam.items())):
        bx.add_patch(Rectangle((x0, y0), w, h, facecolor=PAL[k % len(PAL)],
                               alpha=.92, edgecolor="white", linewidth=1.1))
        if w * h > .026:
            txt = name.replace("_", " ")
            # size by the narrower dimension, and wrap, so nothing overflows
            fs = max(4.0, min(6.6, 190 * min(w, h)))
            est = fs * 0.0062                      # approx char width in axes units
            maxch = max(6, int(w / est))
            if len(txt) > maxch:
                words, lines, cur = txt.split(), [], ""
                for wd in words:
                    if len(cur) + len(wd) + 1 <= maxch:
                        cur = (cur + " " + wd).strip()
                    else:
                        lines.append(cur); cur = wd
                lines.append(cur)
                txt = "\n".join([l for l in lines if l])
                fs = max(4.0, fs - 0.6)
            bx.text(x0 + w / 2, y0 + h / 2, f"{txt}\n{n}",
                    ha="center", va="center", color="white", linespacing=1.2,
                    fontsize=fs, fontweight="bold")
    bx.set_xlim(0, 1); bx.set_ylim(0, 1); bx.set_axis_off()
    bx.text(1, 1.02, f"{len(fam)} families among the {n_mapped} mapped proteins",
            fontsize=6.3, color=MUTE, va="bottom", ha="right")
    bx.set_title("B", loc="left", fontsize=11, fontweight="bold", x=-.03, pad=4)

    # ───────────── C: mosaic of the 2x2 ─────────────────────────────────────
    cx = fig.add_subplot(gs[1, 0])
    fes = d[d.DeJesus_call != "Uncertain"]
    f_e = int(fes.DeJesus_call.isin(["ES", "ESD", "GD"]).sum())
    f_n = len(fes) - f_e
    g_e, g_n = a.gw_esgd, a.gw_total - a.gw_esgd
    orr, p = fisher_exact([[f_e, f_n], [g_e, g_n]])

    groups = [("Fe–S proteins", f_e, f_n), ("Genome-wide", g_e, g_n)]
    raw = np.array([f_e + f_n, g_e + g_n], float)
    widths = .24 + .54 * raw / raw.max()
    gap, xpos = .11, 0.0
    for (lab, e, n), w in zip(groups, widths):
        tot, pe = e + n, e / (e + n)
        cx.add_patch(Rectangle((xpos, 0), w, pe, facecolor=CALL_COL["ES"],
                               alpha=.92, edgecolor="white", linewidth=1.3))
        cx.add_patch(Rectangle((xpos, pe), w, 1 - pe, facecolor=CALL_COL["NE"],
                               alpha=.55, edgecolor="white", linewidth=1.3))
        cx.text(xpos + w / 2, pe / 2, f"{e}\n{100*pe:.1f} %", ha="center",
                va="center", fontsize=6.8, color="white", fontweight="bold",
                linespacing=1.3)
        cx.text(xpos + w / 2, pe + (1 - pe) / 2, f"{n:,}", ha="center",
                va="center", fontsize=6.6, color="white", fontweight="bold")
        cx.text(xpos + w / 2, -.035, f"{lab}\nn = {tot:,}", ha="center",
                va="top", fontsize=6.9, color=INK, linespacing=1.3)
        xpos += w + gap
    cx.set_xlim(-.04, xpos - gap + .04); cx.set_ylim(-.22, 1.18)
    cx.set_axis_off()
    ptxt = "p < 0.001" if p < 1e-3 else f"p = {p:.2g}"
    cx.text((xpos - gap) / 2, 1.16,
            f"odds ratio {orr:.2f}   Fisher exact, {ptxt}",
            ha="center", va="top", fontsize=6.8, color=MUTE)
    cx.legend(handles=[Patch(facecolor=CALL_COL["ES"], alpha=.92,
                             label="essential or growth-defect"),
                       Patch(facecolor=CALL_COL["NE"], alpha=.55,
                             label="other")],
              fontsize=6.2, frameon=False, loc="upper center",
              bbox_to_anchor=(.5, -.04), ncol=2, handlelength=1.1)
    cx.set_title("C", loc="left", fontsize=11, fontweight="bold", x=-.04, pad=4)

    # ───────────── D: selectivity quadrant ──────────────────────────────────
    dx = fig.add_subplot(gs[1, 1])
    q = d[d.DeJesus_call != "Uncertain"].copy()
    q["grp"] = np.where(q.DeJesus_call.isin(["ES", "ESD"]), "essential",
                        np.where(q.DeJesus_call == "GD", "growth defect",
                                 "non-essential"))
    BAND = {"essential": (2.35, 4.00), "growth defect": (1.55, 2.05),
            "non-essential": (0.10, 1.25)}
    rng = np.random.default_rng(4)
    ypos = {}
    for g, (lo, hi) in BAND.items():
        sub = q[q.grp == g]
        nohom = sub[~sub.has_human].sort_values("gene", na_position="last")
        if len(nohom):
            for y, i in zip(np.linspace(hi - .06, lo + .06, len(nohom)),
                            nohom.index):
                ypos[i] = y
        for i in sub[sub.has_human].index:
            ypos[i] = rng.uniform(lo + .08, hi - .08)
    q["ypos"] = [ypos[i] for i in q.index]
    sel = q.DeJesus_call.isin(["ES", "ESD", "GD"]) & (~q.has_human)

    for g, (lo, hi) in BAND.items():
        if g != "non-essential":
            dx.add_patch(Rectangle((-1.6, lo - .06), 3.4, hi - lo + .12,
                                   facecolor=SEL, alpha=.15, lw=0, zorder=0))
    if a.id_cut > 0:
        dx.axvline(a.id_cut, color=MUTE, lw=.7, ls=(0, (3, 3)), zorder=2)
        dx.text(a.id_cut + .9, .10, f"{a.id_cut:g} %", fontsize=5.4,
                color=MUTE, va="bottom")

    pcol = q.DeJesus_call.map(CALL_COL).fillna("#8E9AA3")
    dx.scatter(q.loc[~sel, "pident"], q.loc[~sel, "ypos"], s=18,
               c=pcol[~sel], alpha=.60, edgecolors="white", linewidths=.5,
               zorder=3)
    dx.scatter(q.loc[sel, "pident"], q.loc[sel, "ypos"], s=46,
               c=pcol[sel], edgecolors=INK, linewidths=.9, zorder=4)
    for _, r in q[sel].iterrows():
        lab = (str(r.gene) if pd.notna(r.gene) and str(r.gene) != "nan"
               else str(r.Rv))
        dx.annotate(lab, (r.pident, r.ypos), textcoords="offset points",
                    xytext=(7, 0), ha="left", va="center", fontsize=5.5,
                    color=CALL_COL.get(r.DeJesus_call, INK), fontweight="bold")

    dx.set_yticks([(lo + hi) / 2 for lo, hi in BAND.values()])
    dx.set_yticklabels(["essential", "growth\ndefect", "non-essential"],
                       fontsize=6.8)
    for t, g in zip(dx.get_yticklabels(), BAND):
        t.set_color(CALL_COL["ES"] if g == "essential" else
                    CALL_COL["GD"] if g == "growth defect" else CALL_COL["NE"])
        t.set_fontweight("bold")
    dx.set_ylim(-.05, 4.40); dx.set_xlim(-1.6, 58)
    dx.set_xlabel("Identity to closest human protein (%)", fontsize=8.2)
    dx.text(-1.2, 4.35, f"{int(sel.sum())} essential or growth-defect proteins "
            "with no human homologue", fontsize=6.3, color=SEL,
            fontweight="bold", va="top")
    dx.text(57, .02, "0 % = no detectable human homologue", ha="right",
            fontsize=5.5, color=MUTE, style="italic")
    for s_ in ("top", "right"):
        dx.spines[s_].set_visible(False)
    for s_ in ("left", "bottom"):
        dx.spines[s_].set_linewidth(.7); dx.spines[s_].set_color("#4d4d4d")
    dx.tick_params(labelsize=7.0, colors=MUTE, width=.7)
    dx.tick_params(axis="y", length=0)
    dx.xaxis.grid(True, lw=.5, color=GRID, zorder=0); dx.set_axisbelow(True)
    dx.set_title("D", loc="left", fontsize=11, fontweight="bold", x=-.16, pad=4)

    for ext, kw in [("png", dict(dpi=400)), ("pdf", {}),
                    ("tiff", dict(dpi=600,
                                  pil_kwargs={"compression": "tiff_lzw"}))]:
        fig.savefig(f"Figure7.{ext}", bbox_inches="tight", **kw)

    print(f"\nA: {a.n_total} proteins ({n_mapped} mapped, {n_unmapped} unmapped)")
    for c, n in counts:
        print(f"     {CALL_LAB[c]:20s} {n}")
    print(f"B: {len(fam)} families; largest {fam.index[0]} ({fam.iloc[0]})")
    print(f"C: Fe–S {f_e}/{f_e+f_n} = {100*f_e/(f_e+f_n):.1f}%  vs  "
          f"genome-wide {g_e:,}/{g_e+g_n:,} = {100*g_e/(g_e+g_n):.1f}%")
    print(f"   odds ratio {orr:.2f}, Fisher p = {p:.3g}")
    print(f"D: {int(sel.sum())} selective candidates")
    for _, r in q[sel].iterrows():
        print(f"     {r.Rv:10s} {str(r.gene):12s} {r.family}")
    print("\nwrote Figure7.png / .pdf / .tiff")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-total", type=int, default=72)
    ap.add_argument("--gw-esgd", type=int, default=596)
    ap.add_argument("--gw-total", type=int, default=3943)
    ap.add_argument("--id-cut", type=float, default=0.0,
                    help="panel D reference line only; 0 hides it")
    ap.add_argument("--selectivity", default="selectivity.csv")
    ap.add_argument("--dejesus", default="step11_function/mtb_fes_dejesus.tsv")
    ap.add_argument("--human", default="mtb_vs_human.tsv")
    ap.add_argument("--evalue", type=float, default=1e-3)
    ap.add_argument("--hom-id", type=float, default=0.0)
    ap.add_argument("--hom-cov", type=float, default=0.0)
    main(ap.parse_args())
