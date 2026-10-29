#!/usr/bin/env python3
"""Figure 8 - the Fe-S repertoire across Mycobacteriaceae.

A  Species tree with the 44-family presence/absence matrix aligned to its tips,
   plus a lineage colour strip. The backbone panel: what is universal, what is
   patchy, and whether the patchiness follows the phylogeny.
B  Volcano: per-family enrichment or depletion in the M. tuberculosis complex
   relative to the rest of the family (Fisher exact, Benjamini-Hochberg).
C  Ridgeline: how many Fe-S proteins each lineage carries, as a distribution
   rather than a mean.
D  Principal coordinates of the genomes on Jaccard distance between Fe-S
   repertoires, coloured by lineage - does repertoire structure track lineage?

Run from /media/budheswar-dehury/DATADRIVE0/Fe-S/ :
    python3 build_figure8.py --inspect
    python3 build_figure8.py
"""
import argparse, glob, os, re, sys
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy.stats import fisher_exact, gaussian_kde

try:
    from Bio import Phylo
except ImportError:
    sys.exit("needs biopython:  pip install --break-system-packages biopython")

CAT = "step11_function/fes_v3_categorized.tsv"
INK, MUTE, GRID = "#1a1a1a", "#5c5c5c", "#ebebeb"
PRESENT, ABSENT = "#2F5C8A", "#F2F4F6"
UP, DOWN, NS = "#B5533C", "#2F5C8A", "#C3CBD1"
LPAL = ["#B5533C", "#2F5C8A", "#3E8E6F", "#D9983E", "#6B4E9B", "#386B7A",
        "#8C4F5E", "#B07C3F", "#527B4E", "#9B5B8A", "#4C7FA8", "#7A7A7A"]
FOCUS = None          # set at runtime
N_CLADES = 8          # major clades cut from the species tree
plt.rcParams.update({"font.family": "sans-serif",
                     "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
                     "pdf.fonttype": 42, "ps.fonttype": 42})


def tip_genome(t):
    return str(t).split("|")[0].strip()


def lineage(name):
    """GTDB genus. The set holds one genome per species, so genus is the only
    grouping with enough members to test; GTDB suffixes (_A, _B) are stripped."""
    s = str(name).strip()
    g = s.split()[0] if s.split() else s
    g = re.sub(r"_[A-Z]$", "", g)
    return g


PATHOGENS = {
    "tuberculosis", "bovis", "africanum", "canettii", "microti", "caprae",
    "orygis", "pinnipedii", "leprae", "lepromatosis", "abscessus", "chelonae",
    "fortuitum", "avium", "intracellulare", "chimaera", "kansasii", "marinum",
    "ulcerans", "szulgai", "xenopi", "malmoense", "simiae", "genavense",
    "haemophilum", "mucogenicum", "peregrinum", "immunogenum", "gordonae",
    "scrofulaceum", "celatum", "shimoidei", "asiaticum", "conspicuum",
    "interjectum", "heckeshornense", "bohemicum", "triplex", "lentiflavum",
}


def is_pathogen(name):
    parts = str(name).split()
    if len(parts) < 2:
        return False
    ep = re.sub(r"_[A-Z]$", "", parts[1]).lower()
    return ep in PATHOGENS


def tree_clades(tree, k):
    """Split the species tree into ~k clades by cutting every edge that
    crosses a fixed root-to-node depth. A depth cut gives balanced groups on a
    pectinate tree, where repeatedly splitting the largest subtree does not."""
    depth = {}

    def setd(node, acc):
        depth[id(node)] = acc
        for ch in node.clades:
            setd(ch, acc + (ch.branch_length or 0.0))
    sys.setrecursionlimit(10000)
    setd(tree.root, 0.0)
    tipmax = max(depth[id(t)] for t in tree.get_terminals())

    def cut_at(h):
        """Subtrees whose stem crosses depth h."""
        out, stack = [], [tree.root]
        while stack:
            n = stack.pop()
            for ch in n.clades:
                if depth[id(n)] < h <= depth[id(ch)]:
                    out.append(ch)
                else:
                    stack.append(ch)
        return out

    best, best_h = None, None
    lo, hi = 0.0, tipmax
    for _ in range(60):                      # bisect on the clade count
        h = (lo + hi) / 2
        cl = [c for c in cut_at(h) if len(c.get_terminals()) >= 2]
        n = len(cl)
        if best is None or abs(n - k) < abs(best - k):
            best, best_h = n, h
        if n > k:
            hi = h
        else:
            lo = h
    cl = sorted(cut_at(best_h), key=lambda c: -len(c.get_terminals()))
    out = {}
    for i, c in enumerate(cl, 1):
        for t in c.get_terminals():
            out[t.name] = f"Clade {i}"
    return out


def find(*pats):
    for p in pats:
        h = glob.glob(p)
        if h:
            return h[0]
    return None


def load():
    cat = pd.read_csv(CAT, sep="\t", dtype=str)
    cat.columns = [c.strip() for c in cat.columns]
    pres = (cat.groupby(["genome", "family"]).size().unstack(fill_value=0))
    spf = find("step*/genome_species.tsv", "genome_species.tsv")
    if not spf:
        sys.exit("genome_species.tsv not found")
    sp = pd.read_csv(spf, sep=None, engine="python", header=None, dtype=str)
    if not re.match(r"^GC[AF]_", str(sp.iloc[0, 0])):
        sp.columns = sp.iloc[0]
        sp = sp.iloc[1:]
    sp = sp.iloc[:, :2]
    sp.columns = ["genome", "species"]
    sp["lineage"] = sp.species.map(lineage)
    tre = find("step*/species_tree.treefile", "step*/species_tree.contree")
    tree = None
    if tre:
        tree = Phylo.read(tre, "newick")
        for t in tree.get_terminals():
            t.name = tip_genome(t.name)
        tree.ladderize()
    return pres, sp.set_index("genome"), tree, tre


def tree_coords(tree):
    """x = root-to-node distance, y = tip index; returns dicts and tip order."""
    tips = tree.get_terminals()
    y = {id(t): i for i, t in enumerate(tips)}
    x = {}

    def setx(node, acc):
        x[id(node)] = acc
        for ch in node.clades:
            setx(ch, acc + (ch.branch_length or 0.0))
    sys.setrecursionlimit(10000)
    setx(tree.root, 0.0)

    def sety(node):
        if node.is_terminal():
            return y[id(node)]
        vals = [sety(c) for c in node.clades]
        y[id(node)] = (min(vals) + max(vals)) / 2
        return y[id(node)]
    sety(tree.root)
    return x, y, [t.name for t in tips]


def draw_tree(ax, tree, x, y, colour_of):
    segs_h, segs_v = [], []

    def walk(node):
        for ch in node.clades:
            segs_h.append(((x[id(node)], y[id(ch)]), (x[id(ch)], y[id(ch)])))
            walk(ch)
        if node.clades:
            ys = [y[id(c)] for c in node.clades]
            segs_v.append(((x[id(node)], min(ys)), (x[id(node)], max(ys))))
    walk(tree.root)
    from matplotlib.collections import LineCollection
    ax.add_collection(LineCollection(segs_h + segs_v, colors="#6f6f6f",
                                     linewidths=.28, zorder=2))
    ax.set_xlim(0, max(x.values()) * 1.02)
    ax.set_ylim(-1, len(tree.get_terminals()))
    ax.invert_yaxis()
    ax.set_axis_off()


def pcoa(Dm, k=2):
    """Classical multidimensional scaling on a distance matrix."""
    n = len(Dm)
    J = np.eye(n) - np.ones((n, n)) / n
    B = -.5 * J @ (Dm ** 2) @ J
    w, V = np.linalg.eigh(B)
    idx = np.argsort(w)[::-1][:k]
    L = np.clip(w[idx], 0, None)
    return V[:, idx] * np.sqrt(L), w


def jaccard_matrix(X):
    X = X.astype(bool)
    inter = X.astype(int) @ X.T.astype(int)
    rs = X.sum(1)
    union = rs[:, None] + rs[None, :] - inter
    with np.errstate(divide="ignore", invalid="ignore"):
        J = np.where(union > 0, 1 - inter / union, 0.0)
    np.fill_diagonal(J, 0.0)
    return J


def inspect():
    pres, sp, tree, tre = load()
    print(f"presence matrix: {pres.shape[0]} genomes x {pres.shape[1]} families")
    print(f"species map: {len(sp)} genomes")
    print(sp.lineage.value_counts().to_string())
    if tree is None:
        print("\nno species tree found - panel A will fall back to lineage order")
    else:
        tips = [t.name for t in tree.get_terminals()]
        print(f"\nspecies tree: {tre}  {len(tips)} tips")
        ov = len(set(tips) & set(pres.index))
        print(f"  genomes shared with the presence matrix: {ov}")
        if ov == 0:
            print("  !! no overlap - check the tip label format")


def build():
    pres, sp, tree, tre = load()
    if tree is not None:
        _, _, order = tree_coords(tree)
        order = [g for g in order if g in pres.index]
    else:
        order = list(sp.sort_values("lineage").index)
        order = [g for g in order if g in pres.index]
    pres = pres.loc[order]
    if tree is not None:
        cl = tree_clades(tree, N_CLADES)
        lin = pd.Series([cl.get(g, "unplaced") for g in order], index=order)
    else:
        lin = sp.reindex(order).lineage.fillna("unplaced")
    vc = lin.value_counts()
    small = set(vc[vc < 5].index)
    if small:
        lin = lin.where(~lin.isin(small), "minor clades")
    lins = lin.value_counts().index.tolist()
    lcol = {l: LPAL[i % len(LPAL)] for i, l in enumerate(lins)}

    path = pd.Series([bool(is_pathogen(sp.species.get(g, ""))) for g in order],
                     index=order)
    global FOCUS
    FOCUS = "clinically relevant species"
    if path.sum() < 5 or (~path).sum() < 5:
        FOCUS = lin.value_counts().idxmax()
        path = (lin == FOCUS)
    print(f"panel B focus: {FOCUS}  (n={int(path.sum())} vs "
          f"{int((~path).sum())})")

    binp = (pres > 0).astype(int)
    prevalence = binp.mean() * 100
    fam_order = prevalence.sort_values(ascending=False).index.tolist()
    binp = binp[fam_order]

    fig = plt.figure(figsize=(7.4, 9.3), dpi=200)
    gsA_b, gsA_t = .505, .978
    L, R = .055, .965

    # ───────────────── A: tree + presence matrix ────────────────────────────
    tw, sw = .115, .012
    axT = fig.add_axes([L, gsA_b, tw, gsA_t - gsA_b])
    axS = fig.add_axes([L + tw + .004, gsA_b, sw, gsA_t - gsA_b])
    axM = fig.add_axes([L + tw + sw + .012, gsA_b,
                        R - L - tw - sw - .012, gsA_t - gsA_b])
    if tree is not None:
        x, y, _ = tree_coords(tree)
        draw_tree(axT, tree, x, y, lcol)
    else:
        axT.set_axis_off()
    axS.imshow(np.array([[lins.index(l) for l in lin]]).T, aspect="auto",
               cmap=ListedColormap([lcol[l] for l in lins]),
               vmin=0, vmax=len(lins) - 1, interpolation="nearest")
    axS.set_axis_off()
    axM.imshow(binp.values, aspect="auto", interpolation="nearest",
               cmap=ListedColormap([ABSENT, PRESENT]), vmin=0, vmax=1)
    axM.set_xticks(range(len(fam_order)))
    axM.set_xticklabels([f.replace("_", " ") for f in fam_order],
                        rotation=90, fontsize=4.0, color=INK)
    axM.set_yticks([])
    axM.tick_params(length=0)
    for s in axM.spines.values():
        s.set_visible(False)
    axM.text(0, 1.004, f"{binp.shape[0]} genomes  ×  {binp.shape[1]} families"
             "    ·    families ordered by prevalence    ·    "
             "filled = present",
             transform=axM.transAxes, fontsize=5.6, color=MUTE, va="bottom")
    axT.text(-.02, 1.004, "A", transform=axT.transAxes, fontsize=11,
             fontweight="bold", color=INK, va="bottom")

    # ───────────────── B: volcano, MTB complex vs the rest ──────────────────
    axB = fig.add_axes([L + .015, .232, .40, .140])
    inq = path.values
    rows = []
    for f in fam_order:
        v = binp[f].values.astype(bool)
        a, b = int((v & inq).sum()), int((~v & inq).sum())
        c, d_ = int((v & ~inq).sum()), int((~v & ~inq).sum())
        if (a + b) == 0 or (c + d_) == 0:
            continue
        orr, p = fisher_exact([[a, b], [c, d_]])
        # Haldane-Anscombe correction so a zero cell still plots
        orr = ((a + .5) * (d_ + .5)) / ((b + .5) * (c + .5))
        rows.append(dict(family=f, log2or=np.log2(orr), p=p,
                         pct_focus=100 * a / (a + b),
                         pct_other=100 * c / (c + d_)))
    vo = pd.DataFrame(rows)
    if vo.empty:
        sys.exit("panel B has no testable families - the focus group is "
                 "either empty or the whole set; check the grouping")
    m = len(vo)
    o = vo.p.sort_values()
    q = np.minimum.accumulate((o.values * m / np.arange(1, m + 1))[::-1])[::-1]
    vo.loc[o.index, "q"] = np.clip(q, 0, 1)
    vo["sig"] = vo.q < .05
    col = [UP if (r.sig and r.log2or > 0) else DOWN if (r.sig and r.log2or < 0)
           else NS for _, r in vo.iterrows()]
    axB.axhline(-np.log10(.05), color=MUTE, lw=.7, ls=(0, (3, 3)), zorder=2)
    axB.axvline(0, color=MUTE, lw=.7, zorder=2)
    axB.scatter(vo.log2or, -np.log10(vo.p.clip(lower=1e-300)), s=18, c=col,
                alpha=.85, edgecolors="white", linewidths=.5, zorder=3)
    for k, (_, r) in enumerate(vo[vo.sig].reindex(
            vo[vo.sig].log2or.abs().sort_values(ascending=False).index)
            .head(4).iterrows()):
        axB.annotate(r.family.replace("_", " "),
                     (r.log2or, -np.log10(max(r.p, 1e-300))),
                     textcoords="offset points",
                     xytext=[(7, 7), (-7, 7), (7, -8), (-7, -8)][k % 4],
                     ha="left" if k % 2 == 0 else "right", fontsize=5.0,
                     color=MUTE,
                     arrowprops=dict(arrowstyle="-", color="#c8c8c8", lw=.4,
                                     shrinkA=0, shrinkB=2))
    axB.set_xlabel(f"log$_2$ odds ratio   ({FOCUS} vs the rest)", fontsize=7.4)
    axB.set_ylabel("−log$_{10}$ p", fontsize=7.4)
    axB.text(.985, .16, f"{int((vo.sig & (vo.log2or>0)).sum())} enriched · "
             f"{int((vo.sig & (vo.log2or<0)).sum())} depleted at BH q < 0.05",
             transform=axB.transAxes, va="bottom", ha="right", fontsize=5.4,
             color=MUTE)
    axB.set_title("B", loc="left", fontsize=11, fontweight="bold",
                  x=-.17, pad=4)

    # ───────────────── C: ridgeline of complement size by lineage ───────────
    axC = fig.add_axes([L + .56, .232, .375, .140])
    tot = pres.sum(axis=1)
    med = {l: tot[lin == l].median() for l in lins}
    seq = sorted(lins, key=lambda l: med[l])
    xs = np.linspace(tot.min() * .92, tot.max() * 1.06, 260)
    for i, l in enumerate(seq):
        v = tot[lin == l].values.astype(float)
        if len(v) == 0:
            continue
        if len(v) > 2 and v.std() > 0:
            dens = gaussian_kde(v)(xs)
            dens = dens / dens.max() * .86
        else:
            dens = np.exp(-((xs - v.mean()) ** 2) / (2 * max(v.std(), 3) ** 2)) * .86
        axC.fill_between(xs, i, i + dens, color=lcol[l], alpha=.45, lw=0,
                         zorder=3 + i)
        axC.plot(xs, i + dens, color=lcol[l], lw=.8, zorder=3 + i)
        axC.scatter(v, np.full(len(v), i - .07), s=2.2, color=lcol[l],
                    alpha=.5, lw=0, zorder=20 + i)
        axC.text(xs[0], i + .30, f"{l}  n={len(v)}", fontsize=5.0,
                 color=lcol[l], va="bottom", fontweight="bold")
    axC.set_xlabel("Fe–S proteins per genome", fontsize=7.4)
    axC.set_yticks([])
    axC.set_ylim(-.4, len(seq) + .5)
    axC.set_title("C", loc="left", fontsize=11, fontweight="bold",
                  x=-.06, pad=4)

    # ───────────────── D: PCoA of repertoires ───────────────────────────────
    axD = fig.add_axes([L + .015, .038, .40, .140])
    J = jaccard_matrix(binp.values)
    Y, w = pcoa(J, 2)
    ev = np.clip(w, 0, None)
    ev = 100 * ev[np.argsort(ev)[::-1][:2]] / ev.sum()
    for l in lins:
        s = (lin == l).values
        axD.scatter(Y[s, 0], Y[s, 1], s=11, color=lcol[l], alpha=.78,
                    edgecolors="white", linewidths=.35, zorder=3, label=l)
    axD.set_xlabel(f"PCo1 ({ev[0]:.1f} %)", fontsize=7.4)
    axD.set_ylabel(f"PCo2 ({ev[1]:.1f} %)", fontsize=7.4)
    axD.set_title("D", loc="left", fontsize=11, fontweight="bold",
                  x=-.17, pad=4)
    axD.legend(fontsize=4.8, frameon=False, loc="upper left",
               bbox_to_anchor=(1.01, 1.02), handletextpad=.25,
               labelspacing=.3, borderpad=0)

    for a in (axB, axC, axD):
        for s in ("top", "right"):
            a.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            a.spines[s].set_linewidth(.7); a.spines[s].set_color("#4d4d4d")
        a.tick_params(labelsize=6.6, colors=MUTE, width=.7)
        a.grid(True, lw=.5, color=GRID, zorder=0); a.set_axisbelow(True)
    axC.spines["left"].set_visible(False)

    fig.savefig("Figure8.png", dpi=400, bbox_inches="tight")
    fig.savefig("Figure8.pdf", bbox_inches="tight")
    fig.savefig("Figure8.tiff", dpi=600, bbox_inches="tight",
                pil_kwargs={"compression": "tiff_lzw"})
    vo.sort_values("q").to_csv("fig8_mtb_enrichment.csv", index=False)
    pd.DataFrame({"genome": order, "lineage": lin.values,
                  "fes_total": tot.values,
                  "PCo1": Y[:, 0], "PCo2": Y[:, 1]}).to_csv(
        "fig8_genome_coords.csv", index=False)

    print(f"{binp.shape[0]} genomes x {binp.shape[1]} families")
    print(f"universal (100 % of genomes): {int((prevalence == 100).sum())} families")
    print(f"\nlineage medians (Fe–S proteins per genome):")
    for l in seq[::-1]:
        v = tot[lin == l]
        print(f"  {l:26s} n={len(v):3d}  median {v.median():5.1f}  "
              f"range {v.min():.0f}-{v.max():.0f}")
    sig = vo[vo.sig]
    print(f"\n{FOCUS} vs the rest: {len(sig)} families differ at BH q < 0.05")
    for _, r in sig.reindex(sig.log2or.abs().sort_values(
            ascending=False).index).head(8).iterrows():
        d_ = "enriched" if r.log2or > 0 else "depleted"
        print(f"  {r.family:30s} {d_:8s} log2OR {r.log2or:+5.2f}  "
              f"{r.pct_focus:5.1f}% vs {r.pct_other:5.1f}%  q={r.q:.2g}")
    print(f"\nPCoA: PCo1 {ev[0]:.1f} %, PCo2 {ev[1]:.1f} %")
    print("\nwrote Figure8.png/.pdf/.tiff, fig8_mtb_enrichment.csv, "
          "fig8_genome_coords.csv")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--inspect", action="store_true")
    ap.add_argument("--focus", default=None,
                    help="group to test in panel B (default: the largest)")
    a = ap.parse_args()
    if a.focus:
        FOCUS = a.focus
    inspect() if a.inspect else build()
