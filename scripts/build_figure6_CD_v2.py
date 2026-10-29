#!/usr/bin/env python3
"""Figure 6, panels C and D  (v2 - correct congruence metric).

C  Gene tree vs species tree congruence. For each family, patristic distances
   between genomes in the gene tree are correlated (Spearman) with the
   distances between the same genome pairs in the 138-marker species tree.
   rho near 1 = the family tracks the species phylogeny; low rho = it does not.
   This is rooting-independent, so it is safe on unrooted IQ-TREE output, and
   one copy per genome is used so paralogues cannot distort it.

D  Duplication mode. Copies per genome against the fraction of multi-copy
   genomes whose copies form a single group on the tree, tested in the
   UNROOTED sense (a set is clustered if one edge separates it from the rest).
   Top right = recent lineage-specific duplication; bottom = ancient paralogy.

Run from /media/budheswar-dehury/DATADRIVE0/Fe-S/ :
    python3 build_figure6_CD_v2.py --inspect
    python3 build_figure6_CD_v2.py
"""
import argparse, glob, os, re, sys
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

try:
    from Bio import Phylo
except ImportError:
    sys.exit("needs biopython:  pip install --break-system-packages biopython")

PHY = "step13_phylo"
INK, MUTE, GRID = "#1a1a1a", "#5c5c5c", "#ebebeb"
PAL = ["#2F5C8A", "#B5533C", "#3E8E6F", "#D9983E", "#6B4E9B", "#386B7A",
       "#8C4F5E", "#B07C3F", "#527B4E", "#9B5B8A", "#4C7FA8", "#C47A52"]
MAX_PAIR_GENOMES = 120      # subsample for speed; distances are stable well below this
MIN_SHARED = 25             # a family needs this many shared genomes to be scored
plt.rcParams.update({"font.family": "sans-serif",
                     "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
                     "pdf.fonttype": 42, "ps.fonttype": 42})


def tip_genome(t):
    return str(t).split("|")[0].strip()


def family_of(p):
    b = os.path.basename(p)
    return re.sub(r"(_tree)?(_clean)?\.treefile$", "",
                  re.sub(r"_tree(_clean)?\.treefile$", "", b))


def find_species_tree():
    for pat in ("step*/species_tree.treefile", "step*/species_tree.contree",
                "species_tree.treefile"):
        h = glob.glob(pat)
        if h:
            return h[0]
    sys.exit("species_tree.treefile not found under step*/")


def paths_to_root(tree):
    """tip name -> (list of node ids root..tip, cumulative distance at each)."""
    out = {}
    stack = [(tree.root, [id(tree.root)], [0.0])]
    while stack:
        node, ids, dists = stack.pop()
        if node.is_terminal():
            out[node.name] = (ids, dists)
            continue
        for ch in node.clades:
            bl = ch.branch_length or 0.0
            stack.append((ch, ids + [id(ch)], dists + [dists[-1] + bl]))
    return out


def patristic(tree, want):
    """Symmetric distance matrix over the tips named in `want` (ordered)."""
    P = paths_to_root(tree)
    use = [w for w in want if w in P]
    n = len(use)
    M = np.zeros((n, n))
    for i in range(n):
        ia, da = P[use[i]]
        seta = {v: k for k, v in enumerate(ia)}
        for j in range(i + 1, n):
            ib, db = P[use[j]]
            k = 0
            # walk down from the root while the two paths agree
            while k < len(ib) and ib[k] in seta and seta[ib[k]] == k:
                k += 1
            k -= 1
            d = (da[-1] - da[k]) + (db[-1] - db[k])
            M[i, j] = M[j, i] = d
    return use, M


def unrooted_clustered(tree, members, total):
    """True if `members` is separated from the rest of the tree by one edge."""
    s = set(m.name for m in members)
    for node in tree.get_nonterminals():
        names = set(t.name for t in node.get_terminals())
        if names == s or names == set(
                t.name for t in tree.get_terminals()) - s:
            return True
    return False


def load_species():
    sp_tree = find_species_tree()
    t = Phylo.read(sp_tree, "newick")
    # species-tree tips may be bare accessions or accession|something
    ren = {}
    for tip in t.get_terminals():
        ren[tip.name] = tip_genome(tip.name)
        tip.name = ren[tip.name]
    return sp_tree, t


def collect_trees():
    files = sorted(glob.glob(f"{PHY}/*.treefile"))
    best = {}
    for f in files:
        fam = family_of(f)
        # prefer *_clean over the raw tree for the same family
        if fam not in best or ("clean" in f and "clean" not in best[fam]):
            best[fam] = f
    return best


def inspect():
    sp_path, sp = load_species()
    tips = [t.name for t in sp.get_terminals()]
    print(f"species tree: {sp_path}   {len(tips)} tips")
    print("  e.g.", tips[:3])
    best = collect_trees()
    print(f"\n{len(best)} families after deduplication "
          f"(from {len(glob.glob(f'{PHY}/*.treefile'))} treefiles)")
    t = Phylo.read(list(best.values())[0], "newick")
    gen = {tip_genome(x.name) for x in t.get_terminals()}
    print(f"first family {list(best)[0]}: {len(t.get_terminals())} tips, "
          f"{len(gen)} genomes")
    print(f"  genomes shared with the species tree: {len(gen & set(tips))}")
    if not (gen & set(tips)):
        print("  !! no overlap - the species tree labels its tips differently")


def build():
    sp_path, sp = load_species()
    sp_tips = [t.name for t in sp.get_terminals()]
    best = collect_trees()

    rows = []
    sp_cache = {}
    for fam, f in sorted(best.items()):
        try:
            tr = Phylo.read(f, "newick")
        except Exception as e:
            print(f"  skipped {fam}: {e}")
            continue
        tips = tr.get_terminals()
        by_gen = {}
        for t in tips:
            by_gen.setdefault(tip_genome(t.name), []).append(t)
        n_tips, n_gen = len(tips), len(by_gen)
        copies = n_tips / n_gen if n_gen else np.nan

        shared = [g for g in by_gen if g in set(sp_tips)]
        rho = np.nan
        if len(shared) >= MIN_SHARED:
            if len(shared) > MAX_PAIR_GENOMES:
                shared = list(np.random.default_rng(0).choice(
                    shared, MAX_PAIR_GENOMES, replace=False))
            rep = {g: by_gen[g][0].name for g in shared}
            used_g, Mg = patristic(tr, [rep[g] for g in shared])
            gmap = {v: k for k, v in rep.items()}
            order_g = [gmap[u] for u in used_g]
            key = tuple(sorted(order_g))
            if key not in sp_cache:
                names_s, Ms_full = patristic(sp, sorted(order_g))
                sp_cache[key] = (list(names_s), Ms_full)
            names_s, Ms_full = sp_cache[key]
            pos = {g: i for i, g in enumerate(names_s)}
            if all(g in pos for g in order_g) and len(order_g) >= MIN_SHARED:
                idx = [pos[g] for g in order_g]          # reorder to match Mg
                Ms = Ms_full[np.ix_(idx, idx)]
                iu = np.triu_indices(len(order_g), 1)
                rho = spearmanr(Mg[iu], Ms[iu]).statistic

        multi = {g: v for g, v in by_gen.items() if len(v) > 1}
        clustered = sum(1 for g, v in multi.items()
                        if unrooted_clustered(tr, v, n_tips))
        frac = clustered / len(multi) if multi else np.nan

        rows.append(dict(family=fam, n_tips=n_tips, n_genomes=n_gen,
                         copies_per_genome=copies, n_shared=len(shared),
                         rho=rho, n_multicopy_genomes=len(multi),
                         frac_paralogs_clustered=frac))
        print(f"  {fam:32s} {n_tips:4d} tips  "
              f"rho {'  n/a' if np.isnan(rho) else f'{rho:5.2f}'}  "
              f"copies {copies:.2f}")

    d = pd.DataFrame(rows)
    d.to_csv("fig6_congruence.csv", index=False)
    c = d.dropna(subset=["rho"]).sort_values("rho").reset_index(drop=True)

    fig, (cx, dx) = plt.subplots(1, 2,
                                 figsize=(7.4, max(4.6, .185 * len(c) + 1.4)),
                                 dpi=200,
                                 gridspec_kw=dict(width_ratios=[1.25, 1],
                                                  wspace=.44))

    # ───────────────── C: congruence with the species tree ──────────────────
    y = np.arange(len(c))
    cx.axvspan(-.05, .5, color="#B5533C", alpha=.045, lw=0, zorder=0)
    cx.hlines(y, 0, c.rho, color="#B9C7D1", lw=1.3, zorder=2)
    sc = cx.scatter(c.rho, y, s=16 + 70 * c.n_shared / c.n_shared.max(),
                    c=c.copies_per_genome, cmap="YlGnBu", vmin=1, zorder=4,
                    edgecolors="white", linewidths=.6)
    cx.axvline(0, color=MUTE, lw=.8, zorder=3)
    cx.set_yticks(y)
    cx.set_yticklabels([f.replace("_", " ") for f in c.family], fontsize=5.6,
                       color=INK)
    cx.set_xlabel("Congruence with the species tree\n"
                  "(Spearman ρ, patristic distances)", fontsize=8.0)
    cx.set_xlim(-.08, 1.04)
    cx.set_ylim(-.8, len(c) - .2)
    cx.tick_params(axis="y", length=0)
    cb = fig.colorbar(sc, ax=cx, fraction=.030, pad=.015)
    cb.set_label("copies per genome", fontsize=6.2, color=INK)
    cb.ax.tick_params(labelsize=5.6, colors=MUTE, width=.6, length=2)
    cb.outline.set_visible(False)
    cx.set_title("C", loc="left", fontsize=11, fontweight="bold", x=-.42, pad=8)

    # ───────────────── D: duplication mode ──────────────────────────────────
    m = d.dropna(subset=["frac_paralogs_clustered"])
    m = m[m.n_multicopy_genomes >= 3]
    dx.axhline(.5, color=MUTE, lw=.7, ls=(0, (3, 3)), zorder=2)
    dx.scatter(m.copies_per_genome, m.frac_paralogs_clustered,
               s=16 + 70 * m.n_multicopy_genomes / m.n_multicopy_genomes.max(),
               c=[PAL[i % len(PAL)] for i in range(len(m))], alpha=.80,
               edgecolors="white", linewidths=.6, zorder=3)
    for k, (_, r) in enumerate(m.nlargest(3, "copies_per_genome").iterrows()):
        dx.annotate(r.family.replace("_", " "),
                    (r.copies_per_genome, r.frac_paralogs_clustered),
                    textcoords="offset points",
                    xytext=[(-9, 10), (-9, -11), (10, 9)][k % 3],
                    ha="right" if k < 2 else "left", fontsize=5.4, color=MUTE,
                    arrowprops=dict(arrowstyle="-", color="#c8c8c8", lw=.4,
                                    shrinkA=0, shrinkB=3))
    dx.text(.98, .97, "copies are each other's closest\nrelatives → recent,"
            " lineage-specific\nduplication", transform=dx.transAxes,
            ha="right", va="top", fontsize=5.3, color=MUTE, linespacing=1.35)
    dx.text(.98, .03, "copies scattered across the tree\n→ ancient paralogues",
            transform=dx.transAxes, ha="right", va="bottom", fontsize=5.3,
            color=MUTE, linespacing=1.35)
    dx.set_xlabel("Copies per genome", fontsize=8.0)
    dx.set_ylabel("Multi-copy genomes whose copies cluster", fontsize=8.0)
    dx.set_ylim(-.07, 1.14)
    dx.set_title("D", loc="left", fontsize=11, fontweight="bold", x=-.24, pad=8)

    for a in (cx, dx):
        for s in ("top", "right"):
            a.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            a.spines[s].set_linewidth(.7); a.spines[s].set_color("#4d4d4d")
        a.tick_params(labelsize=7.2, colors=MUTE, width=.7)
        a.grid(True, lw=.5, color=GRID, zorder=0); a.set_axisbelow(True)
    cx.yaxis.grid(False)

    fig.savefig("Figure6_CD.png", dpi=400, bbox_inches="tight")
    fig.savefig("Figure6_CD.pdf", bbox_inches="tight")
    fig.savefig("Figure6_CD.tiff", dpi=600, bbox_inches="tight",
                pil_kwargs={"compression": "tiff_lzw"})

    print(f"\n{len(c)} families scored against the species tree "
          f"({len(d) - len(c)} had fewer than {MIN_SHARED} shared genomes)")
    if len(c):
        print(f"congruence rho: median {c.rho.median():.2f}, "
              f"range {c.rho.min():.2f} to {c.rho.max():.2f}; "
              f"{int((c.rho >= .7).sum())} families at rho >= 0.70")
        print("\nmost congruent (vertical inheritance):")
        for _, r in c.nlargest(5, "rho").iterrows():
            print(f"  {r.family:30s} rho {r.rho:.2f}  "
                  f"({int(r.n_shared)} genomes, {r.copies_per_genome:.1f} copies)")
        print("\nleast congruent:")
        for _, r in c.nsmallest(5, "rho").iterrows():
            print(f"  {r.family:30s} rho {r.rho:.2f}  "
                  f"({int(r.n_shared)} genomes, {r.copies_per_genome:.1f} copies)")
    print("\nwrote Figure6_CD.png/.pdf/.tiff, fig6_congruence.csv")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--inspect", action="store_true")
    a = ap.parse_args()
    inspect() if a.inspect else build()
