#!/usr/bin/env python3
"""Figure 6, panels C and D - what the family trees actually show.

C  Clade retention: for each family, the fraction of species groups that come
   out monophyletic in its gene tree. High = vertical inheritance; low = the
   gene tree disagrees with the species tree (duplication, transfer, or a tree
   too poorly resolved to say).
D  Duplication mode: copies per genome against the fraction of multi-copy
   genomes whose copies are each other's closest relatives. Top right =
   recent, lineage-specific duplication. Bottom right = ancient paralogues
   that predate the genus.

Also writes an iTOL colour strip so panels A and B (your Rieske and WhiB
trees) can have their tips coloured by species group instead of one flat
colour per panel.

Run from /media/budheswar-dehury/DATADRIVE0/Fe-S/ :
    python3 build_figure6_CD.py --inspect     # check the inputs parse
    python3 build_figure6_CD.py
"""
import argparse, glob, os, re, sys
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

try:
    from Bio import Phylo
except ImportError:
    sys.exit("needs biopython:  pip install --break-system-packages biopython")

PHY = "step13_phylo"
INK, MUTE, GRID = "#1a1a1a", "#5c5c5c", "#ebebeb"
PAL = ["#2F5C8A", "#B5533C", "#3E8E6F", "#D9983E", "#6B4E9B", "#386B7A",
       "#8C4F5E", "#B07C3F", "#527B4E", "#9B5B8A", "#4C7FA8", "#C47A52"]
MIN_CLADE = 3            # species groups smaller than this are not testable
plt.rcParams.update({"font.family": "sans-serif",
                     "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
                     "pdf.fonttype": 42, "ps.fonttype": 42})


def find_species_map():
    """genome_species.tsv, wherever it sits under step*/."""
    hits = glob.glob("step*/genome_species.tsv") + glob.glob("genome_species.tsv")
    if not hits:
        sys.exit("genome_species.tsv not found under step*/")
    d = pd.read_csv(hits[0], sep=None, engine="python", header=None,
                    dtype=str).dropna(how="all")
    # first row may or may not be a header; detect by whether col0 looks like an accession
    if not re.match(r"^GC[AF]_", str(d.iloc[0, 0])):
        d.columns = d.iloc[0]
        d = d.iloc[1:]
    d = d.iloc[:, :2]
    d.columns = ["genome", "species"]
    d["genome"] = d.genome.str.strip()
    d["species"] = d.species.str.strip()
    return hits[0], d


def species_group(name):
    """Collapse a species name to the group the figure colours by."""
    s = str(name)
    low = s.lower()
    if "tuberculosis" in low or "bovis" in low or "africanum" in low \
            or "microti" in low or "canetti" in low or "caprae" in low:
        return "MTB complex"
    if "leprae" in low or "lepromatosis" in low:
        return "M. leprae / lepromatosis"
    if "abscessus" in low or "chelonae" in low:
        return "M. abscessus group"
    if "avium" in low or "intracellulare" in low or "chimaera" in low:
        return "MAC"
    if "ulcerans" in low or "marinum" in low:
        return "M. marinum / ulcerans"
    if low.startswith("mycobacteroides"):
        return "Mycobacteroides"
    if low.startswith("mycolicibacter"):
        return "Mycolicibacterium"
    # fall back to the genus + first epithet so nothing is lumped silently
    parts = s.split()
    return " ".join(parts[:2]) if len(parts) > 1 else s


def tip_genome(tip):
    return str(tip).split("|")[0].strip()


def family_of(path):
    b = os.path.basename(path)
    return re.sub(r"(_tree)?(_clean)?\.treefile$", "",
                  re.sub(r"_tree(_clean)?\.treefile$", "", b))


def analyse(tree, g2grp):
    """Clade retention and duplication mode for one gene tree."""
    tips = tree.get_terminals()
    names = [t.name for t in tips]
    gen = [tip_genome(n) for n in names]
    grp = [g2grp.get(g) for g in gen]

    n_tips, n_gen = len(names), len(set(gen))
    copies = n_tips / n_gen if n_gen else np.nan

    # ---- clade retention ---------------------------------------------------
    testable, mono = 0, 0
    by_grp = {}
    for t, g in zip(tips, grp):
        if g:
            by_grp.setdefault(g, []).append(t)
    for g, members in by_grp.items():
        if len(members) < MIN_CLADE:
            continue
        testable += 1
        try:
            anc = tree.common_ancestor(members)
            if len(anc.get_terminals()) == len(members):
                mono += 1
        except Exception:
            pass
    retention = mono / testable if testable else np.nan

    # ---- duplication mode --------------------------------------------------
    by_gen = {}
    for t, g in zip(tips, gen):
        by_gen.setdefault(g, []).append(t)
    multi = {g: v for g, v in by_gen.items() if len(v) > 1}
    clustered = 0
    for g, v in multi.items():
        try:
            anc = tree.common_ancestor(v)
            if len(anc.get_terminals()) == len(v):
                clustered += 1
        except Exception:
            pass
    frac_clustered = clustered / len(multi) if multi else np.nan

    return dict(n_tips=n_tips, n_genomes=n_gen, copies_per_genome=copies,
                clades_testable=testable, clades_monophyletic=mono,
                clade_retention=retention, n_multicopy_genomes=len(multi),
                frac_paralogs_clustered=frac_clustered)


def load_all():
    src, sp = find_species_map()
    sp["group"] = sp.species.map(species_group)
    g2grp = dict(zip(sp.genome, sp.group))
    files = sorted(glob.glob(f"{PHY}/*.treefile"))
    return src, sp, g2grp, files


def inspect():
    src, sp, g2grp, files = load_all()
    print(f"species map: {src}   {len(sp)} genomes")
    print(sp.head(4).to_string(index=False))
    print("\nspecies groups:")
    print(sp.group.value_counts().to_string())
    print(f"\n{len(files)} treefiles under {PHY}/")
    t = Phylo.read(files[0], "newick")
    tips = [x.name for x in t.get_terminals()][:3]
    print(f"first tree {os.path.basename(files[0])}: "
          f"{len(t.get_terminals())} tips, e.g. {tips}")
    hit = sum(1 for x in t.get_terminals() if tip_genome(x.name) in g2grp)
    print(f"  {hit}/{len(t.get_terminals())} tips map to a species group")
    if hit == 0:
        print("  !! tip labels do not match the species map - check the format")


def build():
    src, sp, g2grp, files = load_all()
    rows = []
    for f in files:
        fam = family_of(f)
        try:
            tr = Phylo.read(f, "newick")
        except Exception as e:
            print(f"  skipped {fam}: {e}")
            continue
        r = analyse(tr, g2grp)
        r["family"] = fam
        rows.append(r)
        print(f"  {fam:34s} {r['n_tips']:4d} tips  "
              f"retention {r['clade_retention']:.2f}  "
              f"copies {r['copies_per_genome']:.2f}")
    d = pd.DataFrame(rows).dropna(subset=["clade_retention"])
    d = d.sort_values("clade_retention").reset_index(drop=True)
    d.to_csv("fig6_congruence.csv", index=False)

    fig, (cx, dx) = plt.subplots(1, 2, figsize=(7.4, max(4.6, .185 * len(d) + 1.3)),
                                 dpi=200,
                                 gridspec_kw=dict(width_ratios=[1.25, 1],
                                                  wspace=.42))

    # ───────────────── C: clade retention ───────────────────────────────────
    y = np.arange(len(d))
    cx.axvspan(0, .5, color="#B5533C", alpha=.045, lw=0, zorder=0)
    cx.hlines(y, 0, d.clade_retention, color="#B9C7D1", lw=1.3, zorder=2)
    sz = 16 + 70 * d.clades_testable / d.clades_testable.max()
    sc = cx.scatter(d.clade_retention, y, s=sz, c=d.copies_per_genome,
                    cmap="YlGnBu", vmin=1, zorder=4, edgecolors="white",
                    linewidths=.6)
    cx.set_yticks(y)
    cx.set_yticklabels([f.replace("_", " ") for f in d.family], fontsize=5.6,
                       color=INK)
    cx.set_xlabel("Species groups recovered as monophyletic", fontsize=8.3)
    cx.set_xlim(-.03, 1.06)
    cx.set_ylim(-.8, len(d) - .2)
    cx.tick_params(axis="y", length=0)
    cb = fig.colorbar(sc, ax=cx, fraction=.030, pad=.015)
    cb.set_label("copies per genome", fontsize=6.2, color=INK)
    cb.ax.tick_params(labelsize=5.6, colors=MUTE, width=.6, length=2)
    cb.outline.set_visible(False)
    cx.text(0, len(d) + .25, "gene tree disagrees\nwith species groups",
            fontsize=5.6, color="#B5533C", va="bottom", linespacing=1.3)
    cx.set_title("C", loc="left", fontsize=11, fontweight="bold", x=-.40, pad=12)

    # ───────────────── D: duplication mode ──────────────────────────────────
    m = d.dropna(subset=["frac_paralogs_clustered"])
    m = m[m.n_multicopy_genomes >= 3]
    dx.axhline(.5, color=MUTE, lw=.7, ls=(0, (3, 3)), zorder=2)
    dx.scatter(m.copies_per_genome, m.frac_paralogs_clustered,
               s=16 + 70 * m.n_multicopy_genomes / m.n_multicopy_genomes.max(),
               c=[PAL[i % len(PAL)] for i in range(len(m))], alpha=.80,
               edgecolors="white", linewidths=.6, zorder=3)
    for _, r in m.nlargest(3, "copies_per_genome").iterrows():
        dx.annotate(r.family.replace("_", " "),
                    (r.copies_per_genome, r.frac_paralogs_clustered),
                    textcoords="offset points", xytext=(-9, 9), ha="right",
                    fontsize=5.4, color=MUTE,
                    arrowprops=dict(arrowstyle="-", color="#c8c8c8", lw=.4,
                                    shrinkA=0, shrinkB=3))
    dx.text(.985, .965, "copies are each other's\nclosest relatives\n"
            "→ recent duplication", transform=dx.transAxes, ha="right",
            va="top", fontsize=5.4, color=MUTE, linespacing=1.35)
    dx.text(.985, .035, "copies scattered across the tree\n→ ancient paralogues",
            transform=dx.transAxes, ha="right", va="bottom", fontsize=5.4,
            color=MUTE, linespacing=1.35)
    dx.set_xlabel("Copies per genome", fontsize=8.3)
    dx.set_ylabel("Multi-copy genomes with monophyletic copies", fontsize=8.3)
    dx.set_ylim(-.06, 1.12)
    dx.set_title("D", loc="left", fontsize=11, fontweight="bold", x=-.22, pad=12)

    for a in (cx, dx):
        for s in ("top", "right"):
            a.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            a.spines[s].set_linewidth(.7)
            a.spines[s].set_color("#4d4d4d")
        a.tick_params(labelsize=7.2, colors=MUTE, width=.7)
        a.grid(True, lw=.5, color=GRID, zorder=0)
        a.set_axisbelow(True)
    cx.yaxis.grid(False)

    fig.savefig("Figure6_CD.png", dpi=400, bbox_inches="tight")
    fig.savefig("Figure6_CD.pdf", bbox_inches="tight")
    fig.savefig("Figure6_CD.tiff", dpi=600, bbox_inches="tight",
                pil_kwargs={"compression": "tiff_lzw"})

    # ───────────────── iTOL colour strip for panels A and B ─────────────────
    groups = sorted(sp.group.unique())
    gcol = {g: PAL[i % len(PAL)] for i, g in enumerate(groups)}
    with open("itol_clade_colourstrip.txt", "w") as fh:
        fh.write("DATASET_COLORSTRIP\nSEPARATOR TAB\n"
                 "DATASET_LABEL\tspecies group\nCOLOR\t#2F5C8A\n"
                 "COLOR_BRANCHES\t1\n"
                 "LEGEND_TITLE\tspecies group\n"
                 "LEGEND_SHAPES\t" + "\t".join(["1"] * len(groups)) + "\n"
                 "LEGEND_COLORS\t" + "\t".join(gcol[g] for g in groups) + "\n"
                 "LEGEND_LABELS\t" + "\t".join(groups) + "\nDATA\n")
        for f in files:
            try:
                tr = Phylo.read(f, "newick")
            except Exception:
                continue
            for t in tr.get_terminals():
                g = g2grp.get(tip_genome(t.name))
                if g:
                    fh.write(f"{t.name}\t{gcol[g]}\t{g}\n")

    print(f"\n{len(d)} families scored")
    print(f"clade retention: median {d.clade_retention.median():.2f}, "
          f"{int((d.clade_retention >= .8).sum())} families >= 0.80, "
          f"{int((d.clade_retention < .5).sum())} below 0.50")
    print("\nleast congruent with the species groups:")
    for _, r in d.nsmallest(5, "clade_retention").iterrows():
        print(f"  {r.family:32s} {r.clade_retention:.2f}  "
              f"({int(r.clades_monophyletic)}/{int(r.clades_testable)} clades, "
              f"{r.copies_per_genome:.1f} copies/genome)")
    print("\nmost congruent:")
    for _, r in d.nlargest(4, "clade_retention").iterrows():
        print(f"  {r.family:32s} {r.clade_retention:.2f}")
    print("\nwrote Figure6_CD.png/.pdf/.tiff, fig6_congruence.csv, "
          "itol_clade_colourstrip.txt")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--inspect", action="store_true")
    a = ap.parse_args()
    inspect() if a.inspect else build()
