#!/usr/bin/env Rscript
# Figure 6 - four family phylogenies, rendered identically.
#
# Every panel uses the same support threshold, the same dot size, the same
# line width, the same fan layout and the same panel box, so the four are
# directly comparable. A scale bar is drawn on each.
#
# Run from /media/budheswar-dehury/DATADRIVE0/Fe-S/ :
#     Rscript build_figure6.R
#
# If a panel reports "no support values", that tree was inferred without
# bootstrap replicates - rerun IQ-TREE with -B 1000 for that family, or point
# PANELS at its .contree file.

suppressPackageStartupMessages(library(ape))

## ─────────────────────────── settings ──────────────────────────────────────
PANELS <- list(
  list(tag = "A", name = "DNA glycosylase MutY",        col = "#2F5C8A",
       file = "step13_phylo/DNA_glycosylase_MutY_tree.treefile"),
  list(tag = "B", name = "Ferredoxin 4Fe4S",      col = "#B5533C",
       file = "step13_phylo/Ferredoxin_4Fe4S_tree.treefile"),
  list(tag = "C", name = "Radical SAM", col = "#2E7D5B",
       file = "step13_phylo/Radical_SAM_tree.treefile"),
  list(tag = "D", name = "WhiB",        col = "#6B4E9B",
       file = "step13_phylo/WhiB_tree.treefile")
)

SUPPORT   <- 90      # UFBoot threshold for a filled dot
DOT_CEX   <- 0.45    # dot size, identical in every panel
EDGE_LWD  <- 0.45    # branch width, identical in every panel
OPEN_ANG  <- 8       # gap in the fan, degrees
OUT       <- "Figure6"

## ─────────────────────────── helpers ───────────────────────────────────────
read_tree_safe <- function(p) {
  if (!file.exists(p)) stop("missing file: ", p)
  tr <- read.tree(p)
  if (inherits(tr, "multiPhylo")) tr <- tr[[1]]
  tr
}

support_of <- function(tr) {
  lab <- tr$node.label
  if (is.null(lab) || !length(lab)) return(rep(NA_real_, tr$Nnode))
  s <- suppressWarnings(as.numeric(lab))
  # IQ-TREE sometimes writes "SH-aLRT/UFBoot"; take the last field
  if (all(is.na(s)) && any(grepl("/", lab))) {
    s <- suppressWarnings(as.numeric(sub(".*/", "", lab)))
  }
  # a 0-1 scale is a proportion, not a percentage
  if (sum(!is.na(s)) && max(s, na.rm = TRUE) <= 1) s <- s * 100
  s
}

draw_panel <- function(p) {
  tr <- read_tree_safe(p$file)
  tr <- ladderize(tr)
  sup <- support_of(tr)
  n_ok <- sum(!is.na(sup))
  n_hi <- sum(sup >= SUPPORT, na.rm = TRUE)

  par(mar = c(0.2, 0.2, 1.6, 0.2), xpd = NA)
  plot.phylo(tr, type = "fan", open.angle = OPEN_ANG,
             show.tip.label = FALSE, edge.color = p$col,
             edge.width = EDGE_LWD, no.margin = FALSE)

  if (n_hi > 0) {
    nodelabels(node = which(sup >= SUPPORT) + Ntip(tr),
               pch = 21, cex = DOT_CEX, bg = "black", col = "black")
  }

  title(main = sprintf("%s (n=%d)", p$name, Ntip(tr)),
        cex.main = 0.95, font.main = 2, line = 0.35)
  mtext(p$tag, side = 3, adj = 0, line = 0.1, cex = 1.15, font = 2)

  # scale bar, same style in every panel
  try(add.scale.bar(lwd = 0.8, cex = 0.55, length = NULL), silent = TRUE)

  invisible(list(tag = p$tag, name = p$name, tips = Ntip(tr),
                 nodes = tr$Nnode, with_support = n_ok, above = n_hi))
}

render <- function(dev_open, file) {
  dev_open(file)
  par(mfrow = c(2, 2), oma = c(0.4, 0.4, 0.4, 0.4))
  out <- lapply(PANELS, draw_panel)
  dev.off()
  invisible(out)
}

## ─────────────────────────── render ────────────────────────────────────────
W <- 7.2; H <- 7.0

info <- render(function(f) pdf(f, width = W, height = H, useDingbats = FALSE),
               paste0(OUT, ".pdf"))

render(function(f) tiff(f, width = W, height = H, units = "in", res = 600,
                        compression = "lzw", type = "cairo"),
       paste0(OUT, ".tiff"))

render(function(f) png(f, width = W, height = H, units = "in", res = 400,
                       type = "cairo"),
       paste0(OUT, ".png"))

## ─────────────────────────── report ────────────────────────────────────────
cat(sprintf("\nsupport threshold: UFBoot >= %d\n\n", SUPPORT))
for (x in info) {
  cat(sprintf("  %s  %-12s tips %4d  internal nodes %4d  ",
              x$tag, x$name, x$tips, x$nodes))
  if (x$with_support == 0) {
    cat("NO SUPPORT VALUES IN FILE\n")
  } else {
    cat(sprintf("support on %4d nodes, %4d at >= %d (%.0f%%)\n",
                x$with_support, x$above, SUPPORT,
                100 * x$above / x$with_support))
  }
}
bad <- Filter(function(x) x$with_support == 0, info)
if (length(bad)) {
  cat("\nPanels without support values:",
      paste(sapply(bad, function(x) x$name), collapse = ", "), "\n")
  cat("Those trees were inferred without bootstrap replicates. Rerun:\n")
  cat("   iqtree2 -s <alignment> -m MFP -B 1000 -T AUTO --prefix <family>\n")
  cat("and point PANELS at the new .treefile, or use the .contree if present.\n")
}
cat("\nwrote Figure6.pdf, Figure6.tiff (600 dpi LZW), Figure6.png\n")

