## ============================================================================
##  Figure 9  —  STRING PPI network, 4-panel module view
##  Input : nodes.csv  (gene, module, direction, degree, label)
##          edges.csv  (from, to, score)
##  Output: fig9_ppi_panels.pdf / .png  (2 x 2 panels, intra-module edges only)
##
##  Design rules (see spec):
##   * each panel shows ONE module and ONLY its intra-module edges
##   * node fill = module colour, size proportional to STRING degree
##   * only flagged genes are labelled, with leader lines (ggrepel)
##   * identical layout algorithm + seed for every panel
## ============================================================================
suppressPackageStartupMessages({
  library(data.table); library(igraph); library(ggplot2)
  library(ggrepel); library(patchwork); library(scales); library(systemfonts)
})

FIG <- "results/figures"

## ---------------------------------------------------------------- parameters
CFG <- list(
  modules   = c("M4", "M6", "M3", "M2"),          # top-left, top-right, bottom-left, bottom-right
  title_cn  = c(M4 = "M4 线粒体氧化磷酸化", M6 = "M6 分化代谢丢失",
                M3 = "M3 骨架/黏附重塑",   M2 = "M2 分泌型细胞外基质"),
  title_en  = c(M4 = "M4 Mitochondrial OXPHOS", M6 = "M6 Differentiation-metabolism loss",
                M3 = "M3 Cytoskeleton / adhesion", M2 = "M2 Secreted ECM"),
  arrow     = c(M4 = "↓", M6 = "↓", M3 = "↑", M2 = "↑"),
  use_cn    = TRUE,                    # FALSE -> English panel titles
  palette   = c(M4 = "#08306B", M6 = "#2171B5",
                M3 = "#B2182B", M2 = "#EF8A62",
                other = "#BDBDBD"),
  edge_alpha = c(0.08, 0.42),      # alpha scaled by STRING score
  edge_width = c(0.14, 0.42),
  node_range = c(0.8, 3.2),         # smaller: dense modules stop merging into a blob
  node_alpha = 0.90,
  layout_k   = 6,                   # top-k edges per node used for the LAYOUT only
  layout_method = "kk",             # "kk" (Kamada-Kawai) or "fr"
  quant_clip = c(0.02, 0.98),       # robust rescaling: stops peripheral nodes squeezing the core
  edge_keep_q = 0.40,               # draw only the top 60% of intra-module edges
  label_size = 2.30,
  drop_isolates = TRUE,
  seed = 1,
  out_w_mm = 180, out_h_mm = 168, dpi = 600
)

## Optional manual fine-tuning of individual labels (normalised panel units).
## Leave empty and let ggrepel do the work; add rows only for stubborn clashes.
NUDGE <- data.table(
  gene = c("NDUFS8", "NDUFS2", "HSPA9", "SDHA", "ANXA2", "FLNA", "S100A10", "ECHDC2", "SUCLA2"),
  dx   = c( 0.16,     0.16,    -0.20,  -0.10,   -0.16,   -0.12,  -0.16,     -0.10,   -0.12),
  dy   = c( 0.10,    -0.14,     0.10,  -0.10,    0.14,   -0.16,   0.00,     -0.12,    0.08)
)

## ---------------------------------------------------------------- fonts
pick_font <- function(cands) {
  f <- systemfonts::system_fonts()$family
  hit <- cands[cands %in% f]
  if (length(hit)) hit[1] else "sans"
}
CJK <- pick_font(c("Microsoft YaHei", "SimHei", "Noto Sans CJK SC",
                   "Source Han Sans SC", "WenQuanYi Zen Hei", "Arial Unicode MS"))
FONT <- if (CFG$use_cn) CJK else "Arial"
cat("font in use:", FONT, "\n")

## ---------------------------------------------------------------- data
nodes <- fread(file.path(FIG, "nodes.csv"))
edges <- fread(file.path(FIG, "edges.csv"))
nodes[, module := factor(module, levels = c(CFG$modules, "other"))]

## ---------------------------------------------------------------- one panel
make_panel <- function(mod) {
  e   <- edges[from %in% nodes[module == mod]$gene &
               to   %in% nodes[module == mod]$gene]
  e   <- unique(e, by = c("from", "to"))
  ## keep only nodes that appear in at least one intra-module edge
  keep <- if (CFG$drop_isolates) union(e$from, e$to) else nodes[module == mod]$gene
  nd   <- nodes[gene %in% keep & module == mod]
  e    <- e[from %in% keep & to %in% keep]
  n_drop <- length(unique(nodes[module == mod]$gene)) - nrow(nd)
  if (n_drop > 0) cat(sprintf("  %s: dropped %d isolated node(s)\n", mod, n_drop))

  ## ---- graph used for the LAYOUT: keep each node's top-k strongest edges ----
  ## (all intra-module edges are still DRAWN; this only prevents weak edges
  ##  from stretching the canvas and crushing the core)
  e_lay <- unique(rbind(
    e[order(-score), head(.SD, CFG$layout_k), by = from][, .(from, to, score)],
    e[order(-score), head(.SD, CFG$layout_k), by = to  ][, .(from, to, score)]))
  g <- graph_from_data_frame(e_lay[, .(from, to)], directed = FALSE,
                             vertices = data.table(gene = nd$gene))
  E(g)$weight <- e_lay$score
  set.seed(CFG$seed)
  lay <- if (CFG$layout_method == "kk")
    layout_with_kk(g, weights = E(g)$weight, maxiter = 2000)
  else layout_with_fr(g, niter = 5000, grid = "nogrid", weights = E(g)$weight)

  ## ---- robust rescaling to [0,1] on the 2nd-98th percentile ----
  ## so a handful of peripheral nodes no longer squeeze the core
  qs <- CFG$quant_clip
  for (j in 1:2) {
    q <- quantile(lay[, j], qs, na.rm = TRUE)
    if (diff(q) == 0) q <- range(lay[, j])
    lay[, j] <- (lay[, j] - q[1]) / diff(q)
  }
  lay[lay < -0.08] <- -0.08; lay[lay > 1.08] <- 1.08
  pos <- data.table(gene = V(g)$name, x = lay[, 1], y = lay[, 2])

  ## ---- edges actually drawn: M3 is too dense, keep only its top 60% ----
  ## (other panels keep ALL their intra-module edges; nodes are always kept)
  edge_note <- ""
  if (mod == "M3") {
    e <- e[score >= quantile(score, CFG$edge_keep_q)]
    edge_note <- "，仅示前60%边"
  }
  e2  <- merge(e, pos, by.x = "from", by.y = "gene")
  e2  <- merge(e2, pos, by.x = "to", by.y = "gene", suffixes = c("", ".to"))
  nd  <- merge(nd, pos, by = "gene")
  lab <- merge(nd[label == 1], NUDGE, by = "gene", all.x = TRUE)
  lab[is.na(dx), dx := 0]; lab[is.na(dy), dy := 0]
  ## edge weight shown in the drawing follows the STRING score
  e2[, a := scales::rescale(score, to = CFG$edge_alpha, from = c(0.4, 1))]
  e2[, w := scales::rescale(score, to = CFG$edge_width, from = c(0.4, 1))]

  tt <- if (CFG$use_cn) CFG$title_cn[[mod]] else CFG$title_en[[mod]]
  subtitle <- sprintf("%s  |  %d nodes · %d edges", CFG$arrow[[mod]],
                      nrow(nd), nrow(e2))

  ggplot() +
    geom_segment(data = e2,
                 aes(x = x, y = y, xend = x.to, yend = y.to,
                     alpha = a, linewidth = w),
                 colour = CFG$palette[[mod]], lineend = "round") +
    scale_alpha_identity() + scale_linewidth_identity() +
    geom_point(data = nd,
               aes(x = x, y = y, size = degree, fill = module),
               shape = 21, colour = "white", stroke = 0.45,
               alpha = CFG$node_alpha) +
    geom_text_repel(data = lab,
                    aes(x = x, y = y, label = gene),
                    nudge_x = lab$dx, nudge_y = lab$dy,
                    size = CFG$label_size, family = "Arial", colour = "grey10",
                    bg.colour = "white", bg.r = 0.14,
                    segment.size = 0.24, segment.colour = "grey40",
                    min.segment.length = 0, box.padding = 0.75,
                    point.padding = 0.40, force = 3.0, force_pull = 0.35,
                    max.iter = 50000, max.overlaps = Inf,
                    seed = CFG$seed) +
    scale_size_continuous(range = CFG$node_range, guide = "none") +
    scale_fill_manual(values = CFG$palette, guide = "none") +
    coord_fixed(xlim = c(-0.06, 1.06), ylim = c(-0.06, 1.06), clip = "off") +
    labs(title = tt, subtitle = subtitle) +
    theme_void(base_family = FONT) +
    theme(plot.title    = element_text(size = 8.2, face = "bold", hjust = 0),
          plot.subtitle = element_text(size = 6.2, colour = "grey35", hjust = 0),
          plot.margin   = margin(3, 3, 2, 3))
}

## ---------------------------------------------------------------- assemble
plots <- lapply(CFG$modules, make_panel)
names(plots) <- CFG$modules

## shared legend drawn once (bottom)
legend_df <- data.table(x = seq_along(CFG$modules), y = 1,
                        module = factor(CFG$modules, levels = CFG$modules))
## dummy panel: points are invisible, only the legend keys are drawn
legend <- ggplot(legend_df, aes(x, y, fill = module)) +
  geom_point(size = 3.2, shape = 21, colour = "white", stroke = 0.4, alpha = 0) +
  scale_fill_manual(
    values = CFG$palette[CFG$modules],
    labels = setNames(paste0(CFG$title_cn, "（", CFG$arrow, "）"), CFG$modules),
    name = NULL) +
  guides(fill = guide_legend(nrow = 2, byrow = FALSE,
                             override.aes = list(size = 3.2, alpha = 1))) +
  theme_void(base_family = FONT) +
  theme(legend.position = "bottom",
        legend.text  = element_text(size = 6.4),
        legend.key.height = unit(3.5, "mm"),
        legend.spacing.x = unit(1.2, "mm"))

fig <- wrap_plots(plots, ncol = 2) / legend +
  plot_layout(heights = c(1, 0.12)) +
  plot_annotation(
    title = "肝肾共享一致蛋白的 STRING PPI 网络（按功能模块拆分）",
    caption = paste0(
      "STRING v12, combined score ≥ 0.4；四个主要模块分别展示，仅保留模块内部连边；",
      "其余小模块与模块间连边未显示。\n",
      "节点大小 ∝ 连接度（degree）；蓝色系为下调模块，红色系为上调模块。"),
    theme = theme(
      plot.title   = element_text(size = 9, face = "bold", family = FONT,
                                  margin = margin(0, 0, 2, 0)),
      plot.caption = element_text(size = 5.6, colour = "grey35", family = FONT,
                                  hjust = 0, lineheight = 1.15)))

## ---------------------------------------------------------------- export
w <- CFG$out_w_mm / 25.4; h <- CFG$out_h_mm / 25.4
ggsave(file.path(FIG, "fig9_ppi_panels.pdf"), fig, width = w, height = h,
       device = cairo_pdf, bg = "white")
ggsave(file.path(FIG, "fig9_ppi_panels.png"), fig, width = w, height = h,
       dpi = CFG$dpi, bg = "white")
cat("written: fig9_ppi_panels.pdf / .png\n")
