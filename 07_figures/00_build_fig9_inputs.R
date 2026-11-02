## Build the two generic input tables used by fig9_ppi_panels.R
## (demo version, built from the reproduction outputs; replace with your own tables)
suppressPackageStartupMessages({library(data.table); library(igraph)})
OUT <- "results"
FIG <- "results/figures"
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

mm <- fread(file.path(OUT, "louvain_modules.csv"))
st <- fread(file.path(OUT, "string_network.tsv"))
el <- unique(st[, .(from = preferredName_A, to = preferredName_B, score = score)])

g  <- graph_from_data_frame(el, directed = FALSE)
dg <- degree(g)

## map the four main Louvain communities to the manuscript's M-labels
map <- c(`4` = "M4", `6` = "M6", `3` = "M3", `2` = "M2")
mm[, module := fifelse(as.character(mod) %in% names(map),
                       map[as.character(mod)], "other")]
## single-gene functional-annotation correction: PXDN (peroxidasin, cross-links
## collagen IV) is a secreted ECM protein; Louvain placed it in M3, move to M2.
if (mm[gene == "PXDN", module] == "M3") {
  mm[gene == "PXDN", module := "M2"]
  cat("manual annotation: PXDN M3 -> M2\n")
}
mm[, degree := dg[gene]]; mm[is.na(degree), degree := 0]
mm[, direction := fifelse(dir == "UP", "UP", "DOWN")]

lab_genes <- c("SDHA","HINT2","GPT","SUCLA2","HOGA1","ECHDC2",
               "FLNA","DPYSL3","S100A10","LUM","PXDN","SERPINF1","FBLN5",
               "ANXA2","DBT","NDUFS3","NDUFS8","NDUFS2","HSPA9")
mm[, label := as.integer(gene %in% lab_genes)]

fwrite(mm[, .(gene, module, direction, degree, label)], file.path(FIG, "nodes.csv"))
fwrite(el, file.path(FIG, "edges.csv"))
cat("nodes:", nrow(mm), " edges:", nrow(el), "\n")
print(mm[, .N, by = module])
