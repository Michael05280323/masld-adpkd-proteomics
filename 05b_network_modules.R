## ============================================================================
## STRING PPI network + Louvain modules  (manuscript 1.5 / 2.5)
## Input : results/joint_concordant.csv   (from 03_joint.R)
## Output: results/string_network.tsv
##         results/louvain_modules.csv
## ============================================================================
suppressPackageStartupMessages({ library(data.table); library(igraph) })

OUT <- "results"
conc <- fread(file.path(OUT, "joint_concordant.csv"))
cat("concordant proteins:", nrow(conc), "\n")

## ---------------- STRING network (GET, combined score >= 0.4) ----------------
st_file <- file.path(OUT, "string_network.tsv")
if (!file.exists(st_file) || file.info(st_file)$size < 1000) {
  ids <- paste(conc$gene, collapse = "%0d")           # STRING accepts %0d separators
  url <- paste0("https://string-db.org/api/tsv/network?identifiers=", ids,
                "&species=9606&required_score=400&caller_identity=masld_adpkd")
  cat("downloading STRING network...\n")
  ok <- tryCatch({ download.file(url, st_file, mode = "wb", quiet = TRUE); TRUE },
                 error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); FALSE })
  if (!ok || file.info(st_file)$size < 1000)
    stop("STRING download failed. Download the TSV manually from ",
         "https://string-db.org/ and save it as ", st_file)
}
st <- fread(st_file)
el <- unique(st[, .(from = preferredName_A, to = preferredName_B, score = score)])
cat("edges:", nrow(el), "\n")

G <- graph_from_data_frame(el[, .(from, to)], directed = FALSE,
                           vertices = data.table(gene = unique(c(el$from, el$to))))
E(G)$weight <- el$score
cat("nodes:", vcount(G), "\n")

deg <- sort(degree(G), decreasing = TRUE)
cat("\ntop 16 by degree:\n"); print(head(deg, 16))

## ---------------- Louvain modules ----------------
set.seed(1)
cl <- cluster_louvain(G)
dirn <- setNames(ifelse(conc$lfc > 0, "UP", "DOWN"), conc$gene)
mem <- data.table(gene = V(G)$name, degree = as.integer(degree(G)),
                  mod = as.integer(membership(cl)))
mem[, dir := dirn[gene]]

mod <- mem[, .(n = .N,
               up = sum(dir == "UP",   na.rm = TRUE),
               down = sum(dir == "DOWN", na.rm = TRUE)), by = mod][order(-n)]
cat("\nLouvain modules:\n"); print(mod, row.names = FALSE)

for (m in head(mod$mod, 8)) {
  g <- mem[mod == m]$gene
  cat(sprintf("\n  module %d (n=%d, up=%d, down=%d)\n    %s\n",
              m, length(g), sum(mem[mod == m]$dir == "UP", na.rm = TRUE),
              sum(mem[mod == m]$dir == "DOWN", na.rm = TRUE),
              paste(head(sort(g), 12), collapse = ", ")))
}

fwrite(mem, file.path(OUT, "louvain_modules.csv"))
fwrite(el,  file.path(OUT, "string_edges.csv"))
cat("\nwritten: results/louvain_modules.csv, results/string_edges.csv\n")
cat("NOTE: module ids depend on the STRING version and the Louvain seed;\n",
    "      map them to M2/M3/M4/M6 by their representative genes before plotting.\n")
