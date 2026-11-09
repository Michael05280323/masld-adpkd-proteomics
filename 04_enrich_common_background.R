## ============================================================================
## v1.1 GO/KEGG ORA with the common quantified background (2851 genes;
## 2796 mapped to Entrez). Replaces the default full-genome background used
## in v1.0. Run after 03_joint.R; requires the KEGG tables in results/
## (downloaded by 00_download_data.R).
## ============================================================================
suppressPackageStartupMessages({
  library(clusterProfiler); library(org.Hs.eg.db); library(data.table)
})
OUT <- "results"

if (file.exists(file.path(OUT, "concordant_sets.rds"))) {
  S <- readRDS(file.path(OUT, "concordant_sets.rds"))
} else {
  conc <- fread(file.path(OUT, "joint_concordant.csv"))
  S <- list(up = conc[lfc > 0]$gene, dn = conc[lfc < 0]$gene)
  saveRDS(S, file.path(OUT, "concordant_sets.rds"))
  cat("created results/concordant_sets.rds from joint_concordant.csv\n")
}
J <- fread(file.path(OUT, "joint_all_common.csv"))
ent <- function(s) unique(na.omit(bitr(s, "SYMBOL", "ENTREZID", org.Hs.eg.db)$ENTREZID))
eUp <- ent(S$up); eDn <- ent(S$dn)
universe <- ent(J$gene)
cat("input up/down:", length(eUp), "/", length(eDn),
    " | common universe:", length(universe), "\n")

## ---------------- GO ----------------
go_all <- list()
for (ont in c("BP", "CC", "MF")) {
  for (dir in c("DOWN", "UP")) {
    eg <- if (dir == "DOWN") eDn else eUp
    r <- enrichGO(eg, OrgDb = org.Hs.eg.db, ont = ont, universe = universe,
                  pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
                  minGSSize = 10, readable = TRUE)
    d <- as.data.frame(r)
    if (nrow(d)) {
      d$ontology <- ont; d$direction <- dir
      go_all[[paste(ont, dir)]] <- d[, c("ontology","direction","ID","Description",
                                          "GeneRatio","BgRatio","pvalue","p.adjust","qvalue","Count")]
    }
  }
}
go <- do.call(rbind, go_all); rownames(go) <- NULL
fwrite(go, file.path(OUT, "ora_common_GO.csv"))
cat("GO rows:", nrow(go), "\n")

## ---------------- KEGG ----------------
lk <- fread(file.path(OUT, "kegg_link_hsa.tsv"), header = FALSE)
nm <- fread(file.path(OUT, "kegg_pathway_hsa.tsv"), header = FALSE)
setnames(lk, c("path","gene")); setnames(nm, c("path","desc"))
lk[, `:=`(gene = sub("^hsa:", "", gene), path = sub("^path:", "", path))]
nm[, desc := sub(" - Homo sapiens \\(human\\)$", "", desc)]
t2g <- unique(merge(lk, nm, by = "path")[, .(term = desc, gene)])

kegg_all <- list()
for (dir in c("DOWN", "UP")) {
  eg <- if (dir == "DOWN") eDn else eUp
  r <- enricher(eg, TERM2GENE = t2g, universe = universe,
                pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 10)
  d <- as.data.frame(r)
  if (nrow(d)) {
    d$direction <- dir
    kegg_all[[dir]] <- d[, c("direction","ID","Description","GeneRatio","BgRatio",
                             "pvalue","p.adjust","qvalue","Count")]
  }
}
kegg <- do.call(rbind, kegg_all); rownames(kegg) <- NULL
fwrite(kegg, file.path(OUT, "ora_common_KEGG.csv"))
cat("KEGG rows:", nrow(kegg), "\n")

## ---------------- cited terms ----------------
terms <- c("cellular respiration", "aerobic respiration",
           "mitochondrial electron transport, NADH to ubiquinone",
           "oxidative phosphorylation")
cat("\n--- GO-BP DOWN cited terms (common background) ---\n")
print(go[go$direction == "DOWN" & go$ontology == "BP" & go$Description %in% terms,
         c("Description","GeneRatio","BgRatio","pvalue","p.adjust","Count")],
      row.names = FALSE)
cat("\n--- KEGG DOWN cited pathways (common background) ---\n")
print(kegg[kegg$direction == "DOWN" &
           grepl("Oxidative phosphorylation|fatty liver", kegg$Description, ignore.case = TRUE),
           c("Description","GeneRatio","BgRatio","pvalue","p.adjust","Count")],
      row.names = FALSE)
