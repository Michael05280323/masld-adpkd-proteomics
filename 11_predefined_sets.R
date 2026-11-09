## ============================================================================
## 11_predefined_sets.R  (v1.3)
## FDR-defined core and predefined up-side/M6 KEGG gene-set analyses.
## Run from the repository root, after 01, 02, 06 and 00 (KEGG tables).
## Outputs:
##   results/supp_table_FDR_core.csv
##   results/supp_table_predefined_sets.csv
## ============================================================================
suppressPackageStartupMessages({library(limma); library(org.Hs.eg.db)})

lt <- read.csv("results/liver_limma_full.csv", stringsAsFactors = FALSE)
kt <- read.csv("results/kidney_limma_full.csv", stringsAsFactors = FALSE)
gt <- read.csv("results/gse96804_limma.csv", stringsAsFactors = FALSE)

## ---------------- FDR core ----------------
J <- merge(lt[, c("gene", "logFC", "adj.P.Val")],
           kt[, c("gene", "logFC", "adj.P.Val")],
           by = "gene", suffixes = c(".liver", ".kidney"))
core <- J[J$adj.P.Val.liver < 0.05 & J$adj.P.Val.kidney < 0.05 &
          sign(J$logFC.liver) == sign(J$logFC.kidney) &
          abs(J$logFC.liver) >= log2(1.5) & abs(J$logFC.kidney) >= log2(1.5), ]
core <- core[order(-core$logFC.liver), ]
write.csv(core, "results/supp_table_FDR_core.csv", row.names = FALSE)
cat(sprintf("FDR core: %d (up %d / down %d)\n",
            nrow(core), sum(core$logFC.liver > 0), sum(core$logFC.liver < 0)))

## ---------------- predefined sets ----------------
lk <- read.delim("results/kegg_link_hsa.tsv", header = FALSE, stringsAsFactors = FALSE)
nm <- read.delim("results/kegg_pathway_hsa.tsv", header = FALSE, stringsAsFactors = FALSE)
colnames(lk) <- c("path", "gene"); colnames(nm) <- c("path", "desc")
lk$gene <- sub("^hsa:", "", lk$gene); lk$path <- sub("^path:", "", lk$path)
nm$desc <- sub(" - Homo sapiens \\(human\\)$", "", nm$desc)
t2g <- merge(lk, nm, by = "path")

sets <- data.frame(
  set = c("OXPHOS", "Focal adhesion", "Regulation of actin cytoskeleton",
          "ECM-receptor interaction", "Cytoskeleton in muscle cells",
          "Val/Leu/Ile degradation", "Propanoate metabolism",
          "Glyoxylate/dicarboxylate metabolism", "Butanoate metabolism",
          "TCA cycle", "Fatty acid degradation"),
  kegg = c("Oxidative phosphorylation", "Focal adhesion", "Regulation of actin cytoskeleton",
           "ECM-receptor interaction", "Cytoskeleton in muscle cells",
           "Valine, leucine and isoleucine degradation", "Propanoate metabolism",
           "Glyoxylate and dicarboxylate metabolism", "Butanoate metabolism",
           "Citrate cycle (TCA cycle)", "Fatty acid degradation"),
  direction = c("down", "up", "up", "up", "up", "down", "down", "down", "down", "down", "down"),
  stringsAsFactors = FALSE
)

res <- list()
for (i in seq_len(nrow(sets))) {
  entrez <- t2g$gene[t2g$desc == sets$kegg[i]]
  sym <- unique(na.omit(mapIds(org.Hs.eg.db, keys = entrez, column = "SYMBOL",
                               keytype = "ENTREZID", multiVals = "first")))
  for (dname in c("liver", "kidney", "GSE96804")) {
    tt <- switch(dname, liver = lt, kidney = kt, GSE96804 = gt)
    idx <- which(tt$gene %in% sym)
    if (!length(idx)) next
    a <- tt$logFC[idx]; b <- tt$logFC[-idx]
    alt <- sets$direction[i]
    gs <- geneSetTest(idx, tt$t, alternative = alt)
    wx <- wilcox.test(a, b)$p.value
    res[[length(res) + 1]] <- data.frame(
      set = sets$set[i], dataset = dname, direction = alt,
      n_genes = length(idx),
      median_logFC_set = median(a), median_logFC_other = median(b),
      pct_direction = if (alt == "down") 100 * mean(a < 0) else 100 * mean(a > 0),
      pct_direction_other = if (alt == "down") 100 * mean(b < 0) else 100 * mean(b > 0),
      geneSetTest_P = gs, wilcox_two_sided_P = wx,
      stringsAsFactors = FALSE
    )
  }
}
sup <- do.call(rbind, res)
write.csv(sup, "results/supp_table_predefined_sets.csv", row.names = FALSE)
print(sup, row.names = FALSE)
