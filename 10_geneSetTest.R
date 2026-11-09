## ============================================================================
## 10_geneSetTest.R  (v1.3)
## Correlation-aware geneSetTest for the predefined KEGG OXPHOS set (hsa00190).
## Run from the repository root, after 01, 02, 06 and 00 (KEGG tables).
## ============================================================================
suppressPackageStartupMessages({library(limma); library(org.Hs.eg.db)})

lk <- read.delim("results/kegg_link_hsa.tsv", header = FALSE, stringsAsFactors = FALSE)
nm <- read.delim("results/kegg_pathway_hsa.tsv", header = FALSE, stringsAsFactors = FALSE)
colnames(lk) <- c("path", "gene"); colnames(nm) <- c("path", "desc")
lk$gene <- sub("^hsa:", "", lk$gene); lk$path <- sub("^path:", "", lk$path)
nm$desc <- sub(" - Homo sapiens \\(human\\)$", "", nm$desc)
ox <- merge(lk, nm, by = "path")
ox <- ox[grepl("Oxidative phosphorylation", ox$desc), ]
ox_sym <- unique(na.omit(mapIds(org.Hs.eg.db, keys = ox$gene, column = "SYMBOL",
                                keytype = "ENTREZID", multiVals = "first")))
cat("OXPHOS symbols:", length(ox_sym), "\n")

files <- c(liver = "results/liver_limma_full.csv",
           kidney = "results/kidney_limma_full.csv",
           GSE96804 = "results/gse96804_limma.csv")
res <- list()
for (nmx in names(files)) {
  tt <- read.csv(files[[nmx]], stringsAsFactors = FALSE)
  idx <- which(tt$gene %in% ox_sym)
  gs <- geneSetTest(idx, tt$t, alternative = "down")
  wx <- wilcox.test(tt$logFC[idx], tt$logFC[-idx])$p.value
  res[[nmx]] <- data.frame(
    dataset = nmx,
    n = length(idx),
    median_logFC_set = median(tt$logFC[idx]),
    median_logFC_other = median(tt$logFC[-idx]),
    pct_down = 100 * mean(tt$logFC[idx] < 0),
    pct_down_other = 100 * mean(tt$logFC[-idx] < 0),
    geneSetTest_P = gs,
    wilcox_two_sided_P = wx,
    stringsAsFactors = FALSE
  )
}
out <- do.call(rbind, res)
print(out, row.names = FALSE)
write.csv(out, "results/supp_table_genesettest.csv", row.names = FALSE)
