## ============================================================================
## 12_liver_covariate.R  (v1.3)
## Liver MASH vs control adjusted for age, sex and BMI (29 complete cases).
## Run from the repository root, after 00. Rebuilds the liver matrix from raw.
## Outputs:
##   results/liver_limma_adjusted.csv
##   results/supp_table_liver_covariate.csv
## ============================================================================
suppressPackageStartupMessages({library(data.table); library(limma); library(org.Hs.eg.db)})

Q <- fread("data/PXD051911__liver_protein_quant.txt",
           na.strings = c("NA", "", "NaN"))
M <- fread("data/PXD051911__meta_data.txt",
           na.strings = c("NA", "", "NaN"))
sm <- setdiff(names(Q), c("ProteinAccessions", "Genes", "ProteinDescriptions"))
grp <- M$saf_diagnosis[match(sm, M$liver_proteomics_filename)]

Q[, Genes := trimws(as.character(Genes))]
Q[, Genes := sapply(strsplit(Genes, ";", fixed = TRUE), function(x) trimws(x[1]))]
Q <- Q[!is.na(Genes) & Genes != ""]
ag <- Q[, c("Genes", sm), with = FALSE][, lapply(.SD, function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}), by = Genes, .SDcols = sm]
LM <- as.matrix(ag[, -1, with = FALSE]); rownames(LM) <- ag$Genes
LM[LM <= 0] <- NA; LM <- log2(LM)

keep <- rowSums(!is.na(LM[, grp == "MASH", drop = FALSE])) >= ceiling(.7 * sum(grp == "MASH")) &
        rowSums(!is.na(LM[, grp == "No_MASLD", drop = FALSE])) >= ceiling(.7 * sum(grp == "No_MASLD"))
LM <- LM[keep, , drop = FALSE]
for (gg in unique(grp)) {
  j <- which(grp == gg); s <- LM[, j, drop = FALSE]
  g0 <- apply(s, 1, median, na.rm = TRUE); g0[!is.finite(g0)] <- median(s, na.rm = TRUE)
  for (k in j) { v <- LM[, k]; v[is.na(v)] <- g0[is.na(v)]; LM[, k] <- v }
}
LM <- sweep(LM, 2, apply(LM, 2, median, na.rm = TRUE), "-")

idx <- match(colnames(LM), M$liver_proteomics_filename)
cov <- data.frame(
  group = factor(grp, levels = c("No_MASLD", "MASH")),
  sex = factor(M$gender[idx]),
  age = as.numeric(M$alder[idx]),
  bmi = as.numeric(M$bmi[idx])
)
ok <- complete.cases(cov)
LM <- LM[, ok]; cov <- cov[ok, ]
cat("analysed samples:", ncol(LM), " covariates complete:", sum(ok), "\n")

d0 <- model.matrix(~ group, data = cov)
f0 <- eBayes(lmFit(LM, d0), trend = TRUE)
t0 <- topTable(f0, coef = "groupMASH", number = Inf, sort.by = "none"); t0$gene <- rownames(t0)

d1 <- model.matrix(~ group + sex + age + bmi, data = cov)
f1 <- eBayes(lmFit(LM, d1), trend = TRUE)
t1 <- topTable(f1, coef = "groupMASH", number = Inf, sort.by = "none"); t1$gene <- rownames(t1)

summ <- function(tt, lab) {
  dep <- tt[tt$P.Value < 0.05 & abs(tt$logFC) >= log2(1.5), ]
  fdr <- tt[tt$adj.P.Val < 0.05 & abs(tt$logFC) >= log2(1.5), ]
  cat(sprintf("%-12s DEP %4d (up %3d / down %3d) | FDR+1.5x %4d\n",
              lab, nrow(dep), sum(dep$logFC > 0), sum(dep$logFC < 0), nrow(fdr)))
  c(dep = nrow(dep), up = sum(dep$logFC > 0), down = sum(dep$logFC < 0), fdr = nrow(fdr))
}
s0 <- summ(t0, "unadjusted")
s1 <- summ(t1, "adjusted")

## OXPHOS geneSetTest in the adjusted model
lk <- read.delim("results/kegg_link_hsa.tsv", header = FALSE, stringsAsFactors = FALSE)
nm <- read.delim("results/kegg_pathway_hsa.tsv", header = FALSE, stringsAsFactors = FALSE)
colnames(lk) <- c("path", "gene"); colnames(nm) <- c("path", "desc")
lk$gene <- sub("^hsa:", "", lk$gene); lk$path <- sub("^path:", "", lk$path)
nm$desc <- sub(" - Homo sapiens \\(human\\)$", "", nm$desc)
ox <- merge(lk, nm, by = "path"); ox <- ox[grepl("Oxidative phosphorylation", ox$desc), ]
ox_sym <- unique(na.omit(mapIds(org.Hs.eg.db, keys = ox$gene, column = "SYMBOL",
                                keytype = "ENTREZID", multiVals = "first")))
idx_ox <- which(t1$gene %in% ox_sym)
gs <- geneSetTest(idx_ox, t1$t, alternative = "down")
cat(sprintf("adjusted OXPHOS: n=%d median=%.3f pct_down=%.0f%% geneSetTest P=%.3g\n",
            length(idx_ox), median(t1$logFC[idx_ox]), 100 * mean(t1$logFC[idx_ox] < 0), gs))

out <- data.frame(
  model = c("unadjusted (MASH vs control)", "adjusted for age + sex + BMI"),
  DEP = c(s0["dep"], s1["dep"]), up = c(s0["up"], s1["up"]),
  down = c(s0["down"], s1["down"]), FDR_1.5x = c(s0["fdr"], s1["fdr"]),
  OXPHOS_n = length(idx_ox),
  OXPHOS_median_logFC = median(t1$logFC[idx_ox]),
  OXPHOS_pct_down = 100 * mean(t1$logFC[idx_ox] < 0),
  geneSetTest_P = gs, row.names = NULL
)
write.csv(out, "results/supp_table_liver_covariate.csv", row.names = FALSE)
write.csv(t1, "results/liver_limma_adjusted.csv", row.names = FALSE)
