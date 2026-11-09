## ============================================================================
## v1.1 corrected KEGG + direction-aware cross-cohort AUC
## Fixes relative to v1.0:
##   (1) liver gene names are split by ";" before collapse (same rule as 01_liver.R);
##   (2) kidney duplicate gene symbols are collapsed before filtering/AUC;
##   (3) KEGG ORA uses the common quantified background (2851 genes) and
##       minGSSize = 10, matching the revised Methods.
## Run from the repository root, after 00-03 and 04.
## ============================================================================
suppressPackageStartupMessages({
  library(data.table); library(clusterProfiler); library(org.Hs.eg.db); library(pROC)
})
DATA <- "data"; OUT <- "results"

if (file.exists(file.path(OUT, "concordant_sets.rds"))) {
  S <- readRDS(file.path(OUT, "concordant_sets.rds"))
} else {
  conc <- fread(file.path(OUT, "joint_concordant.csv"))
  S <- list(up = conc[lfc > 0]$gene, dn = conc[lfc < 0]$gene)
  saveRDS(S, file.path(OUT, "concordant_sets.rds"))
}
up <- S$up; dn <- S$dn
J  <- fread(file.path(OUT, "joint_all_common.csv"))
sgn <- setNames(sign(J$kfc), J$gene)

## ---------------- KEGG with common-gene background ----------------
cat("#################  KEGG  #################\n")
lk <- fread(file.path(OUT, "kegg_link_hsa.tsv"), header = FALSE, sep = "\t")
nm <- fread(file.path(OUT, "kegg_pathway_hsa.tsv"), header = FALSE, sep = "\t")
setnames(lk, c("path", "gene")); setnames(nm, c("path", "desc"))
lk[, `:=`(gene = sub("^hsa:", "", gene), path = sub("^path:", "", path))]
nm[, desc := sub(" - Homo sapiens \\(human\\)$", "", desc)]
t2g <- unique(merge(lk, nm, by = "path")[, .(term = desc, gene)])
cat("### TERM2GENE:", nrow(t2g), "pairs /", uniqueN(t2g$term), "pathways\n")

ent <- function(s) unique(na.omit(bitr(s, "SYMBOL", "ENTREZID", org.Hs.eg.db)$ENTREZID))
universe <- ent(J$gene)          # common quantified background (2851 genes)
cat("### universe (Entrez):", length(universe), "\n")

kegg <- function(s, lab) {
  r <- enricher(ent(s), TERM2GENE = t2g, universe = universe,
                pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 10)
  if (is.null(r)) { cat(lab, ": none\n"); return(invisible(NULL)) }
  d <- as.data.frame(r)
  cat(sprintf("\n--- KEGG %s : %d pathways, top 8 ---\n", lab, nrow(d)))
  print(head(d[, c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust","Count")], 8),
        row.names = FALSE)
  h <- d[grepl("Oxidative phosphorylation|fatty liver", d$Description, ignore.case = TRUE),
         c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust","Count")]
  if (nrow(h)) { cat("   >> cited pathways:\n"); print(h, row.names = FALSE) }
  invisible(d)
}
kegg(dn, "DOWN"); kegg(up, "UP")

## ---------------- liver matrix (split multi-symbol first) ----------------
cat("\n\n#################  CROSS-COHORT AUC  #################\n")
Q <- fread(file.path(DATA, "PXD051911__liver_protein_quant.txt"),
           na.strings = c("NA","","NaN"))
M <- fread(file.path(DATA, "PXD051911__meta_data.txt"),
           na.strings = c("NA","","NaN"))
sm <- setdiff(names(Q), c("ProteinAccessions","Genes","ProteinDescriptions"))
gm <- M$saf_diagnosis; names(gm) <- M$liver_proteomics_filename
grp <- gm[sm]; names(grp) <- sm

Q[, Genes := trimws(as.character(Genes))]
## v1.1 FIX: take the first symbol of multi-symbol entries
Q[, Genes := sapply(strsplit(Genes, ";", fixed = TRUE), function(x) trimws(x[1]))]
Q <- Q[!is.na(Genes) & Genes != ""]

ag <- Q[, c("Genes", sm), with = FALSE][, lapply(.SD, function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}), by = Genes, .SDcols = sm]
LM <- as.matrix(ag[, -1, with = FALSE]); rownames(LM) <- ag$Genes
LM[LM <= 0] <- NA; LM <- log2(LM)
k <- rowSums(!is.na(LM[, grp == "MASH", drop = FALSE])) >= ceiling(.7 * sum(grp == "MASH")) &
     rowSums(!is.na(LM[, grp == "No_MASLD", drop = FALSE])) >= ceiling(.7 * sum(grp == "No_MASLD"))
LM <- LM[k, , drop = FALSE]
for (gg in unique(grp)) {
  i <- which(grp == gg); s <- LM[, i, drop = FALSE]
  g0 <- apply(s, 1, median, na.rm = TRUE); g0[!is.finite(g0)] <- median(s, na.rm = TRUE)
  for (j in seq_along(i)) { v <- s[, j]; v[is.na(v)] <- g0[is.na(v)]; s[, j] <- v }
  LM[, i] <- s
}
LMn <- sweep(LM, 2, apply(LM, 2, median, na.rm = TRUE), "-")
lLab <- factor(grp, levels = c("No_MASLD","MASH"))

## ---------------- kidney matrix: collapse duplicates first ----------------
L <- readRDS(file.path(OUT, "kidney_mat.rds")); mi <- L$mi; kgrp <- L$grp
gn <- rownames(mi)
if (anyDuplicated(gn)) {
  num <- rowsum(replace(mi, is.na(mi), 0), gn, reorder = TRUE)
  den <- rowsum((!is.na(mi)) * 1, gn, reorder = TRUE)
  mi <- num / den; mi[den == 0] <- NA
}
A <- colnames(mi)[kgrp == "ADPKD"]; N <- colnames(mi)[kgrp == "normal"]
KM <- mi[rowSums(!is.na(mi[, A, drop = FALSE])) >= ceiling(.5 * length(A)) &
         rowSums(!is.na(mi[, N, drop = FALSE])) >= ceiling(.5 * length(N)), , drop = FALSE]
KM[KM <= 0] <- NA; KM <- log2(KM)
for (gg in unique(kgrp)) {
  i <- which(kgrp == gg); s <- KM[, i, drop = FALSE]
  g0 <- apply(s, 1, median, na.rm = TRUE); g0[!is.finite(g0)] <- median(s, na.rm = TRUE)
  for (j in seq_along(i)) { v <- s[, j]; v[is.na(v)] <- g0[is.na(v)]; s[, j] <- v }
  KM[, i] <- s
}
kLab <- factor(kgrp, levels = c("normal","ADPKD"))
cat("### liver:", dim(LMn), " kidney (un-normalised):", dim(KM), "\n")

## ---------------- AUC ----------------
one <- function(v, lab) as.numeric(auc(roc(lab, v, quiet = TRUE, direction = "<")))
orient <- function(a, s) if (s >= 0) a else 1 - a
allg <- intersect(rownames(LMn), rownames(KM))
cat("### AUC genes:", length(allg), "\n")
AU <- data.table(gene = allg,
  liverAUC  = sapply(allg, function(g) orient(one(LMn[g, ], lLab), sgn[[g]])),
  kidneyAUC = sapply(allg, function(g) orient(one(KM[g, ],  kLab), sgn[[g]])))
AU[, minAUC := pmin(liverAUC, kidneyAUC)]
fwrite(AU, file.path(OUT, "auc_all_shared.csv"))

mol <- c("HINT2","GPT","ECHDC2","DBT","DPYSL3","S100A10","ANXA2")
ms  <- data.table(gene = mol,
                  ms_liver = c(.95,.93,.92,.95,.96,.94,.99),
                  ms_kidney = c(.96,1.00,.89,.78,.87,.81,.65))
cmp <- merge(ms, AU, by = "gene", all.x = TRUE)
cmp[, `:=`(d_liver = round(liverAUC - ms_liver, 2),
           d_kidney = round(kidneyAUC - ms_kidney, 2))]
cat("\n--- Table 3 check ---\n"); print(cmp, row.names = FALSE)

for (nmx in c("KM")) {
  X <- get(nmx)
  f1 <- mean(sapply(rownames(X), function(g) { a <- one(X[g, ], kLab); max(a, 1 - a) >= 0.9999 }))
  cat(sprintf("### %s: AUC==1.00 %.2f%% (%d/%d)\n", nmx, 100 * f1, round(f1 * nrow(X)), nrow(X)))
}
