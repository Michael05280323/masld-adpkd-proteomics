## ============================================================================
## 09_permutation.R  (v1.3)
## Kidney-label permutation null for the 362 direction-concordant proteins.
## Run from the repository root, after 01-03. Requires results/kidney_mat.rds
## (created by 02_kidney.R) and results/liver_limma_full.csv (01_liver.R).
## ============================================================================
suppressPackageStartupMessages(library(limma))

lt <- read.csv("results/liver_limma_full.csv", stringsAsFactors = FALSE)
lt <- lt[, c("gene", "logFC", "P.Value")]

K <- readRDS("results/kidney_mat.rds")
mi <- K$mi; grp <- K$grp

## gene-level collapse of duplicate symbols (same rule as 02_kidney.R)
gn <- rownames(mi)
if (anyDuplicated(gn)) {
  num <- rowsum(replace(mi, is.na(mi), 0), gn, reorder = TRUE)
  den <- rowsum((!is.na(mi)) * 1, gn, reorder = TRUE)
  mi <- num / den; mi[den == 0] <- NA
}

A <- colnames(mi)[grp == "ADPKD"]; N <- colnames(mi)[grp == "normal"]
keep <- rowSums(!is.na(mi[, A, drop = FALSE])) >= ceiling(.5 * length(A)) &
        rowSums(!is.na(mi[, N, drop = FALSE])) >= ceiling(.5 * length(N))
x <- mi[keep, , drop = FALSE]; x[x <= 0] <- NA; x <- log2(x)
for (gg in unique(grp)) {
  j <- which(grp == gg); s <- x[, j, drop = FALSE]
  g0 <- apply(s, 1, median, na.rm = TRUE); g0[!is.finite(g0)] <- median(s, na.rm = TRUE)
  for (k in j) { v <- x[, k]; v[is.na(v)] <- g0[is.na(v)]; x[, k] <- v }
}
x <- sweep(x, 2, apply(x, 2, median, na.rm = TRUE), "-")

common <- intersect(lt$gene, rownames(x))
lt <- lt[lt$gene %in% common, ]
cat("common genes:", length(common), "| kidney matrix:", nrow(x), "x", ncol(x), "\n")

count_conc <- function(labels) {
  d <- model.matrix(~0 + factor(labels, levels = c("normal", "ADPKD")))
  colnames(d) <- c("ctrl", "dis")
  fit <- eBayes(contrasts.fit(lmFit(x, d), makeContrasts(dis - ctrl, levels = d)), trend = TRUE)
  tt <- topTable(fit, number = Inf, sort.by = "none"); tt$gene <- rownames(tt)
  J <- merge(lt, tt[, c("gene", "logFC", "P.Value")], by = "gene",
             suffixes = c(".l", ".k"))
  sum(J$P.Value.l < 0.05 & J$P.Value.k < 0.05 & sign(J$logFC.l) == sign(J$logFC.k))
}

obs <- count_conc(grp)
cat("observed concordant:", obs, "\n")
set.seed(20261104)
B <- 2000
null <- numeric(B)
for (b in 1:B) {
  null[b] <- count_conc(sample(grp))
  if (b %% 500 == 0) cat("perm", b, "\n")
}
summ <- data.frame(
  observed = obs,
  null_median = median(null),
  null_q1 = quantile(null, .25, names = FALSE),
  null_q3 = quantile(null, .75, names = FALSE),
  null_mean = mean(null),
  null_max = max(null),
  n_perm = B,
  empirical_P = (1 + sum(null >= obs)) / (1 + B),
  fold_over_median = obs / median(null)
)
print(summ, row.names = FALSE)
write.csv(summ, "results/supp_table_permutation.csv", row.names = FALSE)
saveRDS(null, "results/concordance_null_2000.rds")
