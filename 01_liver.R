## ============================================================
## Independent reproduction  --  LIVER side (PXD051911)
## Implemented strictly from the manuscript's Methods 1.2-1.3.
## ============================================================
suppressPackageStartupMessages({library(data.table); library(limma)})

DATA <- "data"
OUT  <- "results"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

Q <- fread(file.path(DATA, "PXD051911__liver_protein_quant.txt"), na.strings = c("NA","","NaN"))
M <- fread(file.path(DATA, "PXD051911__meta_data.txt"),        na.strings = c("NA","","NaN"))

cat("### raw quant matrix:", nrow(Q), "protein entries x", ncol(Q), "cols\n")

samples <- setdiff(names(Q), c("ProteinAccessions","Genes","ProteinDescriptions"))
cat("### sample columns:", length(samples), "\n")

gmap <- M$saf_diagnosis; names(gmap) <- M$liver_proteomics_filename
grp  <- gmap[samples]
cat("### grouping over quantified samples:\n"); print(table(grp, useNA="ifany"))

## ---- 1.2 collapse multiple protein entries per gene (mean) ----
Q[, Genes := trimws(as.character(Genes))]
Q[, Genes := sapply(strsplit(Genes, ";", fixed = TRUE), function(x) trimws(x[1]))]
Q <- Q[!is.na(Genes) & Genes != ""]
cat("### entries with a gene symbol:", nrow(Q), " unique genes:", uniqueN(Q$Genes), "\n")

sdt <- Q[, c("Genes", samples), with = FALSE]
agg <- sdt[, lapply(.SD, function(x){x <- suppressWarnings(as.numeric(x))
                                     if (all(is.na(x))) NA_real_ else mean(x, na.rm=TRUE)}),
           by = Genes, .SDcols = samples]
mat <- as.matrix(agg[,-1,with=FALSE]); rownames(mat) <- agg$Genes
cat("### gene-collapsed matrix:", nrow(mat), "x", ncol(mat), "\n")

v <- mat[is.finite(mat)]
cat(sprintf("### value range: min=%.1f median=%.1f max=%.1f  -> linear scale, log2 required\n",
            min(v), median(v), max(v)))

## ---- detection filter scan (manuscript 1.3: liver uses >=70%) ----
det  <- !is.na(mat)
mash <- names(grp)[grp == "MASH"]; nom <- names(grp)[grp == "No_MASLD"]; masl <- names(grp)[grp == "MASL"]
cat("\n=========== DETECTION FILTER SCAN (liver) ===========\n")
cat(sprintf("%6s %10s\n", "thr", "in-model"))
for (thr in c(.5,.6,.7,.8,.9,1)) {
  a <- rowSums(det[,mash,drop=FALSE]) >= ceiling(thr*length(mash))
  b <- rowSums(det[,nom ,drop=FALSE]) >= ceiling(thr*length(nom))
  cat(sprintf("%6.2f %10d\n", thr, sum(a & b)))
}

## ============================================================
run_liver <- function(g1, g2, thr, verbose = FALSE) {
  cols <- names(grp)[grp %in% c(g1, g2)]; g <- grp[cols]
  x <- mat[, cols, drop = FALSE]
  x[x <= 0] <- NA; x <- log2(x)                     # 1.3 log2
  d <- !is.na(x)
  keep <- rowSums(d[, g==g1, drop=FALSE]) >= ceiling(thr*sum(g==g1)) &
          rowSums(d[, g==g2, drop=FALSE]) >= ceiling(thr*sum(g==g2))
  x <- x[keep, , drop = FALSE]
  if (verbose) cat("   in model:", nrow(x), "\n")
  for (gg in unique(g)) {                            # within-group median imputation
    idx <- which(g == gg); sub <- x[, idx, drop = FALSE]
    gm <- apply(sub, 1, median, na.rm = TRUE)
    gm[!is.finite(gm)] <- median(sub, na.rm = TRUE)
    for (j in seq_along(idx)) { v <- sub[, j]; v[is.na(v)] <- gm[is.na(v)]; sub[, j] <- v }
    x[, idx] <- sub
  }
  cm <- apply(x, 2, median, na.rm = TRUE); x <- sweep(x, 2, cm, "-")   # column median norm
  design <- model.matrix(~0 + factor(g, levels = c(g2, g1))); colnames(design) <- c("ref","dis")
  fit <- eBayes(contrasts.fit(lmFit(x, design), makeContrasts(dis-ref, levels=design)), trend = TRUE)
  tt <- topTable(fit, number = Inf, sort.by = "none"); tt$gene <- rownames(tt)
  list(tt = tt, n_model = nrow(x))
}

report <- function(tag, tt, n_model, fc = log2(1.5)) {
  dep <- tt[tt$P.Value < 0.05 & abs(tt$logFC) >= fc, ]
  fdr <- tt[tt$adj.P.Val < 0.05 & abs(tt$logFC) >= fc, ]
  cat(sprintf("%-28s model %5d | DEP(nominal) %4d (up %3d / down %3d) | FDR+1.5x %4d\n",
              tag, n_model, nrow(dep), sum(dep$logFC>0), sum(dep$logFC<0), nrow(fdr)))
  invisible(list(dep=dep, fdr=fdr))
}

cat("\n=========== MAIN: MASH vs No_MASLD ===========\n")
main <- NULL
for (thr in c(.5,.6,.7,.8)) {
  r <- run_liver("MASH","No_MASLD", thr)
  rr <- report(sprintf("MASH vs No_MASLD thr=%.2f", thr), r$tt, r$n_model)
  if (thr == .7) main <- c(r, rr)
}
main_dep <- main$dep
fwrite(main$tt,  file.path(OUT,"liver_limma_full.csv"))
fwrite(main$dep, file.path(OUT,"liver_DEP.csv"))

cat("\n=========== SENSITIVITY: MASLD(MASH+MASL) vs No_MASLD, thr=0.7 ===========\n")
grp2 <- grp; grp2[grp2 == "MASH"] <- "MASLD"; grp2[grp2 == "MASL"] <- "MASLD"
sav <- grp; grp <<- grp2
sens <- run_liver("MASLD","No_MASLD", .7)
grp <<- sav
for (fc in c(log2(1.5), 1)) {
  d <- sens$tt[sens$tt$P.Value < 0.05 & abs(sens$tt$logFC) >= fc, ]
  cat(sprintf("   MASLD vs No_MASLD  fc=%.3f (%.1fx) | model %d | DEP %d\n",
              fc, 2^fc, sens$n_model, nrow(d)))
}

cat("\n=========== manuscript-named liver hits ===========\n")
named <- c("ACSL4","AKR1B10","LAMA5","PLIN1","AGRN","HSPB6","TGFB1I1",
           "NDUFA2","NDUFA7","GCAT","SMIM1","NDUFS6")
s <- main$tt[main$tt$gene %in% named, c("gene","logFC","P.Value","adj.P.Val")]
print(s[order(-s$logFC),], row.names = FALSE)

cat("\n=========== ALL liver DEPs (sorted by logFC) ===========\n")
print(main_dep[order(-main_dep$logFC), c("gene","logFC","P.Value","adj.P.Val")], row.names = FALSE)
