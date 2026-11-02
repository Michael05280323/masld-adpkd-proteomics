## ============================================================
## Independent reproduction  --  KIDNEY side (PXD056281)
## DT = ADPKD tissue (9 individuals x 3 tech reps); NT = normal (6 x 3)
## ============================================================
suppressPackageStartupMessages({library(data.table); library(limma)})

DATA <- "data"
OUT  <- "results"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
P <- file.path(DATA, "PXD056281__proteinGroups.txt")

hdr <- names(fread(P, nrows = 0))
lfq <- hdr[startsWith(hdr, "LFQ intensity ") & (grepl("_DT", hdr) | grepl("_NT", hdr))]
cat("### LFQ DT/NT columns:", length(lfq), "\n")
keep <- c("Gene names","Reverse","Potential contaminant","Only identified by site", lfq)
dt <- fread(P, select = keep, na.strings = c("NA","","NaN"))
cat("### raw protein groups:", nrow(dt), "\n")

for (cc in c("Reverse","Potential contaminant","Only identified by site")) {
  n0 <- nrow(dt); dt <- dt[is.na(get(cc)) | get(cc) != "+"]
  cat(sprintf("### removed '%s': -%d -> %d\n", cc, n0-nrow(dt), nrow(dt)))
}

g <- as.character(dt[["Gene names"]])
g <- sub(";.*$", "", g); g <- sub("[ \t].*$", "", g); g <- trimws(g)
dt[, Gene := g]; dt <- dt[!is.na(Gene) & Gene != ""]
cat("### rows with gene symbol:", nrow(dt), " unique genes:", uniqueN(dt$Gene), "\n")

m <- as.matrix(dt[, ..lfq])
m <- apply(m, 2, function(x){ x <- suppressWarnings(as.numeric(x)); x[x == 0] <- NA; x })
rownames(m) <- dt$Gene
cat(sprintf("### protein-group level LFQ non-missing: %.1f%%\n", 100*mean(!is.na(m))))

indiv <- sub("_[0-9]+$", "", sub("^LFQ intensity [0-9]+_", "", lfq))
cat("### individuals:", paste(sort(unique(indiv)), collapse=" "), "\n")

mi <- sapply(sort(unique(indiv)), function(id)
  apply(m[, indiv == id, drop = FALSE], 1, mean, na.rm = TRUE))
mi[is.nan(mi)] <- NA
cat("### individual-level matrix:", nrow(mi), "x", ncol(mi), "\n")
cat(sprintf("### individual-level LFQ non-missing: %.1f%%\n", 100*mean(!is.na(mi))))

grp <- ifelse(startsWith(colnames(mi), "DT"), "ADPKD", "normal"); names(grp) <- colnames(mi)
print(table(grp))
A <- colnames(mi)[grp=="ADPKD"]; N <- colnames(mi)[grp=="normal"]

## ---- collapse duplicate gene symbols: several plausible variants ----
collapse <- function(mi, how) {
  if (!anyDuplicated(rownames(mi))) return(mi)
  if (how == "mean") {
    ag <- rowsum(mi, rownames(mi), na.rm = TRUE, reorder = TRUE)
    cn <- rowsum(matrix(as.numeric(!is.na(mi)), nrow=nrow(mi), dimnames=dimnames(mi)),
                 rownames(mi), reorder = TRUE)
    o <- ag/cn; o[cn == 0] <- NA; o
  } else if (how == "max") {
    t(sapply(split(seq_len(nrow(mi)), rownames(mi)), function(ii)
      apply(mi[ii,,drop=FALSE], 2, function(v) if (all(is.na(v))) NA_real_ else max(v, na.rm=TRUE))))
  } else if (how == "first") {
    mi[!duplicated(rownames(mi)), , drop = FALSE]
  } else stop("unknown")
}

run <- function(x, tag, thr = 0.5) {
  d <- !is.na(x)
  keep <- rowSums(d[,A,drop=FALSE]) >= ceiling(thr*length(A)) &
          rowSums(d[,N,drop=FALSE]) >= ceiling(thr*length(N))
  x <- x[keep,,drop=FALSE]
  x[x <= 0] <- NA; x <- log2(x)
  for (gg in unique(grp)) {
    idx <- which(grp==gg); sub <- x[,idx,drop=FALSE]
    gm <- apply(sub,1,median,na.rm=TRUE); gm[!is.finite(gm)] <- median(sub,na.rm=TRUE)
    for (j in seq_along(idx)) { v <- sub[,j]; v[is.na(v)] <- gm[is.na(v)]; sub[,j] <- v }
    x[,idx] <- sub
  }
  cm <- apply(x,2,median,na.rm=TRUE); x <- sweep(x,2,cm,"-")
  design <- model.matrix(~0+factor(grp, levels=c("normal","ADPKD"))); colnames(design) <- c("ref","dis")
  fit <- eBayes(contrasts.fit(lmFit(x,design), makeContrasts(dis-ref, levels=design)), trend=TRUE)
  tt <- topTable(fit, number=Inf, sort.by="none"); tt$gene <- rownames(tt)
  dep <- tt[tt$P.Value<0.05 & abs(tt$logFC)>=log2(1.5),]
  fdr <- tt[tt$adj.P.Val<0.05 & abs(tt$logFC)>=log2(1.5),]
  cat(sprintf("%-10s | model %5d | DEP %4d (up %3d / down %3d) | FDR+1.5x %4d (up %3d / down %3d)\n",
              tag, nrow(x), nrow(dep), sum(dep$logFC>0), sum(dep$logFC<0),
              nrow(fdr), sum(fdr$logFC>0), sum(fdr$logFC<0)))
  list(tt=tt, dep=dep, n_model=nrow(x))
}

cat("\n=========== KIDNEY: ADPKD vs normal (thr=0.5) ===========\n")
res <- list()
for (how in c("mean","max","first")) {
  x <- collapse(mi, how)
  if (how == "mean") cat("### gene-level matrix:", nrow(x), "x", ncol(x), "\n")
  res[[how]] <- run(x, how)
}
cat("\n### MANUSCRIPT CLAIMS (mean replicate merge): model 3264 | DEP 1265 (651/614) | FDR+1.5x 1035 (567/468)\n")

r <- res[["mean"]]
fwrite(r$tt,  file.path(OUT,"kidney_limma_full.csv"))
fwrite(r$dep, file.path(OUT,"kidney_DEP.csv"))
saveRDS(list(mi=mi, grp=grp), file.path(OUT,"kidney_mat.rds"))

cat("\n=========== manuscript-named kidney hits ===========\n")
named <- c("LTF","PLIN4","EFEMP1","FLNC","AKR1B10","PDLIM4","PDLIM7","ZYX","LCN2",
           "FABP4","GPT","ANK2","SLC12A1","HMGCS2","NDUFS4","SDHA","HINT2","ECHDC2","DBT")
s <- r$tt[r$tt$gene %in% named, c("gene","logFC","P.Value","adj.P.Val")]
print(s[order(-s$logFC),], row.names=FALSE)
