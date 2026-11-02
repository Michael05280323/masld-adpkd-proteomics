## Scan pipeline variants to locate the manuscript's kidney DEP count (1265; 651up/614down)
suppressPackageStartupMessages({library(data.table); library(limma)})
OUT <- "results"
L <- readRDS(file.path(OUT,"kidney_mat.rds")); mi <- L$mi; grp <- L$grp
A <- colnames(mi)[grp=="ADPKD"]; N <- colnames(mi)[grp=="normal"]

run <- function(x, thr=0.5, log2_first=TRUE, impute=TRUE, norm=TRUE, trend=TRUE,
                collapse_mean=TRUE, use_floor=FALSE) {
  d <- !is.na(x)
  f <- if (use_floor) floor else ceiling
  keep <- rowSums(d[,A,drop=FALSE]) >= f(thr*length(A)) & rowSums(d[,N,drop=FALSE]) >= f(thr*length(N))
  x <- x[keep,,drop=FALSE]
  if (log2_first) { x[x<=0] <- NA; x <- log2(x) }
  if (impute) for (gg in unique(grp)) {
    idx <- which(grp==gg); sub <- x[,idx,drop=FALSE]
    gm <- apply(sub,1,median,na.rm=TRUE); gm[!is.finite(gm)] <- median(sub,na.rm=TRUE)
    for (j in seq_along(idx)) { v <- sub[,j]; v[is.na(v)] <- gm[is.na(v)]; sub[,j] <- v }
    x[,idx] <- sub
  }
  if (!log2_first) { x[x<=0] <- NA; x <- log2(x) }
  if (norm) { cm <- apply(x,2,median,na.rm=TRUE); x <- sweep(x,2,cm,"-") }
  design <- model.matrix(~0+factor(grp, levels=c("normal","ADPKD"))); colnames(design) <- c("ref","dis")
  fit <- eBayes(contrasts.fit(lmFit(x,design), makeContrasts(dis-ref, levels=design)), trend=trend)
  tt <- topTable(fit, number=Inf, sort.by="none")
  dep <- tt[tt$P.Value<0.05 & abs(tt$logFC)>=log2(1.5),]
  c(model=nrow(x), dep=nrow(dep), up=sum(dep$logFC>0), dn=sum(dep$logFC<0))
}

cat(sprintf("%-58s %6s %6s %5s %5s\n","variant","model","DEP","up","dn"))
show <- function(tag, ...) { r <- run(mi, ...)
  cat(sprintf("%-58s %6d %6d %5d %5d\n", tag, r[1], r[2], r[3], r[4])) }

show("baseline (log2,impute,norm,trend=TRUE)")
show("trend = FALSE",                                    trend=FALSE)
show("no column-median normalisation",                   norm=FALSE)
show("log2 AFTER imputation",                            log2_first=FALSE)
show("no imputation",                                    impute=FALSE)
show("floor instead of ceiling on 50% threshold",        use_floor=TRUE)
show("thr=0.4",                                          thr=0.4)
show("thr=0.6",                                          thr=0.6)
cat("\n### MANUSCRIPT: model 3264 | DEP 1265 | up 651 | dn 614\n")

## ---------------------------------------------------------------------------
## Technical-replicate merge: MEAN (manuscript) vs MEDIAN
## Reproduces the kidney FDR discrepancy: mean -> 1035 (567/468),
## median -> 1049 (575/474). Reads the raw proteinGroups and re-merges.
## ---------------------------------------------------------------------------
suppressPackageStartupMessages(library(data.table))
DATA <- "data"
P <- file.path(DATA, "PXD056281__proteinGroups.txt")
if (file.exists(P)) {
  hdr <- names(fread(P, nrows = 0))
  lfq <- hdr[startsWith(hdr, "LFQ intensity ") & (grepl("_DT", hdr) | grepl("_NT", hdr))]
  dt <- fread(P, select = c("Gene names","Reverse","Potential contaminant",
                            "Only identified by site", lfq),
              na.strings = c("NA","","NaN"))
  for (cc in c("Reverse","Potential contaminant","Only identified by site"))
    dt <- dt[is.na(get(cc)) | get(cc) != "+"]
  g <- sub("[ \t].*$", "", sub(";.*$", "", trimws(as.character(dt[["Gene names"]]))))
  dt[, Gene := g]; dt <- dt[!is.na(Gene) & Gene != ""]
  m <- as.matrix(dt[, ..lfq])
  m <- apply(m, 2, function(x){ x <- suppressWarnings(as.numeric(x)); x[x==0] <- NA; x })
  rownames(m) <- dt$Gene
  indiv <- sub("_[0-9]+$", "", sub("^LFQ intensity [0-9]+_", "", lfq))
  grp2 <- ifelse(startsWith(sort(unique(indiv)), "DT"), "ADPKD", "normal")

  build_indiv <- function(how) {
    im <- sapply(sort(unique(indiv)), function(id)
      apply(m[, indiv == id, drop = FALSE], 1, how, na.rm = TRUE))
    im[is.nan(im)] <- NA
    ## collapse duplicate genes by mean
    ag <- rowsum(im, rownames(im), na.rm = TRUE)
    cn <- rowsum(matrix(as.numeric(!is.na(im)), nrow=nrow(im)), rownames(im))
    o <- ag/cn; o[cn == 0] <- NA; o
  }
  de_pipeline <- function(im) {
    grp3 <- ifelse(startsWith(colnames(im), "DT"), "ADPKD", "normal")
    A <- colnames(im)[grp3=="ADPKD"]; N <- colnames(im)[grp3=="normal"]
    d <- !is.na(im)
    keep <- rowSums(d[,A,drop=FALSE]) >= ceiling(.5*length(A)) &
            rowSums(d[,N,drop=FALSE]) >= ceiling(.5*length(N))
    im <- im[keep,,drop=FALSE]; im[im<=0] <- NA; im <- log2(im)
    for (gg in unique(grp3)) { idx <- which(grp3==gg); s <- im[,idx,drop=FALSE]
      gm <- apply(s,1,median,na.rm=TRUE); gm[!is.finite(gm)] <- median(s,na.rm=TRUE)
      for (j in seq_along(idx)) { v <- s[,j]; v[is.na(v)] <- gm[is.na(v)]; s[,j] <- v }
      im[,idx] <- s }
    im <- sweep(im, 2, apply(im,2,median,na.rm=TRUE), "-")
    design <- model.matrix(~0+factor(grp3, levels=c("normal","ADPKD"))); colnames(design) <- c("ref","dis")
    fit <- eBayes(contrasts.fit(lmFit(im,design), makeContrasts(dis-ref, levels=design)), trend=TRUE)
    tt <- topTable(fit, number=Inf, sort.by="none")
    dep <- tt[tt$P.Value<0.05 & abs(tt$logFC)>=log2(1.5),]
    fdr <- tt[tt$adj.P.Val<0.05 & abs(tt$logFC)>=log2(1.5),]
    c(model=nrow(im), dep=nrow(dep), dup=sum(dep$logFC>0), ddn=sum(dep$logFC<0),
      fdr=nrow(fdr), fup=sum(fdr$logFC>0), fdn=sum(fdr$logFC<0))
  }
  cat("\n### technical-replicate merge comparison (3 reps per individual)\n")
  for (how in c("mean","median")) {
    r <- de_pipeline(build_indiv(how))
    cat(sprintf("  %-6s | model %d | DEP %d (%d/%d) | FDR+1.5x %d (%d/%d)\n",
                how, r["model"], r["dep"], r["dup"], r["ddn"],
                r["fdr"], r["fup"], r["fdn"]))
  }
  cat("  MANUSCRIPT uses MEAN -> FDR 1035 (567/468); median would give 1049 (575/474).\n")
} else {
  cat("\n(raw proteinGroups not present; skipped mean-vs-median merge comparison)\n")
}
