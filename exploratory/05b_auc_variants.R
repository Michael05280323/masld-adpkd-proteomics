## Locate the kidney-AUC definition used in Table 3 of the manuscript.
suppressPackageStartupMessages({library(data.table); library(pROC)})
DATA <- "data"
OUT  <- "results"
L <- readRDS(file.path(OUT,"kidney_mat.rds")); mi <- L$mi; grp <- L$grp
A <- colnames(mi)[grp=="ADPKD"]; N <- colnames(mi)[grp=="normal"]
lab <- factor(grp, levels=c("normal","ADPKD"))

## gene-level non-missing counts for the molecules of interest
mol  <- c("HINT2","GPT","ECHDC2","DBT","DPYSL3","S100A10","ANXA2")
TGT  <- c(HINT2=.96, GPT=1.00, ECHDC2=.89, DBT=.78, DPYSL3=.87, S100A10=.81, ANXA2=.65)
cat("### detection of the 7 Table-3 molecules across the 15 individuals:\n")
print(data.table(gene=mol, detected=rowSums(!is.na(mi[mol,,drop=FALSE])), of=15))

dirAUC <- function(v, lfc) {                 # direction-consistent AUC
  ok <- !is.na(v); if (sum(ok) < 3) return(NA_real_)
  r <- as.numeric(pROC::auc(pROC::roc(lab[ok], v[ok], quiet=TRUE, direction="<")))
  if (lfc >= 0) r else 1 - r
}
## per-gene logFC sign taken from the joint table
J <- fread(file.path(OUT,"joint_all_common.csv"))
sgn <- setNames(sign(J$kfc), J$gene)

prep <- function(use_log2=TRUE, impute=TRUE, norm=TRUE, filter=TRUE) {
  x <- mi
  if (filter) x <- x[rowSums(!is.na(x[,A,drop=FALSE]))>=ceiling(.5*length(A)) &
                     rowSums(!is.na(x[,N,drop=FALSE]))>=ceiling(.5*length(N)),,drop=FALSE]
  if (use_log2) { x[x<=0] <- NA; x <- log2(x) }
  if (impute) for (gg in unique(grp)) { i <- which(grp==gg); s <- x[,i,drop=FALSE]
    gm <- apply(s,1,median,na.rm=TRUE); gm[!is.finite(gm)] <- median(s,na.rm=TRUE)
    for (j in seq_along(i)) { v<-s[,j]; v[is.na(v)]<-gm[is.na(v)]; s[,j]<-v }; x[,i] <- s }
  if (norm) x <- sweep(x, 2, apply(x,2,median,na.rm=TRUE), "-")
  x
}
evalv <- function(tag, x) {
  a <- sapply(mol, function(g) if (g %in% rownames(x)) dirAUC(x[g,], sgn[[g]]) else NA_real_)
  d <- round(a - TGT[mol], 3)
  n1 <- if ("ANXA2" %in% rownames(x)) sum(sapply(rownames(x), function(g) {
     v <- x[g,]; if (any(is.na(v))) return(FALSE)
     r <- as.numeric(pROC::auc(pROC::roc(lab, v, quiet=TRUE, direction="<"))); r>=0.9999 })) else NA
  cat(sprintf("%-42s | %s | maxdiff %.2f | AUC==1: %d/%d (%.1f%%)\n", tag,
      paste(sprintf("%.2f", a), collapse=" "), max(abs(d), na.rm=TRUE), n1, nrow(x), 100*n1/nrow(x)))
  invisible(a)
}
cat("\n### target (manuscript Table 3):           ",
    paste(sprintf("%.2f", TGT[mol]), collapse=" "), "\n\n")
cat("### order: HINT2 GPT ECHDC2 DBT DPYSL3 S100A10 ANXA2\n\n")

evalv("log2 + impute + normalise (baseline)",  prep(TRUE , TRUE , TRUE ))
evalv("log2 + impute, NO normalise",          prep(TRUE , TRUE , FALSE))
evalv("log2, NO impute, NO normalise",        prep(TRUE , FALSE, FALSE))
evalv("log2 + impute + normalise, no filter", prep(TRUE , TRUE , TRUE , FALSE))
evalv("linear + impute + normalise",          prep(FALSE, TRUE , TRUE ))
evalv("linear + impute, NO normalise",        prep(FALSE, TRUE , FALSE))

## protein-group level (45 columns, technical replicates NOT merged)
P <- file.path(DATA,"PXD056281__proteinGroups.txt")
hdr <- names(fread(P, nrows=0))
lfq <- hdr[startsWith(hdr,"LFQ intensity ") & (grepl("_DT",hdr)|grepl("_NT",hdr))]
dt <- fread(P, select=c("Gene names","Reverse","Potential contaminant","Only identified by site",lfq),
            na.strings=c("NA","","NaN"))
for (cc in c("Reverse","Potential contaminant","Only identified by site"))
  dt <- dt[is.na(get(cc)) | get(cc)!="+"]
g <- sub("[ \t].*$","",sub(";.*$","",as.character(dt[["Gene names"]]))); dt[, Gene:=trimws(g)]
dt <- dt[!is.na(Gene) & Gene!=""]
m <- as.matrix(dt[,..lfq]); m <- apply(m,2,function(x){x<-suppressWarnings(as.numeric(x));x[x==0]<-NA;x})
rownames(m) <- dt$Gene
ag <- rowsum(m, rownames(m), na.rm=TRUE, reorder=TRUE)
cn <- rowsum(matrix(as.numeric(!is.na(m)),nrow=nrow(m),dimnames=dimnames(m)), rownames(m), reorder=TRUE)
mg <- ag/cn; mg[cn==0] <- NA
gl <- factor(ifelse(startsWith(colnames(mg),"DT"),"ADPKD","normal"), levels=c("normal","ADPKD"))
mg[mg<=0] <- NA; mg <- log2(mg)
for (gg in unique(gl)) { i <- which(gl==gg); s <- mg[,i,drop=FALSE]
  gm <- apply(s,1,median,na.rm=TRUE); gm[!is.finite(gm)] <- median(s,na.rm=TRUE)
  for (j in seq_along(i)) { v<-s[,j]; v[is.na(v)]<-gm[is.na(v)]; s[,j]<-v }; mg[,i] <- s }
mgN <- sweep(mg,2,apply(mg,2,median,na.rm=TRUE),"-")
lab <<- gl
t <- sapply(mol, function(g) dirAUC(mg[g,], sgn[[g]]))
cat(sprintf("%-42s | %s | maxdiff %.2f\n","45 protein-group columns (pseudo-reps), NO norm",
            paste(sprintf("%.2f",t),collapse=" "), max(abs(round(t-TGT[mol],3)))))
t <- sapply(mol, function(g) dirAUC(mgN[g,], sgn[[g]]))
cat(sprintf("%-42s | %s | maxdiff %.2f\n","45 protein-group columns (pseudo-reps), norm",
            paste(sprintf("%.2f",t),collapse=" "), max(abs(round(t-TGT[mol],3)))))
