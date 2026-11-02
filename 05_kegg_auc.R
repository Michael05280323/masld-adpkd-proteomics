## Corrected KEGG + direction-aware cross-cohort AUC
suppressPackageStartupMessages({library(data.table);library(clusterProfiler)
  library(org.Hs.eg.db);library(pROC)})
DATA <- "data"
OUT  <- "results"
S <- readRDS(file.path(OUT,"concordant_sets.rds")); up <- S$up; dn <- S$dn
J <- fread(file.path(OUT,"joint_all_common.csv")); sgn <- setNames(sign(J$kfc), J$gene)

## ================= KEGG (fixed column order) =================
cat("#################  KEGG  #################\n")
lk <- fread(file.path(OUT,"kegg_link_hsa.tsv"), header=FALSE, sep="\t")
nm <- fread(file.path(OUT,"kegg_pathway_hsa.tsv"), header=FALSE, sep="\t")
setnames(lk, c("path","gene")); setnames(nm, c("path","desc"))
lk[, `:=`(gene=sub("^hsa:","",gene), path=sub("^path:","",path))]
nm[, desc := sub(" - Homo sapiens \\(human\\)$","",desc)]
t2g <- unique(merge(lk, nm, by="path")[, .(term=desc, gene)])
cat("### TERM2GENE:", nrow(t2g), "pairs /", uniqueN(t2g$term), "pathways\n")
ent <- function(s) unique(na.omit(bitr(s,"SYMBOL","ENTREZID",org.Hs.eg.db)$ENTREZID))
kegg <- function(s, lab) {
  r <- enricher(ent(s), TERM2GENE=t2g, pvalueCutoff=1, qvalueCutoff=1, minGSSize=5)
  if (is.null(r)) { cat(lab,": none\n"); return(invisible()) }
  d <- as.data.frame(r)
  cat(sprintf("\n--- KEGG %s : %d pathways, top 8 ---\n", lab, nrow(d)))
  print(head(d[,c("ID","Description","GeneRatio","pvalue","p.adjust","Count")],8), row.names=FALSE)
  h <- d[grepl("Oxidative phosphorylation|fatty liver|lipid and atherosclerosis|Carbon metabolism",
               d$Description, ignore.case=TRUE), c("ID","Description","GeneRatio","pvalue","p.adjust","Count")]
  if (nrow(h)) { cat("   >> manuscript-cited:\n"); print(h, row.names=FALSE) }
}
kegg(dn,"DOWN"); kegg(up,"UP")

## ================= AUC =================
cat("\n\n#################  CROSS-COHORT AUC  #################\n")
## liver matrix
Q <- fread(file.path(DATA,"PXD051911__liver_protein_quant.txt"), na.strings=c("NA","","NaN"))
M <- fread(file.path(DATA,"PXD051911__meta_data.txt"), na.strings=c("NA","","NaN"))
sm <- setdiff(names(Q), c("ProteinAccessions","Genes","ProteinDescriptions"))
gm <- M$saf_diagnosis; names(gm) <- M$liver_proteomics_filename; grp <- gm[sm]; names(grp) <- sm
Q[, Genes := trimws(as.character(Genes))]; Q <- Q[!is.na(Genes) & Genes!=""]
ag <- Q[, c("Genes",sm), with=FALSE][, lapply(.SD, function(x){
      x<-suppressWarnings(as.numeric(x)); if(all(is.na(x))) NA_real_ else mean(x,na.rm=TRUE)}),
      by=Genes, .SDcols=sm]
LM <- as.matrix(ag[,-1,with=FALSE]); rownames(LM) <- ag$Genes
LM[LM<=0] <- NA; LM <- log2(LM)
k <- rowSums(!is.na(LM[,grp=="MASH",drop=FALSE]))>=ceiling(.7*sum(grp=="MASH")) &
     rowSums(!is.na(LM[,grp=="No_MASLD",drop=FALSE]))>=ceiling(.7*sum(grp=="No_MASLD"))
LM <- LM[k,,drop=FALSE]
for (gg in unique(grp)) { i<-which(grp==gg); s<-LM[,i,drop=FALSE]
  g0<-apply(s,1,median,na.rm=TRUE); g0[!is.finite(g0)]<-median(s,na.rm=TRUE)
  for (j in seq_along(i)) { v<-s[,j]; v[is.na(v)]<-g0[is.na(v)]; s[,j]<-v }; LM[,i]<-s }
LMn <- sweep(LM,2,apply(LM,2,median,na.rm=TRUE),"-")     # normalised
lLab <- factor(grp, levels=c("No_MASLD","MASH"))
## kidney matrix
L <- readRDS(file.path(OUT,"kidney_mat.rds")); mi <- L$mi; kgrp <- L$grp
A <- colnames(mi)[kgrp=="ADPKD"]; N <- colnames(mi)[kgrp=="normal"]
KM <- mi[rowSums(!is.na(mi[,A,drop=FALSE]))>=ceiling(.5*length(A)) &
         rowSums(!is.na(mi[,N,drop=FALSE]))>=ceiling(.5*length(N)),,drop=FALSE]
KM[KM<=0] <- NA; KM <- log2(KM)
for (gg in unique(kgrp)) { i<-which(kgrp==gg); s<-KM[,i,drop=FALSE]
  g0<-apply(s,1,median,na.rm=TRUE); g0[!is.finite(g0)]<-median(s,na.rm=TRUE)
  for (j in seq_along(i)) { v<-s[,j]; v[is.na(v)]<-g0[is.na(v)]; s[,j]<-v }; KM[,i]<-s }
KMn <- sweep(KM,2,apply(KM,2,median,na.rm=TRUE),"-")
kLab <- factor(kgrp, levels=c("normal","ADPKD"))
cat("### liver:", dim(LM), " kidney:", dim(KM), "\n")

aucd <- function(v, lab, lfc) { r <- as.numeric(pROC::auc(pROC::roc(lab, v, quiet=TRUE, direction="<")))
                                if (lfc >= 0) r else 1-r }
best <- function(v, lab) max(as.numeric(pROC::auc(pROC::roc(lab, v, quiet=TRUE, direction="<"))),
                             1-as.numeric(pROC::auc(pROC::roc(lab, v, quiet=TRUE, direction="<"))))
allg <- intersect(rownames(LM), rownames(KM))
AU <- data.table(gene=allg,
  liverAUC  = sapply(allg, function(g) aucd(LMn[g,], lLab, sgn[[g]])),
  kidneyAUC = sapply(allg, function(g) aucd(KM[g,],  kLab, sgn[[g]])))  # un-normalised, as in Table 3
AU[, minAUC := pmin(liverAUC, kidneyAUC)]
fwrite(AU, file.path(OUT,"auc_all_shared.csv"))

cat("\n--- how common is a perfect kidney AUC at n=9 vs 6? ---\n")
for (nmx in c("KMn","KM")) {
  X <- get(nmx)
  f1 <- mean(sapply(rownames(X), function(g) best(X[g,], kLab) >= 0.9999))
  f2 <- mean(sapply(rownames(X), function(g) best(X[g,], kLab) >= 0.95))
  cat(sprintf("   %-4s (normalised=%s): AUC==1.00 %.2f%% (%d/%d) | >=0.95 %.2f%%\n",
              nmx, nmx=="KMn", 100*f1, round(f1*nrow(X)), nrow(X), 100*f2))
}
cat("   manuscript claims: ~2.9% of modelled proteins reach AUC 1.00\n")

mol <- c("HINT2","GPT","ECHDC2","DBT","DPYSL3","S100A10","ANXA2")
ms  <- data.table(gene=mol, ms_liver=c(.95,.93,.92,.95,.96,.94,.99),
                  ms_kidney=c(.96,1.00,.89,.78,.87,.81,.65))
cmp <- merge(ms, AU, by="gene", all.x=TRUE)
cmp[, `:=`(d_liver=round(liverAUC-ms_liver,2), d_kidney=round(kidneyAUC-ms_kidney,2))]
cat("\n--- Table 3: manuscript vs reproduction ---\n"); print(cmp, row.names=FALSE)
cat(sprintf("\n### liver  AUC  max |diff| = %.2f\n", max(abs(cmp$d_liver))))
cat(sprintf("### kidney AUC  max |diff| = %.2f\n", max(abs(cmp$d_kidney))))

cat("\n--- top 15 by min(AUC) ---\n")
print(head(AU[order(-minAUC)],15), row.names=FALSE)
