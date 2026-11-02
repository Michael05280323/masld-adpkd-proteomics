## ============================================================
## Joint / cross-organ shared-protein analysis  (manuscript 1.4)
## ============================================================
suppressPackageStartupMessages({library(data.table); library(limma)})
DATA <- "data"
OUT  <- "results"

## ---- rebuild liver model matrix ----
Q <- fread(file.path(DATA,"PXD051911__liver_protein_quant.txt"), na.strings=c("NA","","NaN"))
M <- fread(file.path(DATA,"PXD051911__meta_data.txt"), na.strings=c("NA","","NaN"))
sm <- setdiff(names(Q), c("ProteinAccessions","Genes","ProteinDescriptions"))
gm <- M$saf_diagnosis; names(gm) <- M$liver_proteomics_filename; lgrp <- gm[sm]
Q[, Genes := trimws(as.character(Genes))]; Q <- Q[!is.na(Genes) & Genes!=""]
sd <- Q[, c("Genes",sm), with=FALSE]
ag <- sd[, lapply(.SD, function(x){x<-suppressWarnings(as.numeric(x))
                                    if(all(is.na(x))) NA_real_ else mean(x,na.rm=TRUE)}),
         by=Genes, .SDcols=sm]
lm_ <- as.matrix(ag[,-1,with=FALSE]); rownames(lm_) <- ag$Genes
lm_[lm_<=0] <- NA; lm_ <- log2(lm_)
keep <- rowSums(!is.na(lm_[,lgrp=="MASH",drop=FALSE])) >= ceiling(.7*sum(lgrp=="MASH")) &
        rowSums(!is.na(lm_[,lgrp=="No_MASLD",drop=FALSE])) >= ceiling(.7*sum(lgrp=="No_MASLD"))
lm_ <- lm_[keep,,drop=FALSE]
for (gg in unique(lgrp)) { idx <- which(lgrp==gg); s <- lm_[,idx,drop=FALSE]
  gmed <- apply(s,1,median,na.rm=TRUE); gmed[!is.finite(gmed)] <- median(s,na.rm=TRUE)
  for (j in seq_along(idx)) { v <- s[,j]; v[is.na(v)] <- gmed[is.na(v)]; s[,j] <- v }
  lm_[,idx] <- s }
lm_ <- sweep(lm_, 2, apply(lm_,2,median,na.rm=TRUE), "-")
cat("### liver model matrix:", nrow(lm_), "genes x", ncol(lm_), "samples\n")

## ---- rebuild kidney model matrix ----
L <- readRDS(file.path(OUT,"kidney_mat.rds")); mi <- L$mi; kgrp <- L$grp
A <- colnames(mi)[kgrp=="ADPKD"]; N <- colnames(mi)[kgrp=="normal"]
kd <- mi[rowSums(!is.na(mi[,A,drop=FALSE]))>=ceiling(.5*length(A)) &
         rowSums(!is.na(mi[,N,drop=FALSE]))>=ceiling(.5*length(N)),,drop=FALSE]
kd[kd<=0] <- NA; kd <- log2(kd)
for (gg in unique(kgrp)) { idx <- which(kgrp==gg); s <- kd[,idx,drop=FALSE]
  gmed <- apply(s,1,median,na.rm=TRUE); gmed[!is.finite(gmed)] <- median(s,na.rm=TRUE)
  for (j in seq_along(idx)) { v <- s[,j]; v[is.na(v)] <- gmed[is.na(v)]; s[,j] <- v }
  kd[,idx] <- s }
kd <- sweep(kd, 2, apply(kd,2,median,na.rm=TRUE), "-")
cat("### kidney model matrix:", nrow(kd), "genes x", ncol(kd), "samples\n")

## ---- differential results ----
lt <- fread(file.path(OUT,"liver_limma_full.csv"));   setnames(lt, "V1","gene", skip_absent=TRUE)
kt <- fread(file.path(OUT,"kidney_limma_full.csv"));  setnames(kt, "V1","gene", skip_absent=TRUE)
cat("### liver tt rows:", nrow(lt), " kidney tt rows:", nrow(kt), "\n")

common <- intersect(lt$gene, kt$gene)
cat("### COMMON quantified background:", length(common), " (manuscript: 2851)\n")

J <- merge(lt[gene %in% common, .(gene, lfc=logFC, lp=P.Value, ladj=adj.P.Val)],
           kt[gene %in% common, .(gene, kfc=logFC, kp=P.Value, kadj=adj.P.Val)], by="gene")
cat("### merged:", nrow(J), "\n")

conc <- J[lp<0.05 & kp<0.05 & sign(lfc)==sign(kfc)]
cat(sprintf("### CONCORDANT (both nominal P<0.05, same direction): %d  (up %d / down %d)   [manuscript 362 = 220+142]\n",
            nrow(conc), sum(conc$lfc>0), sum(conc$lfc<0)))

core <- J[lp<0.05 & kp<0.05 & abs(lfc)>=log2(1.5) & abs(kfc)>=log2(1.5) & sign(lfc)==sign(kfc)]
cat(sprintf("### STRICT CORE (both also >=1.5x): %d   [manuscript 81]\n", nrow(core)))

liverDEP <- lt[P.Value<0.05 & abs(logFC)>=log2(1.5), gene]
kidneyDEP<- kt[P.Value<0.05 & abs(logFC)>=log2(1.5), gene]
cat(sprintf("### VENN intersection of DEP sets: %d   [manuscript 96]\n", length(intersect(liverDEP,kidneyDEP))))

fwrite(J,    file.path(OUT,"joint_all_common.csv"))
fwrite(conc, file.path(OUT,"joint_concordant.csv"))
fwrite(core, file.path(OUT,"joint_strict_core.csv"))

cat("\n=========== manuscript Table 3 / module representatives ===========\n")
named <- c("HINT2","GPT","ECHDC2","DBT","DPYSL3","S100A10","ANXA2","MARCKS",
           "SDHA","SUCLA2","FLNA","LUM","PXDN","FBLN5","SERPINF1","COBLL1","ZYX","NDUFS3","NDUFS8","NDUFS2","HSPA9")
print(J[gene %in% named][order(-lfc)], row.names=FALSE)

cat("\n=========== complex-I subunits among concordant DOWN ===========\n")
dn <- conc[lfc<0]
nduf <- grep("^NDUF", dn$gene, value=TRUE)
cat("NDUF subunits concordantly down:", length(nduf), "\n"); print(sort(nduf))
