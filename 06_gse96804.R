## Independent validation in GSE96804 (human glomerulus, diabetic nephropathy vs control)
suppressPackageStartupMessages({
  library(data.table); library(limma); library(GEOquery)
  library(hta20transcriptcluster.db); library(AnnotationDbi)
})
DATA <- "data"
OUT  <- "results"

## getGPL=FALSE: probe->symbol mapping uses the local hta20transcriptcluster.db
## package below, so the large NCBI platform annotation (GPL17586.soft.gz) is
## not downloaded; this also avoids a slow/timeout-prone network call.
es <- getGEO(filename = file.path(DATA,"GSE96804_series_matrix.txt.gz"),
             getGPL = FALSE)
cat("### platform:", annotation(es), " | samples:", ncol(es), "\n")
pd <- pData(es)
cat("### pData columns:\n"); print(colnames(pd))

## grouping lives in source_name_ch1:
##   "diabetic nephropathy" vs "unaffected portion of tumor nephrectomies"
key <- "source_name_ch1"
cat("\n### grouping column:", key, "\n"); print(table(pd[[key]]))
g <- as.character(pd[[key]])
g <- ifelse(grepl("diabetic nephropathy", g, ignore.case=TRUE), "DN",
     ifelse(grepl("unaffected|tumor nephrectom", g, ignore.case=TRUE), "CTRL", NA))
cat("### group sizes:\n"); print(table(g, useNA="ifany"))
keep <- !is.na(g); g <- factor(g[keep], levels=c("CTRL","DN"))
cat("### analysed:", ncol(es[keep]), "samples ->", paste(names(table(g)), table(g), collapse=" / "), "\n")

ex <- exprs(es)[, keep]
cat("### probes:", nrow(ex), " range:", paste(round(range(ex, na.rm=TRUE),2), collapse=" .. "), "\n")
if (max(ex, na.rm=TRUE) > 100) { ex <- log2(ex); cat("### values were linear -> log2 applied\n") } else
  cat("### values already log2 (RMA/SST) - left as is\n")

sym <- mapIds(hta20transcriptcluster.db, keys=rownames(ex), column="SYMBOL",
              keytype="PROBEID", multiVals="first")
cat("### probes mapped to a symbol:", sum(!is.na(sym)), "/", nrow(ex), "\n")

## collapse multiple probes per gene: keep the probe with highest mean expression
dt <- data.table(probe=rownames(ex), sym=as.character(sym), mean=rowMeans(ex, na.rm=TRUE))
dt <- dt[!is.na(sym)]
setorder(dt, sym, -mean)
sel <- dt[, .SD[1], by=sym]
M <- ex[sel$probe, , drop=FALSE]; rownames(M) <- sel$sym
cat("### gene-level matrix:", dim(M), "\n")

design <- model.matrix(~0+g); colnames(design) <- levels(g)
fit <- eBayes(contrasts.fit(lmFit(M, design), makeContrasts(DN-CTRL, levels=design)))
tt <- topTable(fit, number=Inf, sort.by="none"); tt$gene <- rownames(tt)

## --- manuscript 2.7 claim: NDUFS3, NDUFS8, NDUFS2, SDHA, HSPA9 all significantly DOWN ---
cat("\n=========== manuscript 2.7 claim (mitochondrial module representatives) ===========\n")
g1 <- c("NDUFS3","NDUFS8","NDUFS2","SDHA","HSPA9")
print(tt[g1, c("gene","logFC","P.Value","adj.P.Val")], row.names=FALSE)
cat("\n### all concordant-down complex-I / OXPHOS genes present in GSE96804:\n")
J <- fread(file.path(OUT,"joint_all_common.csv"))
dn <- J[kfc<0 & lfc<0, gene]
nd <- intersect(grep("^NDUF|^SDH|^COX|^ATP5|^UQCR|^CYC", dn, value=TRUE), rownames(tt))
sub <- tt[nd, c("gene","logFC","P.Value","adj.P.Val")]
print(sub[order(sub$logFC),], row.names=FALSE)
cat(sprintf("\n### of %d shared OXPHOS genes tested: %d down (P<0.05), %d up (P<0.05)\n",
            length(nd), sum(tt[nd,]$logFC<0 & tt[nd,]$P.Value<0.05),
            sum(tt[nd,]$logFC>0 & tt[nd,]$P.Value<0.05)))
fwrite(tt, file.path(OUT,"gse96804_limma.csv"))
cat("\n### done\n")
