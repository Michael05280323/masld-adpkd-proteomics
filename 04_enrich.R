## ============================================================
## GO / KEGG over-representation on the concordant shared set (manuscript 1.5)
## ============================================================
suppressPackageStartupMessages({
  library(data.table); library(clusterProfiler); library(org.Hs.eg.db); library(ggplot2)
})
OUT <- "results"
conc <- fread(file.path(OUT,"joint_concordant.csv"))
up   <- conc[lfc>0]$gene
dn   <- conc[lfc<0]$gene
cat("### concordant up:", length(up), " down:", length(dn), "\n")

map <- function(g) {
  m <- suppressWarnings(bitr(g, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db))
  unique(m$ENTREZID)
}
eUp <- map(up); eDn <- map(dn)
cat("### mapped to Entrez  up:", length(eUp), " down:", length(eDn), "\n")

show <- function(res, lab, n=12) {
  if (is.null(res) || nrow(as.data.frame(res))==0) { cat("  [none]\n"); return(invisible()) }
  d <- as.data.frame(res)
  cat(sprintf("\n--- %s (top %d of %d) ---\n", lab, min(n,nrow(d)), nrow(d)))
  print(head(d[,c("Description","GeneRatio","pvalue","p.adjust","Count")], n), row.names=FALSE)
}
grab <- function(res, terms) {
  if (is.null(res)) return(invisible())
  d <- as.data.frame(res)
  hit <- d[d$Description %in% terms, c("Description","GeneRatio","pvalue","p.adjust","Count")]
  if (nrow(hit)) { cat("\n  >> manuscript-cited terms:\n"); print(hit, row.names=FALSE) }
}

## ---- overall GO (all three ontologies) ----
for (ont in c("BP","CC","MF")) {
  r <- tryCatch(enrichGO(unique(c(up,dn)), OrgDb=org.Hs.eg.db, ont=ont,
                         pAdjustMethod="BH", pvalueCutoff=1, qvalueCutoff=1,
                         readable=TRUE), error=function(e) NULL)
  cat(sprintf("\n========== GO-%s on ALL concordant (%d genes) ==========\n", ont, length(unique(c(up,dn)))))
  show(r, paste0("GO-",ont))
}

## ---- down-regulated only ----
cat("\n\n###############  DOWN-REGULATED SET  ###############\n")
rDn <- tryCatch(enrichGO(eDn, OrgDb=org.Hs.eg.db, ont="BP", pAdjustMethod="BH",
                         pvalueCutoff=1, qvalueCutoff=1, readable=TRUE), error=function(e) NULL)
show(rDn, "GO-BP  DOWN")
grab(rDn, c("cellular respiration","aerobic respiration",
            "mitochondrial electron transport, NADH to ubiquinone",
            "oxidative phosphorylation","electron transport chain"))
rDnC <- tryCatch(enrichGO(eDn, OrgDb=org.Hs.eg.db, ont="CC", pAdjustMethod="BH",
                          pvalueCutoff=1, qvalueCutoff=1, readable=TRUE), error=function(e) NULL)
show(rDnC, "GO-CC  DOWN")
rDnM <- tryCatch(enrichGO(eDn, OrgDb=org.Hs.eg.db, ont="MF", pAdjustMethod="BH",
                          pvalueCutoff=1, qvalueCutoff=1, readable=TRUE), error=function(e) NULL)
show(rDnM, "GO-MF  DOWN")

## ---- up-regulated only ----
cat("\n\n###############  UP-REGULATED SET  ###############\n")
rUp <- tryCatch(enrichGO(eUp, OrgDb=org.Hs.eg.db, ont="BP", pAdjustMethod="BH",
                         pvalueCutoff=1, qvalueCutoff=1, readable=TRUE), error=function(e) NULL)
show(rUp, "GO-BP  UP")
grab(rUp, c("actin cytoskeleton organization","actin filament organization",
            "actin nucleation","cell-substrate adhesion"))
rUpC <- tryCatch(enrichGO(eUp, OrgDb=org.Hs.eg.db, ont="CC", pAdjustMethod="BH",
                          pvalueCutoff=1, qvalueCutoff=1, readable=TRUE), error=function(e) NULL)
show(rUpC, "GO-CC  UP")
grab(rUpC, c("cell-substrate junction","focal adhesion","actin cytoskeleton"))

## ---- KEGG ----
cat("\n\n###############  KEGG  ###############\n")
kk <- function(eg, lab) {
  r <- tryCatch(enrichKEGG(eg, organism="hsa", pvalueCutoff=1, qvalueCutoff=1),
                error=function(e){cat("  KEGG failed:", conditionMessage(e), "\n"); NULL})
  if (is.null(r)) return(invisible())
  d <- as.data.frame(setReadable(r, org.Hs.eg.db, keyType="ENTREZID"))
  cat(sprintf("\n--- KEGG %s (top 12 of %d) ---\n", lab, nrow(d)))
  print(head(d[,c("Description","GeneRatio","pvalue","p.adjust","Count")], 12), row.names=FALSE)
  invisible(d)
}
kkD <- kk(eDn, "DOWN"); kkU <- kk(eUp, "UP")
kall <- kk(unique(c(eUp,eDn)), "ALL concordant")

saveRDS(list(up=up, dn=dn), file.path(OUT,"concordant_sets.rds"))
cat("\n### done\n")
