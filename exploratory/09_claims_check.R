suppressPackageStartupMessages({library(data.table)})
OUT <- "results"

cat("=========== 1. Are the Table-2 'module high-discrimination' molecules in those modules? ===========\n")
mm <- fread(file.path(OUT,"louvain_modules.csv"))
chk <- c("HINT2","GPT","ECHDC2","DBT","DPYSL3","S100A10","ANXA2",
         "SDHA","SUCLA2","HOGA1","FLNA","LUM","PXDN","SERPINF1","FBLN5")
print(mm[gene %in% chk][order(mod)], row.names=FALSE)

cat("\n--- composition of each Louvain module (top-degree + counts) ---\n")
print(mm[, .(n=.N, up=sum(dir=="UP"), down=sum(dir=="DOWN")), by=mod][order(-n)], row.names=FALSE)

cat("\n=========== 2. GSE96804 directions for the Table-3 molecules ===========\n")
g <- fread(file.path(OUT,"gse96804_limma.csv"))
print(g[gene %in% c("HINT2","GPT","ECHDC2","DBT","DPYSL3","S100A10","ANXA2"),
        .(gene, logFC=round(logFC,3), P=signif(P.Value,3), adj=signif(adj.P.Val,3))],
      row.names=FALSE)

cat("\n=========== 3. module-representative molecules in GSE96804 ===========\n")
print(g[gene %in% c("SDHA","SUCLA2","HOGA1","ECHDC2","GPT","FLNA","DPYSL3","S100A10",
                    "LUM","PXDN","SERPINF1","FBLN5","NDUFS2","NDUFS3","NDUFS8","HSPA9"),
        .(gene, logFC=round(logFC,3), P=signif(P.Value,3))][order(logFC)], row.names=FALSE)

cat("\n=========== 4. do the four main modules actually contain the claimed representatives? ===========\n")
for (m in sort(unique(mm$mod))) {
  gg <- mm[mod==m]$gene
  cat(sprintf("module %d (n=%d, up=%d/down=%d): %s\n", m, length(gg),
              sum(mm[mod==m]$dir=="UP"), sum(mm[mod==m]$dir=="DOWN"),
              paste(intersect(chk, gg), collapse=", ")))
}
