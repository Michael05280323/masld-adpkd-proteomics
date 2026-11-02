## ============================================================================
## Download the public datasets used by the analysis.
## Files are written to data/  (not tracked by git).
## ============================================================================
options(timeout = 7200)          # kidney proteinGroups.txt is ~250 MB

dir.create("data",    showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

files <- list(
  list(url  = "https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/03/PXD051911/liver_protein_quant.txt",
       dest = "data/PXD051911__liver_protein_quant.txt",
       mb   = 10),
  list(url  = "https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/03/PXD051911/meta_data.txt",
       dest = "data/PXD051911__meta_data.txt",
       mb   = 1),
  list(url  = "https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/05/PXD056281/proteinGroups.txt",
       dest = "data/PXD056281__proteinGroups.txt",
       mb   = 240),
  list(url  = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE96nnn/GSE96804/matrix/GSE96804_series_matrix.txt.gz",
       dest = "data/GSE96804_series_matrix.txt.gz",
       mb   = 5)
)

for (f in files) {
  if (file.exists(f$dest)) {
    cat(sprintf("[skip] %s already present\n", f$dest))
    next
  }
  cat(sprintf("[get ] %s\n", f$url))
  ok <- tryCatch({
    download.file(f$url, f$dest, mode = "wb", quiet = FALSE)
    TRUE
  }, error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); FALSE })

  if (ok) {
    sz <- file.info(f$dest)$size / 1024^2
    cat(sprintf("[ok  ] %s (%.1f MB)\n", f$dest, sz))
    if (sz < 0.5 * f$mb) warning("file looks smaller than expected: ", f$dest)
  }
}

## ---- KEGG pathway tables (cached for 05_kegg_auc.R) ----
kegg <- list(
  list(url = "https://rest.kegg.jp/link/hsa/pathway",   dest = "results/kegg_link_hsa.tsv"),
  list(url = "https://rest.kegg.jp/list/pathway/hsa",   dest = "results/kegg_pathway_hsa.tsv")
)
for (f in kegg) {
  if (file.exists(f$dest)) { cat(sprintf("[skip] %s\n", f$dest)); next }
  cat(sprintf("[get ] %s\n", f$url))
  tryCatch(download.file(f$url, f$dest, mode = "wb", quiet = TRUE),
           error = function(e) cat("  ERROR:", conditionMessage(e), "\n"))
}

cat("\nDone. Files in data/:\n")
print(file.info(list.files("data", full.names = TRUE))[, "size", drop = FALSE])
cat("\nKEGG cache in results/:\n")
print(file.info(list.files("results", pattern = "^kegg_", full.names = TRUE))[, "size", drop = FALSE])
