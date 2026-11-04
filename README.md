# MASLD–ADPKD cross-organ proteomics — analysis code

Analysis code for the manuscript:

> **Integrative proteomic analysis of shared differential proteins between MASLD and ADPKD: a cross-organ bioinformatic study**

The code reproduces all count-based results, enrichments, the STRING/Louvain
network and the cross-cohort AUC analysis reported in the paper.

---

## Data

All data are public and are **not** redistributed in this repository.
Run `00_download_data.R` to fetch them into `data/`.

| Cohort | Source | Accession | Local file (in `data/`) |
|---|---|---|---|
| Liver (MASLD/MASH) | PRIDE | **PXD051911** | `PXD051911__liver_protein_quant.txt`, `PXD051911__meta_data.txt` |
| Kidney (ADPKD) | PRIDE | **PXD056281** | `PXD056281__proteinGroups.txt` |
| Validation transcriptome | GEO | **GSE96804** | `GSE96804_series_matrix.txt.gz` |

Direct URLs used by `00_download_data.R`:

```
https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/03/PXD051911/liver_protein_quant.txt
https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/03/PXD051911/meta_data.txt
https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/05/PXD056281/proteinGroups.txt
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE96nnn/GSE96804/matrix/GSE96804_series_matrix.txt.gz
```

---

## Requirements

- R 4.4.1 (see `sessionInfo.txt`)
- Packages: `data.table`, `limma`, `clusterProfiler`, `org.Hs.eg.db`,
  `pROC`, `igraph`, `ggplot2`, `ggrepel`, `patchwork`, `scales`,
  `systemfonts`, `GEOquery`, `hta20transcriptcluster.db`
- Internet access for KEGG enrichment (the script caches the KEGG tables) and
  for the STRING API (a pre-downloaded edge list can be supplied instead)

Install the Bioconductor packages once:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("limma", "clusterProfiler", "org.Hs.eg.db", "GEOquery",
                       "hta20transcriptcluster.db"))
install.packages(c("data.table", "pROC", "igraph", "ggplot2", "ggrepel",
                   "patchwork", "scales", "systemfonts"))
```

---

## Run order

Run from the repository root.

| Step | Script | Output |
|---|---|---|
| 0 | `00_download_data.R` | `data/` |
| 1 | `01_liver.R` | `results/liver_limma_full.csv`, `results/liver_DEP.csv` |
| 2 | `02_kidney.R` | `results/kidney_limma_full.csv`, `results/kidney_DEP.csv`, `results/kidney_mat.rds` |
| 3 | `03_joint.R` | `results/joint_*.csv` |
| 4 | `04_enrich.R` | GO enrichment printed to console; `results/concordant_sets.rds` |
| 5 | `05_kegg_auc.R` | KEGG enrichment, AUC table `results/auc_all_shared.csv` |
| 6 | `05b_network_modules.R` | `results/string_network.tsv`, `results/louvain_modules.csv` |
| 7 | `06_gse96804.R` | `results/gse96804_limma.csv` |
| 8 | `07_figures/00_build_fig9_inputs.R` | `results/figures/nodes.csv`, `edges.csv` |
| 9 | `07_figures/fig9_ppi_panels.R` | `results/figures/fig9_ppi_panels.pdf/.png` |

```bash
Rscript 00_download_data.R
Rscript 01_liver.R
Rscript 02_kidney.R
Rscript 03_joint.R
Rscript 04_enrich.R
Rscript 05_kegg_auc.R
Rscript 05b_network_modules.R
Rscript 06_gse96804.R
Rscript 07_figures/00_build_fig9_inputs.R
Rscript 07_figures/fig9_ppi_panels.R
```

`exploratory/` contains the sensitivity analyses referenced in the review
(kidney pipeline variants, AUC definition variants, module-membership checks).
They require the outputs of steps 1–3.

---

## Key parameters

| Step | Parameter |
|---|---|
| Gene symbols | multi-symbol names split on `;`, first symbol kept (liver and kidney alike); duplicates collapsed by mean |
| Detection filter | ≥70% of samples per group (liver); ≥50% (kidney) |
| Missing values | within-group median imputation |
| Normalisation | column-wise median centring |
| Differential analysis | `limma`, `eBayes(trend = TRUE)` |
| DEP threshold | nominal P < 0.05 and \|log2FC\| ≥ log2(1.5) |
| Sensitivity analysis | MASLD vs control re-filtered and re-imputed within that subset (5484 genes) |
| PCA | `prcomp(center = TRUE, scale. = TRUE)` |
| Kidney samples | real tissue only (DT/NT); 3 technical replicates merged by mean |
| STRING | v12, combined score ≥ 0.4 |
| Louvain | `igraph::cluster_louvain`, `set.seed(1)`, resolution 1.0 |
| Kidney AUC | un-normalised individual-level log2 LFQ (see Methods 1.6) |

---

## Notes

- The liver matrix distributed by PRIDE is on a **linear intensity** scale;
  `01_liver.R` applies log2 before modelling.
- `results/kidney_mat.rds` stores both `mi` (pre-collapse, for exploratory
  scripts) and `mi_gene` (gene-level, unique symbols). **Use `mi_gene`** in any
  new script; `03_joint.R` and `05_kegg_auc.R` already prefer it.
- Gene symbols must be treated identically in every script (`;` split, first
  symbol). Mixing rules gives a common background of 2846 instead of 2851.
- `PXD056281/proteinGroups.txt` also contains 3D cyst models (DC/NC). Only the
  real-tissue columns (DT/NT) are used.
- Nominal-P DEP counts vary by ~1% across reasonable implementations; the
  FDR-defined count (kidney: 1035 = 567 up / 468 down, 31.7% of the 3264
  modelled proteins) is stable when technical replicates are merged by mean.
  (Merging by median instead gives 1049 = 575/474; the manuscript uses the
  mean pipeline throughout, see `exploratory/02b_kidney_variants.R`.)
- Raw mass-spectrometry files are large (several GB each) and are **not**
  downloaded by `00_download_data.R`; only the processed matrices are needed.

---

## License

MIT — see `LICENSE`.

## Citation

If you use this code, please cite the manuscript above and the versioned
Zenodo archive (release v1.1):

> Michael05280323. (2026). Integrative proteomic analysis of shared
> differential proteins between MASLD and ADPKD: analysis code (v1.1).
> Zenodo. https://doi.org/10.5281/zenodo.23940679

The DOI for this release is **10.5281/zenodo.23940679**
(https://doi.org/10.5281/zenodo.23940679).
