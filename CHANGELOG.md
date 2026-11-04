# Changelog

## v1.1 — 2026-11-04

Bug fixes that make the published counts reproducible without relying on
undocumented behaviour.

### Fixed

- **Inconsistent gene-symbol handling.** `01_liver.R` split multi-symbol names
  on `;` and kept the first symbol, but `05_kegg_auc.R` did not. The AUC script
  therefore intersected 2846 genes instead of the reported 2851.
  `05_kegg_auc.R` now uses the same rule as `01_liver.R`.
- **Pre-collapse kidney matrix.** `02_kidney.R` saved the individual-level
  matrix before duplicate gene symbols were collapsed, so four duplicated
  symbols (MOCS2, RAB34, RPP14, TMPO) were taken from their first row.
  `kidney_mat.rds` now also stores `mi_gene` (gene-level, unique symbols);
  `03_joint.R` and `05_kegg_auc.R` prefer it, and `05_kegg_auc.R` contains a
  defensive collapse for old `.rds` files.

### Verified after the fix

| Quantity | Expected | Reproduced |
|---|---:|---:|
| Common quantified background | 2851 | 2851 |
| Concordant shared proteins | 362 (220 / 142) | 362 (220 / 142) |
| Strict core | 81 | 81 |
| Venn intersection | 96 | 96 |
| Kidney DEP | 1265 (651 / 614) | 1265 (651 / 614) |
| Kidney FDR + 1.5x | 1035 (567 / 468) | 1035 (567 / 468) |
| Table 3 AUC, liver/kidney | — | max abs. difference 0.00 |

### Documentation

- README now lists the gene-symbol rule, `trend = TRUE`, PCA `scale = TRUE`
  and the sensitivity-analysis re-filtering as explicit parameters.

### Not changed

- No analysis result, figure or conclusion changed. Table 3 molecules were
  unaffected because none had multi-symbol names or duplicated kidney symbols.
