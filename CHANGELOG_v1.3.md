# v1.3 — robustness scripts and Methods clarifications

## Added analysis scripts

| Script | Purpose | Main outputs |
|---|---|---|
| `09_permutation.R` | Kidney-label permutation null (2000 permutations, `set.seed(20261104)`) | `results/supp_table_permutation.csv`, `results/concordance_null_2000.rds` |
| `10_geneSetTest.R` | Correlation-aware `limma::geneSetTest` for the KEGG OXPHOS set (hsa00190) | `results/supp_table_genesettest.csv` |
| `11_predefined_sets.R` | Both-cohort FDR core and predefined up-side/M6 KEGG gene sets | `results/supp_table_FDR_core.csv`, `results/supp_table_predefined_sets.csv` |
| `12_liver_covariate.R` | Liver MASH vs control adjusted for age, sex and BMI (29 complete cases) | `results/supp_table_liver_covariate.csv`, `results/liver_limma_adjusted.csv` |

## Manuscript v3.6 changes

- Methods 1.5: Louvain now explicitly **combined-score weighted, resolution=1, set.seed(1)**; module sizes **M4=35, M6=68, M3=107, M2=40**; PXDN (M3→M2) and HOGA1 (→M6) are manually assigned by functional annotation.
- Methods 1.4 / Results 2.3: permutation uses **set.seed(20261104)**; the null is explicitly labelled **approximate** (filtering/imputation used true labels); empirical P is reported as **P<0.001 (0.0005 for this seed; ~0.001 for independent seeds)**.
- Methods 1.5: `geneSetTest` parameters specified as **moderated t, alternative="down", ranks.only=FALSE**.
- Methods 1.6 / Results 2.6: STRING degree vs min(AUC) reported as **Spearman ρ=−0.043, n=350**.
- Table 1: strict core reported as **81 (63 up / 18 down)**.
- Table 2: added **n** column (35/68/107/40).

## Run order

```text
00_download_data.R
01_liver.R
02_kidney.R
03_joint.R
04_enrich_common_background.R
05_kegg_auc.R
05b_network_modules.R
06_gse96804.R
09_permutation.R
10_geneSetTest.R
11_predefined_sets.R
12_liver_covariate.R
```

## Expected key values

- `09_permutation.R`: observed 362; null median 10 (IQR 3–32); empirical P=0.0005 (seed 20261104).
- `10_geneSetTest.R`: liver 1.44×10⁻⁸, kidney 7.07×10⁻²¹, GSE96804 6.06×10⁻³⁷.
- `11_predefined_sets.R`: FDR core 57 (45 up / 12 down); supplementary predefined-set table.
- `12_liver_covariate.R`: adjusted DEP 493 (252 up / 241 down), FDR+1.5× 375; OXPHOS geneSetTest P=3.14×10⁻¹¹.
