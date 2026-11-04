# Data

No raw data are stored in this repository. Every dataset is public. Processed objects
needed to regenerate the reported statistics are deposited alongside the code; place them
under `data/processed/` in the structure the scripts expect (listed by `scripts/00_setup.R`).

## Sources

| Dataset | Accession / release | Retrieval |
|---|---|---|
| Discovery single-nucleus atlas (Hwang et al. 2022) | GSE202051 | GEO supplementary files |
| Replication single-cell atlas (Steele et al. 2020) | GSE155698 | per-sample filtered matrices from the individual GSM sample records, not the series-level file |
| Spatial transcriptomics (Khaliq et al. 2024) | GSE272362 | GEO; Seurat object with RCTD assignments and precomputed fges assay, used as supplied |
| Resistance transcriptomics (Dilly et al. 2024) | GSE269985 | GEO; the 12 human samples only (Panc 02.03, PANC-1); murine KPC samples not used |
| Bulk tumour expression | TCGA-PAAD | NCI GDC via TCGAbiolinks; STAR gene counts, 178 primary tumours |
| Survival annotations | Raphael et al. 2017 | curated clinical columns from the TCGA PAAD marker paper |
| Tumour purity | GDC PanCanAtlas | TCGA_mastercalls.abs_tables_JSedit.fixed.txt (ABSOLUTE) |
| CRISPR dependency and expression | DepMap 24Q2 | depmap.org: CRISPRGeneEffect, OmicsExpressionProteinCodingGenesTPMLogp1, Model.csv |
| Drug sensitivity | GDSC2 release 8.5 | GDSC2_fitted_dose_response_27Oct23.xlsx |

## Notes

- The Steele matrices were obtained from individual sample records rather than the series
  supplementary file, because the series file does not provide per-sample barcodes in the
  form the scripts expect.
- Two DepMap lines present in an earlier version of the analysis table (KP-3, L3.3) are
  absent from the deposited expression matrix and are excluded; all reported cell-line
  analyses use the 44 lines with both expression and CRISPR data.
- Spatial analyses use the ten primary tumour sections; one section with fewer than ten
  tumour-epithelial anchor spots is excluded from the primary analysis and reported as a
  sensitivity check.
