# Arrest-associated transcriptional programmes in treatment-naive PDAC

R code for the manuscript *Arrest-associated transcriptional programmes predominate over
apoptotic programmes in treatment-naive pancreatic ductal adenocarcinoma*.

The study asks how malignant cells in untreated PDAC are weighted between
apoptosis-associated and cell-cycle-arrest-associated transcriptional programmes, using
public single-nucleus, single-cell, spatial and bulk datasets.

## What is here, and what is not

All analyses are in R. The scripts regenerate the reported statistics from the deposited
processed objects. Two things they do not do:

- They do not download the raw data. Accessions and retrieval steps are in `data/README.md`;
  TCGA and GEO downloads are large and are left as an explicit step.
- A small number of figure panels were produced in a separate session from an intermediate
  object too large to deposit. These are noted in the figure map below.

## Requirements

R 4.3 or later. Package versions are recorded in `renv.lock`:

```r
install.packages('renv')
renv::restore()
```

Bioconductor packages (DESeq2, GSVA, TCGAbiolinks) and GitHub packages (UCell, liana) are
included in the lockfile.

## Running the analysis

```r
source('scripts/00_setup.R')   # checks packages and input data
source('run_all.R')            # runs the analysis in dependency order
```

Scripts use `here::here()` for all paths, so the project runs from any directory.

Script names match those recorded in the manuscript's numbers ledger, which lists the
source script for every reported statistic. They are therefore not renamed into a numbered
sequence; the run order is in `run_all.R`.

## Data

| Dataset | Source | Used for |
|---|---|---|
| Discovery single-nucleus atlas | GEO GSE202051 | Figs 1, 2, 3; Supp Figs 5, 6, 9 |
| Replication single-cell atlas | GEO GSE155698 | Fig 2b, Fig 3a |
| Spatial transcriptomics | GEO GSE272362 | Supp Figs 7, 8a, 8b |
| Resistance transcriptomics | GEO GSE269985 | Supp Figs 3, 4 |
| Bulk tumour expression and survival | GDC TCGA-PAAD; Raphael et al. 2017 clinical | Supp Fig 8c |
| Tumour purity | GDC PanCanAtlas ABSOLUTE calls | purity analysis |
| CRISPR dependency and expression | DepMap 24Q2 | Fig 3c; Supp Fig 1 |
| Drug sensitivity | GDSC2 release 8.5 | Supp Table 3 |

See `data/README.md` for file-level detail.

## Figure and table map

| Output | Script |
|---|---|
| Figure 1 | figures_fix_v1.R (panels from aim2 scoring) |
| Figure 2 | figures_fix_v1.R; steele_patientlevel_and_fig3_v1.R |
| Figure 3 | figures_fix_v1.R; fig4_tgm2_build_v2.R; liana_effect_check_v1.R |
| Supplementary Figure 1 | figures_fix_v1.R; figS2_build_v1.R |
| Supplementary Figure 2 | steele_patientlevel_and_fig3_v1.R |
| Supplementary Figures 3, 4 | resistance transcriptomics scripts |
| Supplementary Figure 5 | figures_fix_v1.R |
| Supplementary Figure 6 | figures_fix_v1.R |
| Supplementary Figure 7 | figures_fix_v1.R (panel a header only; maps generated separately) |
| Supplementary Figure 8 | spatial_bh_recompute_v2.R; km_fig4a_v2.R |
| Supplementary Figure 9 | final_figures_v2.R |
| Supplementary Tables 1, 2 | axis2_robustness_v1.R; permutation_full_v1.R |
| Supplementary Table 3 | gdsc_axis2_response_v1.R; gdsc_bh3_pdac_v1.R |

Check this map against your own files before depositing; the build script prints any
mismatch.

## Reproducing specific results

| Result | Script |
|---|---|
| Spatial gradients and BH adjustment | `scripts/spatial_bh_recompute_v2.R` (asserts ledger values) |
| Bulk Axis 2 score reconstruction | `scripts/tcga_axis2_rebuild_v1.R` |
| Zero-score audit | `scripts/tcga_axis2_rebuild_v1.R`, part A |
| Cell-line contrast and median split | `scripts/methods_gaps_check_v1.R` |

## Licence and citation

Code is released under the MIT licence (`LICENSE`). Cite via `CITATION.cff`.
The source datasets carry their own terms; see `data/README.md`.
