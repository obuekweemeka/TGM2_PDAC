# Numbers ledger, v15

Save in: `~/Desktop/rasolutepaper1/docs/numbers_ledger_v15.md` (supersedes v1 to v14; v15 changes §8 spatial only, 2026-10-23)

**Rule for the rewrite:** every number in the manuscript must appear in this ledger. A number that is not here does not go in. A number marked *Open* must be resolved first, or left out.

**Status key**
- **Verified**: recomputed from the saved data during this audit; the source column names the script and section.
- **Record**: confirmed from an earlier session's printed output, not re-run in this audit.
- **External**: checked against published sources.
- **Open**: not yet verified.

---

## 1. Cohort and classification (Hwang et al. 2022, GSE202051)

| Quantity | Value | Scope | Source | Status |
|---|---|---|---|---|
| Malignant cells, total | 64,538 | 43 patients | aim2_checks_v1 / optionA s2 | Verified |
| Untreated malignant cells | 52,999 | 18 patients (U1-U18) | optionA s1 | Verified |
| Treated malignant cells | 11,539 | 25 patients: CRT 10,268; CRTl 410; Other 861 | optionA s2 | Verified |
| Axis 2 formula | log2((pro-apoptotic + executors + 0.001) / (anti-apoptotic + cytostasis + 0.001)), UCell scores | reconstructs stored values exactly (r = 1.000000, max diff 0) | Axis 2 reconstruction | Verified |
| Cells with no Axis 2 signal (score exactly 0) | 7,825 (14.8%) | untreated | aim2_checks_v1 s1 | Verified |
| Why the score is zero (gene-level check) | 6,217 cells (79.5%) have NO Axis 2 gene detected; 1,608 (20.5%) have a median of 3 detected yet all four modules still score zero | untreated | review2_checks_v1 s1 | Verified |
| Zero-score audit (reviewer) | all 7,825 zero-score cells have all four module scores = 0; 0 cells with equal non-zero sums; 0 cells with 0 < abs(score) < 1e-6; smallest non-zero abs(score) 0.00057 | 52,999 untreated | tcga_axis2_rebuild_v1.R part A | Verified 2026-10-23 |
| Rank check on those 1,608 cells (RESOLVED) | tie handling explains it: detected Axis 2 genes sit in tie blocks of equally expressed genes (median 1,656 genes), so their AVERAGE rank, which UCell uses, exceeds 1500 in 400 of 400 sampled cells, while optimistic (min) tie breaking exceeds it in only 6.2% | | zerocell_diagnostic_v1 | Verified |
| Expression level of those genes | median 1.10 log-normalised (range 0.70 to 1.34), i.e. at or near the lowest detected level | | zerocell_diagnostic_v1 | Verified |
| *Superseded*: "maxRank does NOT explain the zero" | that check used ties.method = min; UCell uses average ranks | | | Do not use |
| Genes detected, zero-score vs informative (median, analysed matrix) | 993 vs 1,549 | untreated | review2_checks_v1 s1 | Verified |
| *Superseded*: "no constituent gene detected in every zero-score cell"; medians 1,053 vs 1,647 | inferred from module scores, not the expression matrix | | | Do not use |
| Informative cells | 45,174 | untreated | optionA s1 | Verified |
| **Cytostatic-prone (Axis 2 < 0), informative cells** | **86.32%** | 45,174 cells | optionA s1; fig2 build | Verified |
| Death-competent | 13.68% | 45,174 cells | optionA s1 | Verified |
| Quadrants: RAS-indep. cytostatic / RAS-dep. cytostatic / RAS-indep. death-competent / RAS-dep. death-competent | 59.30% / 27.02% / 9.80% / 3.88% | 45,174 cells | optionA s1 | Verified |
| Per-patient cytostatic-prone | median 87.7%, mean 84.9%, range 69.0-95.1%; all 18 above 50% | 18 patients | optionA s1; fig2 build | Verified |
| Treated groups, cytostatic-prone (informative) | CRT 81.50% (n = 7,596); CRTl 91.88% (n = 320); Other 84.35% (n = 658); all treated 82.11% (n = 8,574) | treated | optionA s2 | Verified |
| No Axis 2 signal, treated groups | CRT 26.0%; CRTl 22.0%; Other 23.6% | treated | optionA s2 | Verified |
| Threshold sensitivity, median split | threshold -6.077; per-patient 28.5-61.9% | untreated informative | optionA s4 | Verified |
| Threshold sensitivity, 33rd percentile | threshold -7.305; per-patient 14.9-52.2% | untreated informative | optionA s4 | Verified |
| Threshold sweep, POOLED cells > 50% | at any threshold above -6.0 | untreated informative | review2_checks_v1 s2 | Verified |
| Threshold sweep, ALL 18 patients > 50% | only above -1.0; at -6.0 just 7 of 18 patients, lowest patient 28.5% | untreated informative | review2_checks_v1 s2 | Verified |
| maxRank sensitivity (re-scored) | 1000: 40,287 informative (76.0%), pooled 86.90%, patient median 88.4%, range 68.7-95.7 | | maxrank_sensitivity_v1 s2 | Verified |
| | 1500: 45,174 (85.2%), 86.32%, 87.7%, 69.0-95.1 | | | Verified |
| | 2000: 46,396 (87.5%), 86.19%, 87.5%, 69.1-95.2 | | | Verified |
| | 2500: 46,702 (88.1%), 85.98%, 87.4%, 69.4-95.1 | | | Verified |
| | all four settings: 18 of 18 patients above 50% | | | Verified |
| maxRank, whole-cohort denominator | cytostatic-prone as a share of all 52,999 cells: 66.06% (1000), 73.57% (1500), 75.45% (2000), 75.76% (2500); unclassifiable 23.99% to 11.88% | | axis2_robustness_v1 s1 | Verified |
| Difference score vs log-ratio | algebraically equivalent at the zero threshold (both sums positive), so identical classification is guaranteed; NOT an independent robustness check | | algebra; axis2_robustness_v1 s2 | Verified |
| **Gene-set specificity (size-equalised, 1,000 permutations)** | **null median 49.20%, 95% interval 20.27-79.38, max 96.56; observed median across 100 equalisation draws 87.90% (range 85.29-92.79); 4 of 1,000 permutations reached the observed median; empirical p = 0.0050** | | permutation_full_v1 | Verified |
| Cell-subsample sensitivity | 86.34% (seed 42), 85.96% (seed 7), 86.38% (seed 2026) on 8,000 cells; 86.32% on all 45,174 informative cells | | permutation_full_v1 s3 | Verified |
| Permutation stability across subsamples | seed 7: null median 51.21% (max 96.73), observed median 86.69%, 1 of 500, p = 0.0040. seed 2026: null median 50.09% (max 87.42), observed 86.90%, 1 of 500, p = 0.0040. Consistent with seed 42 | | fig4_twocohort_and_permstability_v1 | Verified |
| Figure 4a | rebuilt as a two-cohort paired plot: Hwang 18 patients (12 higher, p = 0.0047) and Steele 16 samples (13 higher, p = 0.0155); Figure4.pdf 823 KB | | same | Verified |
| *Superseded*: 60-permutation null (median 48.3, p = 0.016) | replaced by the 1,000-permutation run | | | Do not use |
| Size-mismatched null (uninformative comparator) | 13-vs-10-gene random sets median 84.1% (range 71.3-95.8), because UCell scores fall as signature size grows; this null is saturated and has little power | | axis2_size_artefact_v1 | Verified |
| Effect of equalising module size on the OBSERVED value | 86.34% to 86.81%: size asymmetry inflates the null, not the estimate | | axis2_size_artefact_v1 | Verified |
| Per-module z-scored difference | 52.58% below zero; expected, since centring each module removes the between-module level difference the score measures | | axis2_size_artefact_v1 | Verified |
| Classifier stability across scoring objects | v2, v3, v4 identical | all cells | Item 1a | Verified |
| UMAP dimensions | Harmony dims 1-20 | untreated | aim2_checks_v1 s7 | Verified |
| *Superseded*: 88.3% / 88.1% / 87.6% / 90.5% | counted zero-signal cells as cytostatic, or used mixed cohorts | | | Do not use |

## 2. Moffitt subtype

| Quantity | Value | Scope | Source | Status |
|---|---|---|---|---|
| Patient-level: cytostatic-prone % vs mean basal score | Spearman rho = -0.595, exact p = 0.0105, 95% bootstrap CI -0.848 to -0.131 (10,000 resamples); basal averaged over ALL untreated cells (as plotted in Fig 3c) | 18 patients | moffitt_confound_v1 | Verified |
| Same, basal over informative cells only | rho = -0.618, p = 0.0074, CI -0.865 to -0.156 (definitional variant; not the reported figure) | 18 patients | moffitt_confound_v1 | Verified |
| Leave-one-out influence | rho ranges -0.730 to -0.551; most influential patient U8 | | moffitt_confound_v1 | Verified |
| **Data-quality confound** | **cytostatic-prone % vs uninformative-cell fraction: rho = -0.591, p = 0.0112. Basal vs that fraction only rho = 0.342, p = 0.165. Partial rho controlling for uninformative fraction -0.519 (p = 0.027); controlling for mean genes -0.536 (p = 0.022). Association attenuates but persists; per-patient prevalence is itself quality-sensitive** | | moffitt_confound_v1 | Verified |
| Other quality checks | cytostatic-prone % vs cells per patient rho = 0.150 (p = 0.55); vs mean genes rho = 0.323 (p = 0.19) | | moffitt_confound_v1 | Verified |
| Same, pre-Option A (for reference only) | rho = -0.602, exact p = 0.0096 (asymptotic 0.0083) | 18 patients | aim2_checks_v1 s3 | Do not use |
| Cell-level: cytostatic-prone among basal-higher vs classical-higher cells | 86.1% (n = 24,864) vs 86.6% (n = 20,310) | untreated informative | optionA s3 | Verified |
| Patient subtype rule (Fig 2b colours) | basal if mean basal score > mean classical score | 18 patients | aim2_checks_v1 s4 | Verified |

## 3. Axis correlations

| Quantity | Value | Scope | Source | Status |
|---|---|---|---|---|
| Axis 3 vs Axis 4 (pre-specified merge test, Pearson; threshold 0.6) | Pearson 0.408 (Spearman 0.411): not met, axes kept separate | untreated, 52,999 | fig2 build; aim2_checks_v1 s2 | Verified |
| Axis 3 vs Axis 4, all cells | Pearson 0.429, Spearman 0.435 | 64,538 | aim2_checks_v1 s2 | Verified |
| Full Spearman matrix (Fig 2c) | as printed in the regenerated panel | untreated, 52,999 | fig2 build | Verified |

## 4. Replication (Steele et al. 2020, GSE155698)

| Quantity | Value | Scope | Source | Status |
|---|---|---|---|---|
| Malignant cells | 14,571 | 16 samples | optionA s5 | Verified |
| No Axis 2 signal | 13.3% | all cells | optionA s5 | Verified |
| **Cytostatic-prone, informative cells, per-sample median** | **71.6%** (range 37.1-90.8%) | 16 samples | optionA s5; fig2 build | Verified |
| Cytostatic-prone, informative cells, pooled | 67.7% | all informative cells | optionA s5 | Verified |
| Hwang per-patient median, same rule | 87.7% | 18 patients | fig2 build | Verified |
| *Superseded*: 58.5% | used Axis 2 < 0 with zero-signal cells counted as not cytostatic (the opposite rule to Hwang) | | optionA s5 | Do not use |
| Steele units | report as 16 samples (11A and 11B share a patient number; describing them as samples is accurate either way) | | | Resolved |

## 5. Framework construction and Axis 1 (DepMap 24Q2)

| Quantity | Value | Scope | Source | Status |
|---|---|---|---|---|
| DepMap PDAC panel | 44 lines: intersection of 51 with expression (after removing 2 duplicate IDs) and 46 with CRISPR | | Item 3a | Verified |
| Gene-level expression vs KRAS Chronos (Pearson) | 18,784 genes; 131 with r < -0.4; 92 with r > 0.4 | 44 lines | figS2 build | Verified |
| DepMap correlates module | 30 most negative genes (r -0.565 to -0.466); module vs KRAS Spearman -0.830; circular by construction | 44 lines | module_provenance s3; Item 3b | Verified |
| MAPK output module (curated ERK targets: DUSP4, DUSP6, ETV4, SPRY4, EGR1, FOSL1, CCND1, SPRED1) | per-gene r -0.058 to -0.302; none beyond -0.4; module vs KRAS Spearman -0.308, p = 0.042 (not circular) | 44 lines | module check; figS2 / rebuild v3 | Verified |
| Bypass potential module (EGFR, AXL, FGFR1, FGFR3, ETV1, ETV5, PIK3CA, MTOR) | all upregulated in the resistance contrast, padj < 0.05, log2FC 0.19 to 2.86 | Panc 02.03 contrast | module check | Verified |
| Axis 3 hypoxia (48) and OXPHOS (87) modules | identical to the Hallmark GSEA leading edges | | module_provenance s2 | Verified |
| Axis 4 full (157) and core (66) modules | union of 5 immune-set leading edges; core = genes in 2 or more sets | | module_provenance s2 | Verified |
| Axis 2 modules | curated BCL-2 family, executor caspase, cytostasis gene sets (gene lists in Methods) | | June 9 record; stored object | Verified |
| DepMap Axis 2 grouping | split at the 44-line median of the log-ratio score: 22 High / 22 Low | | Step 4c | Verified |
| Log-ratio vs ssGSEA Axis 2 score (DepMap) | Spearman 0.96 | 44 lines | Step 4a | Verified |

## 6. Resistance transcriptomics (GSE269985)

| Quantity | Value | Scope | Source | Status |
|---|---|---|---|---|
| Design | ~cell_line + condition + cell_line:condition; Panc 02.03 and PANC-1, 3 parental + 3 resistant each (n = 12) | | Item 2 | Verified |
| Significant interaction | 6,686 of 16,439 genes (40.7%) | | Item 2 | Verified |
| Stored "shared" results (Fig 1) | Panc 02.03 coefficient with shrunken log2FC (identical p-values; mean abs log2FC 81% of unshrunk) | | optionA s8 | Verified |
| Significant and concordant in both lines | 1,501 genes | | Item 2 | Verified |
| Axis 2 genes concordant in both lines | BCL2 up (1.59 / 2.13); PMAIP1 up (0.65 / 0.85); RB1 up (0.64 / 0.28); GADD45A up (0.90 / 0.59); BIRC5 down (-0.43 / -0.15); GADD45B down (-0.25 / -0.32) | log2FC Panc 02.03 / PANC-1 | module_provenance s1 | Verified |
| Line-specific or discordant | MCL1 up in Panc 02.03 only (0.70; PANC-1 0.11, ns); executor caspases and CYCS down in Panc 02.03 only; CDKN1A +1.04 vs -0.61; BBC3 +0.61 vs -0.35 | | module_provenance s1 | Verified |
| Hallmark GSEA (padj < 0.05) | EMT 2.02; IFN-alpha 1.92; IFN-gamma 1.90; TNF-alpha/NF-kB 1.90; inflammatory 1.84; allograft 1.80; IL6-JAK-STAT3 1.74; IL2-STAT5 1.68; hypoxia 1.67; apical junction 1.55; OXPHOS -1.66 (NES) | Panc 02.03 contrast | module_provenance s2 | Verified |

## 7. TGM2 and ligand-receptor inference

| Quantity | Value | Scope | Source | Status |
|---|---|---|---|---|
| LIANA subsample | 16,090 cells; malignant senders 500 per state | untreated | optionA s6 | Verified |
| LIANA methods | NATMI, SingleCellSignalR, Connectome, logFC; rank aggregate | | Item 6 | Verified |
| TGM2-ITGA4 rank-score ratio (death-competent / cytostatic aggregate rank) | term. dysfunctional CD8 14.59; quiescent Tregs 13.83; activated Tregs 13.66; CD8 T 11.97; Tregs 11.71; progenitor dysf. CD8 10.06; cDC2 10.01; dysfunctional CD8 9.56; plasma 8.61; cDC1 7.27; CD4 5.59; NK 5.34; T mito-hi 5.24; B 5.20; macrophage 3.39; mast 1.49 | | liana_effect s1 | Verified |
| Example underlying scores | aggregate rank 9.06e-4 (cytostatic) vs 1.32e-2 (death-competent), term. dysfunctional CD8 | | liana_effect s1 | Verified |
| TGM2-ITGA4 rank position | best 430 of 16,407 (cytostatic senders) vs 1,012 of 14,901 (death-competent senders) | | liana_effect s2 | Verified |
| NATMI specificity ratio (the measure the old Methods described) | 1.22 for every receiver | | optionA s7 | Verified |
| **TGM2 expression, cytostatic-prone vs death-competent** | **pooled detection 28.6% vs 21.0% (patient-level mean difference 6.0 percentage points); mean paired difference 0.076 log-normalised units, 95% CI 0.024 to 0.128, Cohen dz 0.73; median paired difference 0.055, Wilcoxon CI 0.020 to 0.130, p = 0.0047; higher in 12 of 18 patients** | untreated informative | review2_checks_v1 s3 | Verified |
| TGM2 selection | emerged from ligand-receptor ranking, then tested in the same cohort: a follow-up to discovery, not independent confirmation | | | Verified |
| **TGM2 replication in Steele (patient level, REPORTED)** | **13 of 15 patients; mean paired difference 0.114 (95% CI 0.013 to 0.216), paired Wilcoxon p = 0.01245 (report as 0.012; 0.013 was a double-rounding error, corrected in v41). 11A and 11B are one patient and are pooled** | | steele_patientlevel_and_fig3_v1 | Verified 2026-10-23 (figures_fix_v1) |
| Same, by specimen (not reported) | 13 of 16; 0.106 (0.010 to 0.202), p = 0.0155 | | same | Verified |
| Steele matching | 14,571 of 14,571 scored cells matched on sample+barcode (100%); 12,634 informative; detection 36.1% vs 30.8% | | same | Verified |
| Steele replication method | TGM2 counts normalised to CP10K then log1p; classifications NOT recomputed; TGM2 is in none of the gene modules, so the test is independent of the classification | | same | Verified |
| Cells per Steele sample (downloaded matrices) | 48,570 total across 16 samples, of which 14,571 are the scored malignant cells | | same | Verified |
| *Superseded*: means 0.560 / 0.431, p = 0.0040 | re-normalised an already log-normalised matrix | | | Do not use |
| ITGA4 in malignant cells | 0.2% vs 0.3%; p = 0.68 (receptor is immune-side, as expected) | | liana_effect s3 | Verified |
| Single-cell matrix | 64,538 cells x 22,164 genes; stored log-normalised (non-integer, range 0.69-5.99); orientation cells x genes | | matrix check | Verified |
| Harmony integration | malignant cells still group by patient after integration (Supp Fig 3a); UMAP is visualisation only, all prevalence results are summarised per patient | | FigureS3 | Verified |
| LIANA cytostatic senders with no Axis 2 signal | 79 of 659 matched (12.0%) | | fig2 build | Verified |
| Figure 4 checks | panel a 12 of 18, p = 0.0047; panel b ratio range 1.49-14.59 | | fig4_tgm2_build_v2 | Verified |
| *Superseded*: "enriched 14.6-fold" | a rank-score ratio, not a signalling magnitude | | | Do not use |
| TGM2 essentiality (DepMap) | median Chronos +0.192 (High) vs +0.122 (Low); 93.2% of lines > 0; Wilcoxon p = 0.33 | 44 lines, 22/22 | Step 4c | Verified |
| NT5E essentiality | +0.005 vs -0.017; 47.7% > 0; p = 0.44 | 44 lines | Step 4c | Verified |
| KRAS essentiality (positive control) | -2.084 vs -2.177; 0% > 0; p = 0.60 | 44 lines | Step 4c | Verified |

## 8. Spatial transcriptomics (GSE272362)

| Quantity | Value | Scope | Source | Status |
|---|---|---|---|---|
| Spots / samples | 35,458 spots; 10 samples; T9 excluded (9 epithelial anchors), leaving 9 | | Item 0; rebuild v3 | Verified |
| Programme scores | precomputed fges assay supplied with the published Seurat object (27 Bagaev-style signatures), not computed here; RCTD assignments (first_type, spot_class) also supplied with the object (spacexr is not installed) | | June 10-11 records; package check | Verified |
| Median per-sample Spearman, distance vs score | Treg 0.153; checkpoint 0.150; myeloid suppression 0.119; effector 0.185 | 9 samples | Step 3 / rebuild v3 | Verified |
| Wilcoxon signed-rank vs 0 (exact, two-sided) | V = 44, 44, 45, 44; p = 0.0078, 0.0078, 0.0039, 0.0078; positive in 8, 8, 9, 8 of 9; **BH across the 4 programmes = 0.0078 for all four** (BH across 3 suppressive programmes gives the same) | 9 samples (T9 excluded: 9 epithelial anchors vs >= 99 in all others) | spatial_bh_recompute_v1.R (reproduces Step 3b exactly) | Verified 2026-10-23 |
| Primary, T9 included (sensitivity only) | median 0.149, 0.150, 0.110, 0.168; p = 0.0098, 0.0098, 0.0039, 0.0098; BH 0.0098 for all four | 10 samples | spatial_bh_recompute_v1.R | Verified 2026-10-23 |
| Specificity: effector rho vs Treg rho across samples | Spearman 0.983, p = 4.96e-05 | 9 samples | Step 3b | Verified |
| **Singlet-only spots (REPORTED)** | **anchors = all epithelial spots; query = singlet non-epithelial spots. median 0.205, 0.194, 0.052, 0.268; V = 45, 45, 44, 45; positive in 9, 9, 8, 9 of 9; p = 0.0039, 0.0039, 0.0078, 0.0039; BH across 4 = 0.0052, 0.0052, 0.0078, 0.0052; paired vs all-spot p = 0.57, 0.57, 0.73, 0.57** | same 9 samples as primary | Console block 2026-10-23 on spatial_bh_recompute_v1_rho.csv; saved as spatial_singlet_T9excluded_v1.csv; to be folded into spatial_bh_recompute_v2.R | Verified 2026-10-23 |
| Singlet specificity | effector rho vs Treg rho across samples: Spearman 0.933, p = 0.00075 | 9 samples | same | Verified 2026-10-23 |
| *Superseded*: singlet analysis on 10 samples incl. T9 (p = 0.0195, 0.0137, 0.0371, 0.00977; BH 0.026, 0.026, 0.037, 0.026; paired p 0.70-0.85) | inconsistent sample set with primary | | | Do not use |
| *Superseded*: "all 0.016 after BH" (primary) and "all adjusted p < 0.04" (singlet) | BH monotonicity misapplied; no script | | | Do not use |
| RCTD spot classes | singlet 4,320; doublet certain 26,792; doublet uncertain 4,015; reject 331 | | Item 0b | Verified |
| Spot-type panel count | 35,041 = spots whose dominant RCTD type is one of the 7 displayed types (8,260 + 21,770 + 1,113 + 3,002 + 196 + 226 + 474) | | spot-type table | Verified |

## 9. Bulk survival (TCGA PAAD)

| Quantity | Value | Scope | Source | Status |
|---|---|---|---|---|
| Cohort and grouping | n = 146; median split of the stored axis2_death ssGSEA score (reproduces 146 of 146 groups; Axis 1+2 index only 128) | | KM check | Verified |
| Cox reference group and direction | High cytostatic = Axis 2 BELOW median (73/73 cross-tab, exact); it is the reference level and has the longer survival (691 vs 511 days); HR 1.386 is the hazard of Low cytostatic relative to it | | maxrank_sensitivity_v1 s3 | Verified |
| *Superseded*: "patients above the median Axis 2 score had longer survival" | direction reversed | | | Do not use |
| Bulk Axis 2 derivation from component ssGSEA modules | RESOLVED: 0.5 x [(proapoptotic + executors) - (anti_apoptotic + cytostasis)] ssGSEA scores (GSVA, normalize = TRUE) on log2(TPM+1) from GDC STAR files. Rebuild vs stored: Pearson 0.99999, Spearman 0.99996, slope 0.998, max abs diff 0.0029, median-split groups 146/146 identical (so all KM/Cox statistics are unchanged) | 146 | tcga_axis2_rebuild_v1.R | Verified 2026-10-23 |
| TCGA data provenance | expression: GDC STAR gene counts/TPM, 178 primary tumours, via TCGAbiolinks GDCquery/GDCprepare; survival: Raphael et al. 2017 curated paper_ columns; 146 with complete OS (75 events). UCSC Xena NOT used | | build session 2026-06-11; tcga_axis2_rebuild_v1.R | Verified 2026-10-23 |
| Median OS (High vs Low cytostatic) | 691 vs 511 days; events 31 vs 44 | 73 / 73 | km_fig4a | Verified |
| Cox HR (Low vs High) | 1.386 (95% CI 0.874-2.197); Wald p = 0.165 | | km_fig4a; p-value check | Verified |
| Log-rank p | 0.164 | | p-value check | Verified |
| **Axis 2 vs ABSOLUTE purity (RE-DERIVED)** | **Spearman r = 0.153, p = 0.078, n = 133.** PanCanAtlas ABSOLUTE table re-fetched from GDC (881 KB, 10,786 rows); matched on first 15 barcode characters | | purity_refetch / purity_discrepancy | Verified |
| Other ABSOLUTE metrics vs axes | strongest is purity vs Axis 1 r = 0.248 (p = 0.004); ploidy vs Axis 1 -0.249; nothing reaches 0.55 | n = 133 | purity_discrepancy_v1 | Verified |
| *Superseded and UNREPRODUCIBLE*: r = 0.550, p = 7.43e-13, n = 145 | no ABSOLUTE metric, no axis pairing and no stored file reproduces it; the plotted purity source cannot be identified. Claim and Fig 5d removed | | | Do not use |

## 9a. Cell-line apoptosis/arrest contrast (verified 2026-10-23)

| Quantity | Value | Scope | Source | Status |
|---|---|---|---|---|
| Definition of stored axis2_death | exactly 0.5 x [(mean pro-apoptotic + mean caspase/mitochondrial) - (mean anti-apoptotic + mean cytostasis)], each mean log2(TPM+1) over module genes present; lm slope 0.5, intercept 0, max residual 1e-15. Equivalent: mean of the two apoptosis module means minus mean of the two arrest module means | 44 lines | methods_gaps_check_v1 reconstruction, Console check | Verified |
| Median split used in all reported cell-line analyses | re-split at the 44-line median (0.2879): lower-scoring 22 = High cytostatic, higher-scoring 22 = Low cytostatic | 44 lines | fig4_build_v4, fig4_tgm2_build_v2, rebuild_figures_v3 l.101 | Verified |
| Spearman vs ssGSEA alternative | 0.9598 (44 lines); 0.9599 (46-line table) | 44 lines | Console check | Verified |
| Drug-analysis line set | both GDSC scripts restrict to the 44 lines with expression before matching | 44 -> 24 matched | gdsc_axis2_response_v1 l.24-25; gdsc_bh3_pdac_v1 l.21 | Verified |
| *Do not use*: stored cytostatic_group in depmap_pdac_analysis_df.rds | 46-line median split (0.2631); differs from the reported 44-line split for 1 line | | | Do not use |
| *Not reproducible*: axis2_death for KP-3 (ACH-000108) and L3.3 (ACH-000685) | absent from the saved expression object; enter no reported analysis; exclude from deposited object | | | Do not use |

## 9b. Drug response (GDSC2 release 8.5)

| Quantity | Value | Scope | Source | Status |
|---|---|---|---|---|
| Prespecification | primary set = 6 MAPK-pathway inhibitors; metric = ln(IC50); predicted direction NEGATIVE Spearman (lower Axis 2 = less sensitive = higher IC50); BH across the primary set | fixed before computation | gdsc_axis2_response_v1 header | Verified |
| Line matching | 24 of 44 DepMap PDAC lines matched GDSC by stripped cell line name | | gdsc_axis2_response_v1 | Verified |
| **Primary result** | **0 of 6 drugs in the predicted direction; median Spearman r = +0.163; best selumetinib r = 0.412, nominal p = 0.046, BH p = 0.278; none significant after correction** | 20-24 lines per drug | same | Verified |
| Per drug (rho vs ln IC50) | selumetinib 0.412; refametinib 0.351; trametinib 0.194; KRAS G12C inhibitor-12 0.131; SCH772984 0.012; ulixertinib 0.006 | | same | Verified |
| AUC as secondary metric | all positive, 0.203 to 0.555 | | same | Verified |
| Whole-panel context | 272 compounds with 10+ matched lines; 11 reach p < 0.05 vs 13.6 expected by chance | | same | Verified |
| Interpretation limit | viability assays score arrest as sensitivity, biasing cytostatic-prone lines toward LOWER ln(IC50), which is the direction observed; the assay cannot cleanly test the hypothesis | | | Reasoning, stated in text |
| **Second test (BH3 mimetics), run after the first was null** | **9 drugs, 10-24 lines each; 6 of 9 in predicted direction but median rho -0.034; none survive BH (best TW 37 rho -0.643 on 12 lines, BH p = 0.252)** | | gdsc_bh3_pdac_v1 | Verified |
| Cytotoxic chemotherapy (secondary, reported WITHOUT BH adjustment) | 9 drugs; 8 of 9 in predicted direction, median rho -0.186; none individually significant (best irinotecan p = 0.061). NOT independent: 4 of 9 are topoisomerase I agents | | same | Verified |
| Overall conclusion on drug response | two prespecified tests, both null; no evidence the signature predicts cell-line drug sensitivity. Power limit: 24 lines detects only |rho| > ~0.45 | | | Verified |

## 10. Other checks

| Quantity | Value | Source | Status |
|---|---|---|---|
| DUSP6 detection vs mean Axis 1, patient level | Spearman rho = 0.205, p = 0.41, n = 18; lowest Axis 1 is U18 (-4.58), not U2 | aim2_checks_v1 s6 | Verified |
| RASolute 302, overall (ITT) population | OS 13.2 vs 6.7 months (HR 0.40); PFS 7.2 vs 3.6 months; confirmed ORR 31.6% vs 11.2% | NEJM 2026; ASCO 2026; Revolution Medicines, Dana-Farber, CHMP reports | External |
| RASolute 302, RAS G12-mutant population | confirmed ORR 33.2% vs 11.8%; comparator median OS 6.6 months | same sources | External |
| Quemliclustat (ARC-8, all patients) | median OS 15.7 months; ORR 38.5% | Wainberg et al., Nat Med 2026 | External |
| GSE202051 sample count | 44 samples from 43 patients (T13 has two samples); untreated: 18 samples from 18 patients; CRT 15, CRTl 5, Other 6 samples | fig3_build_v1 | Verified |
| Axis 1 formula | log2((MAPK output + 0.001) / (bypass potential + 0.001)); reconstructs exactly (r = 1.000000). The Methods' 0.01 is wrong (max difference 3.3) | fig3_build_v1 | Verified |
| GK921 reference (Lv et al.) | title and PMID confirmed; authors and journal details not | | Open |

## 11. Software versions (for Methods)

| Package | Version |
|---|---|
| R | 4.6.0 (2026-04-24), macOS |
| UCell | 2.16.0 |
| Seurat / SeuratObject | 5.5.1 / 5.4.0 |
| harmony | 2.0.5 |
| DESeq2 | 1.52.0 |
| fgsea | 1.38.0 |
| GSVA | 2.6.3 |
| liana | 0.1.14 |
| survival / survminer | 3.8.9 / 0.5.2 |
| ggplot2 / patchwork / magick | 4.0.3 / 1.3.2 / 2.9.1 |
| Matrix / hdf5r | 1.7.5 / 1.3.12 |
| spacexr (RCTD) | not installed: RCTD output came with the published object |


## 14. Figure regeneration (figures_fix_v1.R, 2026-10-23)

| Item | Value | Source | Status |
|---|---|---|---|
| Supp Fig 1b module rho (44 lines, KRAS Chronos, Spearman) | DepMap correlates -0.830 (p < 0.001); MAPK output -0.308 (0.042); bypass -0.049; pro-apoptotic -0.364 (0.016); caspase/mitochondrial 0.274 (0.072); anti-apoptotic 0.037; cytostasis -0.067; hypoxia -0.022; OXPHOS depleted 0.162; innate core -0.065; innate full -0.063 | figures_fix_v1 part 1 | Verified; figure-only values |
| Supp Fig 5a counts | apoptosis-associated 6,181; arrest-associated 38,993; no signal 7,825 (sum 52,999) | figures_fix_v1 | Verified |
| GSE272362 source | Khaliq AM et al. Nat Genet 2024;56(11):2455-2465 (Masood lab, IU; sample IDs IU_PDA_*). The old figure credit 'Masood et al. 2024' named the senior author and is superseded | web check 2026-10-23 | Verified |
| Supp Fig 6 threshold label | 'Zero crossing (Axis 2 < 0)'; the zero threshold was not specified in advance | Methods | Verified |
| MAPK primary set, BH across the 6 prespecified inhibitors (now saved in results/aim4/gdsc_axis2_drug_response.csv, column padj_BH) | rho all positive (opposite to prediction): KRAS G12C inhibitor-12 0.131 (n 20, p 0.581, adj 0.871); refametinib 0.351 (n 23, p 0.101, adj 0.303); SCH772984 0.012 (n 24, p 0.956, adj 0.979); selumetinib 0.412 (n 24, p 0.0464, adj 0.278); trametinib 0.194 (n 24, p 0.362, adj 0.724); ulixertinib 0.006 (n 24, p 0.979, adj 0.979) | pre-deposit check 2026-10-23 | Verified |
| Pre-deposit file fixes | gdsc_bh3_pdac_v1.R line 14 comment corrected; deposit/ holds the 44-line DepMap table (22/22 at 0.2879); GSE269985 confirmed 12 human + 18 murine (6694C2, 6499C4, 6419C5) samples; Supp Table 3 cytotoxic rows carry no adjusted p | pre-deposit check 2026-10-23 | Verified |
| Supplementary Table 3 file | results/aim4/Supplementary_Table_3.csv: 272 compounds (6 MAPK primary, 9 BH3 primary, 9 cytotoxic secondary, 248 context), all n >= 10; adjusted p only for the two primary sets; source files agree on all 24 shared drugs. BH3: 6/9 negative, median r -0.034, min adjusted p 0.252 (TW 37: r -0.643, n 12, nominal p 0.028) | pre-deposit check 2026-10-23 | Verified |
