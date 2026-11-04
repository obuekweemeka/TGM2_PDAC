## run_all.R -- runs the analysis in dependency order.
## Scripts keep the names used in docs/numbers_ledger_v15.md, which records the source
## script for every statistic in the manuscript. The order below is the dependency order,
## not alphabetical.
source(here::here('scripts', '00_setup.R'))

run <- function(f) { cat('\n==== ', f, ' ====\n'); source(here::here('scripts', f), echo = FALSE) }

## Stage 1: gene modules and resistance transcriptomics (Supplementary Figs 3, 4)
## Stage 2: single-cell scoring, thresholds, replication (Figs 1, 2; Supp Figs 5, 6, 9)
run('aim2_checks_v1.R')
run('axis2_robustness_v1.R')
run('axis2_size_artefact_v1.R')
run('permutation_full_v1.R')
run('maxrank_sensitivity_v1.R')
run('zerocell_diagnostic_v1.R')
run('moffitt_confound_v1.R')
run('moffitt_ci_and_steele_tgm2_v1.R')
run('steele_tgm2_replication_v2.R')
run('steele_patientlevel_and_fig3_v1.R')
## Stage 3: TGM2 and ligand-receptor inference (Fig 3)
run('liana_effect_check_v1.R')
## Stage 4: spatial analysis (Supplementary Figs 7, 8)
run('spatial_bh_recompute_v2.R')
## Stage 5: bulk tumour and cell-line analyses (Supp Fig 8c; Supp Table 3)
run('tcga_axis2_rebuild_v1.R')
run('purity_discrepancy_v1.R')
run('gdsc_axis2_response_v1.R')
run('gdsc_bh3_pdac_v1.R')
run('gdsc_bh3_pancancer_v1.R')
run('methods_gaps_check_v1.R')
run('review2_checks_v1.R')
## Panel builders and verification scripts (add before the figure stage)
run('fig3_twopanel_supp8_v1.R')
run('figS2_build_v1.R')
run('fig4_tgm2_build_v2.R')
run('final_figures_v2.R')
run('km_fig4a_v2.R')
run('purity_refetch_fig5d_v1.R')
run('module_provenance_check_v1.R')
run('optionA_and_fig3_checks_v1.R')
## Stage 6: figures
run('resistance_figures_v1.R')   # Supplementary Figures 3 and 4
run('final_figure_renumber_v2.R') # renames panels to manuscript numbering
run('figures_fix_v1.R')           # regenerates the final figures; must run last
cat('\nDone. Outputs in output/.\n')

## NOTE: verify this order against your own scripts before depositing. A script listed
## here that does not exist, or one that exists but is not listed, will show up when you
## run this file.
