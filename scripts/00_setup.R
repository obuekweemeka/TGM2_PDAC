## 00_setup.R -- run first. Checks the environment and the required input data.
if (!requireNamespace('here', quietly = TRUE)) install.packages('here')
base <- normalizePath('.')  # repository root; run from the project directory
cat('Project root:', base, '\n')

need <- c('Seurat','Matrix','ggplot2','UCell','DESeq2','GSVA','survival','survminer',
          'magick','ggrepel','readxl','httr','jsonlite','TCGAbiolinks','liana','hdf5r')
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) cat('Missing packages:', paste(miss, collapse = ', '),
                      '\nRun renv::restore() to install the recorded versions.\n')

## Inputs the scripts expect under data/ (see data/README.md for how to obtain them)
inputs <- c('data/processed/aim2_scoring/malignant_cells_scored_v4.rds',
            'data/processed/aim2_scoring/malignant_cell_metadata.rds',
            'data/processed/aim2_scoring/malignant_count_matrix.rds',
            'data/processed/aim1_signatures/four_axis_signatures_v1.rds',
            'data/processed/aim1_signatures/DepMap_PDAC_expression.rds',
            'data/processed/aim4_bulk/depmap_pdac_analysis_df.rds',
            'data/processed/aim4_bulk/tcga_paad_analysis_df.rds',
            'data/processed/aim2_validation/steele_cell_scores.csv',
            'data/processed/aim3_spatial/spot_level_primary_PDAC.csv')
present <- file.exists(file.path(base, inputs))
print(data.frame(input = inputs, present = present), row.names = FALSE)
if (!all(present)) cat('\nMissing inputs above must be downloaded or regenerated before the',
                       'scripts that use them will run.\n')
set.seed(1)  # scripts that draw random numbers set their own seed; this is the default
