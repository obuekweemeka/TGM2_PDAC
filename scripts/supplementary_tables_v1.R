## supplementary_tables_v1.R
## Assembles Supplementary Tables 1, 2 and 3 into output/tables/ from the analysis outputs,
## and checks each against docs/numbers_ledger.md. Run with the working directory set to
## the repository root. Reads only; writes only into output/tables/.

base <- normalizePath(".")
out  <- file.path(base, "output", "tables")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
cat(R.version.string, "\nRepository:", base, "\n")

search_dirs <- c(file.path(base, "results"), file.path(base, "output"),
                 file.path(path.expand("~"), "Desktop", "rasolutepaper1", "results"))
search_dirs <- search_dirs[dir.exists(search_dirs)]
all_csv <- unlist(lapply(search_dirs, list.files, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE))
all_rds <- unlist(lapply(search_dirs, list.files, pattern = "\\.rds$", recursive = TRUE, full.names = TRUE))
cat("Searched:", paste(search_dirs, collapse = ", "), "\n")
cat("CSV found:", length(all_csv), "| RDS found:", length(all_rds), "\n")

show <- function(f) { x <- if (grepl("\\.rds$", f)) readRDS(f) else read.csv(f, stringsAsFactors = FALSE)
  x <- as.data.frame(x); cat("\n----", basename(f), "|", nrow(x), "x", ncol(x), "----\n")
  print(names(x)); print(utils::head(x, 3)); invisible(x) }

## ------------------------------------------------------------------ Table 1
cat("\n================ SUPPLEMENTARY TABLE 1: maxRank sensitivity ================\n")
c1 <- c(all_csv, all_rds)
c1 <- c1[grepl("maxrank|max_rank|sensitivity|classifier", c1, ignore.case = TRUE)]
if (!length(c1)) cat("No candidate file found. Re-run scripts/maxrank_sensitivity_v1.R and note what it writes.\n")
for (f in c1) show(f)

## ------------------------------------------------------------------ Table 2
cat("\n================ SUPPLEMENTARY TABLE 2: gene-set specificity ================\n")
c2 <- c(all_csv, all_rds)
c2 <- c2[grepl("permut|null|random|specificity|equalis|equaliz", c2, ignore.case = TRUE)]
if (!length(c2)) cat("No candidate file found. Re-run scripts/permutation_full_v1.R and note what it writes.\n")
for (f in c2) show(f)

## ------------------------------------------------------------------ Table 3
cat("\n================ SUPPLEMENTARY TABLE 3: drug sensitivity ================\n")
t3 <- c(all_csv[grepl("Supplementary_Table_3", all_csv)],
        all_csv[grepl("gdsc_axis2_drug_response", all_csv)])
if (length(t3)) {
  s3 <- read.csv(t3[1], stringsAsFactors = FALSE)
  cat("Using:", t3[1], "|", nrow(s3), "rows\n")
  if (!"drug_set" %in% names(s3)) {
    cat("NOTE: this is the raw per-drug file, not the assembled table. Re-run the\n",
        "      Supplementary Table 3 assembly block from the earlier session first.\n")
  } else {
    stopifnot(all(s3$n_lines >= 10))
    pr_ <- grepl("primary", s3$drug_set)
    ok <- all(!is.na(s3$padj_BH_within_primary_set[pr_])) &&
          all(is.na(s3$padj_BH_within_primary_set[!pr_]))
    cat("  compounds:", nrow(s3), "| all n >= 10: TRUE | adjusted p only in primary sets:", ok, "\n")
    print(table(s3$drug_set))
    stopifnot(ok, nrow(s3) == 272, sum(grepl("MAPK", s3$drug_set)) == 6,
              sum(grepl("BH3", s3$drug_set)) == 9, sum(grepl("Cytotoxic", s3$drug_set)) == 9)
    write.csv(s3, file.path(out, "Supplementary_Table_3.csv"), row.names = FALSE)
    cat("  WRITTEN: output/tables/Supplementary_Table_3.csv (ledger checks passed)\n")
  }
} else cat("Supplementary_Table_3.csv not found.\n")

cat("\n================ LEDGER TARGETS for Tables 1 and 2 ================\n")
cat("Table 1 must show, across maxRank 1000-2500:\n")
cat("  overall arrest-associated percentage between 86.0 and 86.9\n")
cat("  informative-cell share between 76.0 and 88.1 per cent\n")
cat("  every patient above 50 per cent at each setting\n")
cat("Table 2 must show:\n")
cat("  size-equalised null median 49.2, observed median 87.9 across 100 draws\n")
cat("  4 of 1000 permutations reaching the observed value, p = 0.005\n")
cat("\nSend this output and I will write the two assembly blocks.\n")
