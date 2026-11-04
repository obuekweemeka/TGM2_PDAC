## supplementary_tables_v2.R
## Writes Supplementary Tables 1 and 2 into output/tables/, checked against
## docs/numbers_ledger.md. Table 3 is written by supplementary_tables_v1.R.
## Run with the working directory set to the repository root. Reads only; writes only
## into output/tables/.

base <- normalizePath(".")
out  <- file.path(base, "output", "tables"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
res  <- file.path(path.expand("~"), "Desktop", "rasolutepaper1", "results")
if (!dir.exists(res)) res <- file.path(base, "results")
stopifnot(dir.exists(res))
find1 <- function(p) { f <- list.files(res, pattern = p, recursive = TRUE, full.names = TRUE)
  stopifnot(length(f) == 1); f }
chk <- function(what, got, want, tol) {
  ok <- length(got) == length(want) && all(abs(got - want) <= tol)
  cat(sprintf("  %-52s %s\n", what, if (ok) "OK" else paste("MISMATCH:", paste(signif(got, 5), collapse = ", "))))
  if (!ok) stop("Check failed (", what, "). Nothing written.")
}
cat(R.version.string, "\n")

## =============================================== Table 1: maxRank sensitivity
cat("\n=== SUPPLEMENTARY TABLE 1 ===\n")
t1 <- read.csv(find1("^maxrank_table_both_denominators\\.csv$"), stringsAsFactors = FALSE)
print(t1, row.names = FALSE)
chk("maxRank settings tested", nrow(t1), 4, 0)
chk("range of pooled arrest-associated percentage", range(t1$pooled_cytostatic), c(86.0, 86.9), .35)
chk("range of informative-cell share", range(t1$pct_informative), c(76.0, 88.1), .5)
chk("all 18 patients above 50% at every setting", unique(t1$patients_over_50), 18, 0)
chk("informative + unclassifiable = 100% at each setting",
    round(t1$pct_informative + t1$pct_unclassifiable, 0), rep(100, nrow(t1)), 1)

T1 <- data.frame(
  `UCell maxRank` = t1$maxRank,
  `Cells with an informative score` = t1$informative,
  `Informative cells (%)` = round(t1$pct_informative, 1),
  `Unclassifiable cells (%)` = round(t1$pct_unclassifiable, 1),
  `Arrest-associated, of informative cells (%)` = round(t1$pooled_cytostatic, 2),
  `Arrest-associated cells (n)` = t1$n_cytostatic,
  `Arrest-associated, of all 52,999 cells (%)` = round(t1$pct_of_all_cells, 2),
  `Per-patient median (%)` = round(t1$patient_median, 1),
  `Per-patient minimum (%)` = round(t1$patient_min, 1),
  `Per-patient maximum (%)` = round(t1$patient_max, 1),
  `Patients above 50%` = t1$patients_over_50,
  check.names = FALSE)
write.csv(T1, file.path(out, "Supplementary_Table_1.csv"), row.names = FALSE)
cat("  WRITTEN: output/tables/Supplementary_Table_1.csv\n")

## =============================================== Table 2: gene-set specificity
cat("\n=== SUPPLEMENTARY TABLE 2 ===\n")
eq  <- read.csv(find1("^axis2_equalisation_draws\\.csv$"), stringsAsFactors = FALSE)   # observed, 100 draws
n1k <- read.csv(find1("^axis2_null_1000\\.csv$"), stringsAsFactors = FALSE)            # 1000 permutations
sz  <- read.csv(find1("^axis2_size_nulls\\.csv$"), stringsAsFactors = FALSE)           # balanced vs unbalanced
obs <- median(eq$pct)
cat(sprintf("  observed median across %d equalisation draws: %.2f\n", nrow(eq), obs))
cat(sprintf("  null median across %d permutations:           %.2f\n", nrow(n1k), median(n1k$pct)))
cat(sprintf("  size-preserving null median (%d perms):       %.2f\n", nrow(sz), median(sz$unbalanced)))
cat(sprintf("  size-equalised null median (%d perms):        %.2f\n", nrow(sz), median(sz$balanced)))
reach <- sum(n1k$pct >= obs); pval <- (reach + 1) / (nrow(n1k) + 1)
cat(sprintf("  permutations reaching the observed value:     %d of %d, p = %.4f\n",
            reach, nrow(n1k), pval))

chk("observed median across 100 draws (ledger 87.9)", obs, 87.9, .15)
chk("null median (ledger 49.2)", median(n1k$pct), 49.2, .6)
chk("permutations reaching observed (ledger 4 of 1000)", reach, 4, 0)
chk("permutation p-value (ledger 0.005)", pval, 0.005, .0005)

T2 <- data.frame(
  `Gene sets` = c("Observed modules, size-equalised (10 vs 10)",
                  "Random, size-equalised (10 vs 10)",
                  "Random, original module sizes (13 vs 10)",
                  "Observed modules, original sizes (13 vs 10)"),
  `Draws or permutations` = c(nrow(eq), nrow(n1k), nrow(sz), nrow(sz)),
  `Median arrest-associated (%)` = round(c(obs, median(n1k$pct), median(sz$balanced),
                                           median(sz$unbalanced)), 2),
  `Minimum (%)` = round(c(min(eq$pct), min(n1k$pct), min(sz$balanced), min(sz$unbalanced)), 2),
  `Maximum (%)` = round(c(max(eq$pct), max(n1k$pct), max(sz$balanced), max(sz$unbalanced)), 2),
  check.names = FALSE)
write.csv(T2, file.path(out, "Supplementary_Table_2.csv"), row.names = FALSE)

T2b <- data.frame(
  Quantity = c("Observed value (median of 100 size-equalised draws)",
               "Permutations of random gene sets",
               "Permutations reaching or exceeding the observed value",
               "Permutation p-value"),
  Value = c(sprintf("%.1f%%", obs), nrow(n1k), reach, sprintf("%.3f", pval)))
write.csv(T2b, file.path(out, "Supplementary_Table_2_permutation_test.csv"), row.names = FALSE)
cat("  WRITTEN: output/tables/Supplementary_Table_2.csv and _permutation_test.csv\n")

## =============================================== summary
cat("\n=== output/tables/ ===\n")
f <- list.files(out, full.names = TRUE)
print(data.frame(file = basename(f), KB = round(file.size(f) / 1024, 1)), row.names = FALSE)
cat("\nAll ledger checks passed.\n")
