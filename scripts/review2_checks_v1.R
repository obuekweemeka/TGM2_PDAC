# ============================================================================
# review2_checks_v1.R -- READ-ONLY. Save in: scripts/
# Settles the points raised in review 2 that need data, not wording:
#  1. Are the excluded zero-score cells truly devoid of Axis 2 gene detection,
#     or merely below the UCell maxRank cutoff?
#  2. Threshold sweep: is the majority pooled-only, or also per-patient?
#  3. TGM2 patient-level effect size with a confidence interval
#  4. Is gene-level expression available to replicate TGM2 in the Steele atlas?
#  5. TCGA Cox reference group, and how the bulk Axis 2 score was built
# ============================================================================
suppressPackageStartupMessages({library(Matrix); library(survival)})
base <- normalizePath('.'); pr <- file.path(base, "data", "processed")
sec <- function(t, e) { cat("\n===", t, "===\n"); tryCatch(e, error = function(err) cat("FAILED:", conditionMessage(err), "\n")) }

v4  <- readRDS(file.path(pr, "aim2_scoring", "malignant_cells_scored_v4.rds"))
mcm <- readRDS(file.path(pr, "aim2_scoring", "malignant_cell_metadata.rds"))
tr  <- unique(mcm[, c("pid", "new_treatment")]); v4$trt <- tr$new_treatment[match(v4$pid, tr$pid)]
u <- v4[v4$trt == "Untreated", ]; u$tied <- u$axis2_death_competence == 0
sig <- readRDS(file.path(pr, "aim1_signatures", "four_axis_signatures_v1.rds"))
a2genes <- unique(unlist(sig[c("axis2_proapoptotic","axis2_executors","axis2_anti_apoptotic","axis2_cytostasis")]))
cat("Axis 2 genes:", length(a2genes), "|", paste(a2genes, collapse = ", "), "\n")

m <- readRDS(file.path(pr, "aim2_scoring", "malignant_count_matrix.rds"))   # cells x genes, log-normalised
cells_are_rows <- "TGM2" %in% colnames(m)
gi <- intersect(a2genes, if (cells_are_rows) colnames(m) else rownames(m))
cat("Axis 2 genes present in matrix:", length(gi), "of", length(a2genes), "\n")
ci <- match(u$cell_id, if (cells_are_rows) rownames(m) else colnames(m)); ok <- !is.na(ci)
sub <- if (cells_are_rows) m[ci[ok], gi, drop = FALSE] else t(m[gi, ci[ok], drop = FALSE])
ndet_a2  <- Matrix::rowSums(sub > 0)
ndet_all <- if (cells_are_rows) Matrix::rowSums(m[ci[ok], ] > 0) else Matrix::colSums(m[, ci[ok]] > 0)
rm(m); invisible(gc())
uu <- u[ok, ]; tied <- uu$tied

sec("1. Zero-score cells: direct gene-level check", {
  cat(sprintf("Cells checked: %d | zero-score: %d\n", nrow(uu), sum(tied)))
  cat(sprintf("Zero-score cells with NO Axis 2 gene detected: %d (%.2f%%)\n", sum(ndet_a2[tied] == 0), 100 * mean(ndet_a2[tied] == 0)))
  cat(sprintf("Zero-score cells WITH at least one Axis 2 gene detected: %d (%.2f%%)\n", sum(ndet_a2[tied] > 0), 100 * mean(ndet_a2[tied] > 0)))
  if (any(tied & ndet_a2 > 0)) {
    cat("  Among those, Axis 2 genes detected per cell:\n"); print(summary(ndet_a2[tied & ndet_a2 > 0]))
    cat("  Their total genes detected (relevant to the maxRank = 1500 cutoff):\n"); print(summary(ndet_all[tied & ndet_a2 > 0]))
    cat(sprintf("  Of these, cells with > 1500 genes detected (so rank cutoff could explain the zero): %d\n",
                sum(ndet_all[tied & ndet_a2 > 0] > 1500)))
  }
  cat(sprintf("Total genes detected: zero-score median %.0f vs informative median %.0f\n", median(ndet_all[tied]), median(ndet_all[!tied])))
  cat(sprintf("Zero-score cells with > 1500 genes detected overall: %d (%.1f%%)\n", sum(ndet_all[tied] > 1500), 100 * mean(ndet_all[tied] > 1500)))
})

sec("2. Threshold sweep: pooled vs per-patient majority", {
  x <- uu[!tied, ]; ad <- x$axis2_death_competence
  sw <- seq(-10, 0, 0.5)
  r <- data.frame(threshold = sw,
                  pooled_pct = sapply(sw, function(t) 100 * mean(ad < t)),
                  patients_above_50 = sapply(sw, function(t) sum(tapply(ad < t, x$pid, mean) > 0.5)),
                  min_patient_pct = sapply(sw, function(t) round(100 * min(tapply(ad < t, x$pid, mean)), 1)))
  print(r, row.names = FALSE)
  cat("\nLowest threshold with pooled > 50%:", min(sw[r$pooled_pct > 50]),
      "| lowest with ALL 18 patients > 50%:", if (any(r$patients_above_50 == 18)) min(sw[r$patients_above_50 == 18]) else NA, "\n")
})

sec("3. TGM2 patient-level effect size and CI", {
  mm <- readRDS(file.path(pr, "aim2_scoring", "malignant_count_matrix.rds"))
  cr <- "TGM2" %in% colnames(mm); idx <- match(uu$cell_id, if (cr) rownames(mm) else colnames(mm))
  uu$tgm2 <- if (cr) as.numeric(mm[idx, "TGM2"]) else as.numeric(mm["TGM2", idx]); rm(mm); invisible(gc())
  x <- uu[!uu$tied, ]; x$cyto <- grepl("cytostatic", x$cell_state)
  pc <- tapply(x$tgm2[x$cyto], x$pid[x$cyto], mean); pd <- tapply(x$tgm2[!x$cyto], x$pid[!x$cyto], mean)
  p <- intersect(names(pc), names(pd)); d <- pc[p] - pd[p]
  tt <- t.test(d); w <- suppressWarnings(wilcox.test(pc[p], pd[p], paired = TRUE, conf.int = TRUE))
  cat(sprintf("n = %d patients | higher in cytostatic-prone: %d\n", length(d), sum(d > 0)))
  cat(sprintf("Mean paired difference %.3f (95%% CI %.3f to %.3f) log-normalised units\n", mean(d), tt$conf.int[1], tt$conf.int[2]))
  cat(sprintf("Median paired difference %.3f (Wilcoxon 95%% CI %.3f to %.3f), p = %.4f\n", median(d), w$conf.int[1], w$conf.int[2], w$p.value))
  dz <- mean(d) / sd(d); cat(sprintf("Standardised paired effect (Cohen dz) = %.2f\n", dz))
  pct <- tapply(x$tgm2 > 0, list(x$pid, x$cyto), mean)
  cat(sprintf("Detection difference, patient-level mean %.1f percentage points\n", 100 * mean(pct[, "TRUE"] - pct[, "FALSE"], na.rm = TRUE)))
})

sec("4. Steele atlas: is gene-level expression available for TGM2 replication?", {
  d <- file.path(pr, "aim2_validation"); print(list.files(d))
  hits <- list.files(pr, pattern = "steele|Steele|GSE155698", recursive = TRUE, full.names = TRUE)
  cat("\nAll Steele-related files on disk:\n"); print(hits)
  cat("\nIf none contains a gene-by-cell matrix, replicating TGM2 in Steele requires re-downloading GSE155698.\n")
})

sec("5. TCGA Cox reference group and the bulk Axis 2 score", {
  tcga <- readRDS(file.path(pr, "aim4_bulk", "tcga_paad_analysis_df.rds"))
  cat("Columns:", paste(names(tcga), collapse = ", "), "\n")
  tcga$grp <- factor(tcga$cytostatic_group, levels = c("High cytostatic", "Low cytostatic"))
  cx <- coxph(Surv(os_days, os_event) ~ grp, data = tcga)
  cat("\nCox term name:", names(coef(cx)), "| reference level:", levels(tcga$grp)[1], "\n")
  cat(sprintf("HR = %.3f means hazard in '%s' relative to '%s'\n", exp(coef(cx)), levels(tcga$grp)[2], levels(tcga$grp)[1]))
  print(summary(survfit(Surv(os_days, os_event) ~ grp, data = tcga))$table[, c("records","events","median")])
  a2 <- grep("axis2", names(tcga), value = TRUE)
  cat("\nAxis 2 columns in the TCGA object:", paste(a2, collapse = ", "), "\n")
  for (c1 in a2) cat(sprintf("  %-28s range %.3f to %.3f\n", c1, min(tcga[[c1]], na.rm = TRUE), max(tcga[[c1]], na.rm = TRUE)))
  cat("\nIf a single ssGSEA Axis 2 column is stored rather than component modules,\nthe bulk score is that column directly and the Methods should say so.\n")
})
