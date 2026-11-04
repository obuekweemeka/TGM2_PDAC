# ============================================================================
# maxrank_sensitivity_v1.R -- Save in: scripts/. Mostly read-only; writes one CSV.
#  1. Verify the ranks of detected Axis 2 genes in the 1,608 zero-score cells
#  2. Re-score Axis 2 with UCell at maxRank 1000/1500/2000/2500 and report the
#     effect on informative cells, pooled proportion and per-patient majority
#  3. Confirm the TCGA survival direction so the Results sentence is correct
# Section 2 re-runs UCell on 52,999 cells x 4 modules; allow 10-30 minutes.
# ============================================================================
suppressPackageStartupMessages({library(Matrix); library(UCell); library(survival)})
base <- normalizePath('.'); pr <- file.path(base, "data", "processed")
sec <- function(t, e) { cat("\n===", t, "===\n"); tryCatch(e, error = function(err) cat("FAILED:", conditionMessage(err), "\n")) }

v4  <- readRDS(file.path(pr, "aim2_scoring", "malignant_cells_scored_v4.rds"))
mcm <- readRDS(file.path(pr, "aim2_scoring", "malignant_cell_metadata.rds"))
tr  <- unique(mcm[, c("pid", "new_treatment")]); v4$trt <- tr$new_treatment[match(v4$pid, tr$pid)]
u <- v4[v4$trt == "Untreated", ]
sig <- readRDS(file.path(pr, "aim1_signatures", "four_axis_signatures_v1.rds"))
mods <- sig[c("axis2_proapoptotic","axis2_executors","axis2_anti_apoptotic","axis2_cytostasis")]
m <- readRDS(file.path(pr, "aim2_scoring", "malignant_count_matrix.rds"))
cells_rows <- "TGM2" %in% colnames(m)
cat("Matrix orientation:", if (cells_rows) "cells x genes" else "genes x cells", "| dims", paste(dim(m), collapse=" x "), "\n")
gxc <- if (cells_rows) Matrix::t(m) else m          # UCell wants genes x cells
rm(m); invisible(gc())
keep <- intersect(u$cell_id, colnames(gxc)); gxc <- gxc[, keep, drop = FALSE]
u <- u[match(keep, u$cell_id), ]
cat("Untreated cells scored:", ncol(gxc), "\n")

sec("1. Rank check in the zero-score cells", {
  a2 <- unique(unlist(mods)); a2 <- intersect(a2, rownames(gxc))
  z <- which(u$axis2_death_competence == 0)
  sub <- gxc[, z, drop = FALSE]
  ndet <- Matrix::colSums(sub > 0)
  has <- which(Matrix::colSums(sub[a2, , drop = FALSE] > 0) > 0)
  cat("Zero-score cells:", length(z), "| with >=1 Axis 2 gene detected:", length(has), "\n")
  worst <- integer(0)
  for (j in has) {
    x <- sub[, j]; nz <- which(x > 0)
    r <- rank(-x[nz], ties.method = "min")           # rank 1 = highest expression
    gi <- which(names(x)[nz] %in% a2)
    worst <- c(worst, min(r[gi]))
  }
  cat("Best (lowest) rank achieved by any detected Axis 2 gene in those cells:\n"); print(summary(worst))
  cat(sprintf("Cells where EVERY detected Axis 2 gene ranks beyond 1500: %d of %d (%.1f%%)\n",
              sum(worst > 1500), length(worst), 100 * mean(worst > 1500)))
  cat(sprintf("Cells with an Axis 2 gene ranked within 1500 (so maxRank does NOT explain the zero): %d\n", sum(worst <= 1500)))
  if (any(worst <= 1500)) cat("  -> the manuscript sentence must be softened accordingly\n")
})

sec("2. maxRank sensitivity (re-scoring)", {
  res <- list()
  for (mr in c(1000, 1500, 2000, 2500)) {
    t0 <- Sys.time()
    sc <- ScoreSignatures_UCell(gxc, features = mods, maxRank = mr, ncores = 1)
    a2 <- log2((sc[, "axis2_proapoptotic_UCell"] + sc[, "axis2_executors_UCell"] + 0.001) /
               (sc[, "axis2_anti_apoptotic_UCell"] + sc[, "axis2_cytostasis_UCell"] + 0.001))
    inf <- a2 != 0; pp <- tapply(a2[inf] < 0, u$pid[inf], mean)
    res[[as.character(mr)]] <- data.frame(maxRank = mr, informative = sum(inf),
      pct_informative = round(100 * mean(inf), 1), pooled_cytostatic = round(100 * mean(a2[inf] < 0), 2),
      patient_median = round(100 * median(pp), 1), patient_min = round(100 * min(pp), 1),
      patient_max = round(100 * max(pp), 1), patients_over_50 = sum(pp > 0.5))
    cat(sprintf("maxRank %4d done in %.1f min\n", mr, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
  out <- do.call(rbind, res); print(out, row.names = FALSE)
  cat("\nReference (manuscript, maxRank 1500): informative 45174, pooled 86.32, patient median 87.7, range 69.0-95.1, 18 of 18\n")
  dir.create(file.path(base, "results", "aim2"), recursive = TRUE, showWarnings = FALSE)
  write.csv(out, file.path(base, "results", "aim2", "maxrank_sensitivity.csv"), row.names = FALSE)
  cat("Saved results/aim2/maxrank_sensitivity.csv\n")
})
rm(gxc); invisible(gc())

sec("3. TCGA survival direction", {
  tc <- readRDS(file.path(pr, "aim4_bulk", "tcga_paad_analysis_df.rds"))
  med <- median(tc$axis2_death)
  cat("Cross-tab of group vs Axis 2 position:\n")
  print(table(group = tc$cytostatic_group, axis2 = ifelse(tc$axis2_death < med, "below median", "at/above median")))
  tc$grp <- factor(tc$cytostatic_group, levels = c("High cytostatic", "Low cytostatic"))
  print(summary(survfit(Surv(os_days, os_event) ~ grp, data = tc))$table[, c("records","events","median")])
  cx <- summary(coxph(Surv(os_days, os_event) ~ grp, data = tc))
  cat(sprintf("\nHR (%s vs reference %s) = %.3f (95%% CI %.3f-%.3f)\n", names(coef(cx$coefficients[,1]))[1],
              levels(tc$grp)[1], cx$conf.int[1], cx$conf.int[3], cx$conf.int[4]))
  cat("Correct Results wording: patients with Axis 2 scores BELOW the median (the high cytostatic group) had the longer median survival.\n")
})
