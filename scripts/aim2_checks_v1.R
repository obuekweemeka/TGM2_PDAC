# ============================================================================
# aim2_checks_v1.R -- READ-ONLY verification of Figure 2, Supp Fig 3 and 5 claims.
# Save in: scripts/. Writes nothing to disk.
# ============================================================================
suppressPackageStartupMessages(library(Matrix))
base <- normalizePath('.'); a2 <- file.path(base, "data", "processed", "aim2_scoring")
v4  <- readRDS(file.path(a2, "malignant_cells_scored_v4.rds"))
mcm <- readRDS(file.path(a2, "malignant_cell_metadata.rds"))
tr  <- unique(mcm[, c("pid", "new_treatment")]); v4$trt <- tr$new_treatment[match(v4$pid, tr$pid)]
u   <- v4[v4$trt == "Untreated", ]; ad <- u$axis2_death_competence; z <- ad == 0
cyto <- grepl("cytostatic", u$cell_state); pids <- sort(unique(u$pid))

# ---------- 1. Axis 2 ties at exactly zero ----------
cat("=== 1. Axis 2 ties at exactly zero (untreated cells) ===\n")
cat(sprintf("n = %d | Axis 2 < 0: %.2f%% | = 0: %.2f%% | > 0: %.2f%% | classified cytostatic: %.2f%%\n",
            nrow(u), 100 * mean(ad < 0), 100 * mean(z), 100 * mean(ad > 0), 100 * mean(cyto)))
comp <- c("axis2_proapoptotic", "axis2_executors", "axis2_anti_apoptotic", "axis2_cytostasis")
cat(sprintf("Tied cells with all four Axis 2 module scores = 0: %.1f%%\n", 100 * mean(rowSums(u[z, comp] == 0) == 4)))
cat("State assigned to tied cells:\n"); print(table(u$cell_state[z]))
cat(sprintf("Median genes detected: tied cells %.0f vs untied cells %.0f\n", median(u$n_genes[z]), median(u$n_genes[!z])))
cat(sprintf("Cytostatic among informative (untied) cells: %.2f%% of %d cells\n", 100 * mean(ad[!z] < 0), sum(!z)))
pp <- data.frame(pid = pids, n = as.vector(table(u$pid)[pids]),
                 pct_tied = round(100 * tapply(z, u$pid, mean)[pids], 1),
                 pct_cyto_all = round(100 * tapply(cyto, u$pid, mean)[pids], 1),
                 pct_cyto_untied = round(100 * tapply(ad[!z] < 0, u$pid[!z], mean)[pids], 1))
print(pp[order(pp$pct_cyto_all), ], row.names = FALSE)
cat(sprintf("Per-patient median: all cells %.1f%% | informative cells %.1f%% (range %.1f-%.1f)\n",
            median(pp$pct_cyto_all), median(pp$pct_cyto_untied), min(pp$pct_cyto_untied), max(pp$pct_cyto_untied)))

# ---------- 2. Axis correlations (Fig 2c) and the pre-specified merge test ----------
ax <- c("axis1_ras_balance", "axis2_death_competence", "axis3_hypoxia", "axis4_innate_mimicry")
cat("\n=== 2. Spearman correlations, all 64,538 cells (Fig 2c) ===\n"); print(round(cor(v4[, ax], method = "spearman"), 3))
for (lab in c("all", "untreated")) {
  d <- if (lab == "all") v4 else u
  cat(sprintf("Axis 3 vs Axis 4, %s cells (n = %d): Pearson %.3f | Spearman %.3f  (pre-specified rule: merge if Pearson > 0.6)\n",
              lab, nrow(d), cor(d$axis3_hypoxia, d$axis4_innate_mimicry),
              cor(d$axis3_hypoxia, d$axis4_innate_mimicry, method = "spearman")))
}

# ---------- 3. Fig 2d: per-patient cytostatic % vs Moffitt basal ----------
pm <- data.frame(pid = pids, cyto = tapply(cyto, u$pid, mean)[pids],
                 basal = tapply(u$Moffitt_basal, u$pid, mean)[pids],
                 classical = tapply(u$Moffitt_classical, u$pid, mean)[pids],
                 cyto_untied = tapply(ad[!z] < 0, u$pid[!z], mean)[pids])
s1 <- cor.test(pm$cyto, pm$basal, method = "spearman")
s2 <- suppressWarnings(cor.test(pm$cyto, pm$basal, method = "spearman", exact = FALSE))
p1 <- cor.test(pm$cyto, pm$basal)
s3 <- cor.test(pm$cyto_untied, pm$basal, method = "spearman")
cat("\n=== 3. Fig 2d (figure shows r = -0.60, p = 0.008) ===\n")
cat(sprintf("Spearman rho = %.3f | exact p = %.4f | asymptotic p = %.4f | Pearson r = %.3f, p = %.4f\n",
            s1$estimate, s1$p.value, s2$p.value, p1$estimate, p1$p.value))
cat(sprintf("Informative cells only: Spearman rho = %.3f, p = %.4f\n", s3$estimate, s3$p.value))

# ---------- 4. Fig 2b: how might patient subtype have been assigned? ----------
pm$subtype_by_mean <- ifelse(pm$basal > pm$classical, "Basal", "Classical")
pm$pct_cells_basal_gt_classical <- round(100 * tapply(u$Moffitt_basal > u$Moffitt_classical, u$pid, mean)[pids], 1)
cat("\n=== 4. Fig 2b subtype candidates (figure order, highest cytostatic first) ===\n")
print(data.frame(pid = pm$pid, pct_cyto = round(100 * pm$cyto, 1), subtype_by_mean = pm$subtype_by_mean,
                 pct_cells_basal = pm$pct_cells_basal_gt_classical)[order(-pm$cyto), ], row.names = FALSE)

# ---------- 5. Supp Fig 5a thresholds vs stored values ----------
sens <- readRDS(file.path(a2, "classifier_sensitivity_df.rds")); sens$pid <- as.character(sens$pid)
cat("\n=== 5. Supp Fig 5a thresholds (cytostatic = Axis 2 below threshold) ===\n")
for (lab in c("all", "untreated")) {
  vals <- if (lab == "all") v4$axis2_death_competence else ad
  med <- median(vals); q33 <- unname(quantile(vals, 0.33))
  x <- data.frame(pid = pids, med_split = round(100 * tapply(ad < med, u$pid, mean)[pids], 1),
                  p33 = round(100 * tapply(ad < q33, u$pid, mean)[pids], 1))
  x <- merge(x, sens[, c("pid", "atlas_med", "p33")], by = "pid", suffixes = c("", "_stored"))
  cat(sprintf("\n%s-cell thresholds: median = %.3f, 33rd percentile = %.3f\n", lab, med, q33))
  x$atlas_med <- round(x$atlas_med, 1); x$p33_stored <- round(x$p33_stored, 1)
  print(x[order(x$med_split), ], row.names = FALSE)
}

# ---------- 6. Supp Fig 3e: DUSP6 detection vs Axis 1, patient level ----------
cat("\n=== 6. Supp Fig 3e: DUSP6 detection vs mean Axis 1 (patient level) ===\n")
m <- readRDS(file.path(a2, "malignant_count_matrix.rds"))
if (!("DUSP6" %in% rownames(m)) && ("DUSP6" %in% colnames(m))) m <- t(m)
if ("DUSP6" %in% rownames(m)) {
  idx <- match(u$cell_id, colnames(m)); ok <- !is.na(idx)
  cat("Untreated cells matched to the matrix:", sum(ok), "of", nrow(u), "\n")
  det <- as.numeric(m["DUSP6", idx[ok]] > 0)
  dd <- data.frame(pid = pids, pct_dusp6 = round(100 * tapply(det, u$pid[ok], mean)[pids], 1),
                   mean_axis1 = round(tapply(u$axis1_ras_balance[ok], u$pid[ok], mean)[pids], 2))
  print(dd[order(dd$mean_axis1), ], row.names = FALSE)
  ct <- suppressWarnings(cor.test(dd$pct_dusp6, dd$mean_axis1, method = "spearman"))
  cat(sprintf("Spearman rho = %.3f, p = %.4f, n = %d patients | lowest mean Axis 1: %s\n",
              ct$estimate, ct$p.value, nrow(dd), dd$pid[which.min(dd$mean_axis1)]))
} else cat("DUSP6 not found in the matrix dimnames\n")
rm(m); invisible(gc())

# ---------- 7. Supp Fig 3c: PCs actually used ----------
cat("\n=== 7. PCs used for Harmony / neighbours / UMAP ===\n")
tryCatch({
  s <- readRDS(file.path(a2, "seurat_untreated_harmony.rds"))
  for (cn in names(s@commands)) {
    p <- slot(s@commands[[cn]], "params")
    d <- if (!is.null(p$dims)) paste(range(p$dims), collapse = "-") else if (!is.null(p$dims.use)) paste(range(p$dims.use), collapse = "-") else "-"
    cat(sprintf("%-32s dims: %s\n", cn, d))
  }
}, error = function(e) cat("Could not read the Seurat command log:", conditionMessage(e), "\n"))
