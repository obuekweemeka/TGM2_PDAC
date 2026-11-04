# ============================================================================
# zerocell_diagnostic_v1.R -- READ-ONLY, fast. Save in: scripts/
# Why do 1,513 cells with a detected Axis 2 gene still score exactly zero?
# Hypothesis: ties. Low-expression genes sit in very large tie blocks, so their
# AVERAGE rank (what UCell uses) is far worse than their MIN rank (what the
# earlier check used). Compares min / average / max tie handling.
# ============================================================================
suppressPackageStartupMessages(library(Matrix))
base <- normalizePath('.'); pr <- file.path(base, "data", "processed")
v4  <- readRDS(file.path(pr, "aim2_scoring", "malignant_cells_scored_v4.rds"))
mcm <- readRDS(file.path(pr, "aim2_scoring", "malignant_cell_metadata.rds"))
tr  <- unique(mcm[, c("pid", "new_treatment")]); v4$trt <- tr$new_treatment[match(v4$pid, tr$pid)]
u <- v4[v4$trt == "Untreated", ]
sig <- readRDS(file.path(pr, "aim1_signatures", "four_axis_signatures_v1.rds"))
a2 <- unique(unlist(sig[c("axis2_proapoptotic","axis2_executors","axis2_anti_apoptotic","axis2_cytostasis")]))
m <- readRDS(file.path(pr, "aim2_scoring", "malignant_count_matrix.rds"))
gxc <- if ("TGM2" %in% colnames(m)) Matrix::t(m) else m; rm(m); invisible(gc())
a2 <- intersect(a2, rownames(gxc)); ng <- nrow(gxc)
cat("Genes in matrix:", ng, "| Axis 2 genes:", length(a2), "\n")

z <- u$cell_id[u$axis2_death_competence == 0]; z <- intersect(z, colnames(gxc))
sub <- gxc[, z, drop = FALSE]
has <- which(Matrix::colSums(sub[a2, , drop = FALSE] > 0) > 0)
cat("Zero-score cells with >=1 Axis 2 gene detected:", length(has), "\n\n")

set.seed(1); samp <- sample(has, min(400, length(has)))
res <- t(sapply(samp, function(j) {
  x <- as.numeric(sub[, j]); names(x) <- rownames(sub)
  rmin <- rank(-x, ties.method = "min"); ravg <- rank(-x, ties.method = "average")
  g <- a2[x[a2] > 0]
  c(best_min = min(rmin[g]), best_avg = min(ravg[g]), n_detected = sum(x > 0),
    lowest_val = min(x[g]), tie_block = sum(x == min(x[g])))
}))
res <- as.data.frame(res)
cat("Best rank of a detected Axis 2 gene, by tie handling:\n")
cat("  ties = min     :", paste(round(quantile(res$best_min, c(0,.5,1))), collapse = " / "), " (min / median / max)\n")
cat("  ties = average :", paste(round(quantile(res$best_avg, c(0,.5,1))), collapse = " / "), "\n")
cat(sprintf("\nCells where the best AVERAGE rank exceeds 1500: %d of %d (%.1f%%)\n",
            sum(res$best_avg > 1500), nrow(res), 100 * mean(res$best_avg > 1500)))
cat(sprintf("Cells where the best MIN rank exceeds 1500:     %d of %d (%.1f%%)\n",
            sum(res$best_min > 1500), nrow(res), 100 * mean(res$best_min > 1500)))
cat("\nExpression value of the detected Axis 2 gene (log-normalised):\n"); print(summary(res$lowest_val))
cat("\nSize of the tie block it sits in (genes sharing that exact value):\n"); print(summary(res$tie_block))
cat("\nVERDICT: if average ranks exceed 1500 in nearly all cells while min ranks do not,\n",
    "the zero scores are explained by tie handling at low expression, and the\n",
    "manuscript can state the maxRank explanation correctly.\n")
