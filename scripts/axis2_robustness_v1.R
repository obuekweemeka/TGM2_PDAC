# ============================================================================
# axis2_robustness_v1.R -- Save in: scripts/. Writes two small CSVs.
#  1. maxRank table with BOTH denominators (informative cells and whole cohort)
#  2. Alternative Axis 2 construction: difference score vs the log-ratio
#  3. Gene-set specificity: observed contrast vs size- and detection-matched
#     random gene sets, on a cell subsample for speed
# ============================================================================
suppressPackageStartupMessages({library(Matrix); library(UCell)})
base <- normalizePath('.'); pr <- file.path(base, "data", "processed")
sec <- function(t, e) { cat("\n===", t, "===\n"); tryCatch(e, error = function(err) cat("FAILED:", conditionMessage(err), "\n")) }
dir.create(file.path(base, "results", "aim2"), recursive = TRUE, showWarnings = FALSE)

v4  <- readRDS(file.path(pr, "aim2_scoring", "malignant_cells_scored_v4.rds"))
mcm <- readRDS(file.path(pr, "aim2_scoring", "malignant_cell_metadata.rds"))
tr  <- unique(mcm[, c("pid", "new_treatment")]); v4$trt <- tr$new_treatment[match(v4$pid, tr$pid)]
u <- v4[v4$trt == "Untreated", ]
comp <- c("axis2_proapoptotic","axis2_executors","axis2_anti_apoptotic","axis2_cytostasis")
stopifnot(all(comp %in% names(u)))

sec("1. maxRank table with both denominators", {
  f <- file.path(base, "results", "aim2", "maxrank_sensitivity.csv")
  if (!file.exists(f)) { cat("Run maxrank_sensitivity_v1.R first\n") } else {
    x <- read.csv(f); N <- nrow(u)
    x$n_cytostatic <- round(x$informative * x$pooled_cytostatic / 100)
    x$pct_of_all_cells <- round(100 * x$n_cytostatic / N, 2)
    x$pct_unclassifiable <- round(100 * (N - x$informative) / N, 2)
    out <- x[, c("maxRank","informative","pct_informative","pooled_cytostatic","n_cytostatic","pct_of_all_cells","pct_unclassifiable","patient_median","patient_min","patient_max","patients_over_50")]
    print(out, row.names = FALSE)
    cat(sprintf("\nAmong informative cells the proportion spans %.2f to %.2f percentage points; as a share of all %d cells it spans %.2f to %.2f.\n",
                min(x$pooled_cytostatic), max(x$pooled_cytostatic), N, min(x$pct_of_all_cells), max(x$pct_of_all_cells)))
    write.csv(out, file.path(base, "results", "aim2", "maxrank_table_both_denominators.csv"), row.names = FALSE)
    cat("Saved maxrank_table_both_denominators.csv\n")
  }
})

sec("2. Difference score vs log-ratio", {
  P <- u$axis2_proapoptotic + u$axis2_executors; A <- u$axis2_anti_apoptotic + u$axis2_cytostasis
  ratio <- u$axis2_death_competence; diff <- P - A
  inf <- !(P == 0 & A == 0)
  cat(sprintf("Informative cells: %d | agreement on sign: %.2f%% | Spearman between scores: %.3f\n",
              sum(inf), 100 * mean((ratio[inf] < 0) == (diff[inf] < 0)), cor(ratio[inf], diff[inf], method = "spearman")))
  cat(sprintf("Cytostatic-prone, log-ratio: %.2f%% | difference score: %.2f%%\n", 100 * mean(ratio[inf] < 0), 100 * mean(diff[inf] < 0)))
  pr_ <- tapply(ratio[inf] < 0, u$pid[inf], mean); pd_ <- tapply(diff[inf] < 0, u$pid[inf], mean)
  cat(sprintf("Per-patient median: %.1f%% vs %.1f%% | patients above 50%%: %d vs %d of %d | patient-level Spearman %.3f\n",
              100 * median(pr_), 100 * median(pd_), sum(pr_ > .5), sum(pd_ > .5), length(pr_), cor(pr_, pd_, method = "spearman")))
})

sec("3. Gene-set specificity against matched random sets", {
  sig <- readRDS(file.path(pr, "aim1_signatures", "four_axis_signatures_v1.rds"))
  mods <- sig[comp]; sizes <- sapply(mods, length)
  m <- readRDS(file.path(pr, "aim2_scoring", "malignant_count_matrix.rds"))
  gxc <- if ("TGM2" %in% colnames(m)) Matrix::t(m) else m; rm(m); invisible(gc())
  set.seed(42); cells <- sample(intersect(u$cell_id, colnames(gxc)), min(8000, ncol(gxc)))
  gxc <- gxc[, cells, drop = FALSE]
  det <- Matrix::rowMeans(gxc > 0); pool <- names(det)[det > 0]
  obs_genes <- unlist(mods); dec <- cut(det[pool], quantile(det[pool], 0:10/10), include.lowest = TRUE, labels = FALSE)
  names(dec) <- pool
  nperm <- 100
  rand <- list()
  for (i in seq_len(nperm)) for (k in seq_along(mods)) {
    want <- dec[intersect(mods[[k]], pool)]
    pick <- unlist(lapply(want, function(d) sample(setdiff(pool[dec == d], obs_genes), 1)))
    rand[[paste0("p", i, "_", names(mods)[k])]] <- unname(pick)
  }
  cat("Scoring", length(rand), "random signatures on", ncol(gxc), "cells...\n")
  t0 <- Sys.time()
  sc <- ScoreSignatures_UCell(gxc, features = rand, maxRank = 1500, ncores = 1)
  cat(sprintf("done in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  pr_null <- sapply(seq_len(nperm), function(i) {
    g <- function(k) sc[, paste0("p", i, "_", names(mods)[k], "_UCell")]
    r <- log2((g(1) + g(2) + 0.001) / (g(3) + g(4) + 0.001)); r <- r[r != 0]
    100 * mean(r < 0)
  })
  obs <- u$axis2_death_competence[u$cell_id %in% cells]; obs <- obs[obs != 0]
  obs_pct <- 100 * mean(obs < 0)
  cat(sprintf("\nObserved cytostatic-prone (same cells): %.2f%%\n", obs_pct))
  cat(sprintf("Random matched sets: median %.2f%%, 2.5th-97.5th percentile %.2f to %.2f\n",
              median(pr_null), quantile(pr_null, .025), quantile(pr_null, .975)))
  cat(sprintf("Random sets as extreme as observed: %d of %d (empirical p = %.3f)\n",
              sum(pr_null >= obs_pct), nperm, (sum(pr_null >= obs_pct) + 1) / (nperm + 1)))
  write.csv(data.frame(perm = seq_len(nperm), pct_cytostatic = pr_null),
            file.path(base, "results", "aim2", "geneset_specificity_null.csv"), row.names = FALSE)
  cat("Saved geneset_specificity_null.csv\n")
  cat("\nNote: a null centred well below the observed value indicates the result is specific to the chosen genes;\na null centred near it indicates the contrast partly reflects generic detection differences.\n")
})
