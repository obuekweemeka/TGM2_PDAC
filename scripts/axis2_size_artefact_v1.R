# ============================================================================
# axis2_size_artefact_v1.R -- Save in: scripts/. Diagnostic, writes one CSV.
# Question: is the negative skew of Axis 2 driven by unequal module sizes
# (numerator 13 genes vs denominator 10) rather than by apoptosis biology?
#  A. Null with BALANCED sizes: if balanced nulls centre near 50% while the
#     unbalanced null centres near 85%, size asymmetry is the cause.
#  B. Size-robust reconstructions of the observed score, each with its own
#     matched null, to see whether any retains a specific signal.
# ============================================================================
suppressPackageStartupMessages({library(Matrix); library(UCell)})
base <- normalizePath('.'); pr <- file.path(base, "data", "processed")
v4  <- readRDS(file.path(pr, "aim2_scoring", "malignant_cells_scored_v4.rds"))
mcm <- readRDS(file.path(pr, "aim2_scoring", "malignant_cell_metadata.rds"))
tr  <- unique(mcm[, c("pid","new_treatment")]); v4$trt <- tr$new_treatment[match(v4$pid, tr$pid)]
u <- v4[v4$trt == "Untreated", ]
sig <- readRDS(file.path(pr, "aim1_signatures", "four_axis_signatures_v1.rds"))
comp <- c("axis2_proapoptotic","axis2_executors","axis2_anti_apoptotic","axis2_cytostasis")
mods <- sig[comp]; cat("Module sizes:", paste(names(mods), sapply(mods, length), collapse = " | "),
                       "\nNumerator", length(mods[[1]]) + length(mods[[2]]), "vs denominator", length(mods[[3]]) + length(mods[[4]]), "genes\n")
m <- readRDS(file.path(pr, "aim2_scoring", "malignant_count_matrix.rds"))
gxc <- if ("TGM2" %in% colnames(m)) Matrix::t(m) else m; rm(m); invisible(gc())
set.seed(42); cells <- sample(intersect(u$cell_id, colnames(gxc)), min(8000, ncol(gxc)))
gxc <- gxc[, cells, drop = FALSE]
det <- Matrix::rowMeans(gxc > 0); pool <- names(det)[det > 0]
dec <- cut(det[pool], quantile(det[pool], 0:10/10), include.lowest = TRUE, labels = FALSE); names(dec) <- pool
obs_genes <- unlist(mods)
pick_matched <- function(tmpl) unname(unlist(lapply(dec[intersect(tmpl, pool)],
                      function(d) sample(setdiff(pool[dec == d], obs_genes), 1))))
pct_neg <- function(num, den) { r <- log2((num + 0.001)/(den + 0.001)); r <- r[r != 0]; 100*mean(r < 0) }

cat("\n=== A. Null with balanced vs unbalanced module sizes ===\n")
nperm <- 60; feats <- list()
for (i in seq_len(nperm)) {
  for (k in 1:4) feats[[sprintf("unb%d_%d", i, k)]] <- pick_matched(mods[[k]])
  # balanced: all four modules the same size (5 genes), detection-matched to module 3
  for (k in 1:4) feats[[sprintf("bal%d_%d", i, k)]] <- pick_matched(mods[[3]])
}
cat("Scoring", length(feats), "signatures...\n"); t0 <- Sys.time()
sc <- ScoreSignatures_UCell(gxc, features = feats, maxRank = 1500, ncores = 1)
cat(sprintf("done in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units="mins"))))
g <- function(tag, i, k) sc[, sprintf("%s%d_%d_UCell", tag, i, k)]
unb <- sapply(seq_len(nperm), function(i) pct_neg(g("unb",i,1)+g("unb",i,2), g("unb",i,3)+g("unb",i,4)))
bal <- sapply(seq_len(nperm), function(i) pct_neg(g("bal",i,1)+g("bal",i,2), g("bal",i,3)+g("bal",i,4)))
cat(sprintf("Unbalanced random sets (13 vs 10 genes): median %.1f%% negative\n", median(unb)))
cat(sprintf("Balanced random sets   (10 vs 10 genes): median %.1f%% negative\n", median(bal)))
cat("If balanced sits near 50%% and unbalanced near 85%%, module size asymmetry is the driver.\n")

cat("\n=== B. Size-robust reconstructions of the observed score ===\n")
obs_sc <- ScoreSignatures_UCell(gxc, features = mods, maxRank = 1500, ncores = 1)
cn <- paste0(comp, "_UCell")
P <- obs_sc[, cn[1]] + obs_sc[, cn[2]]; A <- obs_sc[, cn[3]] + obs_sc[, cn[4]]
cat(sprintf("1. Original (sum, 13 vs 10):            %.2f%%  | matched null median %.1f%%\n", pct_neg(P, A), median(unb)))
Pm <- (obs_sc[, cn[1]] + obs_sc[, cn[2]])/2; Am <- (obs_sc[, cn[3]] + obs_sc[, cn[4]])/2
cat(sprintf("2. Mean of modules (scale-matched):     %.2f%%\n", pct_neg(Pm, Am)))
z <- function(x) { s <- sd(x); if (s == 0) x*0 else (x - mean(x))/s }
Pz <- z(obs_sc[, cn[1]]) + z(obs_sc[, cn[2]]); Az <- z(obs_sc[, cn[3]]) + z(obs_sc[, cn[4]])
d <- Pz - Az; d <- d[is.finite(d)]
cat(sprintf("3. Per-module z-scores, difference:     %.2f%% below zero\n", 100*mean(d < 0)))
set.seed(7)
eq <- list(pro = sample(mods[[1]], 5), exe = sample(mods[[2]], 5), anti = mods[[3]], cyto = mods[[4]])
eq_sc <- ScoreSignatures_UCell(gxc, features = eq, maxRank = 1500, ncores = 1)
cat(sprintf("4. Equal-size modules (10 vs 10 genes): %.2f%%  | matched null median %.1f%%\n",
            pct_neg(eq_sc[,"pro_UCell"]+eq_sc[,"exe_UCell"], eq_sc[,"anti_UCell"]+eq_sc[,"cyto_UCell"]), median(bal)))
write.csv(data.frame(perm = seq_len(nperm), unbalanced = unb, balanced = bal),
          file.path(base, "results", "aim2", "axis2_size_nulls.csv"), row.names = FALSE)
cat("\nSaved axis2_size_nulls.csv\nThe key comparison is line 4 against its balanced null: that is the test of whether\nany apoptosis-specific signal survives once module size is equalised.\n")
