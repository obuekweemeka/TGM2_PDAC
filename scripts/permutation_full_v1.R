# ============================================================================
# permutation_full_v1.R -- Save in: scripts/. Strengthens the specificity test.
#  1. 1,000 size-equalised permutations (batched), for a stable null
#  2. Variability of the OBSERVED value across 100 random equalisation draws
#  3. Cell-subsample sensitivity: observed value on three seeds and on all cells
#  4. Disk search for any Steele gene-level expression data
# Expect roughly 20-30 minutes for section 1.
# ============================================================================
suppressPackageStartupMessages({library(Matrix); library(UCell)})
base <- normalizePath('.'); pr <- file.path(base, "data", "processed")
dir.create(file.path(base, "results", "aim2"), recursive = TRUE, showWarnings = FALSE)
v4  <- readRDS(file.path(pr, "aim2_scoring", "malignant_cells_scored_v4.rds"))
mcm <- readRDS(file.path(pr, "aim2_scoring", "malignant_cell_metadata.rds"))
tr  <- unique(mcm[, c("pid","new_treatment")]); v4$trt <- tr$new_treatment[match(v4$pid, tr$pid)]
u <- v4[v4$trt == "Untreated", ]
sig <- readRDS(file.path(pr, "aim1_signatures", "four_axis_signatures_v1.rds"))
comp <- c("axis2_proapoptotic","axis2_executors","axis2_anti_apoptotic","axis2_cytostasis")
mods <- sig[comp]
m <- readRDS(file.path(pr, "aim2_scoring", "malignant_count_matrix.rds"))
gxc_all <- if ("TGM2" %in% colnames(m)) Matrix::t(m) else m; rm(m); invisible(gc())
gxc_all <- gxc_all[, intersect(u$cell_id, colnames(gxc_all)), drop = FALSE]
pct_neg <- function(num, den) { r <- log2((num + 0.001)/(den + 0.001)); r <- r[r != 0]; 100*mean(r < 0) }

SEED <- 42; set.seed(SEED)
cells <- sample(colnames(gxc_all), min(8000, ncol(gxc_all)))
gxc <- gxc_all[, cells, drop = FALSE]
det <- Matrix::rowMeans(gxc > 0); pool <- names(det)[det > 0]
dec <- cut(det[pool], quantile(det[pool], 0:10/10), include.lowest = TRUE, labels = FALSE); names(dec) <- pool
obs_genes <- unlist(mods)
pick <- function(tmpl) unname(unlist(lapply(dec[intersect(tmpl, pool)],
            function(d) sample(setdiff(pool[dec == d], obs_genes), 1))))

cat("=== 1. 1,000 size-equalised permutations ===\n")
NPERM <- 1000; BATCH <- 200; null <- numeric(0); t0 <- Sys.time()
for (b in seq_len(ceiling(NPERM/BATCH))) {
  idx <- ((b-1)*BATCH + 1):min(b*BATCH, NPERM); f <- list()
  for (i in idx) for (k in 1:4) f[[sprintf("p%d_%d", i, k)]] <- pick(mods[[3]])   # all 5 genes -> 10 vs 10
  sc <- ScoreSignatures_UCell(gxc, features = f, maxRank = 1500, ncores = 1)
  null <- c(null, sapply(idx, function(i) {
    g <- function(k) sc[, sprintf("p%d_%d_UCell", i, k)]
    pct_neg(g(1)+g(2), g(3)+g(4)) }))
  cat(sprintf("  %d of %d done (%.1f min elapsed)\n", length(null), NPERM, as.numeric(difftime(Sys.time(), t0, units="mins"))))
}
cat(sprintf("Null: median %.2f%%, 95%% interval %.2f to %.2f, max %.2f\n",
            median(null), quantile(null,.025), quantile(null,.975), max(null)))

cat("\n=== 2. Observed value across 100 equalisation draws ===\n")
obs <- numeric(100)
for (j in 1:100) {
  set.seed(1000 + j)
  eq <- list(a = sample(mods[[1]], 5), b = sample(mods[[2]], 5), c = mods[[3]], d = mods[[4]])
  s2 <- ScoreSignatures_UCell(gxc, features = eq, maxRank = 1500, ncores = 1)
  obs[j] <- pct_neg(s2[,"a_UCell"]+s2[,"b_UCell"], s2[,"c_UCell"]+s2[,"d_UCell"])
  if (j %% 25 == 0) cat("  ", j, "draws\n")
}
cat(sprintf("Observed across equalisation draws: median %.2f%%, range %.2f to %.2f\n", median(obs), min(obs), max(obs)))
ex <- sum(null >= median(obs))
cat(sprintf("Permutations reaching the median observed value: %d of %d (empirical p = %.4f)\n", ex, NPERM, (ex+1)/(NPERM+1)))
write.csv(data.frame(perm = seq_along(null), pct = null), file.path(base,"results","aim2","axis2_null_1000.csv"), row.names = FALSE)
write.csv(data.frame(draw = 1:100, pct = obs), file.path(base,"results","aim2","axis2_equalisation_draws.csv"), row.names = FALSE)

cat("\n=== 3. Cell-subsample sensitivity (original, unequalised modules) ===\n")
for (sd in c(42, 7, 2026)) {
  set.seed(sd); cc <- sample(colnames(gxc_all), 8000)
  a <- u$axis2_death_competence[match(cc, u$cell_id)]; a <- a[a != 0]
  cat(sprintf("  seed %4d: %.2f%%\n", sd, 100*mean(a < 0)))
}
a <- u$axis2_death_competence[u$axis2_death_competence != 0]
cat(sprintf("  all %d informative cells: %.2f%%\n", length(a), 100*mean(a < 0)))

cat("\n=== 4. Steele gene-level data on disk ===\n")
hits <- list.files(base, pattern = "steele|Steele|GSE155698|155698", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
if (!length(hits)) cat("No Steele files found anywhere under the project folder.\n") else
  for (h in hits) cat(sprintf("  %-70s %8.1f MB\n", sub(paste0(base,"/"), "", h), file.size(h)/1e6))
big <- list.files(file.path(base,"data"), pattern = "\\.(rds|h5|h5ad|mtx|gz)$", recursive = TRUE, full.names = TRUE)
big <- big[file.size(big) > 5e7]
cat("\nLarge data objects that might contain Steele expression:\n")
if (!length(big)) cat("  none over 50 MB\n") else for (h in big) cat(sprintf("  %-70s %8.1f MB\n", sub(paste0(base,"/"), "", h), file.size(h)/1e6))
