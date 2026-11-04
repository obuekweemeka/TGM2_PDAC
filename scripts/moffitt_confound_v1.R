# ============================================================================
# moffitt_confound_v1.R -- READ-ONLY, fast. Save in: scripts/
# 1. Moffitt association under both basal-score definitions, with CIs
# 2. Is the association confounded by cell quality? (uninformative fraction,
#    genes detected, cells per patient), via correlation and partial correlation
# ============================================================================
base <- normalizePath('.'); pr <- file.path(base, "data", "processed")
v4  <- readRDS(file.path(pr, "aim2_scoring", "malignant_cells_scored_v4.rds"))
mcm <- readRDS(file.path(pr, "aim2_scoring", "malignant_cell_metadata.rds"))
tr  <- unique(mcm[, c("pid","new_treatment")]); v4$trt <- tr$new_treatment[match(v4$pid, tr$pid)]
u <- v4[v4$trt == "Untreated", ]; x <- u[u$axis2_death_competence != 0, ]
x$cyto <- grepl("cytostatic", x$cell_state); p <- sort(unique(x$pid))
g <- function(d, col, idx) as.vector(tapply(d[[col]], d$pid, mean)[idx])
pm <- data.frame(pid = p,
  cyto       = 100 * as.vector(tapply(x$cyto, x$pid, mean)[p]),
  basal_all  = g(u, "Moffitt_basal", p),
  basal_inf  = g(x, "Moffitt_basal", p),
  class_all  = g(u, "Moffitt_classical", p),
  pct_uninf  = 100 * as.vector(tapply(u$axis2_death_competence == 0, u$pid, mean)[p]),
  n_cells    = as.vector(table(u$pid)[p]))
if ("n_genes" %in% names(u)) pm$genes <- g(u, "n_genes", p)
sp <- function(a, b) { ct <- suppressWarnings(cor.test(a, b, method = "spearman"))
  sprintf("rho = %6.3f, p = %.4f", ct$estimate, ct$p.value) }
boot <- function(a, b) { set.seed(1); r <- replicate(10000, { i <- sample(length(a), replace = TRUE)
  if (length(unique(a[i])) < 3) NA else suppressWarnings(cor(a[i], b[i], method = "spearman")) })
  r <- r[is.finite(r)]; sprintf("95%% CI %.3f to %.3f", quantile(r,.025), quantile(r,.975)) }

cat("=== 1. Moffitt association under both definitions ===\n")
cat("  basal over ALL untreated cells (as in Fig 3c):", sp(pm$cyto, pm$basal_all), "|", boot(pm$cyto, pm$basal_all), "\n")
cat("  basal over INFORMATIVE cells only:            ", sp(pm$cyto, pm$basal_inf), "|", boot(pm$cyto, pm$basal_inf), "\n")

cat("\n=== 2. Quality confound ===\n")
cat("  cytostatic-prone % vs uninformative %:", sp(pm$cyto, pm$pct_uninf), "\n")
cat("  basal (all) vs uninformative %:       ", sp(pm$basal_all, pm$pct_uninf), "\n")
cat("  basal (inf) vs uninformative %:       ", sp(pm$basal_inf, pm$pct_uninf), "\n")
cat("  cytostatic-prone % vs cells/patient:  ", sp(pm$cyto, pm$n_cells), "\n")
if ("genes" %in% names(pm)) {
  cat("  cytostatic-prone % vs mean genes:     ", sp(pm$cyto, pm$genes), "\n")
  cat("  basal (all) vs mean genes:            ", sp(pm$basal_all, pm$genes), "\n")
}
pcor <- function(a, b, z) { ra <- resid(lm(rank(a) ~ rank(z))); rb <- resid(lm(rank(b) ~ rank(z)))
  ct <- cor.test(ra, rb); sprintf("partial rho = %6.3f, p = %.4f", ct$estimate, ct$p.value) }
cat("\n  Moffitt association controlling for uninformative %:\n")
cat("    basal (all):", pcor(pm$cyto, pm$basal_all, pm$pct_uninf), "\n")
cat("    basal (inf):", pcor(pm$cyto, pm$basal_inf, pm$pct_uninf), "\n")
if ("genes" %in% names(pm)) cat("    controlling for mean genes:", pcor(pm$cyto, pm$basal_all, pm$genes), "\n")
cat("\n=== Per-patient table ===\n"); print(pm[order(pm$cyto), ], row.names = FALSE, digits = 3)
write.csv(pm, file.path(base, "results", "aim2", "moffitt_confound_table.csv"), row.names = FALSE)
cat("\nIf the partial correlation collapses toward zero, the Moffitt association is largely a\ncell-quality effect and must be reported as such. If it holds, the association stands.\n")
