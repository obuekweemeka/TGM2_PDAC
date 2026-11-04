## spatial_bh_recompute_v2.R
## Single source for every spatial statistic reported in manuscript v36 (ledger v15, section 8).
## Supersedes v1. Base R only. Reads the per-spot table (read-only); writes CSVs to results/aim3/spatial/.
##
## Method
##   anchors : spots with RCTD first_type == "Tumor Epithelial cells" (all spot classes)
##   query   : non-epithelial spots (primary) or non-epithelial singlet spots (sensitivity)
##   per section: Spearman rho between distance to nearest anchor and each programme score
##   inference: exact two-sided Wilcoxon signed-rank test of per-section rho vs 0;
##              Benjamini-Hochberg across the four programmes within each analysis
##   inclusion: sections with >= 10 epithelial anchors (excludes T9, 9 anchors);
##              T9-included primary analysis reported as a sensitivity check

base   <- normalizePath('.')
infile <- file.path(base, "data", "processed", "aim3_spatial", "spot_level_primary_PDAC.csv")
outdir <- file.path(base, "results", "aim3", "spatial")
cat(R.version.string, "\n")

spot  <- read.csv(infile, stringsAsFactors = FALSE)
progs <- c(Treg = "Treg", Checkpoint = "Checkpoint_molecules",
           MyeloidSupp = "Immune_Suppression_by_Myeloid_cells", Effector = "Effector_cells")
stopifnot(all(c("sample_id", "first_type", "spot_class", "imagecol", "imagerow", progs) %in% names(spot)))
spot$is_te <- spot$first_type == "Tumor Epithelial cells"
spot$is_sg <- spot$spot_class == "singlet"
by_s <- split(spot, spot$sample_id)

near_dist <- function(qx, qy, ax, ay) {
  d <- rep(Inf, length(qx))
  for (k in seq_along(ax)) d <- pmin(d, sqrt((qx - ax[k])^2 + (qy - ay[k])^2))
  d
}
rho_one <- function(d, anchor, query) {
  dist <- near_dist(d$imagecol[query], d$imagerow[query], d$imagecol[anchor], d$imagerow[anchor])
  sapply(progs, function(p) cor(dist, d[[p]][query], method = "spearman"))
}
n_anchor <- sapply(by_s, function(d) sum(d$is_te))
rho_all  <- t(sapply(by_s, function(d) rho_one(d, d$is_te, !d$is_te)))
rho_sg   <- t(sapply(by_s, function(d) rho_one(d, d$is_te, !d$is_te & d$is_sg)))

keep9  <- n_anchor >= 10
keep10 <- n_anchor >= 1
cat("\nAnchors per section:\n"); print(n_anchor)
cat("Excluded (< 10 anchors):", names(n_anchor)[!keep9], "\n")

wfam <- function(label, R, keep) {
  R <- R[keep, , drop = FALSE]
  o <- do.call(rbind, lapply(names(progs), function(p) {
    w <- suppressWarnings(wilcox.test(R[, p], mu = 0, exact = TRUE))
    data.frame(analysis = label, programme = p, n = nrow(R), median_rho = median(R[, p]),
               V = unname(w$statistic), n_pos = sum(R[, p] > 0), p = w$p.value)
  }))
  o$p_BH <- p.adjust(o$p, "BH")
  o
}
res <- rbind(wfam("primary",             rho_all, keep9),
             wfam("primary_T9_included", rho_all, keep10),
             wfam("singlet",             rho_sg,  keep9))
res$paired_vs_allspot_p <- NA_real_
i_sg <- res$analysis == "singlet"
res$paired_vs_allspot_p[i_sg] <- sapply(names(progs), function(p)
  suppressWarnings(wilcox.test(rho_sg[keep9, p], rho_all[keep9, p], paired = TRUE, exact = TRUE))$p.value)

spec <- data.frame(
  analysis = c("primary", "singlet"),
  rho = c(cor(rho_all[keep9, "Effector"], rho_all[keep9, "Treg"], method = "spearman"),
          cor(rho_sg[keep9, "Effector"],  rho_sg[keep9, "Treg"],  method = "spearman")),
  p   = c(suppressWarnings(cor.test(rho_all[keep9, "Effector"], rho_all[keep9, "Treg"], method = "spearman", exact = TRUE))$p.value,
          suppressWarnings(cor.test(rho_sg[keep9, "Effector"],  rho_sg[keep9, "Treg"],  method = "spearman", exact = TRUE))$p.value))

cat("\n=== Reported spatial statistics ===\n")
pr <- res; num <- sapply(pr, is.numeric); pr[num] <- lapply(pr[num], signif, 4); print(pr, row.names = FALSE)
cat("\n=== Specificity (effector rho vs Treg rho across sections) ===\n"); print(signif(spec[-1], 4))

## ---- assertions against ledger v15 (fail loudly if anything has drifted)
chk <- function(got, want, tol, what) {
  ok <- all(abs(got - want) <= tol)
  cat(sprintf("%-45s %s\n", what, if (ok) "OK" else paste("MISMATCH:", paste(signif(got, 4), collapse = ", "))))
  ok
}
g <- function(a, col) res[res$analysis == a, col]
ok <- c(
  chk(g("primary", "median_rho"), c(0.153, 0.150, 0.119, 0.185), 5e-4, "primary medians"),
  chk(g("primary", "V"),          c(44, 44, 45, 44), 0,              "primary V"),
  chk(g("primary", "p_BH"),       rep(0.0078125, 4), 1e-6,           "primary BH = 0.0078"),
  chk(g("primary_T9_included", "p_BH"), rep(0.009765625, 4), 1e-6,   "T9-included BH = 0.0098"),
  chk(g("singlet", "median_rho"), c(0.2047, 0.1935, 0.05178, 0.2677), 5e-4, "singlet medians"),
  chk(g("singlet", "p_BH"),       c(0.005208, 0.005208, 0.0078125, 0.005208), 1e-5, "singlet BH"),
  chk(g("singlet", "paired_vs_allspot_p"), c(0.5703, 0.5703, 0.7344, 0.5703), 1e-3, "singlet paired p"),
  chk(spec$rho, c(0.983, 0.933), 5e-4,                                "specificity rho"))
cat(if (all(ok)) "\nALL LEDGER CHECKS PASSED\n" else "\nLEDGER CHECK FAILED: do not use these values\n")

write.csv(res,  file.path(outdir, "spatial_stats_v2.csv"),       row.names = FALSE)
write.csv(spec, file.path(outdir, "spatial_specificity_v2.csv"), row.names = FALSE)
write.csv(data.frame(sample = rownames(rho_all), n_anchor, rho_all = rho_all, rho_singlet = rho_sg),
          file.path(outdir, "spatial_rho_per_section_v2.csv"), row.names = FALSE)
cat("Saved: spatial_stats_v2.csv, spatial_specificity_v2.csv, spatial_rho_per_section_v2.csv\n")

## ---- figure-label check: does any figure script print p-values on Supplementary Fig. 8 panels?
cat("\n=== Figure scripts: lines that could put a p-value or 'BH' on a spatial panel ===\n")
fs <- list.files(file.path(base, "scripts"), pattern = "\\.R$", full.names = TRUE)
fs <- fs[sapply(fs, function(f) any(grepl("spatial|distance|Supp.?Fig.?8|supp8|figS8", readLines(f, warn = FALSE), ignore.case = TRUE)))]
fs <- setdiff(fs, file.path(base, "scripts", c("spatial_bh_recompute_v1.R", "spatial_bh_recompute_v2.R")))
for (f in fs) {
  l <- readLines(f, warn = FALSE)
  i <- grep("p *=|p\\.adjust|BH|wilcox|0\\.016|0\\.0078|annotate|label *=|subtitle|caption", l, ignore.case = TRUE)
  cat("\n----", basename(f), "(", length(i), "matching lines ) ----\n")
  for (j in i) cat(sprintf("%5d: %s\n", j, trimws(l[j])))
}
