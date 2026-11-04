## spatial_bh_recompute_v1.R
## Recompute the spatial distance-to-epithelium statistics and BOTH Benjamini-Hochberg
## families directly from the per-spot source table, because no saved script produced
## the values currently in the manuscript (ledger "Step 3b" and "Item 0b" were interactive).
##
## Reads:  data/processed/aim3_spatial/spot_level_primary_PDAC.csv   (read-only)
## Writes: results/aim3/spatial/spatial_bh_recompute_v1_tests.csv
##         results/aim3/spatial/spatial_bh_recompute_v1_rho.csv
##
## Method (as stated in the previous Methods draft):
##   anchors   = spots with RCTD first_type == "Tumor Epithelial cells"
##   query     = non-epithelial spots; Euclidean distance to the nearest anchor (image pixels;
##               Spearman rho is invariant to the unit)
##   per sample: Spearman rho between distance and each programme score
##   inference : two-sided exact Wilcoxon signed-rank test of the per-sample rho values vs 0
## Base R only.

base   <- here::here()
infile <- file.path(base, "data", "processed", "aim3_spatial", "spot_level_primary_PDAC.csv")
outdir <- file.path(base, "results", "aim3", "spatial")
cat(R.version.string, "\n")

spot  <- read.csv(infile, stringsAsFactors = FALSE)
progs <- c(Treg        = "Treg",
           Checkpoint  = "Checkpoint_molecules",
           MyeloidSupp = "Immune_Suppression_by_Myeloid_cells",
           Effector    = "Effector_cells")
stopifnot(all(c("sample_id", "first_type", "spot_class", "imagecol", "imagerow", progs) %in% names(spot)))
TE <- "Tumor Epithelial cells"
stopifnot(TE %in% spot$first_type)

spot$is_te <- spot$first_type == TE
spot$is_sg <- spot$spot_class == "singlet"
by_s <- split(spot, spot$sample_id)

## ---------------------------------------------------------------- 1. counts
cnt <- do.call(rbind, lapply(by_s, function(d) data.frame(
  n_spots       = nrow(d),
  te_all        = sum(d$is_te),
  te_singlet    = sum(d$is_te & d$is_sg),
  nonte_singlet = sum(!d$is_te & d$is_sg))))
cat("\n=== 1. Per-sample counts ===\n"); print(cnt)

## ---------------------------------------------------------------- helpers
near_dist <- function(qx, qy, ax, ay) {
  d <- rep(Inf, length(qx))
  for (k in seq_along(ax)) d <- pmin(d, sqrt((qx - ax[k])^2 + (qy - ay[k])^2))
  d
}
rho_one <- function(d, anchor, query) {
  if (sum(anchor) < 1 || sum(query) < 3) return(setNames(rep(NA_real_, length(progs)), names(progs)))
  dist <- near_dist(d$imagecol[query], d$imagerow[query], d$imagecol[anchor], d$imagerow[anchor])
  sapply(progs, function(p) suppressWarnings(cor(dist, d[[p]][query], method = "spearman")))
}
wtest <- function(r) {
  r <- r[!is.na(r)]
  w <- suppressWarnings(wilcox.test(r, mu = 0, exact = TRUE))
  c(n = length(r), median_rho = median(r), V = unname(w$statistic), n_pos = sum(r > 0), p = w$p.value)
}
family <- function(label, R, keep) {
  R   <- R[keep, , drop = FALSE]
  out <- as.data.frame(t(apply(R, 2, wtest)))
  out <- cbind(family = label, programme = rownames(out), out)
  out$BH_4tests <- p.adjust(out$p, "BH")                       # all four programmes
  s3 <- out$programme != "Effector"
  out$BH_3tests <- NA_real_
  out$BH_3tests[s3] <- p.adjust(out$p[s3], "BH")               # suppressive programmes only
  out$samples <- paste(rownames(R), collapse = ";")
  rownames(out) <- NULL
  out
}
show <- function(f) {
  g <- f; g$samples <- NULL
  num <- sapply(g, is.numeric); g[num] <- lapply(g[num], signif, 4)
  print(g, row.names = FALSE)
  cat("  samples:", f$samples[1], "\n")
}

## ---------------------------------------------------------------- 2. per-sample rho
rho_all <- t(sapply(by_s, function(d) rho_one(d, d$is_te,          !d$is_te)))
rho_sgA <- t(sapply(by_s, function(d) rho_one(d, d$is_te,          !d$is_te & d$is_sg)))  # anchors: all TE spots
rho_sgB <- t(sapply(by_s, function(d) rho_one(d, d$is_te & d$is_sg, !d$is_te & d$is_sg)))  # anchors: singlet TE only
cat("\n=== 2a. Per-sample rho, all non-epithelial spots ===\n");                print(round(rho_all, 3))
cat("\n=== 2b. Per-sample rho, singlet spots, anchors = all epithelial spots ===\n"); print(round(rho_sgA, 3))
cat("\n=== 2c. Per-sample rho, singlet spots, anchors = singlet epithelial only ===\n"); print(round(rho_sgB, 3))

## ---------------------------------------------------------------- 3. test families
s     <- rownames(cnt)
is_T9 <- grepl("T9$", s)
cat("\nSample treated as T9:", s[is_T9], "\n")

keep_A9  <- !is_T9 & !is.na(rho_all[, 1])                            # primary as reported (T9 excluded)
keep_A10 <- !is.na(rho_all[, 1])                                     # sensitivity: T9 included
keep_SA  <- cnt$te_all     >= 5 & cnt$nonte_singlet >= 10 & !is.na(rho_sgA[, 1])
keep_SB  <- cnt$te_singlet >= 5 & cnt$nonte_singlet >= 10 & !is.na(rho_sgB[, 1])

fams <- list(
  family("A_primary_T9_excluded",            rho_all, keep_A9),
  family("A_primary_T9_included",            rho_all, keep_A10),
  family("B_singlet_anchorsAllTE_rule5_10",  rho_sgA, keep_SA),
  family("B_singlet_anchorsSingletTE_rule5_10", rho_sgB, keep_SB))
cat("\n=== 3. Wilcoxon signed-rank vs 0, with BH over 4 and over 3 programmes ===\n")
for (f in fams) { cat("\n--", f$family[1], "--\n"); show(f) }

## ---------------------------------------------------------------- 4. ledger targets
cat("\n=== 4. Ledger values these should reproduce ===\n")
cat("Primary: median 0.153, 0.150, 0.119, 0.185; V 44, 44, 45, 44; p 0.0078, 0.0078, 0.0039, 0.0078 (n = 9)\n")
cat("Singlet: median 0.20, 0.178, 0.052, 0.251; p 0.0195, 0.0137, 0.0371, 0.00977 (n = 10)\n")

## ---------------------------------------------------------------- 5. paired singlet vs all-spot
paired <- function(label, Rs, keep) {
  sh <- rownames(Rs)[keep & !is.na(rho_all[, 1])]
  t(sapply(names(progs), function(p) {
    w <- suppressWarnings(wilcox.test(Rs[sh, p], rho_all[sh, p], paired = TRUE, exact = TRUE))
    c(n_pairs = length(sh), p = signif(w$p.value, 3))
  }))
}
cat("\n=== 5. Paired singlet vs all-spot rho (ledger: p 0.70-0.85) ===\n")
cat("-- anchors = all TE --\n");     print(paired("A", rho_sgA, keep_SA))
cat("-- anchors = singlet TE --\n"); print(paired("B", rho_sgB, keep_SB))

## ---------------------------------------------------------------- 6. save
tests <- do.call(rbind, fams)
rho   <- data.frame(sample = s, cnt,
                    setNames(as.data.frame(rho_all), paste0("rho_all_", names(progs))),
                    setNames(as.data.frame(rho_sgA), paste0("rho_sgA_", names(progs))),
                    setNames(as.data.frame(rho_sgB), paste0("rho_sgB_", names(progs))))
write.csv(tests, file.path(outdir, "spatial_bh_recompute_v1_tests.csv"), row.names = FALSE)
write.csv(rho,   file.path(outdir, "spatial_bh_recompute_v1_rho.csv"),   row.names = FALSE)
cat("\nSaved: results/aim3/spatial/spatial_bh_recompute_v1_tests.csv and _rho.csv\n")
