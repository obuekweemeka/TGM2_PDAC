# ============================================================================
# moffitt_ci_and_steele_tgm2_v1.R -- Save in: scripts/
# PART A (runs now): confidence interval and leave-one-out sensitivity for the
#        patient-level Moffitt basal-like association.
# PART B (runs once GSE155698 processed matrices are downloaded): replicates the
#        TGM2 comparison in the Steele atlas. Only TGM2 is extracted; the
#        existing cell classifications in steele_cell_scores.csv are reused.
#        Set STEELE_DIR to the folder holding the per-sample 10x directories.
# ============================================================================
suppressPackageStartupMessages({library(Matrix)})
base <- normalizePath('.'); pr <- file.path(base, "data", "processed")
STEELE_DIR <- file.path(base, "data", "raw", "GSE155698")   # <- edit if downloaded elsewhere

cat("=== PART A: Moffitt association, uncertainty and influence ===\n")
v4  <- readRDS(file.path(pr, "aim2_scoring", "malignant_cells_scored_v4.rds"))
mcm <- readRDS(file.path(pr, "aim2_scoring", "malignant_cell_metadata.rds"))
tr  <- unique(mcm[, c("pid","new_treatment")]); v4$trt <- tr$new_treatment[match(v4$pid, tr$pid)]
u <- v4[v4$trt == "Untreated", ]; x <- u[u$axis2_death_competence != 0, ]
x$cyto <- grepl("cytostatic", x$cell_state); pids <- sort(unique(x$pid))
pm <- data.frame(pid = pids,
                 cyto = 100 * as.vector(tapply(x$cyto, x$pid, mean)[pids]),
                 basal = as.vector(tapply(x$Moffitt_basal, x$pid, mean)[pids]),
                 pct_uninf = 100 * as.vector(tapply(u$axis2_death_competence == 0, u$pid, mean)[pids]))
ct <- suppressWarnings(cor.test(pm$cyto, pm$basal, method = "spearman"))
cat(sprintf("Spearman rho = %.3f, exact p = %.4f, n = %d\n", ct$estimate, ct$p.value, nrow(pm)))
set.seed(1); bs <- replicate(10000, { i <- sample(nrow(pm), replace = TRUE)
  if (length(unique(pm$cyto[i])) < 3) NA else suppressWarnings(cor(pm$cyto[i], pm$basal[i], method = "spearman")) })
bs <- bs[is.finite(bs)]
cat(sprintf("Bootstrap 95%% CI (%d resamples, percentile): %.3f to %.3f\n", length(bs), quantile(bs,.025), quantile(bs,.975)))
z <- atanh(ct$estimate); se <- 1.06/sqrt(nrow(pm) - 3)
cat(sprintf("Fisher z 95%% CI (Bonett-Wright SE): %.3f to %.3f\n", tanh(z - 1.96*se), tanh(z + 1.96*se)))
loo <- sapply(seq_len(nrow(pm)), function(i) suppressWarnings(cor(pm$cyto[-i], pm$basal[-i], method = "spearman")))
cat(sprintf("Leave-one-out rho: range %.3f to %.3f | most influential patient: %s (rho without = %.3f)\n",
            min(loo), max(loo), pm$pid[which.max(abs(loo - ct$estimate))], loo[which.max(abs(loo - ct$estimate))]))
cat(sprintf("Correlation with proportion of uninformative cells: rho = %.3f (p = %.3f)\n",
            suppressWarnings(cor(pm$cyto, pm$pct_uninf, method="spearman")),
            suppressWarnings(cor.test(pm$cyto, pm$pct_uninf, method="spearman"))$p.value))
write.csv(pm, file.path(base, "results", "aim2", "moffitt_patient_level.csv"), row.names = FALSE)

cat("\n=== PART B: Steele TGM2 replication ===\n")
sc <- read.csv(file.path(pr, "aim2_validation", "steele_cell_scores.csv"), stringsAsFactors = FALSE)
cat("Scored Steele cells:", nrow(sc), "| samples:", length(unique(sc$sample_id)), "\n")
cat("Barcode examples:", paste(head(sc$cell_barcode, 2), collapse = " | "), "\n")
if (!dir.exists(STEELE_DIR)) {
  cat("\nSTEELE_DIR not found:", STEELE_DIR, "\n")
  cat("To run this part, download the processed per-sample supplementary files for GSE155698\n",
      "from GEO, unpack them into that folder (one directory per sample containing\n",
      "matrix.mtx(.gz), features/genes.tsv(.gz) and barcodes.tsv(.gz)), then re-run.\n")
} else {
  dirs <- unique(dirname(list.files(STEELE_DIR, pattern = "^(matrix\\.mtx|barcodes\\.tsv)", recursive = TRUE, full.names = TRUE)))
  cat("10x directories found:", length(dirs), "\n")
  get_tgm2 <- function(d) {
    f <- function(p) { h <- list.files(d, pattern = p, full.names = TRUE); if (length(h)) h[1] else NA }
    mf <- f("^matrix\\.mtx"); bf <- f("^barcodes\\.tsv"); gf <- f("^(features|genes)\\.tsv")
    if (any(is.na(c(mf, bf, gf)))) return(NULL)
    g <- read.delim(gf, header = FALSE, stringsAsFactors = FALSE)
    sym <- if (ncol(g) >= 2) g[[2]] else g[[1]]
    i <- which(toupper(sym) == "TGM2"); if (!length(i)) return(NULL)
    mm <- Matrix::readMM(mf); bc <- readLines(bf)
    tot <- Matrix::colSums(mm)
    data.frame(dir = basename(d), cell_barcode = bc,
               tgm2 = log1p(as.numeric(mm[i[1], ]) / pmax(tot, 1) * 1e4), stringsAsFactors = FALSE)
  }
  got <- do.call(rbind, Filter(Negate(is.null), lapply(dirs, get_tgm2)))
  if (is.null(got)) cat("No usable matrices found.\n") else {
    cat("Cells with TGM2 recovered:", nrow(got), "\n")
    key <- function(v) sub("-\\d+$", "", sub("^.*_", "", v))
    sc$k <- key(sc$cell_barcode); got$k <- key(got$cell_barcode)
    mg <- merge(sc, got[, c("k","tgm2")], by = "k")
    cat("Matched to scored cells:", nrow(mg), "of", nrow(sc), "\n")
    if (nrow(mg) > 1000) {
      mg <- mg[mg$axis2_death != 0, ]; mg$cyto <- mg$axis2_death < 0
      cat(sprintf("Informative matched cells: %d | TGM2 detected %.1f%% cytostatic-prone vs %.1f%% other\n",
                  nrow(mg), 100*mean(mg$tgm2[mg$cyto] > 0), 100*mean(mg$tgm2[!mg$cyto] > 0)))
      pc <- tapply(mg$tgm2[mg$cyto], mg$sample_id[mg$cyto], mean)
      pd <- tapply(mg$tgm2[!mg$cyto], mg$sample_id[!mg$cyto], mean)
      p <- intersect(names(pc), names(pd)); d <- pc[p] - pd[p]
      w <- suppressWarnings(wilcox.test(pc[p], pd[p], paired = TRUE)); tt <- t.test(d)
      cat(sprintf("Samples: %d | higher in cytostatic-prone: %d | mean paired difference %.3f (95%% CI %.3f to %.3f) | p = %.4f\n",
                  length(d), sum(d > 0), mean(d), tt$conf.int[1], tt$conf.int[2], w$p.value))
      cat("Hwang for comparison: 12 of 18 patients, mean difference 0.076 (0.024 to 0.128), p = 0.0047\n")
      write.csv(data.frame(sample = p, cyto = pc[p], other = pd[p], diff = d),
                file.path(base, "results", "aim2", "steele_tgm2_replication.csv"), row.names = FALSE)
    } else cat("Too few barcodes matched; inspect the barcode formats above before interpreting.\n")
  }
}
