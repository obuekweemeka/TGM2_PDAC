# ============================================================================
# optionA_and_fig3_checks_v1.R -- READ-ONLY. Save in: scripts/
# Option A: cells with no Axis 2 signal (score exactly 0, all four modules zero)
# are unclassifiable and excluded; cytostatic-prone = Axis 2 < 0 among the rest.
# Computes every knock-on number, then checks Figure 3 (LIANA) and Figure 1.
# Each section is wrapped so one failure does not stop the others.
# ============================================================================
suppressPackageStartupMessages(library(Matrix))
base <- normalizePath('.'); pr <- file.path(base, "data", "processed")
v4  <- readRDS(file.path(pr, "aim2_scoring", "malignant_cells_scored_v4.rds"))
mcm <- readRDS(file.path(pr, "aim2_scoring", "malignant_cell_metadata.rds"))
tr  <- unique(mcm[, c("pid", "new_treatment")]); v4$trt <- tr$new_treatment[match(v4$pid, tr$pid)]
v4$tied <- v4$axis2_death_competence == 0; v4$cyto <- grepl("cytostatic", v4$cell_state)
u <- v4$trt == "Untreated"; x <- v4[u & !v4$tied, ]
sec <- function(title, expr) { cat("\n===", title, "===\n"); tryCatch(expr, error = function(e) cat("SECTION FAILED:", conditionMessage(e), "\n")) }

sec("1. Headline under Option A (untreated)", {
  cat(sprintf("Informative cells: %d of %d (%d excluded, %.1f%%)\n", nrow(x), sum(u), sum(u & v4$tied), 100 * mean(v4$tied[u])))
  print(round(100 * prop.table(table(x$cell_state)), 2))
  cat(sprintf("Cytostatic-prone %.2f%% | death-competent %.2f%%\n", 100 * mean(x$cyto), 100 * mean(!x$cyto)))
  pp <- 100 * tapply(x$cyto, x$pid, mean)
  cat(sprintf("Per patient: median %.1f%%, mean %.1f%%, range %.1f-%.1f%%, all above 50%%: %s\n",
              median(pp), mean(pp), min(pp), max(pp), all(pp > 50)))
})

sec("2. Treatment groups", {
  for (g in c("Untreated", "CRT", "CRTl", "Other")) { s <- v4$trt == g
    cat(sprintf("%-9s cells %6d | no Axis 2 signal %.1f%% | cytostatic-prone among informative %.2f%% (n = %d)\n",
                g, sum(s), 100 * mean(v4$tied[s]), 100 * mean(v4$cyto[s & !v4$tied]), sum(s & !v4$tied))) }
  s <- v4$trt != "Untreated"
  cat(sprintf("All treated: %.2f%% (n = %d informative of %d)\n", 100 * mean(v4$cyto[s & !v4$tied]), sum(s & !v4$tied), sum(s)))
})

sec("3. Moffitt, cell level (untreated informative)", {
  b <- x$Moffitt_basal > x$Moffitt_classical
  cat(sprintf("Cytostatic-prone: basal-higher cells %.1f%% (n = %d) | classical-higher cells %.1f%% (n = %d)\n",
              100 * mean(x$cyto[b]), sum(b), 100 * mean(x$cyto[!b]), sum(!b)))
})

sec("4. Supp Fig 5 thresholds (untreated informative)", {
  ad <- x$axis2_death_competence; med <- median(ad); q33 <- unname(quantile(ad, 0.33))
  pm <- 100 * tapply(ad < med, x$pid, mean); p3 <- 100 * tapply(ad < q33, x$pid, mean)
  cat(sprintf("Median threshold %.3f -> per-patient %.1f-%.1f%% | 33rd percentile %.3f -> %.1f-%.1f%%\n",
              med, min(pm), max(pm), q33, min(p3), max(p3)))
  cat(sprintf("Check: %% below 0 = %.2f%% (should equal the headline)\n", 100 * mean(ad < 0)))
  sw <- seq(-10, 10, 0.5); cv <- sapply(sw, function(t) 100 * mean(ad < t))
  cat("Lowest threshold keeping > 50% cytostatic-prone:", min(sw[cv > 50]), "\n")
})

sec("5. Steele replication", {
  sc <- read.csv(file.path(pr, "aim2_validation", "steele_cell_scores.csv"))
  cat("Columns:", paste(names(sc), collapse = ", "), "\n")
  a2c <- grep("axis2", names(sc), value = TRUE)[1]; pc <- grep("pid|patient|sample|orig", names(sc), ignore.case = TRUE, value = TRUE)[1]
  cat("Using Axis 2 column:", a2c, "| patient column:", pc, "\n")
  a <- sc[[a2c]]; t0 <- a == 0
  cat(sprintf("Cells %d | no Axis 2 signal %.1f%% | %% < 0 all cells %.1f%% | %% <= 0 %.1f%% | %% < 0 among informative %.1f%%\n",
              length(a), 100 * mean(t0), 100 * mean(a < 0), 100 * mean(a <= 0), 100 * mean(a[!t0] < 0)))
  p_lt <- 100 * tapply(a < 0, sc[[pc]], mean); p_le <- 100 * tapply(a <= 0, sc[[pc]], mean); p_in <- 100 * tapply(a[!t0] < 0, sc[[pc]][!t0], mean)
  cat(sprintf("Per-patient median (n = %d): < 0 %.1f%% | <= 0 %.1f%% | informative %.1f%% (range %.1f-%.1f)  [manuscript 58.5%%]\n",
              length(p_lt), median(p_lt), median(p_le), median(p_in), min(p_in), max(p_in)))
  sp <- read.csv(file.path(pr, "aim2_validation", "steele_per_patient.csv")); cat("Stored per-patient file:\n"); print(head(sp, 20))
})

sec("6. LIANA subsample and sender groups", {
  for (f in c("liana_subsampled_metadata.rds", "liana_metadata_ordered.rds")) { mm <- readRDS(file.path(pr, "aim3_liana", f))
    cat(f, ":", class(mm)[1], nrow(mm), "rows\n") }
  meta <- readRDS(file.path(pr, "aim3_liana", "liana_subsampled_metadata.rds")); meta <- as.data.frame(meta)
  ov <- sapply(meta, function(c) sum(as.character(c) %in% v4$cell_id)); ovr <- sum(rownames(meta) %in% v4$cell_id)
  ids <- if (max(ov) >= ovr) as.character(meta[[names(which.max(ov))]]) else rownames(meta)
  lab <- names(meta)[sapply(meta, function(c) any(grepl("cytostatic", as.character(c), ignore.case = TRUE)))][1]
  cat("Sender label column:", lab, "\n"); print(head(sort(table(meta[[lab]]), decreasing = TRUE), 25))
  m <- match(ids, v4$cell_id); mal <- !is.na(m)
  cat("Subsampled cells matched to scored malignant cells:", sum(mal), "\n")
  cat("Treatment of matched malignant cells:\n"); print(table(v4$trt[m[mal]]))
  cy <- grepl("cytostatic", meta[[lab]], ignore.case = TRUE) & mal
  cat(sprintf("Cytostatic-labelled senders: %d | with no Axis 2 signal: %d (%.1f%%)\n", sum(cy), sum(v4$tied[m[cy]]), 100 * mean(v4$tied[m[cy]])))
})

sec("7. Fig 3b fold enrichments from saved LIANA results", {
  agg <- as.data.frame(readRDS(file.path(pr, "aim3_liana", "liana_results_aggregated.rds")))
  cat("Sources:", paste(head(unique(agg$source), 30), collapse = " | "), "\n")
  cs <- grepl("cytostatic", agg$source, ignore.case = TRUE); ds <- grepl("death", agg$source, ignore.case = TRUE)
  key <- paste(agg$ligand.complex, agg$receptor.complex, agg$target, sep = " | ")
  mc <- tapply(agg$natmi.edge_specificity[cs], key[cs], mean); md <- tapply(agg$natmi.edge_specificity[ds], key[ds], mean)
  cm <- intersect(names(mc), names(md)); cm <- cm[md[cm] > 0]; fe <- mc[cm] / md[cm]
  cat("Top 20 interactions by cytostatic / death-competent NATMI ratio:\n"); print(round(head(sort(fe, decreasing = TRUE), 20), 2))
  cat("\nAll TGM2 - ITGA4 targets:\n"); print(round(sort(fe[grepl("^TGM2 \\| ITGA4", names(fe))], decreasing = TRUE), 2))
  cat("\nVEGFA - FLT1:\n"); print(round(fe[grepl("^VEGFA \\| FLT1", names(fe))], 2))
})

sec("8. Figure 1: which DE coefficient?", {
  suppressPackageStartupMessages(library(DESeq2))
  dds <- readRDS(file.path(pr, "aim1_signatures", "dds_fitted.rds"))
  r0  <- as.data.frame(results(dds, name = "condition_resistant_vs_parental"))
  rsa <- readRDS(file.path(pr, "aim1_signatures", "res_shared_annot.rds")); mm <- match(rsa$ensembl_id, rownames(r0))
  cat(sprintf("Genes matched: %d | max |padj difference|: %.2e | max |p difference|: %.2e | LFC correlation: %.3f | mean |LFC| stored/unshrunk: %.2f\n",
              sum(!is.na(mm)), max(abs(rsa$padj - r0$padj[mm]), na.rm = TRUE), max(abs(rsa$pvalue - r0$pvalue[mm]), na.rm = TRUE),
              cor(rsa$log2FoldChange, r0$log2FoldChange[mm], use = "complete.obs"),
              mean(abs(rsa$log2FoldChange), na.rm = TRUE) / mean(abs(r0$log2FoldChange[mm]), na.rm = TRUE)))
})
