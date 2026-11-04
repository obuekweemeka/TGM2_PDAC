# ============================================================================
# liana_effect_check_v1.R -- READ-ONLY. Save in: scripts/
# 1. Reproduce the Fig 3b "fold enrichment" from the saved LIANA outputs
# 2. Rank positions of TGM2-ITGA4 for cytostatic vs death-competent senders
# 3. TGM2 / ITGA4 expression in cytostatic vs death-competent malignant cells
#    (untreated, informative cells), with a patient-level paired test
# 4. LIANA subsample ID formats, to fix the cell matching
# ============================================================================
suppressPackageStartupMessages(library(Matrix))
base <- normalizePath('.'); pr <- file.path(base, "data", "processed"); r3 <- file.path(base, "results", "aim3")
sec <- function(title, expr) { cat("\n===", title, "===\n"); tryCatch(expr, error = function(e) cat("SECTION FAILED:", conditionMessage(e), "\n")) }

sec("1. Reproduce the fold enrichment", {
  cat("Files in results/aim3:\n"); print(list.files(r3))
  ms <- read.csv(file.path(r3, "liana_malignant_sender_interactions.csv"), stringsAsFactors = FALSE)
  cat("Columns:", paste(names(ms), collapse = ", "), "\n")
  ms$key <- paste(ms$ligand.complex, ms$receptor.complex, ms$target, sep = " | ")
  cy <- ms$source %in% c("RAS_independent_cytostatic", "RAS_dependent_cytostatic")
  dc <- ms$source %in% c("RAS_independent_death_competent", "RAS_dependent_death_competent")
  rc <- tapply(ms$aggregate_rank[cy], ms$key[cy], min); rd <- tapply(ms$aggregate_rank[dc], ms$key[dc], min)
  k <- intersect(names(rc), names(rd)); f <- rd[k] / rc[k]
  t <- grep("^TGM2 \\| ITGA4", k, value = TRUE)
  out <- data.frame(interaction = t, rank_cytostatic = signif(rc[t], 3), rank_death_comp = signif(rd[t], 3), ratio = round(f[t], 2))
  print(out[order(-out$ratio), ], row.names = FALSE)
  ff <- list.files(r3, pattern = "fold", full.names = TRUE)
  if (length(ff)) { fd <- read.csv(ff[1], stringsAsFactors = FALSE); cat("\nStored file:", basename(ff[1]), "| columns:", paste(names(fd), collapse = ", "), "\n")
    print(head(fd[grepl("TGM2", fd[[1]]) | grepl("TGM2", apply(fd, 1, paste, collapse = " ")), ], 10)) }
})

sec("2. Rank positions among each sender group's interactions", {
  ms <- read.csv(file.path(r3, "liana_malignant_sender_interactions.csv"), stringsAsFactors = FALSE)
  for (grp in list(c("RAS_independent_cytostatic", "RAS_dependent_cytostatic"), c("RAS_independent_death_competent", "RAS_dependent_death_competent"))) {
    s <- ms[ms$source %in% grp, ]; s <- s[order(s$aggregate_rank), ]; s$pos <- seq_len(nrow(s))
    h <- s[s$ligand.complex == "TGM2" & s$receptor.complex == "ITGA4", ]
    cat(sprintf("\n%s: %d interactions | best TGM2-ITGA4 position %s | median position %s\n",
                paste(grp, collapse = " + "), nrow(s), if (nrow(h)) min(h$pos) else "none", if (nrow(h)) median(h$pos) else "none"))
  }
})

sec("3. TGM2 / ITGA4 expression by state (untreated, informative cells)", {
  v4  <- readRDS(file.path(pr, "aim2_scoring", "malignant_cells_scored_v4.rds"))
  mcm <- readRDS(file.path(pr, "aim2_scoring", "malignant_cell_metadata.rds"))
  tr  <- unique(mcm[, c("pid", "new_treatment")]); v4$trt <- tr$new_treatment[match(v4$pid, tr$pid)]
  x <- v4[v4$trt == "Untreated" & v4$axis2_death_competence != 0, ]; x$cyto <- grepl("cytostatic", x$cell_state)
  m <- readRDS(file.path(pr, "aim2_scoring", "malignant_count_matrix.rds"))
  if (!("TGM2" %in% rownames(m)) && ("TGM2" %in% colnames(m))) m <- t(m)
  idx <- match(x$cell_id, colnames(m)); ok <- !is.na(idx); x <- x[ok, ]; idx <- idx[ok]
  lib <- Matrix::colSums(m[, idx])
  for (g in c("TGM2", "ITGA4")) {
    if (!(g %in% rownames(m))) { cat(g, "not in matrix\n"); next }
    cnt <- as.numeric(m[g, idx]); ln <- log1p(cnt / lib * 1e4)
    cat(sprintf("\n%s: detected in %.1f%% of cytostatic vs %.1f%% of death-competent cells | mean log-normalised %.3f vs %.3f\n",
                g, 100 * mean(cnt[x$cyto] > 0), 100 * mean(cnt[!x$cyto] > 0), mean(ln[x$cyto]), mean(ln[!x$cyto])))
    pc <- tapply(ln[x$cyto], x$pid[x$cyto], mean); pd <- tapply(ln[!x$cyto], x$pid[!x$cyto], mean); p <- intersect(names(pc), names(pd))
    w <- suppressWarnings(wilcox.test(pc[p], pd[p], paired = TRUE))
    cat(sprintf("   patient-level paired comparison (n = %d): cytostatic higher in %d | median difference %.3f | Wilcoxon signed-rank p = %.4f\n",
                length(p), sum(pc[p] > pd[p]), median(pc[p] - pd[p]), w$p.value))
  }
  rm(m); invisible(gc())
})

sec("4. LIANA subsample ID formats", {
  meta <- as.data.frame(readRDS(file.path(pr, "aim3_liana", "liana_subsampled_metadata.rds")))
  v4  <- readRDS(file.path(pr, "aim2_scoring", "malignant_cells_scored_v4.rds"))
  cat("Scored cell_id examples:", paste(head(v4$cell_id, 3), collapse = " | "), "\n")
  cat("Subsample rownames:", paste(head(rownames(meta), 3), collapse = " | "), "\n")
  for (cn in names(meta)) cat(sprintf("  %-22s %s\n", cn, paste(head(as.character(meta[[cn]]), 3), collapse = " | ")))
})
