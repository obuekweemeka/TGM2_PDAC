# ============================================================================
# steele_tgm2_replication_v2.R -- Save in: scripts/
# Replicates the TGM2 comparison in the Steele atlas (GSE155698).
# Only TGM2 is extracted; the existing cell classifications in
# steele_cell_scores.csv are reused, so the atlas is NOT reprocessed.
#
# STAGE 1 (DRY_RUN = TRUE) : lists the files GEO holds and their sizes. No download.
# STAGE 2 (DRY_RUN = FALSE): downloads the PDAC tissue samples, unpacks, replicates.
# Downloads are skipped if the file already exists, so it is safe to re-run.
# ============================================================================
DRY_RUN <- TRUE          # <<< set to FALSE for the real run
suppressPackageStartupMessages(library(Matrix))
base <- normalizePath('.'); pr <- file.path(base, "data", "processed")
dest <- file.path(base, "data", "raw", "GSE155698"); dir.create(dest, recursive = TRUE, showWarnings = FALSE)
options(timeout = 7200)

sc <- read.csv(file.path(pr, "aim2_validation", "steele_cell_scores.csv"), stringsAsFactors = FALSE)
want <- sort(unique(sc$sample_id))
cat("Samples needing TGM2:", length(want), "\n  ", paste(want, collapse = ", "), "\n\n")

ftp <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE155nnn/GSE155698/suppl/"
cat("Listing", ftp, "\n")
pg <- tryCatch(paste(readLines(ftp, warn = FALSE), collapse = "\n"), error = function(e) {cat("Listing failed:", conditionMessage(e), "\n"); ""})
files <- unique(regmatches(pg, gregexpr('(?<=href=")[^"]+', pg, perl = TRUE))[[1]])
files <- files[!grepl("^/|\\?|^http", files)]
cat("Files advertised:", length(files), "\n"); for (f in files) cat("   ", f, "\n")

# keep only per-sample PDAC tissue archives
samp <- files[grepl("PDAC_TISSUE", files, ignore.case = TRUE) & grepl("\\.(tar|tar\\.gz|tgz|zip)$", files)]
if (!length(samp)) {
  cat("\nNo per-sample PDAC_TISSUE archives listed at the series level.\n",
      "They are then attached to the individual GSM records instead. In that case use:\n",
      "  if (!requireNamespace('GEOquery', quietly=TRUE)) BiocManager::install('GEOquery')\n",
      "  gse <- GEOquery::getGEO('GSE155698', GSEMatrix = FALSE)\n",
      "  sapply(GEOquery::GSMList(gse), function(g) GEOquery::Meta(g)$supplementary_file_1)\n",
      "and download those URLs into", dest, "\n")
} else cat("\nPDAC tissue archives:", length(samp), "\n")

if (DRY_RUN) {
  cat("\nDRY RUN. Nothing downloaded. Review the list above, then set DRY_RUN <- FALSE and re-run.\n")
} else {
  for (f in samp) { p <- file.path(dest, basename(f))
    if (file.exists(p) && file.size(p) > 1e5) { cat("have", basename(f), "\n"); next }
    cat("downloading", basename(f), "...\n")
    tryCatch(download.file(paste0(ftp, f), p, mode = "wb", quiet = TRUE), error = function(e) cat("  failed:", conditionMessage(e), "\n")) }
  for (p in list.files(dest, pattern = "\\.(tar|tar\\.gz|tgz)$", full.names = TRUE)) {
    od <- file.path(dest, sub("\\.(tar|tar\\.gz|tgz)$", "", basename(p)))
    if (!dir.exists(od)) { dir.create(od, showWarnings = FALSE); tryCatch(untar(p, exdir = od), error = function(e) cat("untar failed:", basename(p), "\n")) } }

  dirs <- unique(dirname(list.files(dest, pattern = "^(matrix\\.mtx|barcodes\\.tsv)", recursive = TRUE, full.names = TRUE)))
  cat("\n10x directories found:", length(dirs), "\n")
  pull <- function(d) {
    sid <- regmatches(d, regexpr("PDAC_TISSUE_[0-9]+[AB]?", d, ignore.case = TRUE))
    if (!length(sid)) return(NULL)
    f <- function(p) { h <- list.files(d, pattern = p, full.names = TRUE); if (length(h)) h[1] else NA }
    mf <- f("^matrix\\.mtx"); bf <- f("^barcodes\\.tsv"); gf <- f("^(features|genes)\\.tsv")
    if (any(is.na(c(mf, bf, gf)))) return(NULL)
    g <- read.delim(gf, header = FALSE, stringsAsFactors = FALSE)
    sym <- if (ncol(g) >= 2) g[[2]] else g[[1]]
    i <- which(toupper(sym) == "TGM2"); if (!length(i)) { cat("  TGM2 absent in", sid, "\n"); return(NULL) }
    mm <- Matrix::readMM(mf); bc <- readLines(bf); tot <- Matrix::colSums(mm)
    data.frame(sample_id = sid, bc = sub("-\\d+$", "", bc),
               tgm2 = log1p(as.numeric(mm[i[1], ]) / pmax(tot, 1) * 1e4), stringsAsFactors = FALSE)
  }
  got <- do.call(rbind, Filter(Negate(is.null), lapply(dirs, pull)))
  if (is.null(got)) stop("No usable matrices recovered.")
  cat("Cells with TGM2:", nrow(got), "across", length(unique(got$sample_id)), "samples\n")

  # match on SAMPLE + BARCODE (barcodes repeat across samples)
  sc$bc <- sub("-\\d+$", "", sub("^PDAC_TISSUE_[0-9]+[AB]?_", "", sc$cell_barcode))
  mg <- merge(sc, got, by = c("sample_id", "bc"))
  cat(sprintf("Matched %d of %d scored cells (%.1f%%)\n", nrow(mg), nrow(sc), 100*nrow(mg)/nrow(sc)))
  if (nrow(mg) < 0.5 * nrow(sc)) {
    cat("Low match rate. Barcode examples:\n  scored:", paste(head(sc$bc,3), collapse=" "), "\n  matrix:", paste(head(got$bc,3), collapse=" "), "\n")
  } else {
    mg <- mg[mg$axis2_death != 0, ]; mg$cyto <- mg$axis2_death < 0
    cat(sprintf("\nInformative matched cells: %d\n", nrow(mg)))
    cat(sprintf("TGM2 detected: %.1f%% cytostatic-prone vs %.1f%% death-competent\n",
                100*mean(mg$tgm2[mg$cyto] > 0), 100*mean(mg$tgm2[!mg$cyto] > 0)))
    pc <- tapply(mg$tgm2[mg$cyto], mg$sample_id[mg$cyto], mean)
    pd <- tapply(mg$tgm2[!mg$cyto], mg$sample_id[!mg$cyto], mean)
    k <- intersect(names(pc), names(pd)); d <- pc[k] - pd[k]
    w <- suppressWarnings(wilcox.test(pc[k], pd[k], paired = TRUE)); tt <- t.test(d)
    cat(sprintf("Samples: %d | higher in cytostatic-prone: %d\n", length(d), sum(d > 0)))
    cat(sprintf("Mean paired difference %.3f (95%% CI %.3f to %.3f), paired Wilcoxon p = %.4f\n",
                mean(d), tt$conf.int[1], tt$conf.int[2], w$p.value))
    cat("Hwang for comparison: 12 of 18 patients, 0.076 (0.024 to 0.128), p = 0.0047\n")
    write.csv(data.frame(sample = k, cytostatic = pc[k], death_competent = pd[k], difference = d),
              file.path(base, "results", "aim2", "steele_tgm2_replication.csv"), row.names = FALSE)
    cat("Saved results/aim2/steele_tgm2_replication.csv\n")
  }
}
