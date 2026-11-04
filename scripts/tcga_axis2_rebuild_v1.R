## tcga_axis2_rebuild_v1.R
## Part A: zero-score audit (reviewer point 1).
## Part B: reconstruct the stored TCGA bulk Axis 2 ssGSEA score from raw GDC files (reviewer point 2).
## Read-only on all inputs; writes results/aim4/tcga_axis2_rebuild_candidates_v1.csv only.

base <- normalizePath('.')
pr   <- file.path(base, "data", "processed")
cat(R.version.string, "\n")

## =================================================================== Part A
cat("\n=================== A. Zero-score audit (untreated cells) ===================\n")
v4   <- readRDS(file.path(pr, "aim2_scoring", "malignant_cells_scored_v4.rds"))
mcm  <- readRDS(file.path(pr, "aim2_scoring", "malignant_cell_metadata.rds"))
tr   <- unique(mcm[, c("pid", "new_treatment")]); v4$trt <- tr$new_treatment[match(v4$pid, tr$pid)]   # as aim2_checks_v1.R line 9
comp <- c("axis2_proapoptotic", "axis2_executors", "axis2_anti_apoptotic", "axis2_cytostasis")
stopifnot(all(c("trt", "axis2_death_competence", comp) %in% names(v4)))
stopifnot(sum(v4$trt == "Untreated", na.rm = TRUE) == 52999)   # ledger: 52,999 untreated cells
u  <- v4[v4$trt == "Untreated", ]
s  <- u$axis2_death_competence
num <- u[, comp[1]] + u[, comp[2]]; den <- u[, comp[3]] + u[, comp[4]]
z  <- s == 0
all4 <- rowSums(u[, comp] == 0) == 4
cat("Untreated cells:", nrow(u), "\n")
cat("Score exactly 0:", sum(z), "\n")
cat("  of which all four module scores = 0:   ", sum(z & all4), "\n")
cat("  of which equal NON-zero composite sums:", sum(z & !all4 & num == den & num > 0), "\n")
cat("  of which anything else:                ", sum(z & !all4 & !(num == den & num > 0)), "\n")
cat("Cells with all four = 0 but score not 0 (should be 0):", sum(all4 & !z), "\n")
cat("Cells with 0 < |score| < 1e-6 (near-ties counted as informative):", sum(!z & abs(s) < 1e-6), "\n")
cat("Smallest non-zero |score|:", signif(min(abs(s[!z])), 4), "\n")

## =================================================================== Part B
cat("\n=================== B. TCGA bulk Axis 2 reconstruction ===================\n")
tc <- readRDS(file.path(pr, "aim4_bulk", "tcga_paad_analysis_df.rds"))
cat("Stored TCGA table:", nrow(tc), "rows x", ncol(tc), "cols\n"); print(names(tc))
cat("\nFirst 3 rows:\n"); print(head(tc, 3))
idcol <- names(tc)[sapply(tc, function(v) is.character(v) && mean(grepl("^TCGA-", v)) > 0.9)][1]
ids   <- if (!is.na(idcol)) tc[[idcol]] else rownames(tc)
cat("\nBarcode source:", ifelse(is.na(idcol), "rownames", idcol), "| example:", ids[1], "| nchar:", paste(unique(nchar(ids)), collapse = ","), "\n")
stopifnot(mean(grepl("^TCGA-", ids)) > 0.9, "axis2_death" %in% names(tc))

## ---- 1. map GDC file UUIDs (folder names) to sample barcodes via the public GDC API
gdir  <- file.path(base, "data", "raw", "TCGA_PAAD", "TCGA-PAAD", "Transcriptome_Profiling", "Gene_Expression_Quantification")
files <- list.files(gdir, pattern = "augmented_star_gene_counts\\.tsv$", recursive = TRUE, full.names = TRUE)
fid   <- basename(dirname(files))
cat("\nGDC expression files:", length(files), "\n")
for (p in c("httr", "jsonlite")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
body <- list(filters = list(op = "in", content = list(field = "file_id", value = as.list(fid))),
             fields = "file_id,cases.samples.submitter_id,cases.samples.sample_type", format = "JSON", size = 1000)
resp <- httr::POST("https://api.gdc.cancer.gov/files", body = body, encode = "json")
stopifnot(httr::status_code(resp) == 200)
hits <- jsonlite::fromJSON(httr::content(resp, "text", encoding = "UTF-8"), simplifyVector = FALSE)$data$hits
map  <- do.call(rbind, lapply(hits, function(h) data.frame(
  file_id = h$file_id,
  barcode = h$cases[[1]]$samples[[1]]$submitter_id,
  sample_type = h$cases[[1]]$samples[[1]]$sample_type)))
cat("Mapped:", nrow(map), "of", length(fid), "\n"); print(table(map$sample_type))

## ---- 2. expression matrices (TPM, FPKM, counts), gene symbols, duplicates resolved by highest mean
rd <- function(f) { x <- read.delim(f, comment.char = "#", stringsAsFactors = FALSE); x[!grepl("^N_", x$gene_id), ] }
x1 <- rd(files[1]); g <- x1$gene_name
mk <- function(col) {
  M <- sapply(files, function(f) rd(f)[[col]]); rownames(M) <- g; colnames(M) <- map$barcode[match(fid, map$file_id)]
  M <- M[order(-rowMeans(M)), ]; M[!duplicated(rownames(M)), ]
}
cat("Reading", length(files), "files x 3 value types (about a minute)...\n")
EX <- list(tpm = mk("tpm_unstranded"), fpkm = mk("fpkm_unstranded"), counts = mk("unstranded"))

## ---- 3. ssGSEA on all axis modules, then candidate combinations
sig <- readRDS(file.path(pr, "aim1_signatures", "four_axis_signatures_v1.rds"))
gs  <- c(sig, list(axis2_APOP_pooled  = unique(c(sig$axis2_proapoptotic, sig$axis2_executors)),
                   axis2_ARREST_pooled = unique(c(sig$axis2_anti_apoptotic, sig$axis2_cytostasis))))
cat("Gene sets scored:", paste(names(gs), collapse = ", "), "\n")

key  <- function(b) substr(b, 1, 15)
tcm  <- data.frame(key = key(ids), stored = tc$axis2_death)
res  <- list(); keep_scores <- NULL
for (et in names(EX)) {
  M  <- log2(EX[[et]] + 1)
  M  <- M[rowSums(M > 0) > 0, ]
  for (nrm in c(TRUE, FALSE)) {
    S  <- GSVA::gsva(GSVA::ssgseaParam(M, gs, normalize = nrm), verbose = FALSE)
    cand <- data.frame(key = key(colnames(S)),
      diff4      = S["axis2_proapoptotic", ] + S["axis2_executors", ] - S["axis2_anti_apoptotic", ] - S["axis2_cytostasis", ],
      half_diff4 = (S["axis2_proapoptotic", ] + S["axis2_executors", ] - S["axis2_anti_apoptotic", ] - S["axis2_cytostasis", ]) / 2,
      pooled     = S["axis2_APOP_pooled", ] - S["axis2_ARREST_pooled", ])
    cand <- cand[!duplicated(cand$key), ]
    mm <- merge(tcm, cand, by = "key")
    for (cn in c("diff4", "half_diff4", "pooled")) {
      fit <- lm(mm$stored ~ mm[[cn]])
      grp_s <- mm$stored < median(mm$stored); grp_c <- mm[[cn]] < median(mm[[cn]])
      res[[length(res) + 1]] <- data.frame(expr = et, normalize = nrm, candidate = cn, n = nrow(mm),
        pearson = cor(mm$stored, mm[[cn]]), spearman = cor(mm$stored, mm[[cn]], method = "spearman"),
        slope = unname(coef(fit)[2]), intercept = unname(coef(fit)[1]), max_resid = max(abs(resid(fit))),
        max_abs_diff = max(abs(mm$stored - mm[[cn]])), split_agree = sum(grp_s == grp_c))
    }
    keep_scores <- rbind(keep_scores, data.frame(expr = et, normalize = nrm, cand))
  }
}
res <- do.call(rbind, res); res <- res[order(-res$pearson), ]
cat("\n=== Candidate reconstructions vs stored axis2_death (sorted by Pearson) ===\n")
num <- sapply(res, is.numeric); r2 <- res; r2[num] <- lapply(r2[num], signif, 5); print(r2, row.names = FALSE)
best <- res[1, ]
cat(sprintf("\nBEST: %s, normalize = %s, %s | Pearson %.6f | max |diff| %.3g | median-split agreement %d of %d\n",
            best$expr, best$normalize, best$candidate, best$pearson, best$max_abs_diff, best$split_agree, best$n))
cat(if (best$max_abs_diff < 1e-6) "EXACT reproduction: state this formula in Methods.\n" else
    if (best$pearson > 0.999 && best$split_agree == best$n) "Reproduced up to scale/rounding with identical groups.\n" else
    "NOT reproduced exactly: regenerate the bulk score from a stated formula and rerun survival/purity.\n")

cat("\n=== How the 146 analysed samples relate to the files ===\n")
cat("Primary tumour files (sample type code 01):", sum(substr(map$barcode, 14, 15) == "01"), "\n")
cat("Stored samples found among files:", sum(key(ids) %in% key(map$barcode)), "of", length(ids), "\n")
cat("Stored sample type codes:", paste(names(table(substr(ids, 14, 15))), table(substr(ids, 14, 15)), sep = "=", collapse = ", "), "\n")

write.csv(keep_scores, file.path(base, "results", "aim4", "tcga_axis2_rebuild_candidates_v1.csv"), row.names = FALSE)
cat("\nSaved results/aim4/tcga_axis2_rebuild_candidates_v1.csv\n")
