## resistance_figures_v1.R
## Regenerates Supplementary Figures 3 and 4 (resistance transcriptomics, GSE269985).
## Written to close the last gap in the figure-to-script map: the original builders for
## these two figures were not retained.
##
## Supp Fig 3  a) volcano, resistant vs parental, Panc 02.03  b) Hallmark GSEA
## Supp Fig 4  a) count density  b) library sizes  c) PCA (both lines)
##             d) dispersion estimates  e) PCA within each line
##
## Part 1 checks every plotted value against docs/numbers_ledger_v15.md and stops if any
## check fails. Part 2 writes. Existing files are backed up to output/figures/_archive/.

suppressPackageStartupMessages({
  library(ggplot2); library(ggrepel); library(magick); library(DESeq2)
})
base <- normalizePath('.')
pr   <- file.path(base, "data", "processed", "aim1_signatures")
out  <- file.path(base, "output", "figures")
arc  <- file.path(out, "_archive")
dir.create(arc, recursive = TRUE, showWarnings = FALSE)
cat(R.version.string, "\n")

W <- 8.27; H <- 11.69; dpi <- 600
ml <- .22; mr <- .06; mt <- .06; mb <- .06; gx <- .22; gy <- .12
wf <- W - ml - mr; cw <- (wf - gx) / 2; xr <- ml + cw + gx
px  <- function(z) as.integer(round(z * dpi))
put <- function(cv, img, x, y, w) image_composite(cv, image_scale(img, paste0(px(w), "x")),
                                                  offset = paste0("+", px(x), "+", px(y)))
hin <- function(img, w) { i <- image_info(img); w * i$height / i$width }
thm <- function(fs) theme_classic(base_size = fs) +
  theme(plot.title = element_text(face = "bold", size = fs + 1),
        plot.subtitle = element_text(size = fs - 1.5, colour = "grey40"))
bk <- function(p) if (file.exists(p) && !file.exists(file.path(arc, basename(p))))
  file.copy(p, file.path(arc, basename(p)))
sv <- function(p, f, w, h) { bk(f); ggsave(f, p, dpi = dpi, width = w, height = h, units = "in")
  stopifnot(file.size(f) > 2e4) }
page <- function(cv, labs, name) {
  pg <- ggplot() + annotation_raster(as.raster(cv), xmin = 0, xmax = W, ymin = 0, ymax = H, interpolate = TRUE) +
    geom_text(data = labs, aes(x = x, y = H - y, label = l), hjust = 0, vjust = 1,
              fontface = "bold", size = 9 / .pt) +
    coord_fixed(xlim = c(0, W), ylim = c(0, H), expand = FALSE) + theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
  f <- file.path(out, name); bk(f); tmp <- sub("\\.pdf$", "_new.pdf", f)
  ggsave(tmp, pg, width = W, height = H, units = "in", device = "pdf")
  stopifnot(file.size(tmp) > 1e5); file.copy(tmp, f, overwrite = TRUE); file.remove(tmp)
  cat(sprintf("  %-14s written, %d KB\n", name, round(file.size(f) / 1024)))
}
chk <- function(what, got, want, tol) {
  ok <- length(got) == length(want) && all(abs(got - want) <= tol)
  cat(sprintf("  %-50s %s\n", what, if (ok) "OK" else paste("MISMATCH:", paste(signif(got, 5), collapse = ", "))))
  if (!ok) stop("Check failed (", what, "). Nothing written.")
}

## =========================================================== PART 1: load and check
cat("\n=== PART 1: values checked against the ledger ===\n")
res <- as.data.frame(readRDS(file.path(pr, "res_shared_annot.rds")))
gse <- as.data.frame(readRDS(file.path(pr, "gsea_hallmark_res.rds")))
cm  <- readRDS(file.path(pr, "count_matrix_filtered.rds"))
md  <- as.data.frame(readRDS(file.path(pr, "sample_metadata.rds")))
dds <- readRDS(file.path(pr, "dds_fitted.rds"))
stopifnot(all(c("log2FoldChange", "padj", "hgnc_symbol") %in% names(res)),
          all(c("pathway_clean", "NES", "padj", "direction") %in% names(gse)))

g <- function(s) res$log2FoldChange[match(s, res$hgnc_symbol)]
chk("Panc 02.03 contrast: CDKN1A log2FC (ledger +1.04)", g("CDKN1A"), 1.04, .01)
chk("module genes up in resistant (BCL2, PMAIP1, RB1, GADD45A)",
    as.numeric(g(c("BCL2", "PMAIP1", "RB1", "GADD45A")) > 0), rep(1, 4), 0)
chk("module genes down in resistant (BIRC5, GADD45B)",
    as.numeric(g(c("BIRC5", "GADD45B")) < 0), rep(1, 2), 0)
chk("genes tested / significant at padj < 0.05", c(nrow(res), sum(res$padj < 0.05, na.rm = TRUE)),
    c(16439, 8144), 0)
sig <- gse[!is.na(gse$padj) & gse$padj < 0.05, ]
chk("Hallmark pathways at padj < 0.05", nrow(sig), 11, 0)
chk("depleted hallmark pathways", sum(sig$NES < 0), 1, 0)
stopifnot(grepl("OXIDATIVE", sig$pathway_clean[sig$NES < 0]))
cat("  depleted pathway is", sig$pathway_clean[sig$NES < 0], "\n")
chk("samples / cell lines / conditions", c(nrow(md), length(unique(md$cell_line)),
    length(unique(md$condition))), c(12, 2, 2), 0)
cat("All checks passed.\n")

## =========================================================== PART 2: Supp Fig 3
cat("\n=== PART 2: figures ===\n")
res$class <- ifelse(is.na(res$padj) | res$padj >= 0.05, "NS",
             ifelse(res$log2FoldChange > 0, "Up in resistant", "Down in resistant"))
res$class <- factor(res$class, levels = c("Down in resistant", "NS", "Up in resistant"))
res$negl  <- -log10(pmax(res$pvalue, .Machine$double.xmin))
cap <- quantile(res$negl[is.finite(res$negl)], .9995, na.rm = TRUE)
res$negl <- pmin(res$negl, cap)
lab <- c("S100A2","MGAT3","LBH","FGF19","ALDH1A1","AQP3","FSCN1","VIM","NT5E","CSF1",
         "SERPINE1","CCN2","CD44","ALDH3A1","DHRS9","ANO1","TMPRSS4","B3GALT5","S100A4")
lb <- res[res$hgnc_symbol %in% lab & !is.na(res$padj) & res$padj < 0.05, ]
h3a <- 6.2
p3a <- ggplot(res, aes(log2FoldChange, negl, colour = class)) +
  geom_point(size = .45, alpha = .55, stroke = 0) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "grey45", linewidth = .3) +
  ggrepel::geom_text_repel(data = lb, aes(label = hgnc_symbol), size = 2.3, fontface = "italic",
                           seed = 1, max.overlaps = 30, min.segment.length = .2, show.legend = FALSE) +
  scale_colour_manual(values = c("Down in resistant" = "#2166AC", "NS" = "grey80",
                                 "Up in resistant" = "#B2182B"), name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
  labs(title = "Differential expression, MRTX1133-resistant versus parental Panc 02.03",
       subtitle = sprintf("DESeq2; %s genes tested, %s significant at adjusted p < 0.05",
                          format(nrow(res), big.mark = ","),
                          format(sum(res$padj < 0.05, na.rm = TRUE), big.mark = ",")),
       x = expression(log[2]*" fold change (resistant versus parental)"),
       y = expression(-log[10]*"(p value)")) +
  thm(9) + theme(legend.position = "top")
f3a <- file.path(out, "SuppFig3a_volcano.png"); sv(p3a, f3a, wf, h3a)

sig <- sig[order(sig$NES), ]; sig$pathway_clean <- factor(sig$pathway_clean, levels = sig$pathway_clean)
h3b <- H - mt - h3a - gy - mb
p3b <- ggplot(sig, aes(NES, pathway_clean, colour = direction, size = -log10(padj))) +
  geom_vline(xintercept = 0, colour = "grey45", linewidth = .3) + geom_point() +
  scale_colour_manual(values = c("Enriched in resistant" = "#B2182B",
                                 "Depleted in resistant" = "#2166AC"), name = NULL) +
  scale_size_continuous(range = c(1.5, 4.5), name = expression(-log[10]*"(adjusted p)")) +
  labs(title = "Hallmark pathways enriched or depleted in resistant cells",
       subtitle = "fgsea on the same contrast; pathways with adjusted p < 0.05",
       x = "Normalised enrichment score", y = NULL) +
  thm(9) + theme(legend.position = "right", legend.box = "vertical")
f3b <- file.path(out, "SuppFig3b_gsea.png"); sv(p3b, f3b, wf, h3b)
cv <- image_blank(px(W), px(H), color = "white")
cv <- put(cv, image_read(f3a), ml, mt, wf); cv <- put(cv, image_read(f3b), ml, mt + h3a + gy, wf)
page(cv, data.frame(l = c("a", "b"), x = c(.03, .03), y = c(mt, mt + h3a + gy)), "FigureS3.pdf")

## =========================================================== Supp Fig 4 (QC)
cpm <- t(t(as.matrix(cm)) / colSums(as.matrix(cm))) * 1e6
d <- do.call(rbind, lapply(colnames(cpm), function(s) data.frame(
  sample = s, cell_line = md$cell_line[match(s, md$sample_id)],
  condition = md$condition[match(s, md$sample_id)], v = log2(cpm[, s] + 0.5))))
p4a <- ggplot(d, aes(v, group = sample, colour = condition)) + geom_density(linewidth = .35) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey45", linewidth = .3) +
  facet_wrap(~cell_line) +
  scale_colour_manual(values = c(parental = "#2166AC", resistant = "#B2182B"), name = NULL) +
  labs(title = "Distribution of normalised counts", subtitle = "One line per sample",
       x = expression(log[2]*"(CPM + 0.5)"), y = "Density") + thm(9) + theme(legend.position = "top")

ls_ <- data.frame(sample = colnames(cm), million = colSums(as.matrix(cm)) / 1e6,
                  cell_line = md$cell_line[match(colnames(cm), md$sample_id)],
                  condition = md$condition[match(colnames(cm), md$sample_id)])
p4b <- ggplot(ls_, aes(sample, million, fill = condition)) + geom_col(width = .7) +
  facet_wrap(~cell_line, scales = "free_x") +
  scale_fill_manual(values = c(parental = "#2166AC", resistant = "#B2182B"), name = NULL) +
  labs(title = "Library size", subtitle = "Total assigned reads per sample",
       x = NULL, y = "Million reads") +
  thm(9) + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")

vs <- assay(vst(dds, blind = TRUE))
pca <- function(m, ids) {
  rv <- apply(m, 1, stats::var); v <- head(order(rv, decreasing = TRUE), 500)   # no extra dependency
  p <- prcomp(t(m[v, ])); pv <- round(100 * p$sdev^2 / sum(p$sdev^2), 1)
  data.frame(PC1 = p$x[, 1], PC2 = p$x[, 2], sample = ids,
             cell_line = md$cell_line[match(ids, md$sample_id)],
             condition = md$condition[match(ids, md$sample_id)], v1 = pv[1], v2 = pv[2])
}
pc <- pca(vs, colnames(vs))
p4c <- ggplot(pc, aes(PC1, PC2, colour = condition, shape = cell_line)) + geom_point(size = 2.4) +
  scale_colour_manual(values = c(parental = "#2166AC", resistant = "#B2182B"), name = NULL) +
  scale_shape_discrete(name = NULL) +
  labs(title = "Principal component analysis, both cell lines",
       subtitle = "Variance-stabilised counts, 500 most variable genes",
       x = sprintf("PC1 (%.1f%%)", pc$v1[1]), y = sprintf("PC2 (%.1f%%)", pc$v2[1])) +
  thm(9) + theme(legend.position = "top")

dp <- data.frame(mean = mcols(dds)$baseMean, gene = mcols(dds)$dispGeneEst,
                 fit = mcols(dds)$dispFit, fin = dispersions(dds))
dp <- dp[dp$mean > 0 & !is.na(dp$gene), ]
p4d <- ggplot(dp, aes(mean, gene)) + geom_point(size = .25, alpha = .3, colour = "grey30") +
  geom_point(aes(y = fin), size = .25, alpha = .3, colour = "#2166AC") +
  geom_line(aes(y = fit), colour = "#B2182B", linewidth = .6) +
  scale_x_log10() + scale_y_log10() +
  labs(title = "DESeq2 dispersion estimates",
       subtitle = "Grey: per-gene estimate. Red: fitted trend. Blue: final shrunken estimate.",
       x = "Mean of normalised counts", y = "Dispersion") + thm(9)

pl <- do.call(rbind, lapply(unique(md$cell_line), function(cl) {
  s <- md$sample_id[md$cell_line == cl]; z <- pca(vs[, s, drop = FALSE], s); z$panel <- cl; z }))
p4e <- ggplot(pl, aes(PC1, PC2, colour = condition)) + geom_point(size = 2.2) +
  ggrepel::geom_text_repel(aes(label = sample), size = 2.1, seed = 1, show.legend = FALSE) +
  facet_wrap(~panel, scales = "free") +
  scale_colour_manual(values = c(parental = "#2166AC", resistant = "#B2182B"), name = NULL) +
  labs(title = "Principal component analysis within each cell line",
       subtitle = "Replicates separate by condition in both lines", x = "PC1", y = "PC2") +
  thm(9) + theme(legend.position = "top")

hr <- (H - mt - mb - 2 * gy) / 3
for (z in list(list(p4a, "SuppFig4a_density.png", cw, hr), list(p4b, "SuppFig4b_libsize.png", cw, hr),
               list(p4c, "SuppFig4c_pca.png", cw, hr), list(p4d, "SuppFig4d_dispersion.png", cw, hr),
               list(p4e, "SuppFig4e_pca_within.png", wf, hr)))
  sv(z[[1]], file.path(out, z[[2]]), z[[3]], z[[4]])
y2 <- mt + hr + gy; y3 <- y2 + hr + gy
cv <- image_blank(px(W), px(H), color = "white")
cv <- put(cv, image_read(file.path(out, "SuppFig4a_density.png")), ml, mt, cw)
cv <- put(cv, image_read(file.path(out, "SuppFig4b_libsize.png")), xr, mt, cw)
cv <- put(cv, image_read(file.path(out, "SuppFig4c_pca.png")), ml, y2, cw)
cv <- put(cv, image_read(file.path(out, "SuppFig4d_dispersion.png")), xr, y2, cw)
cv <- put(cv, image_read(file.path(out, "SuppFig4e_pca_within.png")), ml, y3, wf)
page(cv, data.frame(l = letters[1:5], x = c(.03, xr - .19, .03, xr - .19, .03),
                    y = c(mt, mt, y2, y2, y3)), "FigureS4.pdf")

cat("\nDone. FigureS3.pdf and FigureS4.pdf are in output/figures/.\n")
cat("Compare against the versions in the manuscript before replacing them.\n")
