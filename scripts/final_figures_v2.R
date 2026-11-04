# ============================================================================
# final_figures_v2.R -- final figure set (v2: Supp Fig 3a points enlarged; safe to re-run) for the rewritten manuscript. Save in: scripts/
#   Figure 1   unchanged (Figure1.pdf)
#   Figure 2   state UMAP + 4 module UMAPs, lettered a-e        (was Figure2a.pdf)
#   Figure 3   prevalence, correlations, Moffitt, Steele, a-d   (was Figure2b.pdf)
#   Figure 4   TGM2 candidate mechanism                         (was new Figure3.pdf)
#   Figure 5   tested and not supported: spatial test + TCGA    (new)
#   Supp 1, 2  unchanged; Supp 3 rebuilt; Supp 4 rebuilt; Supp 5 renamed
# Panels are PNG; multipanels are PDF only. Every replaced file goes to _archive/.
# Writes docs/figure_map.md listing every figure and its panel files.
# ============================================================================
suppressPackageStartupMessages({library(ggplot2); library(magick); library(Matrix)})
base <- normalizePath('.'); fig <- file.path(base, "figures"); pr <- file.path(base, "data", "processed")
mp <- file.path(fig, "multipanel"); ar <- file.path(fig, "_archive")
W <- 8.27; H <- 11.69; dpi <- 600; ml <- 0.22; mr <- 0.06; mt <- 0.06; mb <- 0.06; gx <- 0.22; gy <- 0.12; fs <- 6.5
cw <- (W - ml - mr - gx) / 2; wf <- W - ml - mr; xr <- ml + cw + gx
arch_copy <- function(p) { a <- file.path(ar, paste0(tools::file_path_sans_ext(basename(p)), "_superseded.", tools::file_ext(p)))
  if (file.exists(p) && !file.exists(a)) file.copy(p, a) }
arch_move <- function(p) { if (file.exists(p)) { arch_copy(p); file.remove(p); cat("Archived", basename(p), "\n") } }
savep <- function(p, path, w, h) { arch_copy(path); ggsave(path, p, dpi = dpi, width = w, height = h, units = "in"); cat("Saved", basename(path), round(file.size(path) / 1024), "KB\n") }
px  <- function(x) as.integer(round(x * dpi))
put <- function(cv, img, x, y, w) image_composite(cv, image_scale(img, paste0(px(w), "x")), offset = paste0("+", px(x), "+", px(y)))
aspect <- function(f) { i <- image_info(image_read(f)); i$height / i$width }
blank <- function() image_blank(px(W), px(H), color = "white")
page_pdf <- function(cv, labs, name) {
  p <- ggplot() + annotation_raster(as.raster(cv), xmin = 0, xmax = W, ymin = 0, ymax = H, interpolate = TRUE) +
    geom_text(data = labs, aes(x = x, y = H - y, label = l), hjust = 0, vjust = 1, fontface = "bold", size = 9 / ggplot2::.pt) +
    coord_fixed(xlim = c(0, W), ylim = c(0, H), expand = FALSE) + theme_void() + theme(plot.margin = margin(0, 0, 0, 0))
  out <- file.path(mp, name); tmp <- sub("\\.pdf$", "_new.pdf", out); arch_copy(out)
  ggsave(tmp, p, width = W, height = H, units = "in", device = "pdf")
  sz <- round(file.size(tmp) / 1024); stopifnot(sz > 100); file.copy(tmp, out, overwrite = TRUE); file.remove(tmp); cat("Saved", name, sz, "KB\n")
}
thm <- function(fs) theme_classic(base_size = fs) +
  theme(plot.title = element_text(face = "bold", size = fs + 1), plot.subtitle = element_text(size = fs - 1, colour = "grey40"))
L <- function(l, x, y) data.frame(l = l, x = x, y = y)

# ================= Figure 4 <- current TGM2 page (once only) =================
marker <- file.path(mp, ".renumbered_v1")
if (!file.exists(marker)) {
  arch_move(file.path(mp, "Figure4.pdf"))                                   # old KM / purity / DepMap page
  stopifnot(file.exists(file.path(mp, "Figure3.pdf")))
  file.copy(file.path(mp, "Figure3.pdf"), file.path(mp, "Figure4.pdf")); cat("TGM2 page is now Figure4.pdf\n")
  arch_copy(file.path(mp, "Figure3.pdf")); file.remove(file.path(mp, "Figure3.pdf"))
  writeLines(format(Sys.time()), marker)
} else cat("Renumbering already done; skipping\n")

# ================= Figure 2: state + module UMAPs, a-e =================
st <- file.path(fig, "aim2", "Fig2A_UMAP_state.png"); mods <- file.path(fig, "aim2", paste0("Fig2A_UMAP_axis", 1:4, ".png"))
stopifnot(all(file.exists(c(st, mods))))
hm <- cw * max(sapply(mods, aspect)); hs <- wf * aspect(st)
y1 <- mt + hs + gy; y2 <- y1 + hm + gy
cv <- blank(); cv <- put(cv, image_read(st), ml, mt, wf)
cv <- put(cv, image_read(mods[1]), ml, y1, cw); cv <- put(cv, image_read(mods[2]), xr, y1, cw)
cv <- put(cv, image_read(mods[3]), ml, y2, cw); cv <- put(cv, image_read(mods[4]), xr, y2, cw)
page_pdf(cv, L(letters[1:5], c(0.03, 0.03, xr - 0.19, 0.03, xr - 0.19), c(mt, y1, y1, y2, y2)), "Figure2.pdf")
arch_move(file.path(mp, "Figure2a.pdf"))

# ================= Figure 3: prevalence, correlations, Moffitt, Steele, a-d =================
pf <- file.path(fig, "aim2", c("Fig2B_per_patient_cytostatic.png", "Fig2C_axis_correlation_heatmap.png", "Fig2D_moffitt_alignment.png", "Fig2E_Steele_replication.png"))
stopifnot(all(file.exists(pf)))
ch <- (H - mt - mb - gy) / 2; yb <- mt + ch + gy
cv <- blank(); cv <- put(cv, image_read(pf[1]), ml, mt, cw); cv <- put(cv, image_read(pf[2]), xr, mt, cw)
cv <- put(cv, image_read(pf[3]), ml, yb, cw); cv <- put(cv, image_read(pf[4]), xr, yb, cw)
page_pdf(cv, L(letters[1:4], c(0.03, xr - 0.19, 0.03, xr - 0.19), c(mt, mt, yb, yb)), "Figure3.pdf")
arch_move(file.path(mp, "Figure2b.pdf"))

# ================= Figure 5: spatial test + bulk survival =================
spots <- read.csv(file.path(pr, "aim3_spatial", "spot_level_primary_PDAC.csv"))
pm <- c("Treg activation" = "Treg", "Checkpoint molecules" = "Checkpoint_molecules",
        "Myeloid suppression" = "Immune_Suppression_by_Myeloid_cells", "Effector (control)" = "Effector_cells")
dist_to_nearest <- function(pts, anchors) apply(pts, 1, function(r) min(sqrt((r[1] - anchors[, 1])^2 + (r[2] - anchors[, 2])^2)))
heat <- data.frame(); corr <- data.frame()
for (s in unique(spots$sample_id)) {
  df <- spots[spots$sample_id == s, ]; epi <- df[df$first_type == "Tumor Epithelial cells", c("imagecol", "imagerow")]
  if (nrow(epi) < 10) next
  ne <- df[df$first_type != "Tumor Epithelial cells", ]; ne$dist <- dist_to_nearest(ne[, c("imagecol", "imagerow")], as.matrix(epi))
  row <- data.frame(sample = sub("^IU_PDA_", "", s))
  for (lab in names(pm)) { v <- ne[[pm[[lab]]]]
    row[[lab]] <- suppressWarnings(cor(ne$dist, v, method = "spearman", use = "complete.obs"))
    b <- aggregate(z ~ bin, data.frame(z = as.numeric(scale(v)), bin = cut(ne$dist, breaks = 10, labels = FALSE)), mean)
    b$sample <- row$sample; b$programme <- lab; heat <- rbind(heat, b) }
  corr <- rbind(corr, row)
}
cat("Spatial samples used:", nrow(corr), "(expect 9) | medians:", paste(names(pm), round(sapply(corr[, names(pm)], median), 3), collapse = "; "), "\n")
heat$programme <- factor(heat$programme, levels = names(pm)); h1 <- 5.4
p5a <- ggplot(heat, aes(bin, sample, fill = z)) + geom_tile() + facet_wrap(~programme, ncol = 2) +
  scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00", name = "Mean\nz-score") +
  scale_x_continuous(breaks = c(1, 5, 10)) +
  labs(title = "Programme score by distance from tumour epithelium", subtitle = "GSE272362, 9 primary PDAC samples; distance deciles for display",
       x = "Distance decile (1 = nearest epithelium)", y = "Sample") + thm(fs)
f5a <- file.path(fig, "aim3", "Fig5A_distance_heatmap.png"); savep(p5a, f5a, cw, h1)
cl <- data.frame(sample = rep(corr$sample, length(pm)), programme = rep(names(pm), each = nrow(corr)), rho = unlist(corr[, names(pm)]))
cl$programme <- factor(cl$programme, levels = names(pm)); set.seed(3)
p5b <- ggplot(cl, aes(programme, rho, fill = programme)) + geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.6) + geom_jitter(width = 0.12, size = 1.1) +
  scale_fill_manual(values = c("#56B4E9", "#0072B2", "#D55E00", "#009E73"), guide = "none") +
  labs(title = "Distance-score correlation per sample", subtitle = "Spearman rho per sample (n = 9); positive = higher away from epithelium",
       x = NULL, y = "Spearman rho") + thm(fs) + theme(axis.text.x = element_text(angle = 30, hjust = 1))
f5b <- file.path(fig, "aim3", "Fig5B_distance_correlations.png"); savep(p5b, f5b, cw, h1)
km <- file.path(fig, "aim4", "Fig4A_KM_cytostatic.png"); pur <- file.path(fig, "aim4", "Fig4C_axis2_purity_TCGA.png"); stopifnot(file.exists(km), file.exists(pur))
y5 <- mt + h1 + gy
cv <- blank(); cv <- put(cv, image_read(f5a), ml, mt, cw); cv <- put(cv, image_read(f5b), xr, mt, cw)
cv <- put(cv, image_read(km), ml, y5, cw); cv <- put(cv, image_read(pur), xr, y5, cw)
page_pdf(cv, L(letters[1:4], c(0.03, xr - 0.19, 0.03, xr - 0.19), c(mt, mt, y5, y5)), "Figure5.pdf")

# ================= Supp Fig 3: integration QC + DUSP6 at patient level =================
v4  <- readRDS(file.path(pr, "aim2_scoring", "malignant_cells_scored_v4.rds"))
mcm <- readRDS(file.path(pr, "aim2_scoring", "malignant_cell_metadata.rds"))
tr  <- unique(mcm[, c("pid", "new_treatment")]); v4$trt <- tr$new_treatment[match(v4$pid, tr$pid)]; u <- v4[v4$trt == "Untreated", ]
u$pid <- factor(u$pid, levels = paste0("U", 1:18)); hA <- 5.0
s3a <- ggplot(u[sample(nrow(u)), ], aes(UMAP_1_ut, UMAP_2_ut, colour = pid)) + geom_point(size = 0.45, alpha = 0.9, stroke = 0) +
  guides(colour = guide_legend(ncol = 2, override.aes = list(size = 2, alpha = 1))) +
  labs(title = "Untreated malignant cells after Harmony integration", subtitle = "UMAP on Harmony dimensions 1-20; coloured by patient (n = 18)",
       x = "UMAP 1", y = "UMAP 2", colour = "Patient") + thm(fs)
fs3a <- file.path(fig, "aim2", "FigS3A_UMAP_patients.png"); savep(s3a, fs3a, wf, hA)
m <- readRDS(file.path(pr, "aim2_scoring", "malignant_count_matrix.rds"))
if (!("DUSP6" %in% rownames(m)) && ("DUSP6" %in% colnames(m))) m <- t(m)
idx <- match(u$cell_id, colnames(m)); ok <- !is.na(idx); det <- as.numeric(m["DUSP6", idx[ok]] > 0); rm(m); invisible(gc())
pids <- levels(u$pid)
dd <- data.frame(pid = pids, pct = 100 * as.vector(tapply(det, u$pid[ok], mean)[pids]), a1 = as.vector(tapply(u$axis1_ras_balance[ok], u$pid[ok], mean)[pids]))
ct <- suppressWarnings(cor.test(dd$pct, dd$a1, method = "spearman"))
cat(sprintf("DUSP6 vs Axis 1: rho = %.3f, p = %.3f (ledger: 0.205, 0.41)\n", ct$estimate, ct$p.value))
s3c <- ggplot(dd, aes(a1, pct)) + geom_point(size = 1.6, colour = "#0072B2") + geom_text(aes(label = pid), size = 2.1, vjust = -0.8) +
  annotate("text", x = min(dd$a1), y = max(dd$pct) * 1.1, hjust = 0, vjust = 1, size = fs / ggplot2::.pt,
           label = sprintf("Spearman rho = %.2f, p = %.2f\nn = 18 patients", ct$estimate, ct$p.value)) +
  labs(title = "DUSP6 detection vs Axis 1 (patient level)", subtitle = "Untreated malignant cells",
       x = "Mean Axis 1 RAS balance score", y = "Cells with DUSP6 detected (%)") + thm(fs)
fs3c <- file.path(fig, "aim2", "FigS3C_DUSP6_vs_axis1.png"); savep(s3c, fs3c, cw, cw)
elb <- file.path(fig, "aim2", "FigS3C_PCA_elbow.png"); stopifnot(file.exists(elb))
yS3 <- mt + hA + gy
cv <- blank(); cv <- put(cv, image_read(fs3a), ml, mt, wf); cv <- put(cv, image_read(elb), ml, yS3, cw); cv <- put(cv, image_read(fs3c), xr, yS3, cw)
page_pdf(cv, L(letters[1:3], c(0.03, 0.03, xr - 0.19), c(mt, yS3, yS3)), "FigureS3.pdf")
for (o in c("FigS3A_UMAP_harmony_QC.png", "FigS3B_UMAP_patient_QC.png", "FigS3D_axis3_response.png", "FigS3E_DUSP6_validation.png")) arch_move(file.path(fig, "aim2", o))

# ================= Supp Fig 4: spatial maps + spot-type distributions =================
maps <- file.path(fig, "aim3", "Fig3C_spatial_Treg_maps.png"); spt <- file.path(fig, "aim3", "Fig3D_spot_type_scores.png"); stopifnot(file.exists(maps), file.exists(spt))
ha <- wf * aspect(maps); hb <- wf * aspect(spt); sc <- min(1, (H - mt - gy - mb) / (ha + hb)); w4 <- wf * sc
cv <- blank(); cv <- put(cv, image_read(maps), ml, mt, w4); cv <- put(cv, image_read(spt), ml, mt + ha * sc + gy, w4)
page_pdf(cv, L(c("a", "b"), c(0.03, 0.03), c(mt, mt + ha * sc + gy)), "FigureS4.pdf")
arch_move(file.path(fig, "aim3", "SuppFig4_regenerated.png"))

# ================= Supp Fig 5 rename; superseded panels =================
s5o <- file.path(mp, "FigureS5_classifier_sensitivity.pdf")
if (file.exists(s5o)) { file.copy(s5o, file.path(mp, "FigureS5.pdf"), overwrite = TRUE); arch_move(s5o); cat("Supp Fig 5 is now FigureS5.pdf\n") }
arch_move(file.path(fig, "aim4", "Fig4b_regenerated.png"))                 # replaced by aim3/Fig3C_DepMap_essentiality.png

# ================= figure map =================
dir.create(file.path(base, "docs"), showWarnings = FALSE)
writeLines(c("# Figure map (final_figures_v2.R)", "",
  "| Figure | Multipanel (PDF) | Panels (PNG) |", "|---|---|---|",
  "| 1 | Figure1.pdf | aim1/Fig1A_volcano_shared_resistance.png (a); aim1/Fig1B_GSEA_Hallmark_dotplot.png (b) |",
  "| 2 | Figure2.pdf | aim2/Fig2A_UMAP_state.png (a); aim2/Fig2A_UMAP_axis1-4.png (b-e) |",
  "| 3 | Figure3.pdf | aim2/Fig2B_per_patient_cytostatic.png (a); Fig2C_axis_correlation_heatmap.png (b); Fig2D_moffitt_alignment.png (c); Fig2E_Steele_replication.png (d) |",
  "| 4 | Figure4.pdf | aim3/Fig3A_TGM2_paired_expression.png (a); aim3/Fig3B_TGM2_ITGA4_rank_scores.png (b); aim3/Fig3C_DepMap_essentiality.png (c) |",
  "| 5 | Figure5.pdf | aim3/Fig5A_distance_heatmap.png (a); aim3/Fig5B_distance_correlations.png (b); aim4/Fig4A_KM_cytostatic.png (c); aim4/Fig4C_axis2_purity_TCGA.png (d) |",
  "| S1 | FigureS1.pdf | aim1/FigS1A-E (QC) |",
  "| S2 | FigureS2.pdf | aim1/FigS2A_KRAS_gene_correlations.png (a); aim1/FigS2B_DepMap_module_crossval.png (b) |",
  "| S3 | FigureS3.pdf | aim2/FigS3A_UMAP_patients.png (a); aim2/FigS3C_PCA_elbow.png (b); aim2/FigS3C_DUSP6_vs_axis1.png (c) |",
  "| S4 | FigureS4.pdf | aim3/Fig3C_spatial_Treg_maps.png (a); aim3/Fig3D_spot_type_scores.png (b) |",
  "| S5 | FigureS5.pdf | aim2/FigS5_classifier_sensitivity.png (a); aim2/FigS5b_classifier_sweep.png (b) |", "",
  "Panel file names keep their original prefixes; this table is the authority for which figure each panel belongs to."),
  file.path(base, "docs", "figure_map.md"))
cat("\nWrote docs/figure_map.md\nMultipanel folder now:\n"); print(list.files(mp))
