# ============================================================================
# rebuild_figures_v3.R  --  regenerates every corrected figure in one run
# (v3: Figure 4 PDF is written to a temporary name first and only replaces the
#  real file if valid, so a failed write can no longer truncate the existing PDF)
# Save this file in:  ~/Desktop/rasolutepaper1/scripts/
# Safe to run from a fresh R session. Overwrites the files listed below.
#
# Outputs
#   figures/aim3/SuppFig4_regenerated.png        Supp Fig 4 (n = 9 samples, T9 excluded)
#   figures/aim1/SuppFig2a_regenerated.png       Supp Fig 2a (MAPK output module vs KRAS Chronos)
#   figures/aim4/Fig4b_regenerated.png           Fig 4 DepMap block (44 lines, re-split at 44-line median)
#   results/aim4/depmap_essentiality_results_44lines.csv
#   figures/multipanel/Figure4.pdf               Figure 4, A4: (a | b) / c
# ============================================================================
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(patchwork); library(tidyr)})
if (!requireNamespace("magick", quietly = TRUE)) install.packages("magick")
library(magick)

base <- here::here(); fig <- file.path(base, "figures")
cb_skyblue <- "#56B4E9"; cb_vermillion <- "#D55E00"; cb_blue <- "#0072B2"; cb_orange <- "#E69F00"; cb_teal <- "#009E73"
tagtheme <- theme(plot.tag = element_text(size = 14, face = "bold"), plot.tag.position = "topleft")
for (dd in c("aim1", "aim3", "aim4", "multipanel")) dir.create(file.path(fig, dd), recursive = TRUE, showWarnings = FALSE)

# ============================ 1. SUPP FIG 4 ================================
cat("\n=== Supplementary Figure 4 ===\n")
spots <- read.csv(file.path(base, "data", "processed", "aim3_spatial", "spot_level_primary_PDAC.csv"))
programme_map <- c(Treg = "Treg", Checkpoint = "Checkpoint_molecules",
                   ImmuneSupp_Myeloid = "Immune_Suppression_by_Myeloid_cells", Effector = "Effector_cells")
dist_to_nearest <- function(pts, anchors)
  apply(pts, 1, function(r) min(sqrt((r[1] - anchors[, 1])^2 + (r[2] - anchors[, 2])^2)))
heat <- data.frame(); corr <- data.frame(); used <- c()
for (s in unique(spots$sample_id)) {
  df  <- spots[spots$sample_id == s, ]
  epi <- df[df$first_type == "Tumor Epithelial cells", c("imagecol", "imagerow")]
  if (nrow(epi) < 10) { cat("Excluded", s, "- only", nrow(epi), "epithelial anchor spots\n"); next }
  used <- c(used, s)
  ne <- df[df$first_type != "Tumor Epithelial cells", ]
  ne$dist <- dist_to_nearest(ne[, c("imagecol", "imagerow")], as.matrix(epi))
  row <- data.frame(sample = s)
  for (lab in names(programme_map)) {
    col <- programme_map[[lab]]
    row[[lab]] <- suppressWarnings(cor(ne$dist, ne[[col]], method = "spearman", use = "complete.obs"))
    z <- as.numeric(scale(ne[[col]])); bin <- cut(ne$dist, breaks = 10, labels = FALSE)
    b <- aggregate(z ~ bin, FUN = mean); b$sample <- s; b$programme <- lab; heat <- rbind(heat, b)
  }
  corr <- rbind(corr, row)
}
cat("Samples used (expect 9):", length(used), "\n")
print(round(sapply(corr[, names(programme_map)], median), 3))
heat$programme <- factor(heat$programme, levels = names(programme_map))
p_a <- ggplot(heat, aes(x = bin, y = sample, fill = z)) + geom_tile() + facet_wrap(~programme, ncol = 4) +
  scale_fill_gradientn(colours = c("#EEEEEE", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B"), name = "Mean z-score") +
  labs(title = "Programme score by distance decile from tumour epithelium",
       subtitle = "GSE272362, n = 9 samples, per-spot data (deciles for display only)",
       x = "Distance decile (1 = nearest epithelium)", y = NULL) +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12), plot.subtitle = element_text(size = 9, colour = "grey40"),
        axis.text.y = element_text(size = 8))
cl <- pivot_longer(corr, -sample, names_to = "programme", values_to = "rho")
cl$programme <- factor(cl$programme, levels = names(programme_map))
p_b <- ggplot(cl, aes(x = programme, y = rho, fill = programme)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) + geom_jitter(width = 0.15, size = 1.8, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  scale_fill_manual(values = c(cb_skyblue, cb_blue, cb_vermillion, cb_teal), guide = "none") +
  labs(title = "Per-sample Spearman correlation, distance vs programme score",
       subtitle = "n = 9 samples; effector programme included as specificity control", x = NULL, y = "Spearman rho") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12), plot.subtitle = element_text(size = 9, colour = "grey40"))
fig_s4 <- (p_a / p_b) + plot_layout(heights = c(1.1, 1.0)) + plot_annotation(tag_levels = "a", theme = tagtheme)
out1 <- file.path(fig, "aim3", "SuppFig4_regenerated.png")
ggsave(out1, fig_s4, dpi = 600, width = 8.0, height = 9.0, units = "in")
cat("Saved:", out1, "-", round(file.size(out1) / 1024), "KB\n")

# ============================ 2. SUPP FIG 2a ===============================
cat("\n=== Supplementary Figure 2a ===\n")
dep_expr <- readRDS(file.path(base, "data", "processed", "aim1_signatures", "DepMap_PDAC_expression.rds"))
dep4     <- readRDS(file.path(base, "data", "processed", "aim4_bulk", "depmap_pdac_analysis_df.rds"))
dep_expr <- dep_expr[!duplicated(dep_expr$ModelID), ]
valid44  <- intersect(dep_expr$ModelID, dep4$model_id)
genes    <- intersect(c("DUSP4", "DUSP6", "ETV4", "SPRY4", "EGR1", "FOSL1", "CCND1", "SPRED1"), colnames(dep_expr))
d2a <- dep_expr[dep_expr$ModelID %in% valid44, c("ModelID", genes)]
d2a$mapk_score <- rowMeans(d2a[, genes], na.rm = TRUE)
d2a <- merge(d2a, dep4[, c("model_id", "kras")], by.x = "ModelID", by.y = "model_id")
ct <- suppressWarnings(cor.test(d2a$mapk_score, d2a$kras, method = "spearman"))
cat("Plotted: n =", nrow(d2a), "| Spearman r =", round(ct$estimate, 3), "| p =", signif(ct$p.value, 3), "(manuscript: -0.308, 0.042)\n")
p_s2a <- ggplot(d2a, aes(x = mapk_score, y = kras)) +
  geom_point(size = 2.2, colour = cb_skyblue, alpha = 0.85) +
  geom_smooth(method = "lm", se = TRUE, colour = cb_vermillion, linewidth = 0.7, formula = y ~ x) +
  labs(title = "Axis 1 MAPK output module vs KRAS Chronos essentiality",
       subtitle = paste0("n = ", nrow(d2a), " PDAC lines; Spearman r = ", round(ct$estimate, 3), ", p = ", signif(ct$p.value, 3)),
       x = "MAPK output module score (mean expression)", y = "KRAS Chronos effect score") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12), plot.subtitle = element_text(size = 9, colour = "grey40"))
out2 <- file.path(fig, "aim1", "SuppFig2a_regenerated.png")
ggsave(out2, p_s2a, dpi = 600, width = 8.0, height = 5.5, units = "in")
cat("Saved:", out2, "-", round(file.size(out2) / 1024), "KB\n")

# ============================ 3. FIG 4b DepMap block =======================
cat("\n=== Figure 4 DepMap block ===\n")
d <- dep4[dep4$model_id %in% valid44, ]; stopifnot(nrow(d) == 44)
d$grp <- ifelse(d$axis2_death < median(d$axis2_death), "High cytostatic", "Low cytostatic")   # re-split at the 44-line median
hi <- d$grp == "High cytostatic"
cat("Group sizes:", table(d$grp), "| lines whose group changed vs stored grouping:", sum(d$grp != d$cytostatic_group), "\n")
res <- do.call(rbind, lapply(c("tgm2", "nt5e", "kras"), function(g) data.frame(
  gene = toupper(g), n_high = sum(hi), n_low = sum(!hi),
  median_high = round(median(d[[g]][hi]), 3), median_low = round(median(d[[g]][!hi]), 3),
  pct_chronos_positive = round(100 * mean(d[[g]] > 0), 1),
  wilcoxon_p = signif(suppressWarnings(wilcox.test(d[[g]][hi], d[[g]][!hi])$p.value), 3))))
print(res, row.names = FALSE)
write.csv(res, file.path(base, "results", "aim4", "depmap_essentiality_results_44lines.csv"), row.names = FALSE)
lab <- c("High cytostatic" = paste0("High cytostatic (n = ", sum(hi), ")"), "Low cytostatic" = paste0("Low cytostatic (n = ", sum(!hi), ")"))
panel <- function(g, ttl) {
  pct <- res$pct_chronos_positive[res$gene == toupper(g)]
  ggplot(d, aes(x = axis2_death, y = .data[[g]], colour = grp)) +
    geom_point(size = 1.8, alpha = 0.85) +
    geom_smooth(method = "lm", se = TRUE, colour = "grey30", linewidth = 0.5, formula = y ~ x) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    scale_colour_manual(values = c("High cytostatic" = cb_orange, "Low cytostatic" = cb_blue), labels = lab, name = NULL) +
    labs(title = ttl, subtitle = paste0("n = 44 PDAC lines; Chronos > 0: ", pct, "%"),
         x = "Axis 2 death competence score", y = "Chronos effect score") +
    theme_classic(base_size = 9) +
    theme(plot.title = element_text(face = "bold", size = 10), plot.subtitle = element_text(size = 8, colour = "grey40"))
}
block <- ((panel("tgm2", "TGM2") | panel("nt5e", "NT5E")) /
          (panel("kras", "KRAS (positive control)") | guide_area())) + plot_layout(guides = "collect")
out_b <- file.path(fig, "aim4", "Fig4b_regenerated.png")
ggsave(out_b, block, dpi = 600, width = 8.27, height = 6.6, units = "in")
cat("Saved:", out_b, "-", round(file.size(out_b) / 1024), "KB\n")

# ============================ 4. FIGURE 4 multipanel =======================
cat("\n=== Figure 4 multipanel ===\n")
load_panel <- function(path) {
  stopifnot(file.exists(path))
  img <- magick::image_read(path); inf <- magick::image_info(img); w <- inf$width; h <- inf$height
  list(plot = ggplot() + annotation_raster(as.raster(img), xmin = 0, xmax = w, ymin = 0, ymax = h, interpolate = TRUE) +
         coord_fixed(xlim = c(0, w), ylim = c(0, h), expand = FALSE) + theme_void() + theme(plot.margin = margin(2, 2, 2, 2)),
       aspect = h / w, w = w, h = h)
}
pa <- load_panel(file.path(fig, "aim4", "Fig4A_KM_cytostatic.png"))       # a: KM
pb <- load_panel(file.path(fig, "aim4", "Fig4C_axis2_purity_TCGA.png"))   # b: purity
pc <- load_panel(out_b)                                                    # c: DepMap block
W <- 8.27; H <- 11.69
top_in <- max(pa$aspect, pb$aspect) * W / 2; bot_in <- pc$aspect * W
cat(sprintf("Row heights: top %.2f + bottom %.2f = %.2f of %.2f in (%.0f%% of A4)\n", top_in, bot_in, top_in + bot_in, H, 100 * (top_in + bot_in) / H))
fig4 <- ((pa$plot | pb$plot) / pc$plot) + plot_layout(heights = c(top_in, bot_in)) + plot_annotation(tag_levels = "a", theme = tagtheme)
out <- file.path(fig, "multipanel", "Figure4.pdf")
tmp <- file.path(fig, "multipanel", "Figure4_new.pdf")
ggsave(tmp, fig4, width = W, height = H, units = "in", device = "pdf")
sz <- round(file.size(tmp) / 1024); cat("Figure4_new.pdf:", sz, "KB\n"); stopifnot(sz > 100)
file.copy(tmp, out, overwrite = TRUE); file.remove(tmp)
cat("Figure4.pdf replaced:", round(file.size(out) / 1024), "KB\n")
cat("\nAll figures rebuilt.\n")
