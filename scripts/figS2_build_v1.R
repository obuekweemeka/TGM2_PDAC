# ============================================================================
# figS2_build_v1.R -- Supplementary Figure 2 multipanel (A4). Save in: scripts/
#  a: histogram of gene-level Pearson r (expression vs KRAS Chronos) across the
#     44-line DepMap PDAC panel, regenerated from data at its exact page size
#  b: existing module-level bar chart (FigS2B_DepMap_module_crossval.png; its
#     values already match the verified 44-line module results)
# Also checks that stored DepMap-correlate module genes fall beyond r = -0.4.
# Old FigureS2.pdf and the unused MAPK scatter are archived, not deleted.
# ============================================================================
suppressPackageStartupMessages({library(ggplot2); library(magick)})
base <- normalizePath('.'); fig <- file.path(base, "figures")
W <- 8.27; H <- 11.69; dpi <- 600
m_left <- 0.22; m_right <- 0.06; m_top <- 0.06; m_bot <- 0.06; gap_y <- 0.12; lab_pt <- 9
fs_a <- 10.5                     # base font for panel a: tick labels 0.8 * 10.5 = 8.4 pt, matching panel b on the page

# ---------- data: gene-level r across the 44 lines ----------
dep_expr <- readRDS(file.path(base, "data", "processed", "aim1_signatures", "DepMap_PDAC_expression.rds"))
dep4     <- readRDS(file.path(base, "data", "processed", "aim4_bulk", "depmap_pdac_analysis_df.rds"))
dep_expr <- dep_expr[!duplicated(dep_expr$ModelID), ]
valid44  <- intersect(dep_expr$ModelID, dep4$model_id)
e <- as.data.frame(dep_expr[match(valid44, dep_expr$ModelID), ]); k <- dep4$kras[match(valid44, dep4$model_id)]
num <- names(e)[sapply(e, is.numeric)]
X <- as.matrix(e[, num]); X <- X[, apply(X, 2, sd, na.rm = TRUE) > 0, drop = FALSE]
r <- as.numeric(cor(X, k, use = "pairwise.complete.obs")); names(r) <- colnames(X)
cat("Lines:", length(k), "| genes tested:", length(r), "| r < -0.4:", sum(r < -0.4, na.rm = TRUE), "| r > 0.4:", sum(r > 0.4, na.rm = TRUE), "\n")
sig <- readRDS(file.path(base, "data", "processed", "aim1_signatures", "four_axis_signatures_v1.rds"))
for (mod in c("axis1_depmap_correlates", "axis1_mapk_output", "axis1_bypass_potential")) {
  g <- intersect(sig[[mod]], names(r))
  cat(sprintf("%-26s %2d genes | r < -0.4: %2d | r > 0.4: %2d | median r = %.2f\n",
              mod, length(g), sum(r[g] < -0.4), sum(r[g] > 0.4), median(r[g])))
}

# ---------- geometry: b at full width, a fills the rest of the page ----------
imgB <- image_read(file.path(fig, "aim1", "FigS2B_DepMap_module_crossval.png")); ib <- image_info(imgB)
w_full <- W - m_left - m_right
h_b <- w_full * ib$height / ib$width
h_a <- H - m_top - gap_y - h_b - m_bot
cat(sprintf("a: %.2f x %.2f in | b: %.2f x %.2f in (b text scale %.2f)\n", w_full, h_a, w_full, h_b, w_full / (ib$width / dpi)))

# ---------- panel a ----------
dfr <- data.frame(r = r[!is.na(r)])
dfr$class <- factor(ifelse(dfr$r < -0.4, "KRAS-dependent (r < -0.4)", ifelse(dfr$r > 0.4, "KRAS-independent (r > 0.4)", "Not selected")),
                    levels = c("KRAS-dependent (r < -0.4)", "KRAS-independent (r > 0.4)", "Not selected"))
pa <- ggplot(dfr, aes(x = r, fill = class)) +
  geom_histogram(binwidth = 0.01, boundary = 0) +
  geom_vline(xintercept = c(-0.4, 0.4), linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  scale_fill_manual(values = c("KRAS-dependent (r < -0.4)" = "#E07B6A", "KRAS-independent (r > 0.4)" = "#4479B7", "Not selected" = "grey80"), name = NULL) +
  labs(title = "Gene-level correlation with KRAS dependency in DepMap PDAC cell lines",
       subtitle = paste0("Pearson r between gene expression and KRAS Chronos effect; n = ", length(k), " lines, ", nrow(dfr), " genes"),
       x = "Pearson r (expression vs KRAS Chronos effect)", y = "Number of genes") +
  theme_classic(base_size = fs_a) +
  theme(plot.title = element_text(face = "bold", size = fs_a + 1), plot.subtitle = element_text(size = fs_a - 1.5, colour = "grey40"),
        legend.position = "top")
out_a <- file.path(fig, "aim1", "FigS2A_KRAS_gene_correlations.png")
ggsave(out_a, pa, dpi = dpi, width = w_full, height = h_a, units = "in")
cat("Saved:", out_a, "-", round(file.size(out_a) / 1024), "KB\n")

# ---------- compose A4 page ----------
px <- function(x) as.integer(round(x * dpi))
put <- function(cv, img, x_in, y_in, w_in)
  image_composite(cv, image_scale(img, paste0(px(w_in), "x")), offset = paste0("+", px(x_in), "+", px(y_in)))
canvas <- image_blank(px(W), px(H), color = "white")
canvas <- put(canvas, image_read(out_a), m_left, m_top, w_full)
canvas <- put(canvas, imgB, m_left, m_top + h_a + gap_y, w_full)
lab_df <- data.frame(l = c("a", "b"), x = c(0.03, 0.03), y = c(m_top, m_top + h_a + gap_y))
p <- ggplot() +
  annotation_raster(as.raster(canvas), xmin = 0, xmax = W, ymin = 0, ymax = H, interpolate = TRUE) +
  geom_text(data = lab_df, aes(x = x, y = H - y, label = l), hjust = 0, vjust = 1, fontface = "bold", size = lab_pt / ggplot2::.pt) +
  coord_fixed(xlim = c(0, W), ylim = c(0, H), expand = FALSE) + theme_void() + theme(plot.margin = margin(0, 0, 0, 0))

out <- file.path(fig, "multipanel", "FigureS2.pdf"); tmp <- file.path(fig, "multipanel", "FigureS2_new.pdf")
ggsave(tmp, p, width = W, height = H, units = "in", device = "pdf")
sz <- round(file.size(tmp) / 1024); cat("FigureS2_new.pdf:", sz, "KB\n"); stopifnot(sz > 100)
arch <- file.path(fig, "_archive", "FigureS2_superseded.pdf")
if (file.exists(out) && !file.exists(arch)) { file.copy(out, arch); cat("Old FigureS2.pdf archived\n") }
file.copy(tmp, out, overwrite = TRUE); file.remove(tmp)
cat("FigureS2.pdf replaced:", round(file.size(out) / 1024), "KB\n")

# ---------- the MAPK scatter is not part of Figure S2: archive it ----------
sc <- file.path(fig, "aim1", "SuppFig2a_regenerated.png")
if (file.exists(sc)) { file.rename(sc, file.path(fig, "_archive", "SuppFig2a_MAPK_scatter_unused.png")); cat("Unused MAPK scatter archived\n") }
