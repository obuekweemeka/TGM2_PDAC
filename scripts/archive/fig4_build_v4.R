# ============================================================================
# fig4_build_v4.R -- Figure 4 with harmonised text sizes. Save in: scripts/
# Panels a (KM) and b (purity) are existing 600-dpi PNGs, placed at exact sizes;
# panel c (DepMap 2x2 block) is regenerated at the exact size it occupies, with
# font sizes chosen so tick labels match a and b (about 5 pt on the page).
# Safe to run from a fresh session. Writes to a temp PDF first, then replaces.
# ============================================================================
suppressPackageStartupMessages({library(ggplot2); library(patchwork); library(magick)})
base <- here::here(); fig <- file.path(base, "figures")
cb_orange <- "#E69F00"; cb_blue <- "#0072B2"

# ---------- parameters (edit to tune) ----------
W <- 8.27; H <- 11.69; dpi <- 600                  # A4 portrait
m_left <- 0.22; m_right <- 0.06; m_top <- 0.06; m_bot <- 0.06
gap_x <- 0.22; gap_y <- 0.12                       # gap_x leaves room for the letter 'b'
fs <- 6.3                                          # base font (pt) for panel c; tick labels = 0.8 * fs
lab_pt <- 9                                        # panel-letter size (pt)

# ---------- geometry from the real image sizes ----------
imgA <- image_read(file.path(fig, "aim4", "Fig4A_KM_cytostatic.png"))
imgB <- image_read(file.path(fig, "aim4", "Fig4C_axis2_purity_TCGA.png"))
ia <- image_info(imgA); ib <- image_info(imgB)
cw  <- (W - m_left - m_right - gap_x) / 2
h_a <- cw * ia$height / ia$width; h_b <- cw * ib$height / ib$width
y_c <- m_top + max(h_a, h_b) + gap_y
w_c <- W - m_left - m_right; h_c <- H - y_c - m_bot
cat(sprintf("a: %.2f x %.2f in | b: %.2f x %.2f in | c: %.2f x %.2f in\n", cw, h_a, cw, h_b, w_c, h_c))
cat(sprintf("Expected tick-label size on the page: a %.1f pt, b %.1f pt, c %.1f pt\n",
            9.05 * cw / (ia$width / dpi), 9.05 * cw / (ib$width / dpi), 0.8 * fs))

# ---------- panel c: DepMap block, 44 lines, re-split at the 44-line median ----------
dep_expr <- readRDS(file.path(base, "data", "processed", "aim1_signatures", "DepMap_PDAC_expression.rds"))
dep4     <- readRDS(file.path(base, "data", "processed", "aim4_bulk", "depmap_pdac_analysis_df.rds"))
valid44  <- intersect(dep_expr$ModelID[!duplicated(dep_expr$ModelID)], dep4$model_id)
d <- dep4[dep4$model_id %in% valid44, ]; stopifnot(nrow(d) == 44)
d$grp <- ifelse(d$axis2_death < median(d$axis2_death), "High cytostatic", "Low cytostatic")
hi <- d$grp == "High cytostatic"
pct <- sapply(c(tgm2 = "tgm2", nt5e = "nt5e", kras = "kras"), function(g) round(100 * mean(d[[g]] > 0), 1))
cat("Groups:", table(d$grp), "| % Chronos > 0:", paste(names(pct), pct, collapse = ", "), "\n")
leg <- c("High cytostatic" = paste0("High cytostatic (n = ", sum(hi), ")"), "Low cytostatic" = paste0("Low cytostatic (n = ", sum(!hi), ")"))
panel <- function(g, ttl)
  ggplot(d, aes(x = axis2_death, y = .data[[g]], colour = grp)) +
    geom_point(size = 1.3, alpha = 0.85) +
    geom_smooth(method = "lm", se = TRUE, colour = "grey30", linewidth = 0.4, formula = y ~ x) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.3) +
    scale_colour_manual(values = c("High cytostatic" = cb_orange, "Low cytostatic" = cb_blue), labels = leg, name = NULL) +
    labs(title = ttl, subtitle = paste0("n = 44 PDAC lines; Chronos > 0: ", pct[[g]], "%"),
         x = "Axis 2 death competence score", y = "Chronos effect score") +
    theme_classic(base_size = fs) +
    theme(plot.title = element_text(face = "bold", size = fs + 1.2),
          plot.subtitle = element_text(size = fs - 0.8, colour = "grey40"),
          legend.text = element_text(size = fs))
block <- ((panel("tgm2", "TGM2") | panel("nt5e", "NT5E")) /
          (panel("kras", "KRAS (positive control)") | guide_area())) + plot_layout(guides = "collect")
out_c <- file.path(fig, "aim4", "Fig4b_regenerated.png")
ggsave(out_c, block, dpi = dpi, width = w_c, height = h_c, units = "in")
cat("Saved:", out_c, "-", round(file.size(out_c) / 1024), "KB\n")

# ---------- compose the page at exact positions ----------
px <- function(x) as.integer(round(x * dpi))
put <- function(cv, img, x_in, y_in, w_in)
  image_composite(cv, image_scale(img, paste0(px(w_in), "x")), offset = paste0("+", px(x_in), "+", px(y_in)))
canvas <- image_blank(px(W), px(H), color = "white")
canvas <- put(canvas, imgA, m_left, m_top, cw)
canvas <- put(canvas, imgB, m_left + cw + gap_x, m_top, cw)
canvas <- put(canvas, image_read(out_c), m_left, y_c, w_c)

lab_df <- data.frame(l = c("a", "b", "c"),
                     x = c(0.03, m_left + cw + gap_x - 0.19, 0.03),
                     y = c(m_top, m_top, y_c))
p <- ggplot() +
  annotation_raster(as.raster(canvas), xmin = 0, xmax = W, ymin = 0, ymax = H, interpolate = TRUE) +
  geom_text(data = lab_df, aes(x = x, y = H - y, label = l), hjust = 0, vjust = 1,
            fontface = "bold", size = lab_pt / ggplot2::.pt) +
  coord_fixed(xlim = c(0, W), ylim = c(0, H), expand = FALSE) +
  theme_void() + theme(plot.margin = margin(0, 0, 0, 0))

out <- file.path(fig, "multipanel", "Figure4.pdf"); tmp <- file.path(fig, "multipanel", "Figure4_new.pdf")
ggsave(tmp, p, width = W, height = H, units = "in", device = "pdf")
sz <- round(file.size(tmp) / 1024); cat("Figure4_new.pdf:", sz, "KB\n"); stopifnot(sz > 100)
file.copy(tmp, out, overwrite = TRUE); file.remove(tmp)
cat("Figure4.pdf replaced:", round(file.size(out) / 1024), "KB\n")
