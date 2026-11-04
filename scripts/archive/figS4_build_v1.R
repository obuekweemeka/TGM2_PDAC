# ============================================================================
# figS4_build_v1.R -- Supplementary Figure 4 multipanel (A4). Save in: scripts/
# The regenerated Supp Fig 4 PNG already contains panels a (heatmap) and
# b (boxplot), so this places it on an A4 page and writes the PDF safely.
# The old FigureS4.pdf (binned, n = 6) is copied to _archive/ before replacement.
# ============================================================================
suppressPackageStartupMessages({library(ggplot2); library(magick)})
base <- here::here(); fig <- file.path(base, "figures")
W <- 8.27; H <- 11.69; dpi <- 600; m <- 0.15

src <- file.path(fig, "aim3", "SuppFig4_regenerated.png"); stopifnot(file.exists(src))
img <- image_read(src); inf <- image_info(img)
w <- W - 2 * m; h <- w * inf$height / inf$width
cat(sprintf("Supp Fig 4 placed at %.2f x %.2f in on A4 (scale %.2f of native size)\n", w, h, w / (inf$width / dpi)))

px <- function(x) as.integer(round(x * dpi))
canvas <- image_composite(image_blank(px(W), px(H), color = "white"),
                          image_scale(img, paste0(px(w), "x")), offset = paste0("+", px(m), "+", px(m)))
p <- ggplot() +
  annotation_raster(as.raster(canvas), xmin = 0, xmax = W, ymin = 0, ymax = H, interpolate = TRUE) +
  coord_fixed(xlim = c(0, W), ylim = c(0, H), expand = FALSE) +
  theme_void() + theme(plot.margin = margin(0, 0, 0, 0))

out <- file.path(fig, "multipanel", "FigureS4.pdf"); tmp <- file.path(fig, "multipanel", "FigureS4_new.pdf")
ggsave(tmp, p, width = W, height = H, units = "in", device = "pdf")
sz <- round(file.size(tmp) / 1024); cat("FigureS4_new.pdf:", sz, "KB\n"); stopifnot(sz > 100)
arch <- file.path(fig, "_archive", "FigureS4_superseded.pdf")
if (file.exists(out) && !file.exists(arch)) { file.copy(out, arch); cat("Old FigureS4.pdf archived\n") }
file.copy(tmp, out, overwrite = TRUE); file.remove(tmp)
cat("FigureS4.pdf replaced:", round(file.size(out) / 1024), "KB\n")
