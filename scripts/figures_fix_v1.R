## figures_fix_v1.R
## Regenerates the faulty figure panels found in the v39 figure audit and rebuilds the
## affected A4 multipanels: Figure 1, 2, 3 and Supplementary Figures 1, 5, 6, 7.
## PART 1 computes every plotted value and checks it against the ledger; if any check fails
## the script stops before writing anything. PART 2 writes. Every file it replaces is first
## copied to figures/_archive/figfix_v1/.

suppressPackageStartupMessages({ library(ggplot2); library(magick) })
if (!requireNamespace("ggrepel", quietly = TRUE)) install.packages("ggrepel")
base <- normalizePath('.')
pr <- file.path(base, "data", "processed"); fig <- file.path(base, "figures"); mp <- file.path(fig, "multipanel")
ad <- file.path(fig, "_archive", "figfix_v1"); dir.create(ad, recursive = TRUE, showWarnings = FALSE)
W <- 8.27; H <- 11.69; dpi <- 600; ml <- .22; mr <- .06; mt <- .06; mb <- .06; gx <- .22; gy <- .12
wf <- W - ml - mr; cw <- (wf - gx) / 2; xr <- ml + cw + gx
COL_APOP <- "#0072B2"; COL_ARREST <- "#009E73"
chk <- function(what, got, want, tol) {
  ok <- length(got) == length(want) && all(abs(got - want) <= tol)
  cat(sprintf("  %-52s %s\n", what, if (ok) "OK" else paste("MISMATCH:", paste(signif(got, 5), collapse = ", "))))
  if (!ok) stop("Check failed (", what, "). Nothing has been written.")
}
thm <- function(fs) theme_classic(base_size = fs) +
  theme(plot.title = element_text(face = "bold", size = fs + 1), plot.subtitle = element_text(size = fs - 1.5, colour = "grey40"))
bk  <- function(p) if (file.exists(p) && !file.exists(file.path(ad, basename(p)))) file.copy(p, file.path(ad, basename(p)))
sv  <- function(p, f, w, h) { bk(f); ggsave(f, p, dpi = dpi, width = w, height = h, units = "in"); stopifnot(file.size(f) > 2e4) }
px  <- function(z) as.integer(round(z * dpi))
put <- function(cv, img, x, y, w) image_composite(cv, image_scale(img, paste0(px(w), "x")), offset = paste0("+", px(x), "+", px(y)))
hin <- function(img, w) { i <- image_info(img); w * i$height / i$width }
pin <- function(f) { i <- image_info(image_read(f)); c(w = i$width / dpi, h = i$height / dpi) }
page <- function(cv, labs, name) {
  pg <- ggplot() + annotation_raster(as.raster(cv), xmin = 0, xmax = W, ymin = 0, ymax = H, interpolate = TRUE) +
    geom_text(data = labs, aes(x = x, y = H - y, label = l), hjust = 0, vjust = 1, fontface = "bold", size = 9 / .pt) +
    coord_fixed(xlim = c(0, W), ylim = c(0, H), expand = FALSE) + theme_void() + theme(plot.margin = margin(0, 0, 0, 0))
  out <- file.path(mp, name); bk(out); tmp <- sub("\\.pdf$", "_new.pdf", out)
  ggsave(tmp, pg, width = W, height = H, units = "in", device = "pdf"); stopifnot(file.size(tmp) > 1e5)
  file.copy(tmp, out, overwrite = TRUE); file.remove(tmp); cat(sprintf("  %-14s rebuilt, %d KB\n", name, round(file.size(out) / 1024)))
}
L <- function(l, x, y) data.frame(l = l, x = x, y = y)

## =========================================================================== PART 1: compute and check
cat("=== PART 1: values checked against the ledger before any write ===\n")
v4  <- readRDS(file.path(pr, "aim2_scoring", "malignant_cells_scored_v4.rds"))
mcm <- readRDS(file.path(pr, "aim2_scoring", "malignant_cell_metadata.rds"))
tr  <- unique(mcm[, c("pid", "new_treatment")]); v4$trt <- tr$new_treatment[match(v4$pid, tr$pid)]
u   <- v4[v4$trt == "Untreated", ]; x <- u[u$axis2_death_competence != 0, ]
x$arrest <- x$axis2_death_competence < 0
chk("untreated cells / informative cells", c(nrow(u), nrow(x)), c(52999, 45174), 0)
chk("arrest side identical to stored cell_state label", mean(x$arrest == grepl("cytostatic", x$cell_state)), 1, 0)

## Figure 1b and 2a: per-patient proportions and Moffitt association
pid <- sort(unique(x$pid))
pm  <- data.frame(pid = pid, pct = 100 * as.vector(tapply(x$arrest, x$pid, mean)[pid]),
                  basal = as.vector(tapply(u$Moffitt_basal, u$pid, mean)[pid]),
                  classical = as.vector(tapply(u$Moffitt_classical, u$pid, mean)[pid]))
pm$subtype <- ifelse(pm$basal > pm$classical, "Basal-like", "Classical")
chk("per-patient mean / median / min / max (%)", c(mean(pm$pct), median(pm$pct), min(pm$pct), max(pm$pct)), c(84.9, 87.7, 69.0, 95.1), .05)
mc <- suppressWarnings(cor.test(pm$pct, pm$basal, method = "spearman", exact = TRUE))
chk("Moffitt basal Spearman rho / exact p", c(unname(mc$estimate), mc$p.value), c(-0.60, 0.011), c(.005, .0005))

## Figure 3a: TGM2 per patient, both cohorts
m <- readRDS(file.path(pr, "aim2_scoring", "malignant_count_matrix.rds"))
if ("TGM2" %in% colnames(m)) { tg <- as.numeric(m[, "TGM2"]); ids <- rownames(m) } else
if ("TGM2" %in% rownames(m)) { tg <- as.numeric(m["TGM2", ]); ids <- colnames(m) } else stop("TGM2 not found in count matrix")
x$tgm2 <- if (!is.null(ids) && mean(x$cell_id %in% ids) > .99) tg[match(x$cell_id, ids)] else
          if (length(tg) == nrow(v4)) tg[match(x$cell_id, v4$cell_id)] else stop("Cannot align TGM2 values to cells")
rm(m); invisible(gc())
ha_ <- tapply(x$tgm2[x$arrest], x$pid[x$arrest], mean); hp_ <- tapply(x$tgm2[!x$arrest], x$pid[!x$arrest], mean)
k <- intersect(names(ha_), names(hp_)); hw <- data.frame(unit = k, arrest = ha_[k], apop = hp_[k])
hw_p <- suppressWarnings(wilcox.test(hw$arrest, hw$apop, paired = TRUE))$p.value
chk("Hwang TGM2: patients / higher / mean diff / p", c(nrow(hw), sum(hw$arrest > hw$apop), mean(hw$arrest - hw$apop), hw_p),
    c(18, 12, .076, .0047), c(0, 0, .0015, .0002))
st <- read.csv(file.path(base, "results", "aim2", "steele_tgm2_replication_patientlevel.csv"), stringsAsFactors = FALSE)
st <- data.frame(unit = as.character(st$unit), arrest = st$cytostatic, apop = st$death_competent)
st_p <- suppressWarnings(wilcox.test(st$arrest, st$apop, paired = TRUE))$p.value
chk("Steele TGM2: patients / higher / mean diff / p", c(nrow(st), sum(st$arrest > st$apop), mean(st$arrest - st$apop), st_p),
    c(15, 13, .114, .0125), c(0, 0, .0015, .0002))
cat(sprintf("  exact p-values for display: Hwang %.5f -> %.3f | Steele %.5f -> %.3f\n", hw_p, hw_p, st_p, st_p))

## Supplementary Figure 1b: module-level correlations with KRAS Chronos (44 lines)
de  <- readRDS(file.path(pr, "aim1_signatures", "DepMap_PDAC_expression.rds"))
d4  <- readRDS(file.path(pr, "aim4_bulk", "depmap_pdac_analysis_df.rds"))
sig <- readRDS(file.path(pr, "aim1_signatures", "four_axis_signatures_v1.rds"))
e <- as.data.frame(de)[!duplicated(de$ModelID), ]; rownames(e) <- e$ModelID; e$ModelID <- NULL   # tibble -> data.frame keeps row names
d44 <- d4[d4$model_id %in% rownames(e), ]
gm <- function(g) { g <- intersect(g, names(e)); rowMeans(e[d44$model_id, g, drop = FALSE], na.rm = TRUE) }
lab1 <- c(axis1_depmap_correlates = "Axis 1: DepMap correlates", axis1_mapk_output = "Axis 1: MAPK output",
          axis1_bypass_potential = "Axis 1: Bypass potential", axis2_proapoptotic = "Axis 2: Pro-apoptotic",
          axis2_executors = "Axis 2: Caspase/mitochondrial", axis2_anti_apoptotic = "Axis 2: Anti-apoptotic",
          axis2_cytostasis = "Axis 2: Cytostasis", axis3_hypoxia_metabolism = "Axis 3: Hypoxia-metabolic",
          axis3_OXPHOS_depleted = "Axis 3: OXPHOS depleted", axis4_innate_mimicry_core = "Axis 4: Innate mimicry core",
          axis4_innate_mimicry_full = "Axis 4: Innate mimicry full")
stopifnot(all(names(lab1) %in% names(sig)))
chk("Supp 1b: n lines", nrow(d44), 44, 0)
chk("Supp 1b: module means all finite", mean(sapply(names(lab1), function(k) all(is.finite(gm(sig[[k]]))))), 1, 0)
s1 <- do.call(rbind, lapply(names(lab1), function(k) { ct <- suppressWarnings(cor.test(gm(sig[[k]]), d44$kras, method = "spearman"))
  data.frame(module = k, label = lab1[[k]], r = unname(ct$estimate), p = ct$p.value) }))
want <- c(axis1_depmap_correlates = -.830, axis1_mapk_output = -.308, axis1_bypass_potential = -.049, axis2_proapoptotic = -.364,
          axis2_executors = .274, axis2_anti_apoptotic = .037, axis2_cytostasis = -.067, axis3_hypoxia_metabolism = -.022,
          axis3_OXPHOS_depleted = .162, axis4_innate_mimicry_core = -.065, axis4_innate_mimicry_full = -.063)

chk("Supp 1b: 11 module rho values match current figure", s1$r[match(names(want), s1$module)], unname(want), .0015)

## Supplementary Figure 6: threshold sweep
ad2 <- x$axis2_death_competence; tmed <- median(ad2); t33 <- unname(quantile(ad2, .33)); at0 <- 100 * mean(ad2 < 0)
chk("Supp 6: pooled at 0 / median / 33rd percentile", c(at0, tmed, t33), c(86.3, -6.08, -7.30), c(.05, .005, .005))

## Supplementary Figure 5: UMAP columns
uc <- if (all(c("UMAP_1_ut", "UMAP_2_ut") %in% names(u)) && !anyNA(u$UMAP_1_ut)) c("UMAP_1_ut", "UMAP_2_ut") else c("UMAP_1", "UMAP_2")
chk("Supp 5: cells without signal", sum(u$axis2_death_competence == 0), 7825, 0)
cat("  Supp 5 UMAP coordinates:", paste(uc, collapse = ", "), "\n")

## Supplementary Figure 7a: locate the header (title + subtitle) of the map panel
maps_f <- file.path(fig, "aim3", "Fig3C_spatial_Treg_maps.png"); spt_f <- file.path(fig, "aim3", "Fig3D_spot_type_scores.png")
mi <- image_read(maps_f); info <- image_info(mi)
sm <- image_convert(image_scale(mi, "1200x"), colorspace = "gray"); a <- as.integer(sm[[1]])[, , 1]
top <- seq_len(round(nrow(a) * .25)); dk <- rowMeans(a[top, , drop = FALSE] < 128) > .001
rl <- rle(dk); ends <- cumsum(rl$lengths); starts <- ends - rl$lengths + 1; bands <- which(rl$values)
cat("  Supp 7a text bands in top quarter (rows at 1200 px width):", paste(sprintf("%d-%d", starts[bands], ends[bands]), collapse = ", "), "\n")
if (length(bands) < 3) stop("Could not locate the Supp 7a title, subtitle and facet labels. Nothing has been written.")
cut_sm <- floor((ends[bands[2]] + starts[bands[3]]) / 2); cut_px <- round(cut_sm * info$height / nrow(a))
cat("  Supp 7a: header removed above row", cut_px, "of", info$height, "and redrawn with the corrected credit\n")
cat("All checks passed.\n")

## =========================================================================== PART 2: write
cat("\n=== PART 2: regenerate panels and rebuild multipanels ===\n")

## ---- Figure 1: panel b at the width Figure 1 gives it
fa <- image_read(file.path(fig, "aim2", "Fig2A_contrast_distribution.png")); h1a <- hin(fa, wf)
y1b <- mt + h1a + gy; h1b <- H - y1b - mb; stopifnot(h1b > 3)
pb <- pm[order(-pm$pct), ]; pb$pid <- factor(pb$pid, levels = pb$pid)
p1b <- ggplot(pb, aes(pid, pct, fill = subtype)) + geom_col(width = .75) +
  geom_hline(yintercept = mean(pm$pct), linetype = "dashed", colour = "grey40") +
  scale_fill_manual(values = c("Basal-like" = "#D55E00", "Classical" = "#0072B2"), name = "Dominant Moffitt subtype") +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  labs(title = "Cytostatic-prone cells per patient",
       subtitle = sprintf("18 untreated patients, informative cells; dashed line = mean of patients (%.1f%%)", mean(pm$pct)),
       x = "Patient", y = "Cytostatic-prone cells (%)") +
  thm(9) + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
f1b <- file.path(fig, "aim2", "Fig1B_per_patient_cytostatic_fullwidth.png"); sv(p1b, f1b, wf, h1b)
cv <- image_blank(px(W), px(H), color = "white"); cv <- put(cv, fa, ml, mt, wf); cv <- put(cv, image_read(f1b), ml, y1b, wf)
page(cv, L(c("a", "b"), c(.03, .03), c(mt, y1b)), "Figure1.pdf")

## ---- Figure 2: panel a with the statistics in the empty corner and repelled labels
f2a <- file.path(fig, "aim2", "Fig2D_moffitt_alignment.png"); d2a <- pin(f2a)
p2a <- ggplot(pm, aes(pct, basal)) + geom_smooth(method = "lm", formula = y ~ x, colour = "grey35", fill = "grey88", linewidth = .6) +
  geom_point(colour = "#0072B2", size = 1.8) +
  ggrepel::geom_text_repel(aes(label = pid), size = 2.5, seed = 1, min.segment.length = .2, box.padding = .3) +
  annotate("text", x = max(pm$pct) + .5, y = max(pm$basal) + .04, hjust = 1, vjust = 1, size = 2.9,
           label = sprintf("Spearman r = %.2f\np = %.3f (exact)\nn = %d patients", mc$estimate, mc$p.value, nrow(pm))) +
  labs(title = "Cytostatic-prone proportion vs Moffitt basal score", subtitle = "Per-patient analysis; untreated patients",
       x = "Cytostatic-prone cells (%)", y = "Mean Moffitt basal score") + thm(9)
sv(p2a, f2a, d2a[["w"]], d2a[["h"]])
cv <- image_blank(px(W), px(H), color = "white")
cv <- put(cv, image_read(f2a), ml, mt, cw); cv <- put(cv, image_read(file.path(fig, "aim2", "Fig2E_Steele_replication.png")), xr, mt, cw)
page(cv, L(c("a", "b"), c(.03, xr - .19), c(mt, mt)), "Figure2.pdf")

## ---- Figure 3: panel b relabelled from its own build script, panel a rebuilt with the patient-level Steele test
f3b <- file.path(fig, "aim3", "Fig3B_TGM2_ITGA4_rank_scores.png"); bk(f3b); bk(file.path(fig, "aim3", "Fig3A_TGM2_paired_expression.png"))
src <- readLines(file.path(base, "scripts", "fig4_tgm2_build_v2.R"), warn = FALSE)
src <- gsub("Death-competent", "Apoptosis-associated", gsub("death-competent", "apoptosis-associated", src))
endb <- grep("Fig3B_TGM2_ITGA4_rank_scores\\.png", src)[1]
okb <- tryCatch({ eval(parse(text = src[1:endb]), envir = new.env()); TRUE },
                error = function(err) { cat("  NOTE: panel 3b relabel failed (", conditionMessage(err), "); previous 3b kept\n"); FALSE })
if (!okb) file.copy(file.path(ad, basename(f3b)), f3b, overwrite = TRUE)
mb3 <- image_read(f3b); htop <- hin(mb3, cw)
lv <- c("Apoptosis-\nassociated", "Arrest-associated\n(cytostatic-prone)")
long <- function(d, coh) rbind(data.frame(cohort = coh, unit = d$unit, side = lv[1], val = d$apop),
                               data.frame(cohort = coh, unit = d$unit, side = lv[2], val = d$arrest))
c1 <- sprintf("Hwang et al.\n(%d patients)", nrow(hw)); c2 <- sprintf("Steele et al.\n(%d patients)", nrow(st))
d3 <- rbind(long(hw, c1), long(st, c2)); d3$side <- factor(d3$side, levels = lv); d3$cohort <- factor(d3$cohort, levels = c(c1, c2))
ann <- data.frame(cohort = factor(c(c1, c2), levels = c(c1, c2)),
                  lab = c(sprintf("higher in %d of %d\np = %.3f", sum(hw$arrest > hw$apop), nrow(hw), hw_p),
                          sprintf("higher in %d of %d\np = %.3f", sum(st$arrest > st$apop), nrow(st), st_p)))
p3a <- ggplot(d3, aes(side, val)) + geom_boxplot(outlier.shape = NA, width = .5, colour = "grey45", fill = NA) +
  geom_line(aes(group = unit), colour = "grey78", linewidth = .35) + geom_point(aes(colour = side), size = 1.3) +
  geom_text(data = ann, aes(x = 1.5, y = Inf, label = lab), vjust = 1.3, size = 2.7, inherit.aes = FALSE) +
  facet_wrap(~cohort, scales = "free_y") + scale_colour_manual(values = setNames(c(COL_APOP, COL_ARREST), lv), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(.03, .22))) +
  labs(title = "TGM2 by side of the contrast",
       subtitle = "Per-patient mean, cells with Axis 2 signal; lines join each patient", x = NULL, y = "TGM2, mean log-normalised expression") +
  thm(9) + theme(strip.background = element_rect(colour = "grey60"), axis.text.x = element_text(size = 7))
f3a <- file.path(fig, "aim3", "Fig3A_TGM2_paired_expression.png"); sv(p3a, f3a, cw, htop)
mc3 <- image_read(file.path(fig, "aim3", "Fig3C_DepMap_essentiality.png")); yc <- mt + htop + gy
wc <- min(wf, (H - yc - mb) * image_info(mc3)$width / image_info(mc3)$height)
cv <- image_blank(px(W), px(H), color = "white")
cv <- put(cv, image_read(f3a), ml, mt, cw); cv <- put(cv, mb3, xr, mt, cw); cv <- put(cv, mc3, ml, yc, wc)
page(cv, L(c("a", "b", "c"), c(.03, xr - .19, .03), c(mt, mt, yc)), "Figure3.pdf")

## ---- Supplementary Figure 1: panel b relabelled, same pixel size as before
fS1b <- file.path(fig, "aim1", "FigS2B_DepMap_module_crossval.png"); dS1 <- pin(fS1b)
s1 <- s1[order(s1$r), ]; s1$label <- factor(s1$label, levels = s1$label)
s1$txt <- sprintf("r = %.3f  p %s", s1$r, ifelse(s1$p < .001, "< 0.001", sprintf("= %.3f", s1$p)))
pS1b <- ggplot(s1, aes(r, label, fill = r < 0)) + geom_col(width = .7) + geom_vline(xintercept = 0, colour = "grey30", linewidth = .3) +
  geom_text(aes(label = txt, hjust = ifelse(r < 0, 1.05, -.05)), size = 2.4) +
  scale_fill_manual(values = c(`TRUE` = COL_APOP, `FALSE` = "#E69F00"), guide = "none") +
  scale_x_continuous(limits = c(-1.4, .75), breaks = seq(-1, .5, .25)) + coord_cartesian(clip = "off") +
  labs(title = "Module-level KRAS dependency correlations in DepMap PDAC cell lines",
       subtitle = "Spearman r between module mean expression and KRAS Chronos score,\nn = 44 PDAC cell lines; negative r = lower expression in KRAS-dependent lines",
       x = "Spearman r vs KRAS Chronos dependency", y = NULL) + thm(9)
sv(pS1b, fS1b, dS1[["w"]], dS1[["h"]])
ia <- image_read(file.path(fig, "aim1", "FigS2A_KRAS_gene_correlations.png")); ib <- image_read(fS1b)
hb1 <- hin(ib, wf); ha1 <- min(hin(ia, wf), H - mt - gy - hb1 - mb); wa1 <- ha1 * image_info(ia)$width / image_info(ia)$height
cv <- image_blank(px(W), px(H), color = "white"); cv <- put(cv, ia, ml, mt, wa1); cv <- put(cv, ib, ml, mt + ha1 + gy, wf)
page(cv, L(c("a", "b"), c(.03, .03), c(mt, mt + ha1 + gy)), "FigureS1.pdf")

## ---- Supplementary Figure 5: panel a coloured by side of the contrast, legible points
u$side <- ifelse(u$axis2_death_competence == 0, "No Axis 2 signal",
          ifelse(u$axis2_death_competence < 0, "Arrest-associated (cytostatic-prone)", "Apoptosis-associated"))
nn <- table(u$side); lab5 <- setNames(sprintf("%s (n = %s)", names(nn), formatC(as.vector(nn), big.mark = ",", format = "d")), names(nn))
set.seed(1)
u5 <- rbind(u[u$side == "No Axis 2 signal", ], u[u$side != "No Axis 2 signal", ][sample(sum(u$side != "No Axis 2 signal")), ])
ha5 <- 4.6
p5a <- ggplot(u5, aes(.data[[uc[1]]], .data[[uc[2]]], colour = side)) + geom_point(size = .4, alpha = .7, stroke = 0) +
  scale_colour_manual(values = c("No Axis 2 signal" = "grey80", "Arrest-associated (cytostatic-prone)" = COL_ARREST,
                                 "Apoptosis-associated" = COL_APOP), labels = lab5, name = "Side of the contrast") +
  guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
  labs(title = "Malignant cells by side of the apoptosis/arrest contrast",
       subtitle = sprintf("GSE202051, untreated patients (n = %s cells, 18 patients); %s cells (%.1f%%) without Axis 2 signal in grey",
                          format(nrow(u), big.mark = ","), format(sum(u$axis2_death_competence == 0), big.mark = ","),
                          100 * mean(u$axis2_death_competence == 0)), x = "UMAP 1", y = "UMAP 2") + thm(9)
f5a <- file.path(fig, "aim2", "FigS5A_UMAP_contrast_side.png"); sv(p5a, f5a, wf, ha5)
axs <- lapply(1:4, function(i) image_read(file.path(fig, "aim2", paste0("Fig2A_UMAP_axis", i, ".png"))))
hcell <- (H - mt - ha5 - 2 * gy - mb) / 2
wcell <- min(cw, hcell * image_info(axs[[1]])$width / image_info(axs[[1]])$height); hcell <- hin(axs[[1]], wcell)
y5r1 <- mt + ha5 + gy; y5r2 <- y5r1 + hcell + gy
cv <- image_blank(px(W), px(H), color = "white"); cv <- put(cv, image_read(f5a), ml, mt, wf)
cv <- put(cv, axs[[1]], ml, y5r1, wcell); cv <- put(cv, axs[[2]], xr, y5r1, wcell)
cv <- put(cv, axs[[3]], ml, y5r2, wcell); cv <- put(cv, axs[[4]], xr, y5r2, wcell)
page(cv, L(letters[1:5], c(.03, .03, xr - .19, .03, xr - .19), c(mt, y5r1, y5r1, y5r2, y5r2)), "FigureS5.pdf")

## ---- Supplementary Figure 6: threshold labels corrected ("zero crossing", not "pre-specified")
s6 <- do.call(rbind, lapply(pid, function(q) { a <- ad2[x$pid == q]
  data.frame(pid = q, strategy = c("Zero crossing (Axis 2 < 0)", "Median of informative cells", "33rd percentile of informative cells"),
             pct = 100 * c(mean(a < 0), mean(a < tmed), mean(a < t33))) }))
s6$strategy <- factor(s6$strategy, levels = unique(s6$strategy))
ord <- 100 * tapply(ad2 < 0, x$pid, mean)[pid]; s6$pid <- factor(s6$pid, levels = pid[order(ord)])
ha6 <- 5.4; hb6 <- H - mt - gy - ha6 - mb
p6a <- ggplot(s6, aes(pid, pct, fill = strategy)) + geom_col(position = position_dodge(.8), width = .75) +
  geom_hline(yintercept = 50, linetype = "dashed", colour = "grey40") +
  scale_fill_manual(values = c("#0072B2", "#E69F00", "#D55E00"), name = "Threshold") +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  labs(title = "Per-patient cytostatic-prone proportion under three thresholds",
       subtitle = sprintf("Untreated patients, cells with Axis 2 signal; median %.2f, 33rd percentile %.2f", tmed, t33),
       x = "Patient", y = "Cytostatic-prone cells (%)") + thm(10) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
f6a <- file.path(fig, "aim2", "FigS5_classifier_sensitivity.png"); sv(p6a, f6a, wf, ha6)
sw <- data.frame(t = seq(-10, 10, .1)); sw$pct <- sapply(sw$t, function(z) 100 * mean(ad2 < z))
p6b <- ggplot(sw, aes(t, pct)) + geom_line(colour = "#0072B2", linewidth = .8) +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "#0072B2") + geom_hline(yintercept = 50, linetype = "dashed", colour = "grey40") +
  annotate("text", x = .3, y = 25, hjust = 0, size = 3.2, colour = "#0072B2", label = sprintf("Zero crossing (0)\n%.1f%% cytostatic-prone", at0)) +
  annotate("text", x = -9.8, y = 53, hjust = 0, size = 3.2, colour = "grey40", label = "50%") + scale_y_continuous(limits = c(0, 100)) +
  labs(title = "Cytostatic-prone proportion across thresholds", subtitle = sprintf("Untreated cells with Axis 2 signal (n = %s)", format(nrow(x), big.mark = ",")),
       x = "Axis 2 apoptosis/arrest contrast threshold", y = "Cells classified as cytostatic-prone (%)") + thm(10)
f6b <- file.path(fig, "aim2", "FigS5b_classifier_sweep.png"); sv(p6b, f6b, wf, hb6)
cv <- image_blank(px(W), px(H), color = "white"); cv <- put(cv, image_read(f6a), ml, mt, wf); cv <- put(cv, image_read(f6b), ml, mt + ha6 + gy, wf)
page(cv, L(c("a", "b"), c(.03, .03), c(mt, mt + ha6 + gy)), "FigureS6.pdf")

## ---- Supplementary Figure 7: map header redrawn with the correct data credit
bk(maps_f)
body <- image_crop(mi, sprintf("%dx%d+0+%d", info$width, info$height - cut_px, cut_px))
hdr <- ggplot() + theme_void() + thm(9) + theme(axis.line = element_blank(), plot.margin = margin(2, 2, 2, 2)) +
  labs(title = "Spatial distribution of Treg activation score",
       subtitle = "Primary PDAC, n = 35,458 spots across 10 tumours (GSE272362; Khaliq et al. 2024)")
hh <- .42; fhd <- file.path(fig, "aim3", "FigS7A_header.png"); sv(hdr, fhd, wf, hh)
mapsA <- image_append(c(image_scale(image_read(fhd), paste0(info$width, "x")), body), stack = TRUE)
fS7a <- file.path(fig, "aim3", "FigS7A_spatial_Treg_maps_credit.png"); image_write(mapsA, fS7a, format = "png")
ms <- image_read(spt_f); hA7 <- hin(mapsA, wf); hB7 <- hin(ms, wf); sc <- min(1, (H - mt - gy - mb) / (hA7 + hB7))
cv <- image_blank(px(W), px(H), color = "white"); cv <- put(cv, mapsA, ml, mt, wf * sc); cv <- put(cv, ms, ml, mt + hA7 * sc + gy, wf * sc)
page(cv, L(c("a", "b"), c(.03, .03), c(mt, mt + hA7 * sc + gy)), "FigureS7.pdf")

cat("\nDone. Rebuilt: Figure1, Figure2, Figure3, FigureS1, FigureS5, FigureS6, FigureS7 (figures/multipanel/).\n")
cat("Replaced files are backed up in figures/_archive/figfix_v1/.\n")
