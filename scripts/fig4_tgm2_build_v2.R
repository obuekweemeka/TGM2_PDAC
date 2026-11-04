# ============================================================================
# fig4_tgm2_build_v2.R -- Figure 4: TGM2 as a candidate mechanism. Save in: scripts/
# v2: the single-cell matrix is cells x genes and ALREADY log-normalised; panel a now
#     uses the stored values directly (v1 wrongly re-normalised them). Writes Figure4.pdf.
#  a: per-patient TGM2 expression, cytostatic-prone vs death-competent cells (paired)
#  b: TGM2-ITGA4 LIANA aggregate rank scores by sender state, all immune receivers
#  c: DepMap essentiality block (TGM2, NT5E, KRAS), 44 lines, 22/22 split
# Panels -> figures/aim3/*.png ; multipanel -> figures/multipanel/Figure4.pdf (PDF only).
# Superseded old Fig 3 panels are moved to _archive/. Ends with read-only checks
# on two open ledger items (GSE202051 sample count; Axis 1 formula epsilon).
# ============================================================================
suppressPackageStartupMessages({library(ggplot2); library(patchwork); library(magick); library(Matrix)})
base <- normalizePath('.'); fig <- file.path(base, "figures"); pr <- file.path(base, "data", "processed")
W <- 8.27; H <- 11.69; dpi <- 600; ml <- 0.22; mr <- 0.06; mt <- 0.06; mb <- 0.06; gx <- 0.22; gy <- 0.12; fs <- 6.5
col_cy <- "#009E73"; col_dc <- "#0072B2"
arch <- function(p) { a <- file.path(fig, "_archive", paste0(tools::file_path_sans_ext(basename(p)), "_superseded.", tools::file_ext(p)))
  if (file.exists(p) && !file.exists(a)) file.copy(p, a) }
savep <- function(p, path, w, h) { arch(path); ggsave(path, p, dpi = dpi, width = w, height = h, units = "in"); cat("Saved", basename(path), round(file.size(path) / 1024), "KB\n") }
px  <- function(x) as.integer(round(x * dpi))
put <- function(cv, img, x, y, w) image_composite(cv, image_scale(img, paste0(px(w), "x")), offset = paste0("+", px(x), "+", px(y)))
thm <- function(fs) theme_classic(base_size = fs) +
  theme(plot.title = element_text(face = "bold", size = fs + 1), plot.subtitle = element_text(size = fs - 1, colour = "grey40"))
cw <- (W - ml - mr - gx) / 2; h_top <- 4.2; w_c <- W - ml - mr; h_c <- H - mt - h_top - gy - mb

# ---------------- panel a: TGM2 expression, paired by patient ----------------
v4  <- readRDS(file.path(pr, "aim2_scoring", "malignant_cells_scored_v4.rds"))
mcm <- readRDS(file.path(pr, "aim2_scoring", "malignant_cell_metadata.rds"))
tr  <- unique(mcm[, c("pid", "new_treatment")]); v4$trt <- tr$new_treatment[match(v4$pid, tr$pid)]
x <- v4[v4$trt == "Untreated" & v4$axis2_death_competence != 0, ]; x$cyto <- grepl("cytostatic", x$cell_state)
m <- readRDS(file.path(pr, "aim2_scoring", "malignant_count_matrix.rds"))   # cells x genes, stored log-normalised values
if ("TGM2" %in% colnames(m)) { idx <- match(x$cell_id, rownames(m)); ok <- !is.na(idx); x <- x[ok, ]; x$tgm2 <- as.numeric(m[idx[ok], "TGM2"])
} else { idx <- match(x$cell_id, colnames(m)); ok <- !is.na(idx); x <- x[ok, ]; x$tgm2 <- as.numeric(m["TGM2", idx[ok]]) }
rm(m); invisible(gc())
cat(sprintf("TGM2 detected in %.1f%% of cytostatic-prone vs %.1f%% of death-competent cells | mean stored log-normalised %.3f vs %.3f\n",
            100 * mean(x$tgm2[x$cyto] > 0), 100 * mean(x$tgm2[!x$cyto] > 0), mean(x$tgm2[x$cyto]), mean(x$tgm2[!x$cyto])))
pc <- tapply(x$tgm2[x$cyto], x$pid[x$cyto], mean); pdc <- tapply(x$tgm2[!x$cyto], x$pid[!x$cyto], mean); pids <- intersect(names(pc), names(pdc))
w <- suppressWarnings(wilcox.test(pc[pids], pdc[pids], paired = TRUE))
cat(sprintf("Panel a: n = %d patients | cytostatic higher in %d | paired Wilcoxon p = %.4f (v1, mis-normalised: 12 of 18, p = 0.0040)\n",
            length(pids), sum(pc[pids] > pdc[pids]), w$p.value))
pa_d <- data.frame(pid = rep(pids, 2), state = rep(c("Death-competent", "Cytostatic-prone"), each = length(pids)), value = c(pdc[pids], pc[pids]))
pa_d$state <- factor(pa_d$state, levels = c("Death-competent", "Cytostatic-prone"))
p_a <- ggplot(pa_d, aes(state, value)) +
  geom_boxplot(aes(group = state), width = 0.4, outlier.shape = NA, fill = NA, colour = "grey55", linewidth = 0.3) +
  geom_line(aes(group = pid), colour = "grey75", linewidth = 0.3) + geom_point(aes(colour = state), size = 1.4) +
  scale_colour_manual(values = c("Death-competent" = col_dc, "Cytostatic-prone" = col_cy), guide = "none") +
  annotate("text", x = 1.5, y = max(pa_d$value) * 1.06, size = fs / ggplot2::.pt,
           label = sprintf("Higher in cytostatic-prone cells in %d of %d patients\npaired Wilcoxon signed-rank p = %.3f", sum(pc[pids] > pdc[pids]), length(pids), w$p.value)) +
  labs(title = "TGM2 expression by malignant cell state", subtitle = "Per-patient mean, untreated cells with Axis 2 signal; lines join each patient",
       x = NULL, y = "TGM2, mean log-normalised expression (stored values)") + thm(fs)
f_a <- file.path(fig, "aim3", "Fig3A_TGM2_paired_expression.png"); savep(p_a, f_a, cw, h_top)

# ---------------- panel b: TGM2-ITGA4 rank scores by sender state ----------------
ms <- read.csv(file.path(base, "results", "aim3", "liana_malignant_sender_interactions.csv"), stringsAsFactors = FALSE)
t2 <- ms[ms$ligand.complex == "TGM2" & ms$receptor.complex == "ITGA4", ]
t2$grp <- ifelse(grepl("cytostatic", t2$source), "Cytostatic-prone senders", ifelse(grepl("death", t2$source), "Death-competent senders", NA))
t2 <- t2[!is.na(t2$grp), ]; b <- aggregate(aggregate_rank ~ target + grp, t2, min); b$score <- -log10(b$aggregate_rank)
clean <- function(s) { s <- sub("^Immune-", "", s); s <- gsub("CD8pos_Tcells", "CD8+ T cells", s); gsub("_", " ", s) }
b$recv <- clean(b$target); cyo <- b[b$grp == "Cytostatic-prone senders", ]; b$recv <- factor(b$recv, levels = cyo$recv[order(cyo$score)])
rr <- merge(b[b$grp == "Death-competent senders", c("recv", "aggregate_rank")], b[b$grp == "Cytostatic-prone senders", c("recv", "aggregate_rank")], by = "recv")
cat(sprintf("Panel b: %d receivers | rank-score ratio range %.2f-%.2f (ledger: 1.49-14.59)\n",
            nrow(rr), min(rr[, 2] / rr[, 3]), max(rr[, 2] / rr[, 3])))
p_b <- ggplot(b, aes(score, recv)) + geom_line(aes(group = recv), colour = "grey70", linewidth = 0.4) +
  geom_point(aes(colour = grp), size = 1.6) +
  scale_colour_manual(values = c("Cytostatic-prone senders" = col_cy, "Death-competent senders" = col_dc), name = NULL) +
  labs(title = "TGM2-ITGA4 ranking by sender state", subtitle = "LIANA consensus aggregate rank; 16,090-cell untreated subsample",
       x = "-log10(aggregate rank score); higher = ranked higher", y = "Immune receiver") +
  thm(fs) + theme(legend.position = "bottom")
f_b <- file.path(fig, "aim3", "Fig3B_TGM2_ITGA4_rank_scores.png"); savep(p_b, f_b, cw, h_top)

# ---------------- panel c: DepMap essentiality block ----------------
dep_expr <- readRDS(file.path(pr, "aim1_signatures", "DepMap_PDAC_expression.rds"))
dep4     <- readRDS(file.path(pr, "aim4_bulk", "depmap_pdac_analysis_df.rds"))
valid44  <- intersect(dep_expr$ModelID[!duplicated(dep_expr$ModelID)], dep4$model_id)
d <- dep4[dep4$model_id %in% valid44, ]; stopifnot(nrow(d) == 44)
d$grp <- ifelse(d$axis2_death < median(d$axis2_death), "High cytostatic", "Low cytostatic"); hi <- d$grp == "High cytostatic"
pct <- sapply(c(tgm2 = "tgm2", nt5e = "nt5e", kras = "kras"), function(g) round(100 * mean(d[[g]] > 0), 1))
leg <- c("High cytostatic" = paste0("High cytostatic (n = ", sum(hi), ")"), "Low cytostatic" = paste0("Low cytostatic (n = ", sum(!hi), ")"))
pnl <- function(g, ttl) ggplot(d, aes(axis2_death, .data[[g]], colour = grp)) +
  geom_point(size = 1.3, alpha = 0.85) + geom_smooth(method = "lm", se = TRUE, colour = "grey30", linewidth = 0.4, formula = y ~ x) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  scale_colour_manual(values = c("High cytostatic" = "#E69F00", "Low cytostatic" = "#0072B2"), labels = leg, name = NULL) +
  labs(title = ttl, subtitle = paste0("n = 44 PDAC lines; Chronos > 0: ", pct[[g]], "%"), x = "Axis 2 death competence score", y = "Chronos effect score") +
  thm(fs) + theme(legend.text = element_text(size = fs))
blk <- ((pnl("tgm2", "TGM2") | pnl("nt5e", "NT5E")) / (pnl("kras", "KRAS (positive control)") | guide_area())) + plot_layout(guides = "collect")
f_c <- file.path(fig, "aim3", "Fig3C_DepMap_essentiality.png"); savep(blk, f_c, w_c, h_c)

# ---------------- archive superseded old Fig 3 panels; compose page ----------------
for (old in c("Fig3A_LIANA_dotplot.png", "Fig3B_fold_enrichment_bar.png")) {
  p <- file.path(fig, "aim3", old); if (file.exists(p)) { file.rename(p, file.path(fig, "_archive", sub("\\.png$", "_superseded.png", old))); cat("Archived", old, "\n") } }
cv <- image_blank(px(W), px(H), color = "white"); xr <- ml + cw + gx; yc <- mt + h_top + gy
cv <- put(cv, image_read(f_a), ml, mt, cw); cv <- put(cv, image_read(f_b), xr, mt, cw); cv <- put(cv, image_read(f_c), ml, yc, w_c)
labs_df <- data.frame(l = c("a", "b", "c"), x = c(0.03, xr - 0.19, 0.03), y = c(mt, mt, yc))
pg <- ggplot() + annotation_raster(as.raster(cv), xmin = 0, xmax = W, ymin = 0, ymax = H, interpolate = TRUE) +
  geom_text(data = labs_df, aes(x = x, y = H - y, label = l), hjust = 0, vjust = 1, fontface = "bold", size = 9 / ggplot2::.pt) +
  coord_fixed(xlim = c(0, W), ylim = c(0, H), expand = FALSE) + theme_void() + theme(plot.margin = margin(0, 0, 0, 0))
out <- file.path(fig, "multipanel", "Figure4.pdf"); tmp <- file.path(fig, "multipanel", "Figure4_new.pdf"); arch(out)
ggsave(tmp, pg, width = W, height = H, units = "in", device = "pdf")
sz <- round(file.size(tmp) / 1024); stopifnot(sz > 100); file.copy(tmp, out, overwrite = TRUE); file.remove(tmp); cat("Saved Figure4.pdf", sz, "KB\n")

cat("\nDone.\n")
