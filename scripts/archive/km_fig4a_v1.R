# ============================================================================
# km_fig4a_v1.R -- Figure 4a (TCGA PAAD Kaplan-Meier) regenerated from data.
# Save in: scripts/. Corrects the panel title (groups are an Axis 2 median split,
# verified 146/146) and re-derives every number shown in the panel/legend.
# Sized to its exact slot in Figure 4 (3.885 in square) with 6.3 pt base font,
# so its text matches panels b and c. The old PNG is archived before replacing.
# ============================================================================
suppressPackageStartupMessages({library(ggplot2); library(patchwork); library(survival)})
if (!requireNamespace("survminer", quietly = TRUE)) install.packages("survminer")
suppressPackageStartupMessages(library(survminer))
base <- here::here(); fig <- file.path(base, "figures")
fs <- 6.3; side <- 3.885

tcga <- readRDS(file.path(base, "data", "processed", "aim4_bulk", "tcga_paad_analysis_df.rds"))
tcga$grp <- factor(tcga$cytostatic_group, levels = c("High cytostatic", "Low cytostatic"))

# ---------- numbers (compare with the manuscript) ----------
fit <- survfit(Surv(os_days, os_event) ~ grp, data = tcga)
print(summary(fit)$table[, c("records", "events", "median")])
cx  <- coxph(Surv(os_days, os_event) ~ grp, data = tcga)
hr  <- exp(coef(cx)); ci <- exp(confint(cx))
lr  <- survdiff(Surv(os_days, os_event) ~ grp, data = tcga); p_lr <- 1 - pchisq(lr$chisq, df = 1)
cat(sprintf("HR (Low vs High) = %.3f (95%% CI %.3f-%.3f); log-rank p = %.3f\n", hr, ci[1], ci[2], p_lr))
cat("Manuscript: HR 1.386 (0.874-2.197), log-rank p = 0.165, median OS 691 (High) vs 511 (Low) days, n = 73/73\n")

# ---------- plot ----------
g <- ggsurvplot(fit, data = tcga, conf.int = TRUE, risk.table = TRUE, pval = FALSE,
                palette = c("#D55E00", "#0072B2"), legend.title = "Cell state group",
                legend.labs = c("High cytostatic", "Low cytostatic"),
                xlab = "Overall survival (days)", ylab = "Survival probability",
                title = "OS by Axis 2 death competence",
                ggtheme = theme_classic(base_size = fs), tables.theme = theme_cleantable(base_size = fs),
                fontsize = fs / ggplot2::.pt, censor.size = 2, size = 0.4)
pk <- g$plot +
  labs(subtitle = "TCGA PAAD, n = 146, median split on Axis 2 death competence") +
  annotate("text", x = 150, y = 0.22, hjust = 0, size = fs / ggplot2::.pt,
           label = sprintf("Log-rank p = %.3f\nHR = %.2f (95%% CI %.2f-%.2f)", p_lr, hr, ci[1], ci[2])) +
  theme(plot.title = element_text(face = "bold", size = fs + 1.2),
        plot.subtitle = element_text(size = fs - 0.8, colour = "grey40"),
        legend.text = element_text(size = fs), legend.title = element_text(size = fs))
pt <- g$table + labs(title = "Number at risk") + theme(plot.title = element_text(size = fs))
km <- pk / pt + plot_layout(heights = c(3, 1))

# ---------- archive old PNG, then save ----------
out  <- file.path(fig, "aim4", "Fig4A_KM_cytostatic.png")
arch <- file.path(fig, "_archive", "Fig4A_KM_cytostatic_superseded.png")
if (file.exists(out) && !file.exists(arch)) { file.copy(out, arch); cat("Old KM panel archived\n") }
ggsave(out, km, dpi = 600, width = side, height = side, units = "in")
cat("Saved:", out, "-", round(file.size(out) / 1024), "KB\n")
