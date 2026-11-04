# ============================================================================
# module_provenance_check_v1.R -- READ-ONLY. Save in: scripts/
# Verifies how every stored four-axis module was built, so the Methods and
# Results can describe the modules exactly. Nothing is written to disk.
#  1. Axis 2 curated genes: direction and significance in BOTH resistance lines
#  2. Axis 3 / Axis 4 modules vs the saved Hallmark GSEA leading edges
#  3. DepMap correlates module vs the 44-line gene-level correlations
# ============================================================================
suppressPackageStartupMessages(library(DESeq2))
base <- normalizePath('.'); ap <- file.path(base, "data", "processed", "aim1_signatures")
sig <- readRDS(file.path(ap, "four_axis_signatures_v1.rds"))
rsa <- readRDS(file.path(ap, "res_shared_annot.rds")); map <- unique(rsa[, c("ensembl_id", "hgnc_symbol")])

# ---------- 1. Axis 2 genes in both line-specific contrasts ----------
dds <- readRDS(file.path(ap, "dds_fitted.rds"))
r0 <- as.data.frame(results(dds, name = "condition_resistant_vs_parental")); r0$ensembl_id <- rownames(r0)
r1 <- as.data.frame(results(dds, contrast = list(c("condition_resistant_vs_parental", "cell_linePANC1.conditionresistant")))); r1$ensembl_id <- rownames(r1)
genes <- c("BCL2", "MCL1", "BCL2L1", "BCL2L2", "BCL2A1", "BIRC5", "BBC3", "PMAIP1", "BCL2L11", "BID", "BAD", "BMF", "BIK",
           "CASP3", "CASP7", "CASP8", "CASP9", "CYCS", "DIABLO", "APAF1", "CDKN1A", "CDKN2A", "RB1", "GADD45A", "GADD45B")
t2 <- merge(map[map$hgnc_symbol %in% genes, ], r0[, c("ensembl_id", "log2FoldChange", "padj")], by = "ensembl_id")
t2 <- merge(t2, r1[, c("ensembl_id", "log2FoldChange", "padj")], by = "ensembl_id", suffixes = c("_Panc0203", "_PANC1"))
t2 <- t2[order(match(t2$hgnc_symbol, genes)), ]
out2 <- data.frame(gene = t2$hgnc_symbol,
                   LFC_Panc0203 = round(t2$log2FoldChange_Panc0203, 2), padj_Panc0203 = signif(t2$padj_Panc0203, 2),
                   LFC_PANC1 = round(t2$log2FoldChange_PANC1, 2), padj_PANC1 = signif(t2$padj_PANC1, 2))
out2$both_sig_same_direction <- with(t2, padj_Panc0203 < 0.05 & padj_PANC1 < 0.05 & sign(log2FoldChange_Panc0203) == sign(log2FoldChange_PANC1))
cat("=== 1. Axis 2 curated genes (resistant vs parental, per line) ===\n"); print(out2, row.names = FALSE)

# ---------- 2. Axis 3 / 4 modules vs GSEA leading edges ----------
gs <- as.data.frame(readRDS(file.path(ap, "gsea_hallmark_res.rds")))
cat("\n=== 2. Significant Hallmark sets (padj < 0.05) ===\n")
s <- gs[gs$padj < 0.05, c("pathway", "NES", "padj", "size")]; s$NES <- round(s$NES, 2); s$padj <- signif(s$padj, 2)
print(s[order(-s$NES), ], row.names = FALSE)
le <- setNames(gs$leadingEdge, gs$pathway)
tosym <- function(v) { v <- unlist(v); if (mean(grepl("^ENSG", v)) > 0.5) map$hgnc_symbol[match(v, map$ensembl_id)] else v }
chk <- function(label, genes_le, stored) cat(sprintf("%-34s leading edge %3d | stored %3d | overlap %3d | identical: %s\n",
                                                    label, length(genes_le), length(stored), length(intersect(genes_le, stored)), setequal(genes_le, stored)))
chk("Hypoxia -> axis3_hypoxia_metabolism", tosym(le[["HALLMARK_HYPOXIA"]]), sig$axis3_hypoxia_metabolism)
chk("OXPHOS -> axis3_OXPHOS_depleted", tosym(le[["HALLMARK_OXIDATIVE_PHOSPHORYLATION"]]), sig$axis3_OXPHOS_depleted)
imm <- c("HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
         "HALLMARK_IL6_JAK_STAT3_SIGNALING", "HALLMARK_INFLAMMATORY_RESPONSE")
imm_le <- lapply(imm[imm %in% names(le)], function(p) tosym(le[[p]]))
cat("Immune sets found:", sum(imm %in% names(le)), "of 5\n")
cnt <- table(unlist(lapply(imm_le, unique)))
chk("Immune union -> innate_mimicry_full", names(cnt), sig$axis4_innate_mimicry_full)
chk("Immune in >=2 sets -> innate_mimicry_core", names(cnt)[cnt >= 2], sig$axis4_innate_mimicry_core)

# ---------- 3. DepMap correlates vs 44-line correlations ----------
dep_expr <- readRDS(file.path(ap, "DepMap_PDAC_expression.rds"))
dep4 <- readRDS(file.path(base, "data", "processed", "aim4_bulk", "depmap_pdac_analysis_df.rds"))
dep_expr <- dep_expr[!duplicated(dep_expr$ModelID), ]; valid44 <- intersect(dep_expr$ModelID, dep4$model_id)
e <- as.data.frame(dep_expr[match(valid44, dep_expr$ModelID), ]); k <- dep4$kras[match(valid44, dep4$model_id)]
num <- names(e)[sapply(e, is.numeric)]; X <- as.matrix(e[, num]); X <- X[, apply(X, 2, sd, na.rm = TRUE) > 0, drop = FALSE]
r <- as.numeric(cor(X, k, use = "pairwise.complete.obs")); names(r) <- colnames(X); r <- r[!is.na(r)]
top30 <- names(sort(r))[1:30]
cat("\n=== 3. DepMap correlates ===\n")
cat("Genes with r < -0.4:", sum(r < -0.4), "| stored module:", length(sig$axis1_depmap_correlates),
    "| stored genes among the 30 most negative:", length(intersect(top30, sig$axis1_depmap_correlates)), "\n")
cat("Stored module r range:", paste(round(range(r[sig$axis1_depmap_correlates]), 3), collapse = " to "), "\n")
