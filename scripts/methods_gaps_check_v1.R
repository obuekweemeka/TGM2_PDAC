# ============================================================================
# methods_gaps_check_v1.R -- READ-ONLY. Save in: scripts/
# Two Methods gaps flagged in review that must be answered from the data,
# not from memory:
#  1. How was the Axis 2 contrast computed for the 44 DepMap cell lines?
#  2. Were the eight bypass-potential genes significant in one resistant line
#     or in both?
# ============================================================================
base <- normalizePath('.'); pr <- file.path(base,"data","processed")
cat("=== 1. Cell-line Axis 2 score ===\n")
d4 <- readRDS(file.path(pr,"aim4_bulk","depmap_pdac_analysis_df.rds"))
cat("columns:", paste(names(d4), collapse=", "), "\n")
cat("axis2_death  range:", paste(round(range(d4$axis2_death),3), collapse=" to "), "\n")
cat("axis2_ssgsea range:", paste(round(range(d4$axis2_ssgsea),3), collapse=" to "), "\n")
cat("Spearman between them:", round(cor(d4$axis2_death, d4$axis2_ssgsea, method="spearman"),3), "\n")
de <- readRDS(file.path(pr,"aim1_signatures","DepMap_PDAC_expression.rds"))
cat("\nexpression object: ", nrow(de), "rows x", ncol(de), "cols; values look like",
    ifelse(max(de[[3]], na.rm=TRUE) < 20, "log2(TPM+1)", "linear TPM"), "\n")
sig <- readRDS(file.path(pr,"aim1_signatures","four_axis_signatures_v1.rds"))
comp <- c("axis2_proapoptotic","axis2_executors","axis2_anti_apoptotic","axis2_cytostasis")
e <- de[!duplicated(de$ModelID), ]; rn <- e$ModelID; e$ModelID <- NULL
gm <- function(g){ g <- intersect(g, names(e)); rowMeans(e[, g, drop=FALSE], na.rm=TRUE) }
P <- gm(sig[[comp[1]]]) + gm(sig[[comp[2]]]); A <- gm(sig[[comp[3]]]) + gm(sig[[comp[4]]])
cand <- data.frame(ModelID = rn,
                   mean_logratio = log2((P + 0.001)/(A + 0.001)),
                   mean_difference = P - A,
                   z_difference = as.numeric(scale(P) - scale(A)))
m <- merge(d4[, c("model_id","axis2_death","axis2_ssgsea")], cand, by.x="model_id", by.y="ModelID")
cat("\nCorrelation of stored scores with candidate reconstructions (n =", nrow(m), "):\n")
for (col in c("mean_logratio","mean_difference","z_difference"))
  cat(sprintf("  axis2_death vs %-16s Spearman %.3f | Pearson %.3f | identical: %s\n", col,
      cor(m$axis2_death, m[[col]], method="spearman"), cor(m$axis2_death, m[[col]]),
      isTRUE(all.equal(m$axis2_death, m[[col]], tolerance=1e-6))))
cat("\nIf none is an exact match, the cell-line score was produced by a method not\nreconstructable from saved objects, and the Methods must say exactly that.\n")

cat("\n=== 2. Bypass-potential gene evidence ===\n")
suppressPackageStartupMessages(library(DESeq2))
dds <- readRDS(file.path(pr,"aim1_signatures","dds_fitted.rds"))
r0 <- as.data.frame(results(dds, name="condition_resistant_vs_parental")); r0$ensembl_id <- rownames(r0)
r1 <- as.data.frame(results(dds, contrast=list(c("condition_resistant_vs_parental","cell_linePANC1.conditionresistant")))); r1$ensembl_id <- rownames(r1)
rsa <- readRDS(file.path(pr,"aim1_signatures","res_shared_annot.rds")); map <- unique(rsa[,c("ensembl_id","hgnc_symbol")])
g <- sig$axis1_bypass_potential
tb <- merge(map[map$hgnc_symbol %in% g, ], r0[,c("ensembl_id","log2FoldChange","padj")], by="ensembl_id")
tb <- merge(tb, r1[,c("ensembl_id","log2FoldChange","padj")], by="ensembl_id", suffixes=c("_Panc0203","_PANC1"))
tb$both_sig_same_dir <- with(tb, padj_Panc0203 < .05 & padj_PANC1 < .05 & sign(log2FoldChange_Panc0203)==sign(log2FoldChange_PANC1))
print(tb[order(match(tb$hgnc_symbol,g)), c("hgnc_symbol","log2FoldChange_Panc0203","padj_Panc0203","log2FoldChange_PANC1","padj_PANC1","both_sig_same_dir")], row.names=FALSE, digits=3)
cat(sprintf("\nSignificant and concordant in BOTH lines: %d of %d\n", sum(tb$both_sig_same_dir), nrow(tb)))
