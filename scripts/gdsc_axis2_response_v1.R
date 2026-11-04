# ============================================================================
# gdsc_axis2_response_v1.R -- Save in: scripts/. Writes one CSV.
# Does the Axis 2 apoptosis/arrest contrast predict measured drug sensitivity
# in PDAC cell lines?
#
# PRESPECIFIED (fixed before any association was computed):
#   Primary drugs : MAPK-pathway inhibitors (MEK, ERK, KRAS-G12C)
#   Primary metric: LN_IC50 from GDSC2 release 8.5
#   Direction     : cytostatic-prone = LOWER Axis 2 = LESS sensitive = HIGHER
#                   ln(IC50), so the expected Spearman correlation between
#                   Axis 2 and ln(IC50) is NEGATIVE.
#   Multiplicity  : Benjamini-Hochberg across the primary drug set.
#   Caveat        : viability assays cannot separate growth arrest from death,
#                   so an arrested line can appear resistant without the drug
#                   having failed. Reported, not adjusted for.
# ============================================================================
base <- normalizePath('.'); pr <- file.path(base,"data","processed")
raw <- file.path(base,"data","raw","GDSC")
if (!requireNamespace("readxl", quietly=TRUE)) install.packages("readxl")
dir.create(file.path(base,"results","aim4"), recursive=TRUE, showWarnings=FALSE)

de <- readRDS(file.path(pr,"aim1_signatures","DepMap_PDAC_expression.rds"))
d4 <- readRDS(file.path(pr,"aim4_bulk","depmap_pdac_analysis_df.rds"))
v44 <- intersect(de$ModelID[!duplicated(de$ModelID)], d4$model_id)
d <- d4[d4$model_id %in% v44, ]; cat("PDAC lines with Axis 2:", nrow(d), "\n")

clean <- function(x) toupper(gsub("[^A-Za-z0-9]","",as.character(x)))
d$key <- clean(d$cell_line)
mod <- tryCatch(read.csv(file.path(base,"data","raw","DepMap","Model.csv"), stringsAsFactors=FALSE), error=function(e) NULL)
if (!is.null(mod)) { nmcol <- intersect(c("StrippedCellLineName","CellLineName","CCLEName"), names(mod))[1]
  cat("Model.csv name column:", nmcol, "| COSMIC column:", intersect(c("COSMICID","COSMICId"), names(mod))[1], "\n")
  mp <- mod[match(d$model_id, mod$ModelID), ]
  d$key2 <- clean(mp[[nmcol]])
  cosc <- intersect(c("COSMICID","COSMICId"), names(mod))[1]
  d$cosmic <- if (!is.na(cosc)) mp[[cosc]] else NA
} else { d$key2 <- d$key; d$cosmic <- NA }

g <- readxl::read_excel(file.path(raw,"GDSC2_fitted_dose_response_27Oct23.xlsx"))
g <- as.data.frame(g); g$key <- clean(g$CELL_LINE_NAME)
PRIMARY <- c("Trametinib","Selumetinib","Refametinib","Ulixertinib","SCH772984","KRAS (G12C) Inhibitor-12")
cat("\nMatching (name, then COSMIC):\n")
m1 <- sum(d$key %in% g$key); m2 <- sum(d$key2 %in% g$key, na.rm=TRUE)
d$match <- ifelse(d$key2 %in% g$key, d$key2, ifelse(d$key %in% g$key, d$key, NA))
if (!all(is.na(d$cosmic))) { cm <- unique(g[,c("COSMIC_ID","key")]); extra <- is.na(d$match) & d$cosmic %in% cm$COSMIC_ID
  d$match[extra] <- cm$key[match(d$cosmic[extra], cm$COSMIC_ID)] }
cat("  by dep4 cell_line:", m1, "| by Model.csv name:", m2, "| total matched:", sum(!is.na(d$match)), "of", nrow(d), "\n")
if (sum(is.na(d$match))) cat("  unmatched:", paste(head(d$cell_line[is.na(d$match)],12), collapse=", "), "\n")
dm <- d[!is.na(d$match), ]

res <- list()
for (drug in sort(unique(g$DRUG_NAME))) {
  sub <- g[g$DRUG_NAME == drug & g$key %in% dm$match, c("key","LN_IC50","AUC")]
  if (!nrow(sub)) next
  sub <- aggregate(cbind(LN_IC50, AUC) ~ key, sub, mean)          # average duplicate screens
  mm <- merge(dm[,c("match","axis2_death","axis1_ras","cell_line")], sub, by.x="match", by.y="key")
  if (nrow(mm) < 10) next
  c1 <- suppressWarnings(cor.test(mm$axis2_death, mm$LN_IC50, method="spearman"))
  c2 <- suppressWarnings(cor.test(mm$axis2_death, mm$AUC, method="spearman"))
  c3 <- suppressWarnings(cor.test(mm$axis1_ras,  mm$LN_IC50, method="spearman"))
  res[[drug]] <- data.frame(drug=drug, n=nrow(mm), primary=drug %in% PRIMARY,
    rho_axis2_lnIC50=round(c1$estimate,3), p_axis2_lnIC50=signif(c1$p.value,3),
    rho_axis2_AUC=round(c2$estimate,3), p_axis2_AUC=signif(c2$p.value,3),
    rho_axis1_lnIC50=round(c3$estimate,3), p_axis1_lnIC50=signif(c3$p.value,3))
}
out <- do.call(rbind, res)
cat("\nDrugs with at least 10 matched PDAC lines:", nrow(out), "\n")
pri <- out[out$primary, ]
if (nrow(pri)) {
  pri$padj_BH <- signif(p.adjust(pri$p_axis2_lnIC50, "BH"), 3)
  cat("\n=== PRIMARY (MAPK-pathway) drugs: prespecified direction is NEGATIVE rho ===\n")
  print(pri[order(pri$p_axis2_lnIC50), c("drug","n","rho_axis2_lnIC50","p_axis2_lnIC50","padj_BH","rho_axis2_AUC","rho_axis1_lnIC50")], row.names=FALSE)
  cat(sprintf("\nPrimary drugs with rho in the predicted (negative) direction: %d of %d\n", sum(pri$rho_axis2_lnIC50 < 0), nrow(pri)))
  cat(sprintf("Median rho across primary drugs: %.3f\n", median(pri$rho_axis2_lnIC50)))
} else cat("\nNo primary drug reached 10 matched lines.\n")
cat("\n=== Strongest associations across ALL screened drugs (context only, not prespecified) ===\n")
print(head(out[order(out$p_axis2_lnIC50), c("drug","n","rho_axis2_lnIC50","p_axis2_lnIC50")], 12), row.names=FALSE)
cat("\nAcross all", nrow(out), "drugs,", sum(out$p_axis2_lnIC50 < 0.05), "reach p < 0.05;",
    round(0.05*nrow(out),1), "expected by chance.\n")
write.csv(out, file.path(base,"results","aim4","gdsc_axis2_drug_response.csv"), row.names=FALSE)
cat("Saved results/aim4/gdsc_axis2_drug_response.csv\n")
