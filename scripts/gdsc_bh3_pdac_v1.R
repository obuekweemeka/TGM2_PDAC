# ============================================================================
# gdsc_bh3_pdac_v1.R -- Save in: scripts/. Uses data already on disk.
# Axis 2 is a BCL-2 family balance score, so the pharmacological readout that
# matches the construct is BH3 mimetic sensitivity, not MEK inhibition.
#
# PRESPECIFIED (fixed before any association is computed):
#   Primary   : BH3 mimetics and direct apoptosis modulators
#   Secondary : cytotoxic chemotherapy (apoptosis-dependent killing)
#   Reference : MAPK inhibitors, already tested and null
#   Metric    : LN_IC50 (primary), AUC (secondary)
#   Direction : HIGHER Axis 2 = more death-competent = MORE sensitive
#               = LOWER ln(IC50), so expected Spearman rho is NEGATIVE
#   Multiplicity: Benjamini-Hochberg within the primary set
#   Axis 2 here is the cell-line apoptosis/arrest contrast (axis2_death): the mean of the pro-apoptotic and caspase/mitochondrial module means minus the mean of the anti-apoptotic and cytostasis module means, on log2(TPM+1). It is not the single-cell UCell score.
# ============================================================================
base <- normalizePath('.'); pr <- file.path(base,"data","processed")
raw <- file.path(base,"data","raw","GDSC")
if (!requireNamespace("readxl", quietly=TRUE)) install.packages("readxl")
de <- readRDS(file.path(pr,"aim1_signatures","DepMap_PDAC_expression.rds"))
d4 <- readRDS(file.path(pr,"aim4_bulk","depmap_pdac_analysis_df.rds"))
d <- d4[d4$model_id %in% intersect(de$ModelID[!duplicated(de$ModelID)], d4$model_id), ]
clean <- function(x) toupper(gsub("[^A-Za-z0-9]","",as.character(x)))
mod <- read.csv(file.path(base,"data","raw","DepMap","Model.csv"), stringsAsFactors=FALSE)
mp <- mod[match(d$model_id, mod$ModelID), ]
d$key <- clean(mp$StrippedCellLineName); d$key[is.na(d$key)] <- clean(d$cell_line[is.na(d$key)])
d$cosmic <- suppressWarnings(as.integer(mp$COSMICID))

g <- as.data.frame(readxl::read_excel(file.path(raw,"GDSC2_fitted_dose_response_27Oct23.xlsx")))
g$key <- clean(g$CELL_LINE_NAME)
d$match <- ifelse(d$key %in% g$key, d$key, NA)
cm <- unique(g[,c("COSMIC_ID","key")]); nx <- is.na(d$match) & d$cosmic %in% cm$COSMIC_ID
d$match[nx] <- cm$key[match(d$cosmic[nx], cm$COSMIC_ID)]
dm <- d[!is.na(d$match), ]; cat("PDAC lines matched to GDSC:", nrow(dm), "of", nrow(d), "\n")

BH3 <- c("Navitoclax","Venetoclax","ABT737","Obatoclax Mesylate","Sabutoclax","AZD5991","MIM1","TW 37","UMI-77","A-1155463","AMG-176","S63845","ABT-199","Venotoclax")
CHEMO <- c("Gemcitabine","Cisplatin","Oxaliplatin","Paclitaxel","Docetaxel","Doxorubicin","Irinotecan","SN-38","Fluorouracil","Etoposide","Camptothecin","Topotecan","Vincristine")
MAPK <- c("Trametinib","Selumetinib","Refametinib","Ulixertinib","SCH772984","KRAS (G12C) Inhibitor-12")
avail <- function(v) intersect(v, unique(g$DRUG_NAME))
cat("\nBH3 mimetics in GDSC2:", paste(avail(BH3), collapse=", "), "\n")
cat("Chemotherapy in GDSC2:", paste(avail(CHEMO), collapse=", "), "\n")

tst <- function(drug) {
  s <- g[g$DRUG_NAME==drug & g$key %in% dm$match, c("key","LN_IC50","AUC")]
  if (!nrow(s)) return(NULL)
  s <- aggregate(cbind(LN_IC50,AUC) ~ key, s, mean)
  m <- merge(dm[,c("match","axis2_death")], s, by.x="match", by.y="key")
  if (nrow(m) < 10) return(NULL)
  a <- suppressWarnings(cor.test(m$axis2_death, m$LN_IC50, method="spearman"))
  b <- suppressWarnings(cor.test(m$axis2_death, m$AUC, method="spearman"))
  data.frame(drug=drug, n=nrow(m), rho_lnIC50=round(a$estimate,3), p_lnIC50=signif(a$p.value,3),
             rho_AUC=round(b$estimate,3), p_AUC=signif(b$p.value,3)) }

out <- list()
for (set in list(c("PRIMARY: BH3 mimetic", list(avail(BH3))), c("SECONDARY: chemotherapy", list(avail(CHEMO))), c("REFERENCE: MAPK inhibitor", list(avail(MAPK))))) {
  lab <- set[[1]]; dr <- set[[2]]
  r <- do.call(rbind, lapply(dr, tst))
  cat("\n---", lab, "| predicted direction NEGATIVE rho ---\n")
  if (is.null(r)) { cat("  no drug reached 10 matched lines\n"); next }
  if (grepl("PRIMARY", lab)) r$padj_BH <- signif(p.adjust(r$p_lnIC50,"BH"),3)
  print(r[order(r$p_lnIC50), ], row.names=FALSE)
  cat(sprintf("  in predicted direction: %d of %d | median rho %.3f\n", sum(r$rho_lnIC50<0), nrow(r), median(r$rho_lnIC50)))
  r$set <- lab; out[[lab]] <- r }
res <- do.call(rbind, out)
dir.create(file.path(base,"results","aim4"), recursive=TRUE, showWarnings=FALSE)
write.csv(res, file.path(base,"results","aim4","gdsc_axis2_bh3_pdac.csv"), row.names=FALSE)
cat("\nSaved results/aim4/gdsc_axis2_bh3_pdac.csv\n")
cat("\nPower note: with about 20-24 lines only a strong association (|rho| > 0.45) reaches p < 0.05.\n")
