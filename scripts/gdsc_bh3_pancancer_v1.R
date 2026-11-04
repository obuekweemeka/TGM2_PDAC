# ============================================================================
# gdsc_bh3_pancancer_v1.R -- Save in: scripts/
# A better-powered, mechanistically direct test of Axis 2.
#
# RATIONALE: Axis 2 contrasts pro-apoptotic + executor-caspase genes against
# anti-apoptotic + arrest genes. The pharmacological readout matching that
# construct is BH3 mimetic sensitivity, not MEK inhibition. Axis 2 uses generic
# apoptosis genes, so it can be computed for every cell line, not just PDAC.
#
# PRESPECIFIED (fixed before any association is computed):
#   Primary drugs : BH3 mimetics / direct apoptosis modulators
#   Secondary     : cytotoxic chemotherapy (apoptosis-dependent killing)
#   Tertiary      : MAPK-pathway inhibitors (already tested in PDAC, null)
#   Metric        : LN_IC50, with AUC as secondary
#   Direction     : HIGHER Axis 2 = more death-competent = MORE sensitive
#                   = LOWER ln(IC50), so expected Spearman rho is NEGATIVE
#   Multiplicity  : Benjamini-Hochberg within each drug set
#   Reported      : pan-cancer primary analysis, PDAC subset separately
# ============================================================================
base <- normalizePath('.'); pr <- file.path(base,"data","processed")
raw <- file.path(base,"data","raw","GDSC"); dir.create(raw, recursive=TRUE, showWarnings=FALSE)
options(timeout=7200)
if (!requireNamespace("readxl", quietly=TRUE)) install.packages("readxl")

cat("=== 1. GDSC basal expression (about 250 MB, one-off) ===\n")
zp <- file.path(raw,"Cell_line_RMA_proc_basalExp.txt.zip"); ex <- file.path(raw,"Cell_line_RMA_proc_basalExp.txt")
if (!file.exists(ex)) {
  if (!file.exists(zp) || file.size(zp) < 1e6) {
    cat("downloading...\n")
    tryCatch(download.file("https://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources//Data/preprocessed/Cell_line_RMA_proc_basalExp.txt.zip",
             zp, mode="wb", quiet=TRUE), error=function(e) cat("  failed:", conditionMessage(e), "\n")) }
  if (file.exists(zp) && file.size(zp) > 1e6) { cat("  ",round(file.size(zp)/1e6,1),"MB; unzipping\n"); unzip(zp, exdir=raw) }
}
if (!file.exists(ex)) stop("Expression file not available; download manually from https://www.cancerrxgene.org/downloads/bulk_download")
e <- read.delim(ex, check.names=FALSE, stringsAsFactors=FALSE)
cat("Expression matrix:", nrow(e), "genes x", ncol(e)-2, "lines\n")
sym <- e[[1]]; mat <- as.matrix(e[, grep("^DATA\\.", names(e))])
colnames(mat) <- sub("^DATA\\.", "", colnames(mat))          # COSMIC IDs
rownames(mat) <- sym; mat <- mat[!is.na(sym) & sym != "", , drop=FALSE]
cat("Usable:", nrow(mat), "genes x", ncol(mat), "lines (columns are COSMIC IDs)\n")

cat("\n=== 2. Axis 2 across all lines ===\n")
sig <- readRDS(file.path(pr,"aim1_signatures","four_axis_signatures_v1.rds"))
mods <- sig[c("axis2_proapoptotic","axis2_executors","axis2_anti_apoptotic","axis2_cytostasis")]
for (m in names(mods)) cat(sprintf("  %-24s %d of %d genes present\n", m, length(intersect(mods[[m]], rownames(mat))), length(mods[[m]])))
zs <- function(v) { s <- sd(v, na.rm=TRUE); if (is.na(s) || s == 0) v*0 else (v - mean(v, na.rm=TRUE))/s }
modscore <- function(gs) { g <- intersect(gs, rownames(mat)); if (!length(g)) return(rep(NA, ncol(mat)))
  colMeans(t(apply(mat[g, , drop=FALSE], 1, zs)), na.rm=TRUE) }   # gene-wise z, then mean across genes
A2 <- (modscore(mods[[1]]) + modscore(mods[[2]])) - (modscore(mods[[3]]) + modscore(mods[[4]]))
cat("Axis 2 (z-based difference) computed for", sum(!is.na(A2)), "lines | range",
    paste(round(range(A2, na.rm=TRUE),2), collapse=" to "), "\n")
ax <- data.frame(COSMIC_ID=as.integer(colnames(mat)), axis2=as.numeric(A2))

cat("\n=== 3. Drug response ===\n")
g <- as.data.frame(readxl::read_excel(file.path(raw,"GDSC2_fitted_dose_response_27Oct23.xlsx")))
BH3 <- c("Navitoclax","Venetoclax","ABT737","Obatoclax Mesylate","Sabutoclax","AZD5991","MIM1","TW 37","UMI-77","A-1155463","AMG-176","S63845")
CHEMO <- c("Gemcitabine","Cisplatin","Oxaliplatin","Paclitaxel","Docetaxel","Doxorubicin","Irinotecan","SN-38","Fluorouracil","Etoposide","Camptothecin")
MAPK <- c("Trametinib","Selumetinib","Refametinib","Ulixertinib","SCH772984","KRAS (G12C) Inhibitor-12")
have <- function(v) intersect(v, unique(g$DRUG_NAME))
cat("BH3 mimetics present:", paste(have(BH3), collapse=", "), "\n")
cat("Chemotherapy present:", paste(have(CHEMO), collapse=", "), "\n")
sets <- list(`1 BH3 mimetic`=have(BH3), `2 Chemotherapy`=have(CHEMO), `3 MAPK inhibitor`=have(MAPK))

test <- function(drug, cos_keep=NULL, label="") {
  s <- g[g$DRUG_NAME==drug, c("COSMIC_ID","LN_IC50","AUC")]
  if (!is.null(cos_keep)) s <- s[s$COSMIC_ID %in% cos_keep, ]
  if (!nrow(s)) return(NULL)
  s <- aggregate(cbind(LN_IC50,AUC) ~ COSMIC_ID, s, mean)
  m <- merge(s, ax, by="COSMIC_ID"); m <- m[!is.na(m$axis2), ]
  if (nrow(m) < 15) return(NULL)
  c1 <- suppressWarnings(cor.test(m$axis2, m$LN_IC50, method="spearman"))
  c2 <- suppressWarnings(cor.test(m$axis2, m$AUC, method="spearman"))
  data.frame(cohort=label, drug=drug, n=nrow(m), rho_lnIC50=round(c1$estimate,3), p_lnIC50=signif(c1$p.value,3),
             rho_AUC=round(c2$estimate,3), p_AUC=signif(c2$p.value,3)) }

paad <- unique(g$COSMIC_ID[g$TCGA_DESC=="PAAD"])
for (lab in names(sets)) {
  dr <- sets[[lab]]; if (!length(dr)) next
  pan <- do.call(rbind, lapply(dr, test, cos_keep=NULL, label="pan-cancer"))
  if (is.null(pan)) next
  pan$padj_BH <- signif(p.adjust(pan$p_lnIC50,"BH"),3)
  cat("\n---", lab, "(pan-cancer; predicted direction NEGATIVE rho) ---\n")
  print(pan[order(pan$p_lnIC50), c("drug","n","rho_lnIC50","p_lnIC50","padj_BH","rho_AUC")], row.names=FALSE)
  cat(sprintf("  in predicted direction: %d of %d | median rho %.3f | significant after BH: %d\n",
              sum(pan$rho_lnIC50<0), nrow(pan), median(pan$rho_lnIC50), sum(pan$padj_BH<0.05)))
  pd <- do.call(rbind, lapply(dr, test, cos_keep=paad, label="PDAC"))
  if (!is.null(pd)) { cat("  PDAC subset:\n"); print(pd[, c("drug","n","rho_lnIC50","p_lnIC50")], row.names=FALSE) }
  assign(paste0("res_", substr(lab,1,1)), list(pan=pan, pdac=pd))
}
all <- do.call(rbind, lapply(ls(pattern="^res_"), function(x) do.call(rbind, get(x))))
dir.create(file.path(base,"results","aim4"), recursive=TRUE, showWarnings=FALSE)
write.csv(all, file.path(base,"results","aim4","gdsc_axis2_bh3_pancancer.csv"), row.names=FALSE)
cat("\nSaved results/aim4/gdsc_axis2_bh3_pancancer.csv\n")
cat("\nNOTE: Axis 2 here is a z-score difference on microarray data, not the UCell log-ratio\n",
    "used in the single-cell work. It is the same gene sets and contrast, computed in the way\n",
    "the data allow. That difference must be stated in the Methods.\n")
