# ============================================================================
# gdsc_feasibility_v1.R -- READ-ONLY feasibility check. Save in: scripts/
# Can the Axis 2 score be tested against measured drug sensitivity?
#  1. What cell-line identifiers do we hold for the 44 DepMap PDAC lines?
#  2. Download GDSC2 fitted dose response + cell line details (about 21 MB)
#  3. How many of the 44 lines match, and which MAPK-pathway drugs cover them?
# PRESPECIFIED DIRECTION (recorded before any association is computed):
#   cytostatic-prone lines (LOWER Axis 2 death competence) are expected to be
#   LESS sensitive, i.e. HIGHER ln(IC50), to MAPK-pathway inhibitors.
# ============================================================================
base <- here::here(); pr <- file.path(base,"data","processed")
raw <- file.path(base,"data","raw","GDSC"); dir.create(raw, recursive=TRUE, showWarnings=FALSE)
options(timeout=3600)
if (!requireNamespace("readxl", quietly=TRUE)) install.packages("readxl")

cat("=== 1. Identifiers available for the 44 PDAC lines ===\n")
de <- readRDS(file.path(pr,"aim1_signatures","DepMap_PDAC_expression.rds"))
d4 <- readRDS(file.path(pr,"aim4_bulk","depmap_pdac_analysis_df.rds"))
v44 <- intersect(de$ModelID[!duplicated(de$ModelID)], d4$model_id)
cat("dep_expr columns (first 12):", paste(head(names(de),12), collapse=", "), "\n")
cat("dep4 columns:", paste(names(d4), collapse=", "), "\n")
idcols <- grep("name|Name|CCLE|stripped|cosmic|COSMIC|sanger|Sanger", c(names(de),names(d4)), value=TRUE)
cat("identifier-like columns:", if(length(idcols)) paste(unique(idcols),collapse=", ") else "NONE beyond ModelID", "\n")
mm <- list.files(base, pattern="Model|sample_info|cell_line", recursive=TRUE, full.names=TRUE, ignore.case=TRUE)
mm <- mm[grepl("\\.csv$|\\.tsv$|\\.rds$", mm)]
cat("local DepMap metadata candidates:", if(length(mm)) paste(sub(paste0(base,"/"),"",mm),collapse=" | ") else "none", "\n")

cat("\n=== 2. GDSC2 download ===\n")
f1 <- file.path(raw,"GDSC2_fitted_dose_response_27Oct23.xlsx")
f2 <- file.path(raw,"Cell_Lines_Details.xlsx")
get <- function(url,dest,label){ if (file.exists(dest) && file.size(dest)>1e5) { cat("have",label,"\n"); return(TRUE) }
  cat("downloading",label,"...\n")
  ok <- tryCatch({download.file(url,dest,mode="wb",quiet=TRUE); TRUE}, error=function(e){cat("  failed:",conditionMessage(e),"\n"); FALSE})
  if (ok) cat("  ",round(file.size(dest)/1e6,1),"MB\n"); ok }
ok1 <- get("https://cog.sanger.ac.uk/cancerrxgene/GDSC_release8.5/GDSC2_fitted_dose_response_27Oct23.xlsx", f1, "GDSC2 fitted dose response")
ok2 <- get("https://cog.sanger.ac.uk/cancerrxgene/GDSC_release8.5/Cell_Lines_Details.xlsx", f2, "cell line details")
if (!ok1) { cat("\nIf the download failed, get the file manually from\n  https://www.cancerrxgene.org/downloads/bulk_download\n  and save it to", raw, "\n") } else {
  g <- readxl::read_excel(f1)
  cat("GDSC2 rows:", nrow(g), "| columns:", paste(names(g), collapse=", "), "\n")

  cat("\n=== 3. Drug coverage and line matching ===\n")
  dn <- unique(g$DRUG_NAME)
  mapk <- grep("tramet|selumet|cobimet|binimet|refamet|mirdamet|ulixert|ravoxert|SCH772984|LY3214996|erlotin|gefitin|afatin|sotorasib|adagrasib|MRTX|AMG 510|RMC-|KRAS", dn, ignore.case=TRUE, value=TRUE)
  cat("MAPK/RAS-pathway drugs present:", if(length(mapk)) paste(sort(mapk),collapse=", ") else "none found", "\n")
  panc <- grep("panc|pancrea", unique(g$TCGA_DESC), ignore.case=TRUE, value=TRUE)
  cat("pancreatic TCGA_DESC codes:", paste(unique(c(panc,"PAAD")), collapse=", "), "\n")
  gp <- g[g$TCGA_DESC %in% c("PAAD", panc), ]
  cat("GDSC2 pancreatic cell lines:", length(unique(gp$CELL_LINE_NAME)), "| rows:", nrow(gp), "\n")
  clean <- function(x) toupper(gsub("[^A-Za-z0-9]","",x))
  cat("\nExample GDSC pancreatic lines:", paste(head(sort(unique(gp$CELL_LINE_NAME)),12), collapse=", "), "\n")
  if (length(mapk)) {
    sub <- gp[gp$DRUG_NAME %in% mapk, ]
    tb <- table(sub$DRUG_NAME)
    cat("\nPancreatic lines screened per MAPK drug:\n"); print(sort(tb, decreasing=TRUE))
  }
  cat("\nNEXT: the 44 DepMap lines are identified only by ModelID above. To match GDSC we need\n",
      "DepMap Model.csv (ModelID -> StrippedCellLineName / COSMICID). If it is not listed in\n",
      "section 1, it can be downloaded from https://depmap.org/portal/data_page/?tab=allData\n")
}
