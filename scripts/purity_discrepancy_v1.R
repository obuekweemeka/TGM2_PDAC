# ============================================================================
# purity_discrepancy_v1.R -- READ-ONLY diagnostic. Save in: scripts/
# The manuscript reports Spearman r = 0.550, p = 7.43e-13, n = 145 between the
# Axis 2 ssGSEA score and ABSOLUTE tumour purity. Recomputation from the
# PanCanAtlas ABSOLUTE table gives r = 0.153, n = 133. This finds out why.
# ============================================================================
base <- normalizePath('.'); pr <- file.path(base,"data","processed")
ab <- read.delim(file.path(base,"data","raw","TCGA_purity","TCGA_mastercalls.abs_tables_JSedit.fixed.txt"), stringsAsFactors=FALSE)
tc <- readRDS(file.path(pr,"aim4_bulk","tcga_paad_analysis_df.rds"))
cat("TCGA object: n =", nrow(tc), "| columns:", paste(names(tc), collapse=", "), "\n")
cat("moffitt column class:", class(tc$moffitt), "| values:", paste(head(unique(tc$moffitt),6), collapse=", "), "\n\n")

cat("=== 1. Matching: which key and how many ===\n")
for (col in c("array","sample")) {
  ab$k <- substr(gsub("\\.","-",ab[[col]]),1,15); tc$k <- substr(gsub("\\.","-",tc$sample_id),1,15)
  m <- merge(tc, unique(ab[,c("k","purity")]), by="k")
  cat(sprintf("  key from ab$%-7s: matched %d, with purity %d\n", col, nrow(m), sum(!is.na(m$purity))))
}
ab$k <- substr(gsub("\\.","-",ab$array),1,15); tc$k <- substr(gsub("\\.","-",tc$sample_id),1,15)
cat("  unmatched TCGA samples:", sum(!tc$k %in% ab$k), "\n")
cat("  TCGA sample types present:", paste(sort(unique(substr(tc$sample_id,14,15))), collapse=", "), "\n")

cat("\n=== 2. Every ABSOLUTE metric against every axis ===\n")
mg <- merge(tc, unique(ab[,c("k","purity","ploidy","Cancer.DNA.fraction","Subclonal.genome.fraction","Genome.doublings")]), by="k")
ax <- grep("^axis", names(mg), value=TRUE)
met <- c("purity","ploidy","Cancer.DNA.fraction","Subclonal.genome.fraction","Genome.doublings")
for (mt in met) { v <- suppressWarnings(as.numeric(mg[[mt]])); if (all(is.na(v))) next
  out <- sapply(ax, function(a){ ct <- suppressWarnings(cor.test(v, mg[[a]], method="spearman")); sprintf("%.3f (p=%.3g)", ct$estimate, ct$p.value) })
  cat(sprintf("  %-26s n=%3d | %s\n", mt, sum(!is.na(v)), paste(sprintf("%s %s", ax, out), collapse=" | "))) }

cat("\n=== 3. Could 0.550 come from a different pairing? ===\n")
num <- names(mg)[sapply(mg, is.numeric)]
hits <- list()
for (i in seq_along(num)) for (j in seq_along(num)) if (i < j) {
  r <- suppressWarnings(cor(mg[[num[i]]], mg[[num[j]]], method="spearman", use="complete.obs"))
  if (!is.na(r) && abs(abs(r) - 0.550) < 0.03) hits[[length(hits)+1]] <- sprintf("  %s vs %s: r = %.3f", num[i], num[j], r) }
if (length(hits)) cat(paste(unlist(hits), collapse="\n"), "\n") else cat("  no numeric pair in this object gives |r| near 0.550\n")

cat("\n=== 4. Is an alternative purity estimate stored anywhere? ===\n")
f <- list.files(base, pattern="\\.(rds|csv|tsv|txt)$", recursive=TRUE, full.names=TRUE)
f <- f[grepl("aim4|tcga|TCGA", f) & file.size(f) < 5e7]
for (p in f) { d <- tryCatch(if (grepl("\\.rds$",p)) readRDS(p) else read.delim(p, nrows=2), error=function(e) NULL)
  if (is.data.frame(d)) { hit <- grep("purity|estimate|stromal|immune|CPE|tumor_frac", names(d), ignore.case=TRUE, value=TRUE)
    if (length(hit)) cat("  ", sub(paste0(base,"/"),"",p), "->", paste(hit, collapse=", "), "\n") } }

cat("\n=== 5. What the figure plotted ===\n")
cat("  Figure 5d shows roughly 145 points spanning purity 0.1 to 0.9 with a clear positive slope.\n")
cat("  ABSOLUTE purity here matches only 133 PAAD samples and gives r = 0.153.\n")
cat("  If no stored alternative is found above, the plotted purity source cannot be identified\n")
cat("  from saved files, and the reported r = 0.550 cannot be reproduced.\n")
