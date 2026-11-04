# ============================================================================
# final_figure_renumber_v2.R -- Save in: scripts/
# Final figure set to match manuscript v30: 3 main figures, 9 supplementary,
# supplementary numbered in order of first citation.
#
#   Figure 1 <- Figure2.pdf    contrast distribution per patient
#   Figure 2 <- Figure3.pdf    Moffitt subtype, Steele replication
#   Figure 3 <- Figure4.pdf    TGM2
#   Supp 1   <- FigureS2.pdf   Axis 1 and KRAS dependency (DepMap)
#   Supp 2   <- FigureS8.pdf   pairwise axis correlations
#   Supp 3   <- Figure1.pdf    resistance volcano and GSEA
#   Supp 4   <- FigureS1.pdf   resistance sequencing QC
#   Supp 5   <- FigureS6.pdf   UMAP landscape
#   Supp 6   <- FigureS5.pdf   threshold sensitivity
#   Supp 7   <- FigureS4.pdf   spatial data characteristics
#   Supp 8   <- FigureS7.pdf   spatial distance analysis and TCGA survival
#   Supp 9   <- FigureS3.pdf   Harmony integration QC and DUSP6
# Renames via temporary names; backs everything up first; verifies the result.
# ============================================================================
base <- normalizePath('.')
mp <- file.path(base,"figures","multipanel"); dir.create(file.path(base,"docs"), showWarnings=FALSE)
bk <- file.path(base,"figures","_archive","preFinalRenumber"); dir.create(bk, recursive=TRUE, showWarnings=FALSE)

map <- c("Figure2.pdf"="Figure1.pdf",  "Figure3.pdf"="Figure2.pdf",  "Figure4.pdf"="Figure3.pdf",
         "FigureS2.pdf"="FigureS1.pdf","FigureS8.pdf"="FigureS2.pdf","Figure1.pdf"="FigureS3.pdf",
         "FigureS1.pdf"="FigureS4.pdf","FigureS6.pdf"="FigureS5.pdf","FigureS5.pdf"="FigureS6.pdf",
         "FigureS4.pdf"="FigureS7.pdf","FigureS7.pdf"="FigureS8.pdf","FigureS3.pdf"="FigureS9.pdf")
cur <- list.files(mp, pattern="\\.pdf$")
cat("Before:", paste(sort(cur), collapse=", "), "\n\n")
miss <- setdiff(names(map), cur)
if (length(miss)) { cat("ABORTED, missing:", paste(miss, collapse=", "), "\n") } else {
  for (f in cur) file.copy(file.path(mp,f), file.path(bk,f), overwrite=TRUE)
  cat("Backed up", length(cur), "files to figures/_archive/preFinalRenumber/\n")
  for (f in names(map)) file.rename(file.path(mp,f), file.path(mp, paste0("TMP_", map[[f]])))
  for (f in list.files(mp, pattern="^TMP_")) file.rename(file.path(mp,f), file.path(mp, sub("^TMP_","",f)))
  after <- sort(list.files(mp, pattern="\\.pdf$"))
  want <- c(paste0("Figure",1:3,".pdf"), paste0("FigureS",1:9,".pdf"))
  cat("\nAfter:", paste(after, collapse=", "), "\n")
  cat("Exactly the expected set:", setequal(after, want), "\n")
  if (!setequal(after, want)) cat("  unexpected:", paste(setdiff(after,want),collapse=", "),
                                  "| missing:", paste(setdiff(want,after),collapse=", "), "\n")
  for (f in want) cat(sprintf("  %-15s %6.0f KB\n", f, file.size(file.path(mp,f))/1024))
  writeLines(c("# Final figure map (manuscript v30)","",
    "Main figures and supplementary figures, supplementary numbered in order of first citation.","",
    "| Figure | Content |","|---|---|",
    "| 1 | Contrast distribution per patient; per-patient proportions |",
    "| 2 | Moffitt subtype association; Steele replication |",
    "| 3 | TGM2 in both cohorts; LIANA ranks; DepMap essentiality |",
    "| S1 | Axis 1 and KRAS dependency in DepMap |",
    "| S2 | Pairwise axis correlations |",
    "| S3 | Resistance transcriptomics: volcano and GSEA |",
    "| S4 | Resistance sequencing quality control |",
    "| S5 | UMAP landscape of malignant cells |",
    "| S6 | Threshold sensitivity |",
    "| S7 | Spatial data characteristics |",
    "| S8 | Spatial distance analysis and TCGA survival |",
    "| S9 | Harmony integration QC and DUSP6 |","",
    "Supplementary tables 1-3 are CSV files in results/aim2 and results/aim4.",
    "Previous file names are preserved in figures/_archive/preFinalRenumber/."),
    file.path(base,"docs","figure_map.md"))
  write.csv(data.frame(new=unname(map), old=names(map)), file.path(base,"docs","figure_renumber_map.csv"), row.names=FALSE)
  cat("\nWrote docs/figure_map.md and docs/figure_renumber_map.csv\n")
}
