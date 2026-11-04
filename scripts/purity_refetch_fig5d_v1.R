# ============================================================================
# purity_refetch_fig5d_v1.R -- Save in: scripts/
# The ABSOLUTE purity values behind Figure 5d are not on disk, so the reported
# correlation (r = 0.550) has never been re-verified. This script re-fetches the
# PanCanAtlas ABSOLUTE purity table from the GDC, recomputes the correlation,
# and regenerates panel 5d with the corrected axis label and no unsupported
# subtitle. Figure5.pdf is rebuilt only if the data are recovered.
# ============================================================================
suppressPackageStartupMessages({library(ggplot2); library(magick)})
if (!requireNamespace("jsonlite", quietly = TRUE)) install.packages("jsonlite")
base <- normalizePath('.'); pr <- file.path(base,"data","processed"); fig <- file.path(base,"figures")
raw <- file.path(base,"data","raw","TCGA_purity"); dir.create(raw, recursive=TRUE, showWarnings=FALSE)
options(timeout = 3600)
AX <- "Axis 2 apoptosis/arrest contrast"
FN <- "TCGA_mastercalls.abs_tables_JSedit.fixed.txt"
dest <- file.path(raw, FN)

cat("=== 1. Obtain the ABSOLUTE purity table ===\n")
if (file.exists(dest) && file.size(dest) > 1e4) cat("Already present:", dest, "\n") else {
  q <- paste0('https://api.gdc.cancer.gov/files?filters=',
    URLencode('{"op":"=","content":{"field":"file_name","value":"TCGA_mastercalls.abs_tables_JSedit.fixed.txt"}}', reserved=TRUE),
    '&fields=file_id,file_name,file_size&format=JSON&size=5')
  j <- tryCatch(jsonlite::fromJSON(q), error=function(e){cat("GDC query failed:",conditionMessage(e),"\n"); NULL})
  hits <- if (!is.null(j)) j$data$hits else NULL
  if (!is.null(hits) && length(hits) && nrow(hits) > 0) {
    print(hits)
    url <- paste0("https://api.gdc.cancer.gov/data/", hits$file_id[1])
    cat("Downloading", url, "\n")
    tryCatch(download.file(url, dest, mode="wb", quiet=TRUE), error=function(e) cat("download failed:", conditionMessage(e), "\n"))
  } else cat("GDC lookup returned nothing.\n")
}
if (!file.exists(dest) || file.size(dest) < 1e4) {
  cat("\nAutomatic download did not work. Fetch it manually:\n",
      " 1. Open https://gdc.cancer.gov/about-data/publications/pancanatlas\n",
      " 2. Download 'ABSOLUTE purity/ploidy file -", FN, "'\n",
      " 3. Save it to", raw, "and re-run this script.\n")
} else {
  ab <- read.delim(dest, stringsAsFactors = FALSE)
  cat("Rows:", nrow(ab), "| columns:", paste(head(names(ab),8), collapse=", "), "\n")
  pcol <- grep("^purity$", names(ab), ignore.case=TRUE, value=TRUE)[1]
  scol <- grep("array|sample", names(ab), ignore.case=TRUE, value=TRUE)[1]
  cat("Using purity column:", pcol, "| sample column:", scol, "\n")
  ab$short <- substr(gsub("\\.", "-", ab[[scol]]), 1, 15)

  tc <- readRDS(file.path(pr,"aim4_bulk","tcga_paad_analysis_df.rds"))
  tc$short <- substr(gsub("\\.", "-", tc$sample_id), 1, 15)
  cat("TCGA object sample id example:", tc$sample_id[1], "-> key", tc$short[1], "\n")
  mg <- merge(tc, unique(ab[, c("short", pcol)]), by="short")
  names(mg)[names(mg)==pcol] <- "purity"; mg <- mg[!is.na(mg$purity), ]
  cat(sprintf("Matched %d of %d TCGA PAAD samples with a purity estimate\n", nrow(mg), nrow(tc)))

  if (nrow(mg) < 100) {
    cat("Too few matched; inspect the key formats above before interpreting.\n")
  } else {
    ct <- suppressWarnings(cor.test(mg$purity, mg$axis2_death_z, method="spearman"))
    cat(sprintf("\nSpearman r = %.3f, p = %.3g, n = %d\n", ct$estimate, ct$p.value, nrow(mg)))
    cat("Manuscript currently reports: r = 0.550, p = 7.43e-13, n = 145\n")
    sub <- if ("moffitt" %in% names(mg)) mg$moffitt else NA
    p5d <- ggplot(mg, aes(purity, axis2_death_z)) +
      {if (!all(is.na(sub))) geom_point(aes(colour = moffitt), size=1.3, alpha=.85) else geom_point(size=1.3, alpha=.85, colour="#0072B2")} +
      geom_smooth(method="lm", se=TRUE, colour="grey30", linewidth=.4, formula=y~x) +
      {if (!all(is.na(sub))) scale_colour_manual(values=c("#D55E00","#0072B2"), name="Moffitt subtype")} +
      labs(title=paste(AX, "tracks tumour purity"),
           subtitle=sprintf("TCGA PAAD, n = %d; Spearman r = %.3f, p = %.2g", nrow(mg), ct$estimate, ct$p.value),
           x="Tumour purity (ABSOLUTE)", y=paste(AX, "(scaled)")) +
      theme_classic(base_size=6.5) +
      theme(plot.title=element_text(face="bold",size=7.5), plot.subtitle=element_text(size=5.5,colour="grey40"), legend.position="bottom")
    out <- file.path(fig,"aim4","Fig4C_axis2_purity_TCGA.png")
    a <- file.path(fig,"_archive","Fig4C_axis2_purity_TCGA_superseded.png"); if (file.exists(out) && !file.exists(a)) file.copy(out,a)
    ggsave(out, p5d, dpi=600, width=(8.27-.28-.22)/2, height=5.4, units="in")
    cat("Saved panel 5d:", round(file.size(out)/1024), "KB\n")

    W<-8.27; H<-11.69; dpi<-600; ml<-.22; mr<-.06; mt<-.06; mb<-.06; gx<-.22; gy<-.12
    cw <- (W-ml-mr-gx)/2; xr <- ml+cw+gx; h1 <- 5.4; y5 <- mt+h1+gy
    px <- function(z) as.integer(round(z*dpi))
    put <- function(cv,img,x,y,w) image_composite(cv, image_scale(img,paste0(px(w),"x")), offset=paste0("+",px(x),"+",px(y)))
    cv <- image_blank(px(W),px(H),color="white")
    cv <- put(cv, image_read(file.path(fig,"aim3","Fig5A_distance_heatmap.png")), ml, mt, cw)
    cv <- put(cv, image_read(file.path(fig,"aim3","Fig5B_distance_correlations.png")), xr, mt, cw)
    cv <- put(cv, image_read(file.path(fig,"aim4","Fig4A_KM_cytostatic.png")), ml, y5, cw)
    cv <- put(cv, image_read(out), xr, y5, cw)
    labs <- data.frame(l=letters[1:4], x=c(.03,xr-.19,.03,xr-.19), y=c(mt,mt,y5,y5))
    pg <- ggplot() + annotation_raster(as.raster(cv),xmin=0,xmax=W,ymin=0,ymax=H,interpolate=TRUE) +
      geom_text(data=labs, aes(x=x,y=H-y,label=l), hjust=0,vjust=1,fontface="bold",size=9/ggplot2::.pt) +
      coord_fixed(xlim=c(0,W),ylim=c(0,H),expand=FALSE) + theme_void() + theme(plot.margin=margin(0,0,0,0))
    o <- file.path(fig,"multipanel","Figure5.pdf"); tmp <- sub("\\.pdf$","_new.pdf",o)
    a2 <- file.path(fig,"_archive","Figure5_superseded.pdf"); if (file.exists(o) && !file.exists(a2)) file.copy(o,a2)
    ggsave(tmp,pg,width=W,height=H,units="in",device="pdf"); stopifnot(file.size(tmp)>1e5)
    file.copy(tmp,o,overwrite=TRUE); file.remove(tmp); cat("Figure5.pdf rebuilt:", round(file.size(o)/1024),"KB\n")
    write.csv(mg[,c("sample_id","purity","axis2_death_z")], file.path(base,"results","aim4","tcga_purity_axis2.csv"), row.names=FALSE)
    cat("Saved results/aim4/tcga_purity_axis2.csv\n")
  }
}
