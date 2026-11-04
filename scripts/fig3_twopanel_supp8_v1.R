# ============================================================================
# fig3_twopanel_supp8_v1.R -- Save in: scripts/
# Figure 3 becomes two panels (a = Moffitt subtype, b = Steele replication).
# The four-axis correlation heatmap moves to Supplementary Figure 8, since the
# main narrative now concerns Axis 2 only.
# ============================================================================
suppressPackageStartupMessages({library(ggplot2); library(magick)})
base <- normalizePath('.'); fig <- file.path(base,"figures"); mp <- file.path(fig,"multipanel")
W<-8.27; H<-11.69; dpi<-600; ml<-.22; mr<-.06; mt<-.06; mb<-.06; gx<-.22; gy<-.12
cw <- (W-ml-mr-gx)/2; wf <- W-ml-mr; xr <- ml+cw+gx
arch <- function(p){a<-file.path(fig,"_archive",paste0(tools::file_path_sans_ext(basename(p)),"_superseded.",tools::file_ext(p))); if(file.exists(p)&&!file.exists(a)) file.copy(p,a)}
px <- function(z) as.integer(round(z*dpi))
put <- function(cv,img,x,y,w) image_composite(cv, image_scale(img,paste0(px(w),"x")), offset=paste0("+",px(x),"+",px(y)))
page <- function(cv,labs,name){ p <- ggplot() + annotation_raster(as.raster(cv),xmin=0,xmax=W,ymin=0,ymax=H,interpolate=TRUE) +
    geom_text(data=labs,aes(x=x,y=H-y,label=l),hjust=0,vjust=1,fontface="bold",size=9/ggplot2::.pt) +
    coord_fixed(xlim=c(0,W),ylim=c(0,H),expand=FALSE) + theme_void() + theme(plot.margin=margin(0,0,0,0))
  out <- file.path(mp,name); tmp <- sub("\\.pdf$","_new.pdf",out); arch(out)
  ggsave(tmp,p,width=W,height=H,units="in",device="pdf"); stopifnot(file.size(tmp)>1e5)
  file.copy(tmp,out,overwrite=TRUE); file.remove(tmp); cat("  wrote",name,round(file.size(out)/1024),"KB\n") }

fm <- file.path(fig,"aim2","Fig2D_moffitt_alignment.png")      # Moffitt
fs <- file.path(fig,"aim2","Fig2E_Steele_replication.png")     # Steele
fc <- file.path(fig,"aim2","Fig2C_axis_correlation_heatmap.png") # axis correlations
stopifnot(all(file.exists(c(fm,fs,fc))))

cat("=== Figure 3: two panels ===\n")
ch <- (H-mt-mb-gy)/2; yb <- mt+ch+gy
cv <- image_blank(px(W),px(H),color="white")
cv <- put(cv, image_read(fm), ml, mt, cw)
cv <- put(cv, image_read(fs), xr, mt, cw)
page(cv, data.frame(l=c("a","b"), x=c(.03,xr-.19), y=c(mt,mt)), "Figure3.pdf")

cat("\n=== Supplementary Figure 8: axis correlations ===\n")
inf <- image_info(image_read(fc)); hh <- wf * inf$height/inf$width
cv <- image_blank(px(W),px(H),color="white")
cv <- put(cv, image_read(fc), ml, mt, wf)
page(cv, data.frame(l=character(0), x=numeric(0), y=numeric(0)), "FigureS8.pdf")

cat("\nMultipanel folder now:\n"); print(sort(list.files(mp)))
cat("Expected: Figure1-4 and FigureS1-S8\n")
