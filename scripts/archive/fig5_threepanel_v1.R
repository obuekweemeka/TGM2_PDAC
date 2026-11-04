# ============================================================================
# fig5_threepanel_v1.R -- Save in: scripts/
# The purity panel is removed: ABSOLUTE purity does not track Axis 2
# (r = 0.153, p = 0.078, n = 133), so the panel's premise no longer holds.
# Figure 5 becomes three panels: spatial heatmap, per-sample correlations, KM.
# ============================================================================
suppressPackageStartupMessages({library(ggplot2); library(magick)})
base <- here::here(); fig <- file.path(base,"figures")
W<-8.27; H<-11.69; dpi<-600; ml<-.22; mr<-.06; mt<-.06; mb<-.06; gx<-.22; gy<-.12
cw <- (W-ml-mr-gx)/2; xr <- ml+cw+gx; h1 <- 5.4
px <- function(z) as.integer(round(z*dpi))
put <- function(cv,img,x,y,w) image_composite(cv, image_scale(img,paste0(px(w),"x")), offset=paste0("+",px(x),"+",px(y)))
cv <- image_blank(px(W),px(H),color="white"); y5 <- mt+h1+gy
cv <- put(cv, image_read(file.path(fig,"aim3","Fig5A_distance_heatmap.png")), ml, mt, cw)
cv <- put(cv, image_read(file.path(fig,"aim3","Fig5B_distance_correlations.png")), xr, mt, cw)
cv <- put(cv, image_read(file.path(fig,"aim4","Fig4A_KM_cytostatic.png")), ml, y5, cw)
labs <- data.frame(l=c("a","b","c"), x=c(.03,xr-.19,.03), y=c(mt,mt,y5))
pg <- ggplot() + annotation_raster(as.raster(cv),xmin=0,xmax=W,ymin=0,ymax=H,interpolate=TRUE) +
  geom_text(data=labs, aes(x=x,y=H-y,label=l), hjust=0,vjust=1,fontface="bold",size=9/ggplot2::.pt) +
  coord_fixed(xlim=c(0,W),ylim=c(0,H),expand=FALSE) + theme_void() + theme(plot.margin=margin(0,0,0,0))
o <- file.path(fig,"multipanel","Figure5.pdf"); tmp <- sub("\\.pdf$","_new.pdf",o)
a <- file.path(fig,"_archive","Figure5_superseded_4panel.pdf"); if (file.exists(o) && !file.exists(a)) file.copy(o,a)
ggsave(tmp,pg,width=W,height=H,units="in",device="pdf"); stopifnot(file.size(tmp)>1e5)
file.copy(tmp,o,overwrite=TRUE); file.remove(tmp)
cat("Figure5.pdf rebuilt with three panels:", round(file.size(o)/1024),"KB\n")
pp <- file.path(fig,"aim4","Fig4C_axis2_purity_TCGA.png")
ap <- file.path(fig,"_archive","Fig4C_axis2_purity_TCGA_unreproducible.png")
if (file.exists(pp)) { file.copy(pp,ap); file.remove(pp); cat("Purity panel moved to _archive (its r = 0.550 could not be reproduced)\n") }
