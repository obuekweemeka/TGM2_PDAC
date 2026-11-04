# ============================================================================
# fig2_distributions_and_renumber_v1.R -- Save in: scripts/
# Restructured figure set:
#   NEW Figure 2 = distribution of the continuous contrast per patient (a)
#                  plus per-patient proportions (b)
#   Figure 3     = axis correlations, Moffitt, Steele replication (a-c)
#   old Figure 2 (UMAP landscape) -> Supplementary Figure 6
#   old Figure 5 (spatial + KM)   -> Supplementary Figure 7
# Figures 1 and 4 are unchanged. Replaced files go to _archive/.
# ============================================================================
suppressPackageStartupMessages({library(ggplot2); library(magick)})
base <- here::here(); pr <- file.path(base,"data","processed"); fig <- file.path(base,"figures")
mp <- file.path(fig,"multipanel")
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
thm <- function(fs) theme_classic(base_size=fs) + theme(plot.title=element_text(face="bold",size=fs+1), plot.subtitle=element_text(size=fs-1.5,colour="grey40"))

cat("=== 1. New Figure 2a: continuous contrast per patient ===\n")
v4 <- readRDS(file.path(pr,"aim2_scoring","malignant_cells_scored_v4.rds"))
mcm <- readRDS(file.path(pr,"aim2_scoring","malignant_cell_metadata.rds"))
tr <- unique(mcm[,c("pid","new_treatment")]); v4$trt <- tr$new_treatment[match(v4$pid,tr$pid)]
x <- v4[v4$trt=="Untreated" & v4$axis2_death_competence!=0, ]
x$cyto <- grepl("cytostatic", x$cell_state)
ordp <- names(sort(tapply(x$axis2_death_competence, x$pid, median)))
x$pid <- factor(x$pid, levels=ordp)
cat("cells:", nrow(x), "| patients:", nlevels(x$pid), "| median overall:", round(median(x$axis2_death_competence),3), "\n")
p2a <- ggplot(x, aes(pid, axis2_death_competence)) +
  geom_hline(yintercept=0, linetype="dashed", colour="grey40") +
  geom_violin(fill="#0072B2", colour=NA, alpha=.55, scale="width") +
  geom_boxplot(width=.12, outlier.shape=NA, fill="white", colour="grey30", linewidth=.3) +
  coord_flip() +
  labs(title="Apoptosis/arrest contrast per patient",
       subtitle="Untreated malignant cells with signal (n = 45,174); negative = anti-apoptotic and arrest dominant",
       x="Patient", y="Contrast score (log ratio)") + thm(9)
f2a <- file.path(fig,"aim2","Fig2A_contrast_distribution.png"); arch(f2a)
ggsave(f2a, p2a, dpi=dpi, width=wf, height=5.6, units="in"); cat("  saved panel a\n")

cat("\n=== 2. Figure 2 page (a = distributions, b = per-patient proportions) ===\n")
f2b <- file.path(fig,"aim2","Fig2B_per_patient_cytostatic.png"); stopifnot(file.exists(f2b))
ib <- image_info(image_read(f2b)); hb <- wf * ib$height/ib$width
cv <- image_blank(px(W),px(H),color="white"); y2 <- mt+5.6+gy
cv <- put(cv, image_read(f2a), ml, mt, wf)
cv <- put(cv, image_read(f2b), ml, y2, wf)
page(cv, data.frame(l=c("a","b"), x=c(.03,.03), y=c(mt,y2)), "Figure2.pdf")

cat("\n=== 3. Figure 3 page (a = axis correlations, b = Moffitt, c = Steele) ===\n")
fa <- file.path(fig,"aim2","Fig2C_axis_correlation_heatmap.png")
fb <- file.path(fig,"aim2","Fig2D_moffitt_alignment.png")
fc <- file.path(fig,"aim2","Fig2E_Steele_replication.png")
stopifnot(all(file.exists(c(fa,fb,fc))))
ch <- (H-mt-mb-gy)/2; yb <- mt+ch+gy
cv <- image_blank(px(W),px(H),color="white")
cv <- put(cv, image_read(fa), ml, mt, cw); cv <- put(cv, image_read(fb), xr, mt, cw)
cv <- put(cv, image_read(fc), ml, yb, cw)
page(cv, data.frame(l=c("a","b","c"), x=c(.03,xr-.19,.03), y=c(mt,mt,yb)), "Figure3.pdf")

cat("\n=== 4. Demote old Figure 2 and Figure 5 to supplementary ===\n")
mv <- function(from,to){ f <- file.path(mp,from); if (!file.exists(f)) { cat("  missing",from,"\n"); return(invisible()) }
  arch(file.path(mp,to)); file.copy(f, file.path(mp,to), overwrite=TRUE); cat("  ",from,"->",to,"\n") }
# old Figure2 content was overwritten above, so take the archived copy
oldfig2 <- file.path(fig,"_archive","Figure2_superseded.pdf")
if (file.exists(oldfig2)) { file.copy(oldfig2, file.path(mp,"FigureS6.pdf"), overwrite=TRUE); cat("   archived Figure2 -> FigureS6.pdf\n") } else cat("   NOTE: archived old Figure2 not found; rebuild S6 from the UMAP panels\n")
mv("Figure5.pdf","FigureS7.pdf")
f5 <- file.path(mp,"Figure5.pdf"); if (file.exists(f5)) { arch(f5); file.remove(f5); cat("   removed Figure5.pdf (now S7)\n") }
cat("\nMultipanel folder now:\n"); print(sort(list.files(mp)))
cat("\nExpected: Figure1-4 and FigureS1-S7\n")
