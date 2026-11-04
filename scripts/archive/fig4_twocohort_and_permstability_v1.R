# ============================================================================
# fig4_twocohort_and_permstability_v1.R -- Save in: scripts/
#  PART 1: rebuild Figure 4 panel a as a two-cohort paired plot (Hwang + Steele)
#          and re-assemble Figure4.pdf
#  PART 2: repeat the size-equalised gene-set permutation on two further
#          independent 8,000-cell subsamples (500 permutations each)
# Part 2 takes roughly 15 minutes.
# ============================================================================
suppressPackageStartupMessages({library(ggplot2); library(patchwork); library(magick); library(Matrix); library(UCell)})
base <- here::here(); fig <- file.path(base,"figures"); pr <- file.path(base,"data","processed")
W <- 8.27; H <- 11.69; dpi <- 600; ml <- .22; mr <- .06; mt <- .06; mb <- .06; gx <- .22; gy <- .12; fs <- 6.5
cw <- (W-ml-mr-gx)/2; h_top <- 4.2; w_c <- W-ml-mr; h_c <- H-mt-h_top-gy-mb
thm <- function(fs) theme_classic(base_size=fs) + theme(plot.title=element_text(face="bold",size=fs+1), plot.subtitle=element_text(size=fs-1,colour="grey40"))
arch <- function(p){a<-file.path(fig,"_archive",paste0(tools::file_path_sans_ext(basename(p)),"_superseded.",tools::file_ext(p))); if(file.exists(p)&&!file.exists(a)) file.copy(p,a)}

cat("=== PART 1: two-cohort Figure 4a ===\n")
v4 <- readRDS(file.path(pr,"aim2_scoring","malignant_cells_scored_v4.rds"))
mcm <- readRDS(file.path(pr,"aim2_scoring","malignant_cell_metadata.rds"))
tr <- unique(mcm[,c("pid","new_treatment")]); v4$trt <- tr$new_treatment[match(v4$pid,tr$pid)]
x <- v4[v4$trt=="Untreated" & v4$axis2_death_competence!=0,]; x$cyto <- grepl("cytostatic",x$cell_state)
m <- readRDS(file.path(pr,"aim2_scoring","malignant_count_matrix.rds"))
idx <- match(x$cell_id, rownames(m)); x$tgm2 <- as.numeric(m[idx,"TGM2"]); rm(m); invisible(gc())
hc <- tapply(x$tgm2[x$cyto],x$pid[x$cyto],mean); hd <- tapply(x$tgm2[!x$cyto],x$pid[!x$cyto],mean)
hk <- intersect(names(hc),names(hd))
st <- read.csv(file.path(base,"results","aim2","steele_tgm2_replication.csv"), stringsAsFactors=FALSE)
d1 <- data.frame(cohort=sprintf("Hwang et al.\n(%d patients)",length(hk)), id=hk, dc=hd[hk], cy=hc[hk])
d2 <- data.frame(cohort=sprintf("Steele et al.\n(%d samples)",nrow(st)), id=st$sample, dc=st$death_competent, cy=st$cytostatic)
dd <- rbind(d1,d2)
lg <- rbind(data.frame(cohort=dd$cohort,id=dd$id,state="Death-competent",value=dd$dc),
            data.frame(cohort=dd$cohort,id=dd$id,state="Cytostatic-prone",value=dd$cy))
lg$state <- factor(lg$state, levels=c("Death-competent","Cytostatic-prone"))
pv <- sapply(split(dd,dd$cohort), function(z) suppressWarnings(wilcox.test(z$cy,z$dc,paired=TRUE))$p.value)
ann <- data.frame(cohort=names(pv), lab=sprintf("higher in %d of %d\np = %.3f",
        sapply(split(dd,dd$cohort),function(z) sum(z$cy>z$dc)), sapply(split(dd,dd$cohort),nrow), pv))
for (n in names(pv)) cat(sprintf("  %-28s n=%d, higher in %d, p=%.4f\n", gsub("\n"," ",n), sum(dd$cohort==n), sum(dd$cy[dd$cohort==n]>dd$dc[dd$cohort==n]), pv[[n]]))
p_a <- ggplot(lg, aes(state,value)) + facet_wrap(~cohort, scales="free_y") +
  geom_boxplot(aes(group=state), width=.4, outlier.shape=NA, fill=NA, colour="grey55", linewidth=.3) +
  geom_line(aes(group=id), colour="grey75", linewidth=.3) +
  geom_point(aes(colour=state), size=1.3) +
  geom_text(data=ann, aes(x=1.5, y=Inf, label=lab), vjust=1.3, size=fs/ggplot2::.pt, inherit.aes=FALSE) +
  scale_colour_manual(values=c("Death-competent"="#0072B2","Cytostatic-prone"="#009E73"), guide="none") +
  scale_y_continuous(expand=expansion(mult=c(.05,.22))) +
  labs(title="TGM2 expression by malignant cell state, both cohorts",
       subtitle="Per-sample mean, cells with Axis 2 signal; lines join each sample", x=NULL, y="TGM2, mean log-normalised expression") + thm(fs)
f_a <- file.path(fig,"aim3","Fig3A_TGM2_paired_expression.png"); arch(f_a)
ggsave(f_a, p_a, dpi=dpi, width=cw, height=h_top, units="in"); cat("Saved panel a\n")
px <- function(z) as.integer(round(z*dpi))
put <- function(cv,img,x,y,w) image_composite(cv, image_scale(img,paste0(px(w),"x")), offset=paste0("+",px(x),"+",px(y)))
cv <- image_blank(px(W),px(H),color="white"); xr <- ml+cw+gx; yc <- mt+h_top+gy
cv <- put(cv, image_read(f_a), ml, mt, cw)
cv <- put(cv, image_read(file.path(fig,"aim3","Fig3B_TGM2_ITGA4_rank_scores.png")), xr, mt, cw)
cv <- put(cv, image_read(file.path(fig,"aim3","Fig3C_DepMap_essentiality.png")), ml, yc, w_c)
labs <- data.frame(l=c("a","b","c"), x=c(.03,xr-.19,.03), y=c(mt,mt,yc))
pg <- ggplot() + annotation_raster(as.raster(cv),xmin=0,xmax=W,ymin=0,ymax=H,interpolate=TRUE) +
  geom_text(data=labs, aes(x=x,y=H-y,label=l), hjust=0,vjust=1,fontface="bold",size=9/ggplot2::.pt) +
  coord_fixed(xlim=c(0,W),ylim=c(0,H),expand=FALSE) + theme_void() + theme(plot.margin=margin(0,0,0,0))
out <- file.path(fig,"multipanel","Figure4.pdf"); tmp <- sub("\\.pdf$","_new.pdf",out); arch(out)
ggsave(tmp,pg,width=W,height=H,units="in",device="pdf"); stopifnot(file.size(tmp)>1e5)
file.copy(tmp,out,overwrite=TRUE); file.remove(tmp); cat("Figure4.pdf rebuilt:", round(file.size(out)/1024),"KB\n")

cat("\n=== PART 2: permutation stability across cell subsamples ===\n")
sig <- readRDS(file.path(pr,"aim1_signatures","four_axis_signatures_v1.rds"))
mods <- sig[c("axis2_proapoptotic","axis2_executors","axis2_anti_apoptotic","axis2_cytostasis")]
m <- readRDS(file.path(pr,"aim2_scoring","malignant_count_matrix.rds"))
gxa <- Matrix::t(m); rm(m); invisible(gc()); gxa <- gxa[, intersect(x$cell_id, colnames(gxa)), drop=FALSE]
pn <- function(a,b){r<-log2((a+.001)/(b+.001)); r<-r[r!=0]; 100*mean(r<0)}
res <- list()
for (sd in c(7, 2026)) {
  set.seed(sd); cc <- sample(colnames(gxa), 8000); g <- gxa[,cc,drop=FALSE]
  det <- Matrix::rowMeans(g>0); pool <- names(det)[det>0]
  dec <- cut(det[pool], quantile(det[pool],0:10/10), include.lowest=TRUE, labels=FALSE); names(dec) <- pool
  og <- unlist(mods)
  pk <- function() unname(unlist(lapply(dec[intersect(mods[[3]],pool)], function(dd) sample(setdiff(pool[dec==dd],og),1))))
  nul <- numeric(0)
  for (b in 1:5) { f <- list()
    for (i in ((b-1)*100+1):(b*100)) for (k in 1:4) f[[sprintf("p%d_%d",i,k)]] <- pk()
    s2 <- ScoreSignatures_UCell(g, features=f, maxRank=1500, ncores=1)
    nul <- c(nul, sapply(((b-1)*100+1):(b*100), function(i){gg<-function(k) s2[,sprintf("p%d_%d_UCell",i,k)]; pn(gg(1)+gg(2),gg(3)+gg(4))}))
  }
  ob <- sapply(1:20, function(j){ set.seed(5000+j)
    eq <- list(a=sample(mods[[1]],5), b=sample(mods[[2]],5), c=mods[[3]], d=mods[[4]])
    s3 <- ScoreSignatures_UCell(g, features=eq, maxRank=1500, ncores=1)
    pn(s3[,"a_UCell"]+s3[,"b_UCell"], s3[,"c_UCell"]+s3[,"d_UCell"]) })
  ex <- sum(nul >= median(ob))
  cat(sprintf("seed %4d: null median %.2f%% (max %.2f) | observed median %.2f%% | %d of %d exceed | p = %.4f\n",
              sd, median(nul), max(nul), median(ob), ex, length(nul), (ex+1)/(length(nul)+1)))
  res[[as.character(sd)]] <- data.frame(seed=sd, null_median=median(nul), null_max=max(nul), obs_median=median(ob), exceed=ex, n=length(nul))
}
cat("\nSeed 42 (reported): null median 49.20%, observed median 87.90%, 4 of 1000, p = 0.0050\n")
write.csv(do.call(rbind,res), file.path(base,"results","aim2","perm_stability_subsamples.csv"), row.names=FALSE)
