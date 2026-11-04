# ============================================================================
# axis2_rename_all_v1.R -- Save in: scripts/
# Renames Axis 2 on every panel that shows it ("death competence" ->
# "apoptosis/arrest contrast"), removes the unsupported subtitle from the
# purity panel, and rebuilds Figures 3, 4, 5 and Supplementary Figure 5.
# Panels are PNG; multipanels are PDF only. Replaced files go to _archive/.
# ============================================================================
suppressPackageStartupMessages({library(Matrix); library(ggplot2); library(patchwork); library(magick); library(survival); library(survminer)})
base <- here::here(); pr <- file.path(base,"data","processed"); fig <- file.path(base,"figures")
W<-8.27; H<-11.69; dpi<-600; ml<-.22; mr<-.06; mt<-.06; mb<-.06; gx<-.22; gy<-.12
cw <- (W-ml-mr-gx)/2; wf <- W-ml-mr; xr <- ml+cw+gx
AX <- "Axis 2 apoptosis/arrest contrast"
arch <- function(p){a<-file.path(fig,"_archive",paste0(tools::file_path_sans_ext(basename(p)),"_superseded.",tools::file_ext(p))); if(file.exists(p)&&!file.exists(a)) file.copy(p,a)}
sv <- function(p,path,w,h){arch(path); ggsave(path,p,dpi=dpi,width=w,height=h,units="in"); cat("  saved",basename(path),round(file.size(path)/1024),"KB\n")}
px <- function(z) as.integer(round(z*dpi))
put <- function(cv,img,x,y,w) image_composite(cv, image_scale(img,paste0(px(w),"x")), offset=paste0("+",px(x),"+",px(y)))
page <- function(cv,labs,name){ p <- ggplot() + annotation_raster(as.raster(cv),xmin=0,xmax=W,ymin=0,ymax=H,interpolate=TRUE) +
    geom_text(data=labs,aes(x=x,y=H-y,label=l),hjust=0,vjust=1,fontface="bold",size=9/ggplot2::.pt) +
    coord_fixed(xlim=c(0,W),ylim=c(0,H),expand=FALSE) + theme_void() + theme(plot.margin=margin(0,0,0,0))
  out <- file.path(fig,"multipanel",name); tmp <- sub("\\.pdf$","_new.pdf",out); arch(out)
  ggsave(tmp,p,width=W,height=H,units="in",device="pdf"); stopifnot(file.size(tmp)>1e5)
  file.copy(tmp,out,overwrite=TRUE); file.remove(tmp); cat("  rebuilt",name,round(file.size(out)/1024),"KB\n") }
thm <- function(fs) theme_classic(base_size=fs) + theme(plot.title=element_text(face="bold",size=fs+1), plot.subtitle=element_text(size=fs-1.5,colour="grey40"))

cat("=== 0. Locate ABSOLUTE purity values ===\n")
tc <- readRDS(file.path(pr,"aim4_bulk","tcga_paad_analysis_df.rds"))
pur <- grep("purity|ABSOLUTE|absolute", names(tc), ignore.case=TRUE, value=TRUE)
cat("purity columns in the TCGA object:", if(length(pur)) paste(pur,collapse=", ") else "none", "\n")
cand <- list.files(base, pattern="purity|ABSOLUTE", recursive=TRUE, full.names=TRUE, ignore.case=TRUE)
cat("files mentioning purity:", if(length(cand)) paste(sub(paste0(base,"/"),"",cand),collapse=" | ") else "none", "\n")
HAVE_PURITY <- length(pur)>0
if (!HAVE_PURITY && length(cand)) { for (f in cand[grepl("\\.(csv|tsv|txt|rds)$",cand)]) {
  d <- tryCatch(if(grepl("\\.rds$",f)) readRDS(f) else read.delim(f, nrows=5), error=function(e) NULL)
  if (!is.null(d)) cat("  ", basename(f), "columns:", paste(head(names(d),12),collapse=", "), "\n") } }

cat("\n=== 1. Figure 3b: axis correlation heatmap ===\n")
v4 <- readRDS(file.path(pr,"aim2_scoring","malignant_cells_scored_v4.rds"))
mcm <- readRDS(file.path(pr,"aim2_scoring","malignant_cell_metadata.rds"))
tr <- unique(mcm[,c("pid","new_treatment")]); v4$trt <- tr$new_treatment[match(v4$pid,tr$pid)]
u <- v4[v4$trt=="Untreated",]
axn <- c(axis1_ras_balance="Axis 1\nRAS balance", axis2_death_competence="Axis 2\napoptosis/arrest",
         axis3_hypoxia="Axis 3\nHypoxia-metabolic", axis4_innate_mimicry="Axis 4\nInnate immune mimicry")
cm <- cor(u[,names(axn)],method="spearman"); pr34 <- cor(u$axis3_hypoxia,u$axis4_innate_mimicry)
cd <- expand.grid(a=names(axn),b=names(axn)); cd$r <- cm[cbind(as.character(cd$a),as.character(cd$b))]
cd$a <- factor(axn[as.character(cd$a)],levels=axn); cd$b <- factor(axn[as.character(cd$b)],levels=axn)
p3b <- ggplot(cd,aes(a,b,fill=r)) + geom_tile(colour="white") + geom_text(aes(label=sprintf("%.3f",r)),size=2.6) +
  scale_fill_gradient2(low="#0072B2",mid="white",high="#D55E00",limits=c(-1,1),name="Spearman r") +
  labs(title="Pairwise axis correlations",
       subtitle=sprintf("Untreated cells (n = 52,999); Axis 3 vs 4 Pearson r = %.3f",pr34),x=NULL,y=NULL) +
  thm(9) + theme(axis.text.x=element_text(angle=45,hjust=1), legend.position="bottom")
ch <- (H-mt-mb-gy)/2; sv(p3b, file.path(fig,"aim2","Fig2C_axis_correlation_heatmap.png"), cw, ch)
cv <- image_blank(px(W),px(H),color="white"); yb <- mt+ch+gy
cv <- put(cv,image_read(file.path(fig,"aim2","Fig2B_per_patient_cytostatic.png")),ml,mt,cw)
cv <- put(cv,image_read(file.path(fig,"aim2","Fig2C_axis_correlation_heatmap.png")),xr,mt,cw)
cv <- put(cv,image_read(file.path(fig,"aim2","Fig2D_moffitt_alignment.png")),ml,yb,cw)
cv <- put(cv,image_read(file.path(fig,"aim2","Fig2E_Steele_replication.png")),xr,yb,cw)
page(cv, data.frame(l=letters[1:4],x=c(.03,xr-.19,.03,xr-.19),y=c(mt,mt,yb,yb)), "Figure3.pdf")

cat("\n=== 2. Figure 4c: DepMap block ===\n")
de <- readRDS(file.path(pr,"aim1_signatures","DepMap_PDAC_expression.rds"))
d4 <- readRDS(file.path(pr,"aim4_bulk","depmap_pdac_analysis_df.rds"))
v44 <- intersect(de$ModelID[!duplicated(de$ModelID)], d4$model_id); d <- d4[d4$model_id %in% v44,]
d$grp <- ifelse(d$axis2_death < median(d$axis2_death),"High cytostatic","Low cytostatic"); hi <- d$grp=="High cytostatic"
pct <- sapply(c(tgm2="tgm2",nt5e="nt5e",kras="kras"), function(g) round(100*mean(d[[g]]>0),1))
leg <- c("High cytostatic"=paste0("High cytostatic (n = ",sum(hi),")"),"Low cytostatic"=paste0("Low cytostatic (n = ",sum(!hi),")"))
pnl <- function(g,ttl) ggplot(d,aes(axis2_death,.data[[g]],colour=grp)) + geom_point(size=1.3,alpha=.85) +
  geom_smooth(method="lm",se=TRUE,colour="grey30",linewidth=.4,formula=y~x) +
  geom_hline(yintercept=0,linetype="dashed",colour="grey50",linewidth=.3) +
  scale_colour_manual(values=c("High cytostatic"="#E69F00","Low cytostatic"="#0072B2"),labels=leg,name=NULL) +
  labs(title=ttl, subtitle=paste0("n = 44 PDAC lines; Chronos > 0: ",pct[[g]],"%"), x=AX, y="Chronos effect score") +
  thm(6.5) + theme(legend.text=element_text(size=6.5))
h_top <- 4.2; h_c <- H-mt-h_top-gy-mb
blk <- ((pnl("tgm2","TGM2") | pnl("nt5e","NT5E")) / (pnl("kras","KRAS (positive control)") | guide_area())) + plot_layout(guides="collect")
sv(blk, file.path(fig,"aim3","Fig3C_DepMap_essentiality.png"), wf, h_c)
cv <- image_blank(px(W),px(H),color="white"); yc <- mt+h_top+gy
cv <- put(cv,image_read(file.path(fig,"aim3","Fig3A_TGM2_paired_expression.png")),ml,mt,cw)
cv <- put(cv,image_read(file.path(fig,"aim3","Fig3B_TGM2_ITGA4_rank_scores.png")),xr,mt,cw)
cv <- put(cv,image_read(file.path(fig,"aim3","Fig3C_DepMap_essentiality.png")),ml,yc,wf)
page(cv, data.frame(l=c("a","b","c"),x=c(.03,xr-.19,.03),y=c(mt,mt,yc)), "Figure4.pdf")

cat("\n=== 3. Figure 5c (KM) and 5d (purity) ===\n")
tc$grp <- factor(tc$cytostatic_group, levels=c("High cytostatic","Low cytostatic"))
fitk <- survfit(Surv(os_days,os_event)~grp, data=tc); cx <- coxph(Surv(os_days,os_event)~grp,data=tc)
hr <- exp(coef(cx)); ci <- exp(confint(cx)); lr <- survdiff(Surv(os_days,os_event)~grp,data=tc); plr <- 1-pchisq(lr$chisq,1)
g <- ggsurvplot(fitk, data=tc, conf.int=TRUE, risk.table=TRUE, pval=FALSE,
                palette=c("#D55E00","#0072B2"), legend.title="Cell state group",
                legend.labs=c("High cytostatic","Low cytostatic"),
                xlab="Overall survival (days)", ylab="Survival probability",
                title=paste("OS by", AX), ggtheme=theme_classic(base_size=6.5),
                tables.theme=theme_cleantable(base_size=6.5), fontsize=6.5/ggplot2::.pt, censor.size=2, size=.4)
pk <- g$plot + labs(subtitle=sprintf("TCGA PAAD, n = %d, median split on Axis 2", nrow(tc))) +
  annotate("text", x=1300, y=.97, hjust=0, vjust=1, size=6.5/ggplot2::.pt,
           label=sprintf("Log-rank p = %.3f\nHR = %.2f (95%% CI %.2f-%.2f)", plr, hr, ci[1], ci[2])) +
  theme(plot.title=element_text(face="bold",size=7.7), plot.subtitle=element_text(size=5.5,colour="grey40"),
        legend.text=element_text(size=6.5), legend.title=element_text(size=6.5))
km <- pk / (g$table + labs(title="Number at risk") + theme(plot.title=element_text(size=6.5))) + plot_layout(heights=c(3,1))
h1 <- 5.4; sv(km, file.path(fig,"aim4","Fig4A_KM_cytostatic.png"), cw, h1)
if (HAVE_PURITY) {
  pc <- pur[1]
  dp <- tc[!is.na(tc[[pc]]),]
  ct <- suppressWarnings(cor.test(dp[[pc]], dp$axis2_death_z, method="spearman"))
  p5d <- ggplot(dp, aes(.data[[pc]], axis2_death_z, colour=moffitt)) + geom_point(size=1.3,alpha=.85) +
    geom_smooth(method="lm",se=TRUE,colour="grey30",linewidth=.4,formula=y~x) +
    scale_colour_manual(values=c("#D55E00","#0072B2"), name="Moffitt subtype") +
    labs(title=paste(AX, "tracks tumour purity"),
         subtitle=sprintf("TCGA PAAD, n = %d; Spearman r = %.3f, p = %.2e", nrow(dp), ct$estimate, ct$p.value),
         x="Tumour purity (ABSOLUTE)", y=paste(AX,"(scaled)")) + thm(6.5) + theme(legend.position="bottom")
  sv(p5d, file.path(fig,"aim4","Fig4C_axis2_purity_TCGA.png"), cw, h1)
  cat("  purity panel regenerated\n")
} else cat("  PURITY VALUES NOT FOUND: panel 5d left as is; its baked-in subtitle still reads\n  'explains null Moffitt association in bulk data', which the text does not claim.\n")
cv <- image_blank(px(W),px(H),color="white"); y5 <- mt+h1+gy
cv <- put(cv,image_read(file.path(fig,"aim3","Fig5A_distance_heatmap.png")),ml,mt,cw)
cv <- put(cv,image_read(file.path(fig,"aim3","Fig5B_distance_correlations.png")),xr,mt,cw)
cv <- put(cv,image_read(file.path(fig,"aim4","Fig4A_KM_cytostatic.png")),ml,y5,cw)
cv <- put(cv,image_read(file.path(fig,"aim4","Fig4C_axis2_purity_TCGA.png")),xr,y5,cw)
page(cv, data.frame(l=letters[1:4],x=c(.03,xr-.19,.03,xr-.19),y=c(mt,mt,y5,y5)), "Figure5.pdf")

cat("\n=== 4. Supplementary Figure 5 ===\n")
x <- u[u$axis2_death_competence!=0,]; x$cyto <- grepl("cytostatic",x$cell_state)
p <- sort(unique(x$pid)); ad <- x$axis2_death_competence; tmed <- median(ad); t33 <- unname(quantile(ad,.33))
s5 <- do.call(rbind, lapply(p, function(q){ a <- ad[x$pid==q]
  data.frame(pid=q, strategy=c("Pre-specified (Axis 2 < 0)","Median of informative cells","33rd percentile of informative cells"),
             pct=100*c(mean(a<0),mean(a<tmed),mean(a<t33))) }))
s5$strategy <- factor(s5$strategy, levels=unique(s5$strategy))
ord <- 100*tapply(ad<0, x$pid, mean)[p]; s5$pid <- factor(s5$pid, levels=p[order(ord)])
ha <- 5.4; hb <- H-mt-gy-ha-mb
p5a <- ggplot(s5,aes(pid,pct,fill=strategy)) + geom_col(position=position_dodge(.8),width=.75) +
  geom_hline(yintercept=50,linetype="dashed",colour="grey40") +
  scale_fill_manual(values=c("#0072B2","#E69F00","#D55E00"),name="Threshold strategy") +
  scale_y_continuous(limits=c(0,100),expand=c(0,0)) +
  labs(title="Per-patient cytostatic-prone proportion under three thresholds",
       subtitle=sprintf("Untreated patients, cells with Axis 2 signal; median %.2f, 33rd percentile %.2f",tmed,t33),
       x="Patient", y="Cytostatic-prone cells (%)") + thm(10) + theme(axis.text.x=element_text(angle=45,hjust=1))
sv(p5a, file.path(fig,"aim2","FigS5_classifier_sensitivity.png"), wf, ha)
sw <- data.frame(t=seq(-10,10,.1)); sw$pct <- sapply(sw$t, function(z) 100*mean(ad<z)); at0 <- 100*mean(ad<0)
p5b <- ggplot(sw,aes(t,pct)) + geom_line(colour="#0072B2",linewidth=.8) +
  geom_vline(xintercept=0,linetype="dotted",colour="#0072B2") + geom_hline(yintercept=50,linetype="dashed",colour="grey40") +
  annotate("text",x=.3,y=25,hjust=0,size=3.2,colour="#0072B2",label=sprintf("Pre-specified threshold (0)\n%.1f%% cytostatic-prone",at0)) +
  annotate("text",x=-9.8,y=53,hjust=0,size=3.2,colour="grey40",label="50%") + scale_y_continuous(limits=c(0,100)) +
  labs(title="Cytostatic-prone proportion across thresholds",
       subtitle="Untreated cells with Axis 2 signal (n = 45,174)",
       x=paste(AX,"threshold"), y="Cells classified as cytostatic-prone (%)") + thm(10)
sv(p5b, file.path(fig,"aim2","FigS5b_classifier_sweep.png"), wf, hb)
cv <- image_blank(px(W),px(H),color="white")
cv <- put(cv,image_read(file.path(fig,"aim2","FigS5_classifier_sensitivity.png")),ml,mt,wf)
cv <- put(cv,image_read(file.path(fig,"aim2","FigS5b_classifier_sweep.png")),ml,mt+ha+gy,wf)
page(cv, data.frame(l=c("a","b"),x=c(.03,.03),y=c(mt,mt+ha+gy)), "FigureS5.pdf")
cat("\nDone. Figures 3, 4, 5 and S5 rebuilt.\n")
