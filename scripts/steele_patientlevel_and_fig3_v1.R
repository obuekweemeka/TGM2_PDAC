# ============================================================================
# steele_patientlevel_and_fig3_v1.R -- Save in: scripts/
# PART A: recompute the Steele TGM2 test at PATIENT level (11A and 11B pooled)
# PART B: regenerate Figure 3 panels a, b, d with subtitles that fit, and
#         rebuild Figure3.pdf
# ============================================================================
suppressPackageStartupMessages({library(Matrix); library(ggplot2); library(magick)})
base <- normalizePath('.'); pr <- file.path(base,"data","processed"); fig <- file.path(base,"figures")
dest <- file.path(base,"data","raw","GSE155698")
arch <- function(p){a<-file.path(fig,"_archive",paste0(tools::file_path_sans_ext(basename(p)),"_superseded.",tools::file_ext(p))); if(file.exists(p)&&!file.exists(a)) file.copy(p,a)}

cat("=== PART A: Steele TGM2 at patient level ===\n")
sc <- read.csv(file.path(pr,"aim2_validation","steele_cell_scores.csv"), stringsAsFactors=FALSE)
dirs <- unique(dirname(list.files(dest, pattern="^(matrix\\.mtx|barcodes\\.tsv)", recursive=TRUE, full.names=TRUE)))
pull <- function(d){ sid <- regmatches(d, regexpr("PDAC_TISSUE_[0-9]+[AB]?", d, ignore.case=TRUE)); if(!length(sid)) return(NULL)
  f <- function(p){h<-list.files(d,pattern=p,full.names=TRUE); if(length(h)) h[1] else NA}
  mf<-f("^matrix\\.mtx"); bf<-f("^barcodes\\.tsv"); gf<-f("^(features|genes)\\.tsv"); if(any(is.na(c(mf,bf,gf)))) return(NULL)
  g <- read.delim(gf, header=FALSE, stringsAsFactors=FALSE); sym <- if(ncol(g)>=2) g[[2]] else g[[1]]
  i <- which(toupper(sym)=="TGM2"); if(!length(i)) return(NULL)
  mm <- Matrix::readMM(if(grepl("\\.gz$",mf)) gzfile(mf) else mf); bc <- readLines(bf); tot <- Matrix::colSums(mm)
  data.frame(sample_id=sid, bc=sub("-\\d+$","",bc), tgm2=log1p(as.numeric(mm[i[1],])/pmax(tot,1)*1e4), stringsAsFactors=FALSE) }
got <- do.call(rbind, Filter(Negate(is.null), lapply(dirs, pull)))
sc$bc <- sub("-\\d+$","",sub("^PDAC_TISSUE_[0-9]+[AB]?_","",sc$cell_barcode))
mg <- merge(sc, got, by=c("sample_id","bc")); mg <- mg[mg$axis2_death!=0,]; mg$cyto <- mg$axis2_death<0
mg$patient <- sub("([0-9]+)[AB]$","\\1", mg$sample_id)          # 11A and 11B -> 11
cat("Samples:", length(unique(mg$sample_id)), "| patients after pooling:", length(unique(mg$patient)), "\n")
cat("Samples sharing a patient:\n"); print(table(mg$patient)[table(unique(mg[,c("sample_id","patient")])$patient)>1])
run <- function(unit){
  a <- tapply(mg$tgm2[mg$cyto], mg[[unit]][mg$cyto], mean); b <- tapply(mg$tgm2[!mg$cyto], mg[[unit]][!mg$cyto], mean)
  k <- intersect(names(a),names(b)); d <- a[k]-b[k]
  w <- suppressWarnings(wilcox.test(a[k],b[k],paired=TRUE)); tt <- t.test(d)
  cat(sprintf("%-8s n = %2d | higher in cytostatic-prone: %2d | mean diff %.3f (95%% CI %.3f to %.3f) | p = %.4f\n",
              unit, length(k), sum(d>0), mean(d), tt$conf.int[1], tt$conf.int[2], w$p.value))
  data.frame(unit=k, cytostatic=a[k], death_competent=b[k], difference=d) }
cat("\n"); s_samp <- run("sample_id"); s_pat <- run("patient")
cat("\nHwang for comparison: n = 18, 12 higher, 0.076 (0.024 to 0.128), p = 0.0047\n")
write.csv(s_pat, file.path(base,"results","aim2","steele_tgm2_replication_patientlevel.csv"), row.names=FALSE)
cat("Saved steele_tgm2_replication_patientlevel.csv\n")

cat("\n=== PART B: Figure 3 panels with subtitles that fit ===\n")
W<-8.27; H<-11.69; dpi<-600; ml<-.22; mr<-.06; mt<-.06; mb<-.06; gx<-.22; gy<-.12; fs<-9
cw <- (W-ml-mr-gx)/2; ch <- (H-mt-mb-gy)/2
thm <- function(fs) theme_classic(base_size=fs) + theme(plot.title=element_text(face="bold",size=fs+1), plot.subtitle=element_text(size=fs-1.5,colour="grey40"))
v4 <- readRDS(file.path(pr,"aim2_scoring","malignant_cells_scored_v4.rds"))
mcm <- readRDS(file.path(pr,"aim2_scoring","malignant_cell_metadata.rds"))
tr <- unique(mcm[,c("pid","new_treatment")]); v4$trt <- tr$new_treatment[match(v4$pid,tr$pid)]
u <- v4[v4$trt=="Untreated",]; x <- u[u$axis2_death_competence!=0,]; x$cyto <- grepl("cytostatic",x$cell_state)
p <- sort(unique(x$pid))
pm <- data.frame(pid=p, pct=100*as.vector(tapply(x$cyto,x$pid,mean)[p]), n=as.vector(table(x$pid)[p]),
                 basal=as.vector(tapply(u$Moffitt_basal,u$pid,mean)[p]), classical=as.vector(tapply(u$Moffitt_classical,u$pid,mean)[p]))
pm$subtype <- ifelse(pm$basal>pm$classical,"Basal-like","Classical")
pb <- pm[order(-pm$pct),]; pb$pid <- factor(pb$pid, levels=pb$pid)
p_a <- ggplot(pb, aes(pid,pct,fill=subtype)) + geom_col(width=.75) +
  geom_hline(yintercept=mean(pm$pct), linetype="dashed", colour="grey40") +
  scale_fill_manual(values=c("Basal-like"="#D55E00","Classical"="#0072B2"), name="Moffitt subtype") +
  scale_y_continuous(limits=c(0,100), expand=c(0,0)) +
  labs(title="Cytostatic-prone cells per patient",
       subtitle=sprintf("18 untreated patients; dashed = mean of patients (%.1f%%)", mean(pm$pct)),
       x="Patient", y="Cytostatic-prone cells (%)") + thm(fs) +
  theme(axis.text.x=element_text(angle=45,hjust=1), legend.position="bottom")
f_a <- file.path(fig,"aim2","Fig2B_per_patient_cytostatic.png"); arch(f_a); ggsave(f_a,p_a,dpi=dpi,width=cw,height=ch,units="in")
axn <- c(axis1_ras_balance="Axis 1\nRAS balance", axis2_death_competence="Axis 2\nDeath competence",
         axis3_hypoxia="Axis 3\nHypoxia-metabolic", axis4_innate_mimicry="Axis 4\nInnate immune mimicry")
cm <- cor(u[,names(axn)], method="spearman"); pr34 <- cor(u$axis3_hypoxia,u$axis4_innate_mimicry)
cd <- expand.grid(a=names(axn), b=names(axn)); cd$r <- cm[cbind(as.character(cd$a),as.character(cd$b))]
cd$a <- factor(axn[as.character(cd$a)],levels=axn); cd$b <- factor(axn[as.character(cd$b)],levels=axn)
p_b <- ggplot(cd, aes(a,b,fill=r)) + geom_tile(colour="white") + geom_text(aes(label=sprintf("%.3f",r)),size=2.6) +
  scale_fill_gradient2(low="#0072B2",mid="white",high="#D55E00",limits=c(-1,1),name="Spearman r") +
  labs(title="Pairwise axis correlations",
       subtitle=sprintf("Untreated cells (n = 52,999); Axis 3 vs 4 Pearson r = %.3f", pr34), x=NULL,y=NULL) +
  thm(fs) + theme(axis.text.x=element_text(angle=45,hjust=1), legend.position="bottom")
f_b <- file.path(fig,"aim2","Fig2C_axis_correlation_heatmap.png"); arch(f_b); ggsave(f_b,p_b,dpi=dpi,width=cw,height=ch,units="in")
stc <- read.csv(file.path(pr,"aim2_validation","steele_cell_scores.csv")); stc <- stc[stc$axis2_death!=0,]
se <- data.frame(id=names(tapply(stc$axis2_death<0,stc$sample_id,mean)),
                 pct=100*as.vector(tapply(stc$axis2_death<0,stc$sample_id,mean)),
                 n=as.vector(table(stc$sample_id)))
ed <- rbind(data.frame(atlas="Hwang et al. 2022\n(18 patients)", id=pm$pid, pct=pm$pct, n=pm$n),
            data.frame(atlas=sprintf("Steele et al. 2020\n(%d samples)",nrow(se)), id=se$id, pct=se$pct, n=se$n))
med <- aggregate(pct~atlas, ed, median); set.seed(7)
p_d <- ggplot(ed, aes(atlas,pct)) + geom_hline(yintercept=50,linetype="dashed",colour="grey60") +
  geom_jitter(aes(size=n,colour=atlas), width=.15, alpha=.8) +
  geom_crossbar(data=med, aes(ymin=pct,ymax=pct), width=.5, colour="black", linewidth=.4) +
  scale_colour_manual(values=c("#D55E00","#0072B2"), guide="none") + scale_size_continuous(range=c(.8,4), name="Cells") +
  scale_y_continuous(limits=c(0,100)) +
  labs(title="Replication in an independent PDAC atlas", subtitle="Crossbar = median; cells with Axis 2 signal",
       x=NULL, y="Cytostatic-prone cells (%)") + thm(fs) + theme(legend.position="bottom")
f_d <- file.path(fig,"aim2","Fig2E_Steele_replication.png"); arch(f_d); ggsave(f_d,p_d,dpi=dpi,width=cw,height=ch,units="in")
cat("Panels a, b, d regenerated\n")
px <- function(z) as.integer(round(z*dpi))
put <- function(cv,img,x,y,w) image_composite(cv, image_scale(img,paste0(px(w),"x")), offset=paste0("+",px(x),"+",px(y)))
cv <- image_blank(px(W),px(H),color="white"); xr <- ml+cw+gx; yb <- mt+ch+gy
cv <- put(cv, image_read(f_a), ml, mt, cw); cv <- put(cv, image_read(f_b), xr, mt, cw)
cv <- put(cv, image_read(file.path(fig,"aim2","Fig2D_moffitt_alignment.png")), ml, yb, cw)
cv <- put(cv, image_read(f_d), xr, yb, cw)
labs <- data.frame(l=letters[1:4], x=c(.03,xr-.19,.03,xr-.19), y=c(mt,mt,yb,yb))
pg <- ggplot() + annotation_raster(as.raster(cv),xmin=0,xmax=W,ymin=0,ymax=H,interpolate=TRUE) +
  geom_text(data=labs, aes(x=x,y=H-y,label=l), hjust=0,vjust=1,fontface="bold",size=9/ggplot2::.pt) +
  coord_fixed(xlim=c(0,W),ylim=c(0,H),expand=FALSE) + theme_void() + theme(plot.margin=margin(0,0,0,0))
out <- file.path(fig,"multipanel","Figure3.pdf"); tmp <- sub("\\.pdf$","_new.pdf",out); arch(out)
ggsave(tmp,pg,width=W,height=H,units="in",device="pdf"); stopifnot(file.size(tmp)>1e5)
file.copy(tmp,out,overwrite=TRUE); file.remove(tmp); cat("Figure3.pdf rebuilt:", round(file.size(out)/1024),"KB\n")
