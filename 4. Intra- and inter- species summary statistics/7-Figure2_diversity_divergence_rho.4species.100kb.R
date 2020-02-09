library(circlize)
library(data.table)
library(RColorBrewer)
library(dplyr)
library(vioplot)
#install.packages("beanplot")
library(beanplot)

##########################################################################
###This scripts have two main goals:(P.tremula,P.davidiana,P. tremuloides,P. trichocarpa)
#1. Construct the circlize distribution for gene density, pi, rho, fst, dxy
#2. Make pairwise correlations of various parameters:
###2.1 Between species: pi; rho; D
###2.2 Within species: pi vs. rho; pi vs. gene_density; pi vs. mu
###2.3 inter- and intra- species parameters correlations: pi vs.Fst; rho vs. Fst; D vs. Fst; pi vs. PBS; rho vs. PBS; D vs. PBS; pi vs dxy; rho vs. dxy; D vs dxy
###2.4 interspecies divergence correlations: Fst vs Fst; dxy vs. dxy; PBS vs PBS; Fst vs. dxy

##########################################################################
##read in the datasets
setwd("~/Dropbox/davidiana_paper/data/pairwise_diversity_divergence_rho")

##make a function to add a column called "WIN" for all windows
win<-function(data){
  data$Win=paste(data$Chr,data$Pos,sep=":")
  return(data)
}
win2<-function(data){ ####for divergence estimates, the Chr is replaced by Pop
  data$Win=paste(data$Pop,data$Pos,sep=":")
  return(data)
}

##all summary file and pairwise divergence comparison
summary_100kb=fread("old_100kb/100kb.summary.txt",header=T)
summary_100kb_w=win(summary_100kb)
tra_dav=read.table("data/100kb/tremula_davidiana.all.window1e+05_step1e+05.fst.dxy.txt",header=T)
tra_dav_w=win2(tra_dav)
tra_trs=read.table("data/100kb/tremula_tremuloides.all.window1e+05_step1e+05.fst.dxy.txt",header=T)
tra_trs_w=win2(tra_trs)
tra_tri=read.table("data/100kb/tremula_trichocarpa.all.window1e+05_step1e+05.fst.dxy.txt",header=T)
tra_tri_w=win2(tra_tri)
trs_tri=read.table("data/100kb/tremuloides_trichocarpa.all.window1e+05_step1e+05.fst.dxy.txt",header=T)
trs_tri_w=win2(trs_tri)
dav_tri=read.table("data/100kb/davidiana_trichocarpa.all.window1e+05_step1e+05.fst.dxy.txt",header=T)
dav_tri_w=win2(dav_tri)
dav_trs=read.table("data/100kb/tremuloides_davidiana.all.window1e+05_step1e+05.fst.dxy.txt",header=T)
dav_trs_w=win2(dav_trs)
##diversity estimates
tra=read.table("data/100kb/tremula_all.1e+05bp1e+05bp.thetas.txt",header=T)
tra_w=win(tra)
dav=read.table("data/100kb/davidiana_all.1e+05bp1e+05bp.thetas.txt",header=T)
dav_w=win(dav)
trs=read.table("data/100kb/tremuloides_all.1e+05bp1e+05bp.thetas.txt",header=T)
trs_w=win(trs)
tri=read.table("data/100kb/trichocarpa_all.1e+05bp1e+05bp.thetas.txt",header=T)
tri_w=win(tri)
##recombination rates estimated by ldhelmet
tra_r=read.table("ldhelmet/tremula.window1e+05.ldhelmet.summary.txt",header=T)
dav_r=read.table("ldhelmet/davidiana.window1e+05.ldhelmet.summary.txt",header=T)
trs_r=read.table("ldhelmet/tremuloides.window1e+05.ldhelmet.summary.txt",header=T)
tri_r=read.table("ldhelmet/trichocarpa.window1e+05.ldhelmet.summary.txt",header=T)
#########needs to remove those regions with too few SNPs
tra_rho=filter(tra_r,N_SNP>50)
tra_rho_w=win(tra_rho)
dav_rho=filter(dav_r,N_SNP>50)
dav_rho_w=win(dav_rho)
trs_rho=filter(trs_r,N_SNP>50)
trs_rho_w=win(trs_rho)
tri_rho=filter(tri_r,N_SNP>50)
tri_rho_w=win(tri_rho)

##########################################################################
###extract the common windows for all datasets
common_win=Reduce(intersect,list(summary_100kb_w$Win,
                                 tra_w$Win,
                                 dav_w$Win,
                                 trs_w$Win,
                                 tri_w$Win,
                                 tra_rho_w$Win,
                                 dav_rho_w$Win,
                                 trs_rho_w$Win,
                                 tri_rho_w$Win,
                                 tra_dav_w$Win,
                                 tra_trs_w$Win,
                                 tra_tri_w$Win,
                                 trs_tri_w$Win,
                                 dav_tri_w$Win,
                                 dav_trs_w$Win))
                  
summary_100kb_w_n=summary_100kb_w[which(summary_100kb_w$Win %in% common_win),]
tra_w_n=tra_w[which(tra_w$Win %in% common_win),]
dav_w_n=dav_w[which(dav_w$Win %in% common_win),]
trs_w_n=trs_w[which(trs_w$Win %in% common_win),]
tri_w_n=tri_w[which(tri_w$Win %in% common_win),]
tra_rho_w_n=tra_rho_w[which(tra_rho_w$Win %in% common_win),]
dav_rho_w_n=dav_rho_w[which(dav_rho_w$Win %in% common_win),]
trs_rho_w_n=trs_rho_w[which(trs_rho_w$Win %in% common_win),]
tri_rho_w_n=tri_rho_w[which(tri_rho_w$Win %in% common_win),]
tra_dav_w_n=tra_dav_w[which(tra_dav_w$Win %in% common_win),]
tra_trs_w_n=tra_trs_w[which(tra_trs_w$Win %in% common_win),]
tra_tri_w_n=tra_tri_w[which(tra_tri_w$Win %in% common_win),]
trs_tri_w_n=trs_tri_w[which(trs_tri_w$Win %in% common_win),]
dav_tri_w_n=dav_tri_w[which(dav_tri_w$Win %in% common_win),]
dav_trs_w_n=dav_trs_w[which(dav_trs_w$Win %in% common_win),]


##########################################################################
###output the summary table for the parameters that are most interested

summary_out=data.frame(cbind(summary_100kb_w_n$Chr,summary_100kb_w_n$Pos,summary_100kb_w_n$GC,summary_100kb_w_n$Gene_num,summary_100kb_w_n$Gene_prop,summary_100kb_w_n$Coding_prop,
                             tra_w_n$numSites,tra_w_n$tW.norm,tra_w_n$tP.norm,tra_w_n$tajD,
                             dav_w_n$numSites,dav_w_n$tW.norm,dav_w_n$tP.norm,dav_w_n$tajD,
                             trs_w_n$numSites,trs_w_n$tW.norm,trs_w_n$tP.norm,trs_w_n$tajD,
                             tri_w_n$numSites,tri_w_n$tW.norm,tri_w_n$tP.norm,tri_w_n$tajD,
                             tra_rho_w_n$N_SNP,tra_rho_w_n$rho_mean,tra_rho_w_n$rho_p0.5,
                             dav_rho_w_n$N_SNP,dav_rho_w_n$rho_mean,dav_rho_w_n$rho_p0.5,
                             trs_rho_w_n$N_SNP,trs_rho_w_n$rho_mean,trs_rho_w_n$rho_p0.5,
                             tri_rho_w_n$N_SNP,tri_rho_w_n$rho_mean,tri_rho_w_n$rho_p0.5,
                             tra_dav_w_n$fst,tra_dav_w_n$fixed,tra_dav_w_n$dxy,
                             tra_trs_w_n$fst,tra_trs_w_n$fixed,tra_trs_w_n$dxy,
                             dav_trs_w_n$fst,dav_trs_w_n$fixed,dav_trs_w_n$dxy,
                             tra_tri_w_n$fst,tra_tri_w_n$fixed,tra_tri_w_n$dxy,
                             dav_tri_w_n$fst,dav_tri_w_n$fixed,dav_tri_w_n$dxy,
                             trs_tri_w_n$fst,trs_tri_w_n$fixed,trs_tri_w_n$dxy
))

names(summary_out)=c("Chr","Pos","GC","Gene_num","Gene_prop","Coding_prop",
                     "tra_nsites","tra_tW","tra_tP","tra_tajD",
                     "dav_nsites","dav_tW","dav_tP","dav_tajD",
                     "trs_nsites","trs_tW","trs_tP","trs_tajD",
                     "tri_nsites","tri_tW","tri_tP","tri_tajD",
                     "tra_rho_nSNP","tra_rho_mean","tra_rho_median",
                     "dav_rho_nSNP","dav_rho_mean","dav_rho_median",
                     "trs_rho_nSNP","trs_rho_mean","trs_rho_median",
                     "tri_rho_nSNP","tri_rho_mean","tri_rho_median",
                     "tra_dav_fst","tra_dav_df","tra_dav_dxy",
                     "tra_trs_fst","tra_trs_df","tra_trs_dxy",
                     "dav_trs_fst","dav_trs_df","dav_trs_dxy",
                     "tra_tri_fst","tra_tri_df","tra_tri_dxy",
                     "dav_tri_fst","dav_tri_df","dav_tri_dxy",
                     "trs_tri_fst","trs_tri_df","trs_tri_dxy")


write.table(summary_out,file="/Users/Jing/Dropbox/davidiana_paper/data/linked_selection/summary.div.diverg.rho.100kb.txt",sep="\t", quote=F, row.names=F, col.names=T)


##########################################################################
####make the circlize plots

###group fst
tra_dav_fst<-tra_dav_w_n %>% select(Pop,Pos,fst) %>% mutate(start=Pos-1) %>% select(Pop,start,Pos,fst)
names(tra_dav_fst)=c("chr","start","end","value")
tra_trs_fst=tra_trs_w_n %>% select(Pop,Pos,fst) %>% mutate(start=Pos-1) %>% select(Pop,start,Pos,fst)
names(tra_trs_fst)=c("chr","start","end","value")
tra_tri_fst=tra_tri_w_n %>% select(Pop,Pos,fst) %>% mutate(start=Pos-1) %>% select(Pop,start,Pos,fst)
names(tra_tri_fst)=c("chr","start","end","value")
trs_tri_fst=trs_tri_w_n %>% select(Pop,Pos,fst) %>% mutate(start=Pos-1) %>% select(Pop,start,Pos,fst)
names(trs_tri_fst)=c("chr","start","end","value")
dav_tri_fst=dav_tri_w_n %>% select(Pop,Pos,fst) %>% mutate(start=Pos-1) %>% select(Pop,start,Pos,fst)
names(dav_tri_fst)=c("chr","start","end","value")
dav_trs_fst=dav_trs_w_n %>% select(Pop,Pos,fst) %>% mutate(start=Pos-1) %>% select(Pop,start,Pos,fst)
names(dav_trs_fst)=c("chr","start","end","value")

fst_list=list(tra_trs_fst[!is.na(tra_trs_fst$value),],tra_tri_fst[!is.na(tra_tri_fst$value),],trs_tri_fst[!is.na(trs_tri_fst$value),],dav_trs_fst[!is.na(dav_trs_fst$value),],tra_dav_fst[!is.na(tra_dav_fst$value),],dav_tri_fst[!is.na(dav_tri_fst$value),])

###group dxy
tra_dav_dxy=tra_dav_w_n %>% select(Pop,Pos,dxy) %>% mutate(start=Pos-1) %>% select(Pop,start,Pos,dxy)
names(tra_dav_dxy)=c("chr","start","end","value")
tra_trs_dxy=tra_trs_w_n %>% select(Pop,Pos,dxy) %>% mutate(start=Pos-1) %>% select(Pop,start,Pos,dxy)
names(tra_trs_dxy)=c("chr","start","end","value")
tra_tri_dxy=tra_tri_w_n %>% select(Pop,Pos,dxy) %>% mutate(start=Pos-1) %>% select(Pop,start,Pos,dxy)
names(tra_tri_dxy)=c("chr","start","end","value")
trs_tri_dxy=trs_tri_w_n %>% select(Pop,Pos,dxy) %>% mutate(start=Pos-1) %>% select(Pop,start,Pos,dxy)
names(trs_tri_dxy)=c("chr","start","end","value")
dav_tri_dxy=dav_tri_w_n %>% select(Pop,Pos,dxy) %>% mutate(start=Pos-1) %>% select(Pop,start,Pos,dxy)
names(dav_tri_dxy)=c("chr","start","end","value")
dav_trs_dxy=dav_trs_w_n %>% select(Pop,Pos,dxy) %>% mutate(start=Pos-1) %>% select(Pop,start,Pos,dxy)
names(dav_trs_dxy)=c("chr","start","end","value")

dxy_list=list(tra_trs_dxy[!is.na(tra_trs_dxy$value),],tra_tri_dxy[!is.na(tra_tri_dxy$value),],trs_tri_dxy[!is.na(trs_tri_dxy$value),],dav_trs_dxy[!is.na(dav_trs_dxy$value),],tra_dav_dxy[!is.na(tra_dav_dxy$value),],dav_tri_dxy[!is.na(dav_tri_dxy$value),])

### group diversity
tra_pi=tra_w_n %>% select(Chr,Pos,tP.norm) %>% mutate(start=Pos-1) %>% select(Chr,start,Pos,tP.norm)
names(tra_pi)=c("chr","start","end","value")
tra_pi=tra_pi[which(tra_pi$value<0.04),]
dav_pi=dav_w_n %>% select(Chr,Pos,tP.norm) %>% mutate(start=Pos-1) %>% select(Chr,start,Pos,tP.norm)
names(dav_pi)=c("chr","start","end","value")
dav_pi=dav_pi[which(dav_pi$value<0.04),]
trs_pi=trs_w_n %>% select(Chr,Pos,tP.norm) %>% mutate(start=Pos-1) %>% select(Chr,start,Pos,tP.norm)
names(trs_pi)=c("chr","start","end","value")
trs_pi=trs_pi[which(trs_pi$value<0.04),]
tri_pi=tri_w_n %>% select(Chr,Pos,tP.norm) %>% mutate(start=Pos-1) %>% select(Chr,start,Pos,tP.norm)
names(tri_pi)=c("chr","start","end","value")
tri_pi=tri_pi[which(tri_pi$value<0.04),]

pi_list=list(tra_pi[!is.na(tra_pi$value),],dav_pi[!is.na(dav_pi$value),],trs_pi[!is.na(trs_pi$value),],tri_pi[!is.na(tri_pi$value),])

###group da=dxy-(pi_a+pi_b)/2
#tra_dav_da=tra_dav_dxy
#tra_dav_da$da=tra_dav_dxy$value-((tra_pi$value+dav_pi$value)/2)


### group recombination rate
tra_recom=tra_rho_w_n %>% select(Chr,Pos,rho_mean) %>% mutate(start=Pos-1) %>% select(Chr,start,Pos,rho_mean)
names(tra_recom)=c("chr","start","end","value")
tra_recom=tra_recom[which(tra_recom$value<0.1),]

dav_recom=dav_rho_w_n %>% select(Chr,Pos,rho_mean) %>% mutate(start=Pos-1) %>% select(Chr,start,Pos,rho_mean)
names(dav_recom)=c("chr","start","end","value")
dav_recom=dav_recom[which(dav_recom$value<0.1),]

trs_recom=trs_rho_w_n %>% select(Chr,Pos,rho_mean) %>% mutate(start=Pos-1) %>% select(Chr,start,Pos,rho_mean)
names(trs_recom)=c("chr","start","end","value")
trs_recom=tra_recom[which(trs_recom$value<0.1),]

tri_recom=tri_rho_w_n %>% select(Chr,Pos,rho_mean) %>% mutate(start=Pos-1) %>% select(Chr,start,Pos,rho_mean)
names(tri_recom)=c("chr","start","end","value")
tri_recom=tra_recom[which(tri_recom$value<0.1),]

rho_list=list(tra_recom[!is.na(tra_recom$value),],dav_recom[!is.na(dav_recom$value),],trs_recom[!is.na(trs_recom$value),],tri_recom[!is.na(tri_recom$value),])

##coding density
gene_density=summary_100kb_w_n %>% select(Chr,Pos,Coding_prop) %>% mutate(start=Pos-50000,end=Pos+50000,coding=Coding_prop*1000) %>% select(Chr,start,end,coding)
gene_density=gene_density[!is.na(gene_density$coding),]
gene_density=gene_density[which(gene_density$coding!=0),]

#color
colfunc=colorRampPalette(c("blue","red"))
f=colfunc(329)
gene_density$col=f[round(gene_density$coding,0)]
names(gene_density)=c("chr","start","end","value","col")


###circlize plots
karyo=read.table("Potri_ideogram.txt",header=T)
karyo=karyo %>% mutate(start=0,end=chromEnd-1) %>% select(chrom,start,end)


#########################
##making the plot
colors <- brewer.pal(8,"Set2")

##ggplot default color values
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
n = 4
cols = gg_color_hue(n)

png("4species.rho_pi_fst_dxy.circlize.png",width=8,height=8,units='in',res=400)

par(mfrow=c(1,1), mar = c(1, 1, 1, 1))
circos.par(start.degree = 85,"gap.degree"=c(rep(1,18),6))
circos.initializeWithIdeogram(data.frame(karyo), track.height=0.03, ideogram.height = 0.01, chromosome.index = sort(unique(karyo$chrom)))

#gene density
#Coding density
circos.genomicTrackPlotRegion(gene_density,track.height=0.05,panel.fun = function(region, value,...){
  circos.genomicRect(region,value,col=value$col,border=NA,...)
})

#pi
colors_pi <- cols
circos.genomicTrackPlotRegion(pi_list,track.height=0.1,panel.fun = function(region, value,...){
  i=getI(...)
  circos.genomicLines(region,value,type="l",col=colors_pi[i],border=NA,ylim=c(0,0.035),...)
})

#rho
colors_rho <- cols
circos.genomicTrackPlotRegion(rho_list,track.height=0.1,panel.fun = function(region, value,...){
  i=getI(...)
  circos.genomicLines(region,value,type="l",col=colors_rho[i],border=NA,ylim=c(0,0.1),...)
})

##fst
colors <- brewer.pal(8,"Set2")
circos.genomicTrackPlotRegion(fst_list,track.height=0.1,panel.fun = function(region, value,...){
  i=getI(col=colors,...)
  circos.genomicLines(region,value,type="l",col=colors[i],border=NA,ylim=c(0,1),...)
})

#dxy
circos.genomicTrackPlotRegion(dxy_list,track.height=0.1,panel.fun = function(region, value,...){
  i=getI(col=colors,...)
  circos.genomicLines(region,value,type="l",col=colors[i],border=NA,ylim=c(0,0.1),...)
})

add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0),
              mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}

add_legend("top", legend=c("P.tra-P.trs", "P.tra-P.tri","P.trs-P.tri","P.dav-P.trs","P.tra-P.dav","P.dav-P.tri"), lwd=2,col=c(colors[1:6]),horiz=TRUE, bty='n', cex=0.6)
add_legend("bottom",legend=c("P.tra","P.dav","P.trs","P.tri"),lwd=2,col=c(colors_pi[1:4]),horiz=TRUE, bty='n', cex=0.6)

dev.off()

 ##########################################################################
##########################################################################

#2. Make pairwise correlations of various parameters:
###2.1 Between species: pi; rho; D
###2.2 Within species: pi vs. rho; pi vs. gene_density; pi vs. mu; rho vs. gene_density; rho vs. mu
###2.3 inter- and intra- species parameters correlations: pi vs.Fst; rho vs. Fst; D vs. Fst; pi vs dxy; rho vs. dxy; D vs dxy
###2.4 interspecies divergence correlations: Fst vs Fst; dxy vs. dxy; Fst vs. dxy


##define functions to calculate spearman's rho test
corrFunc <- function(var1, var2, data) {
  result = cor.test(data[,var1], data[,var2],method="spearman",na.rm=T)
  data.frame(var1, var2, result[c("estimate","p.value","statistic")], 
             stringsAsFactors=FALSE)
}

##define the function to summarize the correlation coefficients into a table
corrs_coef<-function(data){
corrs_new=data.frame()
vars=names(data)
comparision_n=ncol(data)  ####the number of comparision

for (i in 1:(comparision_n-1))
{
  if (i<(comparision_n-1)) {
    z=i+1
    for (j in z:comparision_n)
    {
      corrs= do.call(rbind, mapply(corrFunc, vars[i], vars[j], MoreArgs=list(data=data),SIMPLIFY=F))
      corrs_new=rbind(corrs_new,corrs)
    }
  }
  else{
    j=i+1
    corrs= do.call(rbind, mapply(corrFunc, vars[i], vars[j], MoreArgs=list(data=data),SIMPLIFY=F))
    corrs_new=rbind(corrs_new,corrs)
  }
}
return(corrs_new)
}

###function to summarize the correlations between two parameters:
#mean,median,min,max,sig,n_para
summary_cor=function(x){
  return(c(mean(x$estimate),median(x$estimate),min(x$estimate),max(x$estimate),length(which(x$p.value<0.001)),nrow(x)))
}

##########################################################################
###2.1 Between species: pi; rho; D

##pi
pi_table=data.frame(cbind(tra_w_n$tP.norm,dav_w_n$tP.norm,trs_w_n$tP.norm,tri_w_n$tP.norm))
names(pi_table)=c("tra_tP","dav_tP","trs_tP","tri_tP")
corrs_pi=corrs_coef(pi_table)
summary_cor(corrs_pi)

###rho
rho_table=data.frame(cbind(tra_rho_w_n$rho_mean,dav_rho_w_n$rho_mean,trs_rho_w_n$rho_mean,tri_rho_w_n$rho_mean))
names(rho_table)=c("tra_rho_mean","dav_rho_mean","trs_rho_mean","tri_rho_mean")
corrs_rho=corrs_coef(rho_table)
summary_cor(corrs_rho)

###tajD
tajD_table=data.frame(cbind(tra_w_n$tajD,dav_w_n$tajD,trs_w_n$tajD,tri_w_n$tajD))
names(tajD_table)=c("tra_tajD","dav_tajD","trs_tajD","tri_tajD")
corrs_tajD=corrs_coef(tajD_table)
summary_cor(corrs_tajD)

##########################################################################
###2.2 Within species: pi vs. rho; pi vs. gene_density; pi vs. mu; rho vs. gene_density; rho vs. mu
data_total=data.frame(cbind(as.character(summary_100kb_w_n$Chr),
                            summary_100kb_w_n$Pos,
                            summary_100kb_w_n$GC,
                            summary_100kb_w_n$Gene_num,
                            summary_100kb_w_n$Gene_prop,
                            summary_100kb_w_n$Coding_prop,
                            summary_100kb_w_n$tremula_trichocarpa_four_fold_fixed,
                            tra_w_n$tP.norm,
                            tra_w_n$tajD,
                            dav_w_n$tP.norm,
                            dav_w_n$tajD,
                            trs_w_n$tP.norm,
                            trs_w_n$tajD,
                            tri_w_n$tP.norm,
                            tri_w_n$tajD,
                            tra_rho_w_n$rho_mean,
                            dav_rho_w_n$rho_mean,
                            trs_rho_w_n$rho_mean,
                            tri_rho_w_n$rho_mean,
                            tra_dav_w_n$fst,
                            tra_dav_w_n$dxy,
                            tra_trs_w_n$fst,
                            tra_trs_w_n$dxy,
                            tra_tri_w_n$fst,
                            tra_tri_w_n$dxy,
                            dav_trs_w_n$fst,
                            dav_trs_w_n$dxy,
                            dav_tri_w_n$fst,
                            dav_tri_w_n$dxy,
                            trs_tri_w_n$fst,
                            trs_tri_w_n$dxy))

names(data_total)=c("Chr","Pos","GC","Gene_num","Gene_prop","Coding_prop",
                    "four_fold_fixed","tra_tP","tra_tajD","dav_tP","dav_tajD",
                    "trs_tP","trs_tajD","tri_tP","tri_tajD","tra_rho","dav_rho",
                    "trs_rho","tri_rho","tra_dav_fst","tra_dav_dxy","tra_trs_fst",
                    "tra_trs_dxy","tra_tri_fst","tra_tri_dxy","dav_trs_fst","dav_trs_dxy",
                    "dav_tri_fst","dav_tri_dxy","trs_tri_fst","trs_tri_dxy")
data_total$Chr=as.character(data_total$Chr)
###transfer from factor to numeric
indx=sapply(data_total,is.factor)
data_total[indx]=lapply(data_total[indx],function(x) as.numeric(as.character(x)))

##############################################
####add PBS estimates for P.tremula, P. tremuloides and P. trichocarpa
pbs<-function(fst1,fst2,fst3){   ###here assuming there are three populations, PBS of population1
  pbs1=((-log(1-fst1))+(-log(1-fst2))-(-log(1-fst3)))  ###fst1 is the fst(pop1 vs. pop2); fst2 is the fst(pop1 vs. pop3); fst3 is fst(pop2 vs. pop3)
  return(pbs1)
}

data_total$tra_pbs=pbs(data_total$tra_dav_fst,data_total$tra_trs_fst,data_total$dav_trs_fst)
data_total$dav_pbs=pbs(data_total$tra_dav_fst,data_total$dav_trs_fst,data_total$tra_trs_fst)
data_total$trs_pbs=pbs(data_total$dav_trs_fst,data_total$tra_trs_fst,data_total$tra_dav_fst)


#pi vs. rho; pi vs. gene_density; pi vs. mu; rho vs. gene_density; rho vs. mu
###pi vs. rho
pi_rho=data.frame(rbind(corrFunc("tra_tP","tra_rho",data_total),
                        corrFunc("dav_tP","dav_rho",data_total),
                        corrFunc("trs_tP","trs_rho",data_total),
                        corrFunc("tri_tP","tri_rho",data_total)))
summary_cor(pi_rho)


###pi vs. coding prop
pi_coding=data.frame(rbind(corrFunc("tra_tP","Coding_prop",data_total),
                        corrFunc("dav_tP","Coding_prop",data_total),
                        corrFunc("trs_tP","Coding_prop",data_total),
                        corrFunc("tri_tP","Coding_prop",data_total)))
summary_cor(pi_coding)

###rho vs mu
rho_mu=data.frame(rbind(corrFunc("tra_rho","four_fold_fixed",data_total),
                        corrFunc("dav_rho","four_fold_fixed",data_total),
                        corrFunc("trs_rho","four_fold_fixed",data_total),
                        corrFunc("tri_rho","four_fold_fixed",data_total)))
summary_cor(rho_mu)

###pi vs mu
pi_mu=data.frame(rbind(corrFunc("tra_tP","four_fold_fixed",data_total),
                       corrFunc("dav_tP","four_fold_fixed",data_total),
                       corrFunc("trs_tP","four_fold_fixed",data_total),
                       corrFunc("tri_tP","four_fold_fixed",data_total)))
summary_cor(pi_mu)


###pi vs tajD
pi_tajD=data.frame(rbind(corrFunc("tra_tP","tra_tajD",data_total),
                        corrFunc("dav_tP","dav_tajD",data_total),
                        corrFunc("trs_tP","trs_tajD",data_total),
                        corrFunc("tri_tP","tri_tajD",data_total)))

###rho vs. coding prop                        
rho_coding=data.frame(rbind(corrFunc("tra_rho","Coding_prop",data_total),
                        corrFunc("dav_rho","Coding_prop",data_total),
                        corrFunc("trs_rho","Coding_prop",data_total),
                        corrFunc("tri_rho","Coding_prop",data_total)))
summary_cor(rho_coding)

##########################################################################
###2.3 inter- and intra- species parameters correlations: pi vs.Fst; rho vs. Fst; D vs. Fst; pi vs dxy; rho vs. dxy; D vs dxy
##pi vs. Fst
pi_fst=data.frame(rbind(corrFunc("tra_tP","tra_dav_fst",data_total),
                        corrFunc("tra_tP","tra_trs_fst",data_total),
                        corrFunc("tra_tP","tra_tri_fst",data_total),
                        corrFunc("dav_tP","tra_dav_fst",data_total),
                        corrFunc("dav_tP","dav_trs_fst",data_total),
                        corrFunc("dav_tP","dav_tri_fst",data_total),
                        corrFunc("trs_tP","tra_trs_fst",data_total),
                        corrFunc("trs_tP","dav_trs_fst",data_total),
                        corrFunc("trs_tP","trs_tri_fst",data_total),
                        corrFunc("tri_tP","tra_tri_fst",data_total),
                        corrFunc("tri_tP","dav_tri_fst",data_total),
                        corrFunc("tri_tP","trs_tri_fst",data_total)))
summary_cor(pi_fst)

##rho vs. Fst
rho_fst=data.frame(rbind(corrFunc("tra_rho","tra_dav_fst",data_total),
                        corrFunc("tra_rho","tra_trs_fst",data_total),
                        corrFunc("tra_rho","tra_tri_fst",data_total),
                        corrFunc("dav_rho","tra_dav_fst",data_total),
                        corrFunc("dav_rho","dav_trs_fst",data_total),
                        corrFunc("dav_rho","dav_tri_fst",data_total),
                        corrFunc("trs_rho","tra_trs_fst",data_total),
                        corrFunc("trs_rho","dav_trs_fst",data_total),
                        corrFunc("trs_rho","trs_tri_fst",data_total),
                        corrFunc("tri_rho","tra_tri_fst",data_total),
                        corrFunc("tri_rho","dav_tri_fst",data_total),
                        corrFunc("tri_rho","trs_tri_fst",data_total)))
summary_cor(rho_fst)

##Coding prop vs. Fst
coding_fst=data.frame(rbind(corrFunc("Coding_prop","tra_dav_fst",data_total),
                        corrFunc("Coding_prop","tra_trs_fst",data_total),
                        corrFunc("Coding_prop","tra_tri_fst",data_total),
                        corrFunc("Coding_prop","tra_dav_fst",data_total),
                        corrFunc("Coding_prop","dav_trs_fst",data_total),
                        corrFunc("Coding_prop","dav_tri_fst",data_total),
                        corrFunc("Coding_prop","tra_trs_fst",data_total),
                        corrFunc("Coding_prop","dav_trs_fst",data_total),
                        corrFunc("Coding_prop","trs_tri_fst",data_total),
                        corrFunc("Coding_prop","tra_tri_fst",data_total),
                        corrFunc("Coding_prop","dav_tri_fst",data_total),
                        corrFunc("Coding_prop","trs_tri_fst",data_total)))
summary_cor(coding_fst)

##tajD vs. Fst
tajD_fst=data.frame(rbind(corrFunc("tra_tajD","tra_dav_fst",data_total),
                        corrFunc("tra_tajD","tra_trs_fst",data_total),
                        corrFunc("tra_tajD","tra_tri_fst",data_total),
                        corrFunc("dav_tajD","tra_dav_fst",data_total),
                        corrFunc("dav_tajD","dav_trs_fst",data_total),
                        corrFunc("dav_tajD","dav_tri_fst",data_total),
                        corrFunc("trs_tajD","tra_trs_fst",data_total),
                        corrFunc("trs_tajD","dav_trs_fst",data_total),
                        corrFunc("trs_tajD","trs_tri_fst",data_total),
                        corrFunc("tri_tajD","tra_tri_fst",data_total),
                        corrFunc("tri_tajD","dav_tri_fst",data_total),
                        corrFunc("tri_tajD","trs_tri_fst",data_total)))

##pi vs. dxy
pi_dxy=data.frame(rbind(corrFunc("tra_tP","tra_dav_dxy",data_total),
                        corrFunc("tra_tP","tra_trs_dxy",data_total),
                        corrFunc("tra_tP","tra_tri_dxy",data_total),
                        corrFunc("dav_tP","tra_dav_dxy",data_total),
                        corrFunc("dav_tP","dav_trs_dxy",data_total),
                        corrFunc("dav_tP","dav_tri_dxy",data_total),
                        corrFunc("trs_tP","tra_trs_dxy",data_total),
                        corrFunc("trs_tP","dav_trs_dxy",data_total),
                        corrFunc("trs_tP","trs_tri_dxy",data_total),
                        corrFunc("tri_tP","tra_tri_dxy",data_total),
                        corrFunc("tri_tP","dav_tri_dxy",data_total),
                        corrFunc("tri_tP","trs_tri_dxy",data_total)))
summary_cor(pi_dxy)

##rho vs. dxy
rho_dxy=data.frame(rbind(corrFunc("tra_rho","tra_dav_dxy",data_total),
                        corrFunc("tra_rho","tra_trs_dxy",data_total),
                        corrFunc("tra_rho","tra_tri_dxy",data_total),
                        corrFunc("dav_rho","tra_dav_dxy",data_total),
                        corrFunc("dav_rho","dav_trs_dxy",data_total),
                        corrFunc("dav_rho","dav_tri_dxy",data_total),
                        corrFunc("trs_rho","tra_trs_dxy",data_total),
                        corrFunc("trs_rho","dav_trs_dxy",data_total),
                        corrFunc("trs_rho","trs_tri_dxy",data_total),
                        corrFunc("tri_rho","tra_tri_dxy",data_total),
                        corrFunc("tri_rho","dav_tri_dxy",data_total),
                        corrFunc("tri_rho","trs_tri_dxy",data_total)))
summary_cor(rho_dxy)

##Coding prop vs. Fst
coding_dxy=data.frame(rbind(corrFunc("Coding_prop","tra_dav_dxy",data_total),
                        corrFunc("Coding_prop","tra_trs_dxy",data_total),
                        corrFunc("Coding_prop","tra_tri_dxy",data_total),
                        corrFunc("Coding_prop","tra_dav_dxy",data_total),
                        corrFunc("Coding_prop","dav_trs_dxy",data_total),
                        corrFunc("Coding_prop","dav_tri_dxy",data_total),
                        corrFunc("Coding_prop","tra_trs_dxy",data_total),
                        corrFunc("Coding_prop","dav_trs_dxy",data_total),
                        corrFunc("Coding_prop","trs_tri_dxy",data_total),
                        corrFunc("Coding_prop","tra_tri_dxy",data_total),
                        corrFunc("Coding_prop","dav_tri_dxy",data_total),
                        corrFunc("Coding_prop","trs_tri_dxy",data_total)))
summary_cor(coding_dxy)

##tajD vs. dxy
tajD_dxy=data.frame(rbind(corrFunc("tra_tajD","tra_dav_dxy",data_total),
                        corrFunc("tra_tajD","tra_trs_dxy",data_total),
                        corrFunc("tra_tajD","tra_tri_dxy",data_total),
                        corrFunc("dav_tajD","tra_dav_dxy",data_total),
                        corrFunc("dav_tajD","dav_trs_dxy",data_total),
                        corrFunc("dav_tajD","dav_tri_dxy",data_total),
                        corrFunc("trs_tajD","tra_trs_dxy",data_total),
                        corrFunc("trs_tajD","dav_trs_dxy",data_total),
                        corrFunc("trs_tajD","trs_tri_dxy",data_total),
                        corrFunc("tri_tajD","tra_tri_dxy",data_total),
                        corrFunc("tri_tajD","dav_tri_dxy",data_total),
                        corrFunc("tri_tajD","trs_tri_dxy",data_total)))

##########################################################################
###2.4 interspecies divergence correlations: Fst vs Fst; dxy vs. dxy; Fst vs. dxy

##fst
fst_table=data.frame(cbind(tra_dav_w_n$fst,
                          tra_trs_w_n$fst,
                          tra_tri_w_n$fst,
                          dav_trs_w_n$fst,
                          dav_tri_w_n$fst,
                          trs_tri_w_n$fst))

names(fst_table)=c("tra_dav_fst","tra_trs_fst","tra_tri_fst","dav_trs_fst","dav_tri_fst","trs_tri_fst")
corrs_fst=corrs_coef(fst_table)
summary_cor(corrs_fst)

##dxy
dxy_table=data.frame(cbind(tra_dav_w_n$dxy,
                          tra_trs_w_n$dxy,
                          tra_tri_w_n$dxy,
                          dav_trs_w_n$dxy,
                          dav_tri_w_n$dxy,
                          trs_tri_w_n$dxy))

names(dxy_table)=c("tra_dav_dxy","tra_trs_dxy","tra_tri_dxy","dav_trs_dxy","dav_tri_dxy","trs_tri_dxy")
corrs_dxy=corrs_coef(dxy_table)
summary_cor(corrs_dxy)

##fst vs. dxy
fst_dxy=data.frame(rbind(corrFunc("tra_dav_fst","tra_dav_dxy",data_total),
                         corrFunc("tra_trs_fst","tra_trs_dxy",data_total),
                         corrFunc("tra_tri_fst","tra_tri_dxy",data_total),
                         corrFunc("dav_trs_fst","dav_trs_dxy",data_total),
                         corrFunc("dav_tri_fst","dav_tri_dxy",data_total),
                         corrFunc("trs_tri_fst","trs_tri_dxy",data_total)))
summary_cor(fst_dxy)


############################################################################
####making the plots
par(mar=c(3,6,1,1))
par(las=1)
beanplot(corrs_pi$estimate,corrs_rho$estimate,names=c(expression(paste(pi[i],"-",pi[j],sep="")),expression(paste(rho[i],"-",rho[j],sep=""))),col=cols[1],horizontal = T,xaxt='n',frame.plot=F,ylim=c(-1,1),log="")


beanplot(corrs_fixed$estimate,corrs_dxy$estimate,corrs_fst$estimate,names=c(expression(paste(d[f][I],"-",d[f][J],sep="")),expression(paste(d[xy][I],"-",d[xy][J],sep="")),expression(paste(F[ST][I],"-",F[ST][J],sep=""))),col="light blue",horizontal = T,xaxt='n',frame.plot=F,ylim=c(-1,1),log="")
#vioplot(corrs_fixed$estimate,corrs_dxy$estimate,corrs_fst$estimate,names=c(expression(paste(d[f][I],"-",d[f][J],sep="")),expression(paste(d[xy][I],"-",d[xy][J],sep="")),expression(paste(F[ST][I],"-",F[ST][J],sep=""))),col="light blue",horizontal = T,xaxt='n',ylim=c(-1,1),log="")
beanplot(fst_dxy$estimate,fst_df$estimate,dxy_df$estimate,names=c(expression(paste(F[ST][I],"-",d[xy][I],sep="")),expression(paste(F[ST][I],"-",d[f][I],sep="")),expression(paste(d[xy][I],"-",d[f][I],sep=""))),col="orange",horizontal = T,frame.plot=F,ylim=c(-1,1),log="")



png("4species.pairwise_divergence.cor.png",width=7,height=6,units='in',res=400)
pdf("4species.pairwise_pi_rho_fst_dxy.cor.100kb.pdf",width=7,height=6)


##19. corrs_pi
##18. corrs_rho
##17. pi_mu
##16. rho_mu
##15. 
##14. pi_rho
##13. pi_coding
##12. rho_coding
##11
##10. corrs_fst
##9. corrs_dxy
##8. fst_dxy
#7
##6. pi_fst
##5. rho_fst
##4. coding_fst
##3. pi_dxy
##2. rho_dxy
##1. coding_dxy



par(mar=c(3,6,2,2))
par(las=1)
plot(0:1,0:1,type="n",ylim=c(1,19),xlim=c(-1,1),axes=FALSE,ann=FALSE)
##the first group
vioplot(corrs_pi$estimate,corrs_rho$estimate,pi_mu$estimate,rho_mu$estimate,add=T,at=19:16,horizontal=T,col=cols[1],xlim=c(-1,1))
axis(side=2,at=19:16,labels=c(expression(paste(pi[i],"-",pi[j],sep="")),expression(paste(rho[i],"-",rho[j],sep="")),expression(paste(mu,"-",pi[i],sep="")),expression(paste(mu,"-",rho[i],sep=""))),xaxt='n',axes=F)
axis(side=1,at=seq(-1,1,0.2),labels=seq(-1,1,0.2))  
abline(v=0,lty=2)

#the second group
vioplot(pi_rho$estimate,pi_coding$estimate,rho_coding$estimate,at=14:12,col=cols[2],horizontal = T,frame.plot=F,xlim=c(-1,1),add=T)
axis(side=2,at=14:12,labels=c(expression(paste(rho[i],"-",pi[i],sep="")),expression(paste(CD,"-",pi[i],sep="")),expression(paste(CD,"-",rho[i],sep=""))))

#the thrid group
vioplot(corrs_fst$estimate,corrs_dxy$estimate,fst_dxy$estimate,at=10:8,col=cols[4],horizontal = T,frame.plot=F,xlim=c(-1,1),add=T)
axis(side=2,at=10:8,labels=c(expression(paste(F[ST][I],"-",F[ST][J],sep="")),expression(paste(d[xy][I],"-",d[xy][J],sep="")),expression(paste(F[ST][I],"-",d[xy][I],sep=""))))

#the fouth group
vioplot(pi_fst$estimate,rho_fst$estimate,coding_fst$estimate,pi_dxy$estimate,rho_dxy$estimate,coding_fst$estimate,at=6:1,col=cols[3],horizontal = T,frame.plot=F,xlim=c(-1,1),add=T)
axis(side=2,at=6:1,labels=c(expression(paste(pi[i],"-",F[ST][I],sep="")),expression(paste(rho[i],"-",F[ST][I],sep="")),expression(paste("CD","-",F[ST][I],sep="")),expression(paste(pi[i],"-",d[xy][I],sep="")),expression(paste(rho[i],"-",d[xy][I],sep="")),expression(paste("CD","-",d[xy][I],sep=""))))


dev.off()

############################
##In the second version, we removed the relationship between divergence and pi, rho, CD
png("4species.pairwise_divergence.cor.2.png",width=7,height=6,units='in',res=400)
pdf("4species.pairwise_pi_rho_fst_dxy.cor.100kb.2.pdf",width=6,height=5)


##14. corrs_pi
##13. corrs_rho
##12. pi_mu
##11. rho_mu
##10. 
##9. pi_rho
##8. pi_coding
##7. rho_coding
##6
##5. corrs_fst
##4. corrs_dxy
##3. fst_pi
##2. dxy_pi
##1. fst_dxy



par(mar=c(3,6,2,2))
par(las=1)
plot(0:1,0:1,type="n",ylim=c(1,14),xlim=c(-1,1),axes=FALSE,ann=FALSE)
##the first group
vioplot(corrs_pi$estimate,corrs_rho$estimate,pi_mu$estimate,rho_mu$estimate,add=T,at=14:11,horizontal=T,col=cols[1],xlim=c(-1,1))
axis(side=2,at=14:11,labels=c(expression(paste(pi[i],"-",pi[j],sep="")),expression(paste(rho[i],"-",rho[j],sep="")),expression(paste(mu,"-",pi[i],sep="")),expression(paste(mu,"-",rho[i],sep=""))),xaxt='n',axes=F)
axis(side=1,at=seq(-1,1,0.2),labels=seq(-1,1,0.2))  
abline(v=0,lty=2)

#the second group
vioplot(pi_rho$estimate,pi_coding$estimate,rho_coding$estimate,at=9:7,col=cols[2],horizontal = T,frame.plot=F,xlim=c(-1,1),add=T)
axis(side=2,at=9:7,labels=c(expression(paste(rho[i],"-",pi[i],sep="")),expression(paste(CD,"-",pi[i],sep="")),expression(paste(CD,"-",rho[i],sep=""))))

#the thrid group
vioplot(corrs_fst$estimate,corrs_dxy$estimate,pi_fst$estimate,pi_dxy$estimate,fst_dxy$estimate,at=5:1,col=cols[4],horizontal = T,frame.plot=F,xlim=c(-1,1),add=T)
axis(side=2,at=5:1,labels=c(expression(paste(F[ST][I],"-",F[ST][J],sep="")),expression(paste(d[xy][I],"-",d[xy][J],sep="")),expression(paste(pi[i],"-",F[ST][I],sep="")),expression(paste(pi[i],"-",d[xy][I],sep="")),expression(paste(F[ST][I],"-",d[xy][I],sep=""))))

dev.off()


############################
##In the second version, we removed the relationship between divergence and pi, rho, CD
png("4species.pairwise_divergence.cor.2.png",width=7,height=6,units='in',res=400)
pdf("4species.pairwise_pi_rho_fst_dxy.cor.100kb.2.pdf",width=6,height=5)


##10. corrs_pi
##9. corrs_rho
##8. pi_mu
##7. rho_mu
##6.
##5. corrs_fst
##4. corrs_dxy
##3. fst_pi
##2. dxy_pi
##1. fst_dxy



par(mar=c(3,6,2,2))
par(las=1)
plot(0:1,0:1,type="n",ylim=c(1,10),xlim=c(-1,1),axes=FALSE,ann=FALSE)
##the first group
vioplot(corrs_pi$estimate,corrs_rho$estimate,pi_mu$estimate,rho_mu$estimate,add=T,at=10:7,horizontal=T,col=cols[1],xlim=c(-1,1))
axis(side=2,at=10:7,labels=c(expression(paste(pi[i],"-",pi[j],sep="")),expression(paste(rho[i],"-",rho[j],sep="")),expression(paste(mu,"-",pi[i],sep="")),expression(paste(mu,"-",rho[i],sep=""))),xaxt='n',axes=F)
axis(side=1,at=seq(-1,1,0.2),labels=seq(-1,1,0.2))  
abline(v=0,lty=2)

#the second group
vioplot(corrs_fst$estimate,corrs_dxy$estimate,pi_fst$estimate,pi_dxy$estimate,fst_dxy$estimate,at=5:1,col=cols[1],horizontal = T,frame.plot=F,xlim=c(-1,1),add=T)
axis(side=2,at=5:1,labels=c(expression(paste(F[ST][I],"-",F[ST][J],sep="")),expression(paste(d[xy][I],"-",d[xy][J],sep="")),expression(paste(pi[i],"-",F[ST][I],sep="")),expression(paste(pi[i],"-",d[xy][I],sep="")),expression(paste(F[ST][I],"-",d[xy][I],sep=""))))

dev.off()




