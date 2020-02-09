library(data.table)
library(RColorBrewer)
library(dplyr)
library(reshape2)
library(ggplot2)
library("VennDiagram")


##folded
setwd("~/Dropbox/davidiana_paper/data/selection/BetaScan/1kb/unfolded/")
###############
#read in betascan results
beta_tra=fread("Betas_tremula.all.betascan.txt",header=T)
beta_trs=fread("Betas_tremuloides.all.betascan.txt",header=T)
beta_dav=fread("Betas_davidiana.all.betascan.txt",header=T)


#summary
summary_10kb=fread("~/Dropbox/davidiana_paper/data/SnIPRE/scripts/out/summary.24235win.div.fd.gene.txt",header=T)
#ILS
nls_10kb=fread("~/Dropbox/davidiana_paper/data/snp_allele_freq/allele_frq_aspen_10kb.nls.txt",header=T)
nls_10kb_new=nls_10kb %>% filter(n_sites>50)
#twisst
twisst_10kb=read.table("~/Dropbox/davidiana_paper/data/twisst/data/4species.gatk.hap.snp.rm_indel.bed.biallelic.DP5.GQ10.rm_missing.recode.gt.beagle.twisst.weights.10kb.csv")
twisst_10kb_new=twisst_10kb[-1,]
names(twisst_10kb_new)=c("top1","top2","top3")

indx <- sapply(twisst_10kb_new, is.factor)
twisst_10kb_new[indx] <- lapply(twisst_10kb_new[indx], function(x) as.numeric(as.character(x)))
colSums(twisst_10kb_new)/sum(colSums(twisst_10kb_new))

twisst_10kb_prop=prop.table(as.matrix(twisst_10kb_new),1)
###twisst sites
twisst_10kb_sites=read.table("~/Dropbox/davidiana_paper/data/twisst/data/4species.gatk.hap.snp.rm_indel.bed.biallelic.DP5.GQ10.rm_missing.recode.gt.beagle.raxml.w10kb.data.tsv",header=T)
twisst_10kb_sites_new=twisst_10kb_sites[which(twisst_10kb_sites$sites>50),]
twisst_10kb_prop_sites=cbind(twisst_10kb_sites_new,twisst_10kb_prop)

###########################################
###common windows
twisst_10kb_prop_sites$win=paste(twisst_10kb_prop_sites$scaffold,twisst_10kb_prop_sites$end+5000,sep=":")
nls_10kb_new$win=paste(nls_10kb_new$chro,nls_10kb_new$end+5000,sep=":")

common_win=Reduce(intersect,list(summary_10kb$win,
                                 twisst_10kb_prop_sites$win,
                                 nls_10kb_new$win
))

summary_10kb_w=summary_10kb[which(summary_10kb$win %in% common_win),]
twisst_10kb_prop_sites_w=twisst_10kb_prop_sites[which(twisst_10kb_prop_sites$win %in% common_win),]
nls_10kb_new_w=nls_10kb_new[which(nls_10kb_new$win %in% common_win),]

summary_10kb_w$top1=twisst_10kb_prop_sites_w$top1
summary_10kb_w$top2=twisst_10kb_prop_sites_w$top2
summary_10kb_w$top3=twisst_10kb_prop_sites_w$top3
summary_10kb_w$nls_n=nls_10kb_new_w$nls_n

summary=summary_10kb_w

##############################################################################
##############################################################################
#The script has three aims:
#1. For selective regions detected by Sweepfinder2, choose the top 1% windows from all three aspen species, and compare the 
#2. SnIPRE, positive selection windows 
#3. fd outliers 

###for the above three groups, compare: diveristy(pi), divergence (fst,dxy),ILS,topology(%),fd


###making the venndiagram plot
colors <- brewer.pal(12,"Paired")
cols=colors[c(1,2,3,4,5,6,7,8)]


##############################################################################
source("../script/beta_manhanttan_plot.R")
##############################################################################
#making plot
tra_chr1=gsub("Chr0","",beta_tra$Chromo)
tra_chr2=gsub("Chr","",tra_chr1)

###tra
png("beta_tra.unfolded.png",width=16,height=3.5,units='in',res=400)
par(mar=c(4,5,1,1))
tra_beta=data.frame(cbind(CHR=as.numeric(tra_chr2),BP=as.numeric(as.character(beta_tra$Pos)),P=as.numeric(as.character(beta_tra$Score))))
tra_beta$SNP=paste(tra_beta$CHR,tra_beta$BP,sep=":")
manhattan(tra_beta,chr="CHR",bp="BP",p="P",highligt=NA,xlabel="Chromosome",ylabel="Beta Score",col=c(cols[4],cols[3]),mar_value=c(4,5,1,1)+0.1,cex_value=1.5)
tra_beta_cutoff=quantile(tra_beta$P,0.99)
abline(h=tra_beta_cutoff,lty=1,col="red",lwd=2)
dev.off()

###dav
dav_chr1=gsub("Chr0","",beta_dav$Chromo)
dav_chr2=gsub("Chr","",dav_chr1)

png("beta_dav.unfolded.png",width=16,height=3.5,units='in',res=400)
par(mar=c(4,5,1,1))
dav_beta=data.frame(cbind(CHR=as.numeric(dav_chr2),BP=as.numeric(as.character(beta_dav$Pos)),P=as.numeric(as.character(beta_dav$Score))))
manhattan(dav_beta,chr="CHR",bp="BP",p="P",xlabel="Chromosome",ylabel="Beta Score",col=c(cols[2],cols[1]),mar_value=c(4,5,1,1)+0.1,cex_value=1.5)
dav_beta_cutoff=quantile(dav_beta$P,0.99)
abline(h=dav_beta_cutoff,lty=1,col="red",lwd=2)
dev.off()

###trs
trs_chr1=gsub("Chr0","",beta_trs$Chromo)
trs_chr2=gsub("Chr","",trs_chr1)

png("beta_trs.unfolded.png",width=16,height=3.5,units='in',res=400)
par(mar=c(4,5,1,1))
trs_beta=data.frame(cbind(CHR=as.numeric(trs_chr2),BP=as.numeric(as.character(beta_trs$Pos)),P=as.numeric(as.character(beta_trs$Score))))
manhattan(trs_beta,chr="CHR",bp="BP",p="P",xlabel="Chromosome",ylabel="Beta Score",col=c(cols[6],cols[5]),mar_value=c(4,5,1,1)+0.1,cex_value=1.5)
trs_beta_cutoff=quantile(trs_beta$P,0.99)
abline(h=trs_beta_cutoff,lty=1,col="red",lwd=2)
dev.off()


##############################################################################
###VennDiagram plot of the balancing selected regions among the four species
tra_beta$win=paste(tra_beta$CHR,tra_beta$BP,sep=":")
tra_beta_outlier=tra_beta[which(tra_beta$P>tra_beta_cutoff),]

dav_beta$win=paste(dav_beta$CHR,dav_beta$BP,sep=":")
dav_beta_outlier=dav_beta[which(dav_beta$P>dav_beta_cutoff),]

trs_beta$win=paste(trs_beta$CHR,trs_beta$BP,sep=":")
trs_beta_outlier=trs_beta[which(trs_beta$P>trs_beta_cutoff),]

#outlier SNPs with extremely high beta for each species
three_beta=venn.diagram(list(a=tra_beta_outlier$win,b=dav_beta_outlier$win,c=trs_beta_outlier$win),fill=c(cols[3],cols[1],cols[5]),category.names=c("P.tra","P.dav","P.trs"),filename=NULL,cat.cex = 2.5,cex = 2.5,margin = 0.05,cat.fontface = "bold",fontface = "bold")

png("VennDiagram.3species.unfolded.png",height=6,width=6,units='in',res=600)
grid.draw(three_beta)
dev.off()


##############################################################################
setwd("/Users/jiwa0011/Dropbox/davidiana_paper/data/selection/BetaScan/1kb/unfolded/sig_snps_three_species")

###2089 outlier SNPs detected in at least two species
common_two=unique(c(intersect(tra_beta_outlier$win,dav_beta_outlier$win),
                  intersect(tra_beta_outlier$win,trs_beta_outlier$win),
                  intersect(trs_beta_outlier$win,dav_beta_outlier$win)))

##############################################################################
###making plot with common SNPs showing as red 
#making plot
tra_chr1=gsub("Chr0","",beta_tra$Chromo)
tra_chr2=gsub("Chr","",tra_chr1)

###tra
png("beta_tra.mark_common.unfolded.png",width=16,height=3.5,units='in',res=400)
par(mar=c(4,5,1,1))
tra_beta=data.frame(cbind(CHR=as.numeric(tra_chr2),BP=as.numeric(as.character(beta_tra$Pos)),P=as.numeric(as.character(beta_tra$Score))))
tra_beta$SNP=paste(tra_beta$CHR,tra_beta$BP,sep=":")
manhattan(tra_beta,chr="CHR",bp="BP",p="P",highlight=common_two,col=c(cols[4],cols[3]),xlabel="Chromosome",ylabel="Beta Score",mar_value=c(4,5,1,1)+0.1,cex_value=1.5)
tra_beta_cutoff=quantile(tra_beta$P,0.99)
abline(h=tra_beta_cutoff,lty=1,col="red",lwd=2)
dev.off()

###dav
dav_chr1=gsub("Chr0","",beta_dav$Chromo)
dav_chr2=gsub("Chr","",dav_chr1)

png("beta_dav.mark_common.unfolded.png",width=16,height=3.5,units='in',res=400)
par(mar=c(4,5,1,1))
dav_beta=data.frame(cbind(CHR=as.numeric(dav_chr2),BP=as.numeric(as.character(beta_dav$Pos)),P=as.numeric(as.character(beta_dav$Score))))
dav_beta$SNP=paste(dav_beta$CHR,dav_beta$BP,sep=":")
manhattan(dav_beta,chr="CHR",bp="BP",p="P",highlight=common_two,col=c(cols[2],cols[1]),xlabel="Chromosome",ylabel="Beta Score",mar_value=c(4,5,1,1)+0.1,cex_value=1.5)
dav_beta_cutoff=quantile(dav_beta$P,0.99)
abline(h=dav_beta_cutoff,lty=1,col="red",lwd=2)
dev.off()

###trs
trs_chr1=gsub("Chr0","",beta_trs$Chromo)
trs_chr2=gsub("Chr","",trs_chr1)

png("beta_trs.mark_common.unfolded.png",width=16,height=3.5,units='in',res=400)
par(mar=c(4,5,1,1))
trs_beta=data.frame(cbind(CHR=as.numeric(trs_chr2),BP=as.numeric(as.character(beta_trs$Pos)),P=as.numeric(as.character(beta_trs$Score))))
trs_beta$SNP=paste(trs_beta$CHR,trs_beta$BP,sep=":")
manhattan(trs_beta,chr="CHR",bp="BP",p="P",highlight=common_two,col=c(cols[6],cols[5]),xlabel="Chromosome",ylabel="Beta Score",mar_value=c(4,5,1,1)+0.1,cex_value=1.5)
trs_beta_cutoff=quantile(trs_beta$P,0.99)
abline(h=trs_beta_cutoff,lty=1,col="red",lwd=2)
dev.off()


##############################################################################
###Compare all population genetic parameters 

summary$window=gsub("Chr0","",summary$win)
summary$window=gsub("Chr","",summary$window)

summary$chromo=gsub("Chr0","",summary$Chr)
summary$chromo=gsub("Chr","",summary$chromo)


###Define a function to summarize how many common SNPs (balancing selected SNPs detected in all three species) in each window
bal_win=function(data,snp){
  snp_data=data.frame(do.call(rbind, strsplit(snp, ":", fixed=TRUE)))
  names(snp_data)=c("chromo","SNP")
  n=nrow(snp_data)
  n_data=nrow(data)
  data$bal_win=0
  for (i in 1:n_data) {
    data[i,]$bal_win=length(which(as.numeric(as.character(snp_data$chromo))==data[i,]$chromo & as.numeric(as.character(snp_data$SNP))>data[i,]$start & as.numeric(as.character(snp_data$SNP))<data[i,]$end))
  }
  return(data)
}

##summarize the number of SNPs detected as under balancing selection in at least two of the three species
summary_bal_two=bal_win(summary,common_two)

##127 windows for selection detected in at least two species
bal_outlier=summary_bal_two[which(summary_bal_two$bal_win>0),]
##21438 windows with no sselection signals
bal_non_outlier=summary_bal_two[which(summary_bal_two$bal_win==0),]

########################
#write the windows and genes into txt file
write.table(bal_outlier,file="balancing_selected.at_least_two_species.win.txt",sep="\t", quote=F, row.names=F, col.names=T)


##############################################################################
##############################################################################
# * P<0.05
# ** P<0.01
# *** P<0.001

###1. Diversity
pi_bal=c(bal_non_outlier$tra_tP,bal_outlier$tra_tP,
         bal_non_outlier$dav_tP,bal_outlier$dav_tP,
         bal_non_outlier$trs_tP,bal_outlier$trs_tP)

wilcox.test(bal_non_outlier$tra_tP,bal_outlier$tra_tP)
wilcox.test(bal_non_outlier$dav_tP,bal_outlier$dav_tP)
wilcox.test(bal_non_outlier$trs_tP,bal_outlier$trs_tP)


pi_bal_group=c(rep("P.tra",21565),rep("P.dav",21565),rep("P.trs",21565))
pi_bal_fill=c(rep("bal_non_outlier",nrow(bal_non_outlier)),
              rep("bal_outlier",nrow(bal_outlier)),
              rep("bal_non_outlier",nrow(bal_non_outlier)),
              rep("bal_outlier",nrow(bal_outlier)),
              rep("bal_non_outlier",nrow(bal_non_outlier)),
              rep("bal_outlier",nrow(bal_outlier)))

pi_bal_table=as.data.frame(cbind(pi_bal,pi_bal_group,pi_bal_fill))

pi_bal_table$pi_bal=as.numeric(as.character(pi_bal_table$pi_bal))
pi_bal_table$pi_bal_group=factor(pi_bal_table$pi_bal_group,levels=c("P.tra","P.dav","P.trs"))
pi_bal_table$pi_bal_fill=factor(pi_bal_table$pi_bal_fill,levels=c("bal_non_outlier","bal_outlier"))

ggplot(pi_bal_table, aes(x=pi_bal_group, y=pi_bal,fill=pi_bal_fill))+ 
  geom_boxplot(outlier.shape=NA)+scale_fill_manual(values=c("grey60","red"))+
  scale_y_continuous(limits=c(0,0.03))+
  #xlab("")+ylab(expression(paste(bal," (P.tra)",sep="")))+
  xlab("")+ylab("")+
  annotate(geom="text",x=c(1,2,3),y=c(0.029,0.029,0.029),label=c(rep("***",3)),size=8)+
  theme(legend.position="none",axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),axis.text.x = element_blank(),axis.text.y = element_text(size=16))

ggsave("bal.pi.boxplot.png",width=3.5,height=4.5,units='in',dpi=500)

##############################################################################
###2. recombination
rho_bal=c(bal_non_outlier$tra_rho_mean,bal_outlier$tra_rho_mean,
          bal_non_outlier$dav_rho_mean,bal_outlier$dav_rho_mean,
          bal_non_outlier$trs_rho_mean,bal_outlier$trs_rho_mean)

wilcox.test(bal_outlier$tra_rho_mean,bal_non_outlier$tra_rho_mean)
wilcox.test(bal_outlier$dav_rho_mean,bal_non_outlier$dav_rho_mean)
wilcox.test(bal_outlier$trs_rho_mean,bal_non_outlier$trs_rho_mean)


rho_bal_group=c(rep("P.tra",21565),rep("P.dav",21565),rep("P.trs",21565))
rho_bal_fill=c(rep("bal_non_outlier",nrow(bal_non_outlier)),
               rep("bal_outlier",nrow(bal_outlier)),
               rep("bal_non_outlier",nrow(bal_non_outlier)),
               rep("bal_outlier",nrow(bal_outlier)),
               rep("bal_non_outlier",nrow(bal_non_outlier)),
               rep("bal_outlier",nrow(bal_outlier)))

rho_bal_table=as.data.frame(cbind(rho_bal,rho_bal_group,rho_bal_fill))

rho_bal_table$rho_bal=as.numeric(as.character(rho_bal_table$rho_bal))
rho_bal_table$rho_bal_group=factor(rho_bal_table$rho_bal_group,levels=c("P.tra","P.dav","P.trs"))
rho_bal_table$rho_bal_fill=factor(rho_bal_table$rho_bal_fill,levels=c("bal_non_outlier","bal_outlier"))

ggplot(rho_bal_table, aes(x=rho_bal_group, y=rho_bal,fill=rho_bal_fill))+ 
  geom_boxplot(outlier.shape=NA)+scale_fill_manual(values=c("grey60","red"))+
  scale_y_continuous(limits=c(0,0.07))+
  #xlab("")+ylab(expression(paste(bal," (P.tra)",sep="")))+
  xlab("")+ylab("")+
  annotate(geom="text",x=c(1,2,3),y=c(0.068,0.068,0.068),label=c(c("n.s.","**","n.s.")),size=8)+
  theme(legend.position="none",axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),axis.text.x = element_blank(),axis.text.y = element_text(size=16))

ggsave("bal.rho.boxplot.png",width=3.5,height=4.5,units='in',dpi=500)


##############################################################################
###3. Divergence, Fst
fst_bal=c(bal_non_outlier$tra_dav_fst,bal_outlier$tra_dav_fst,
          bal_non_outlier$tra_trs_fst,bal_outlier$tra_trs_fst,
          bal_non_outlier$dav_trs_fst,bal_outlier$dav_trs_fst)

wilcox.test(bal_non_outlier$tra_dav_fst,bal_outlier$tra_dav_fst)
wilcox.test(bal_non_outlier$tra_trs_fst,bal_outlier$tra_trs_fst)
wilcox.test(bal_non_outlier$dav_trs_fst,bal_outlier$dav_trs_fst)


fst_bal_group=c(rep("P.tra-P.dav",21565),rep("P.tra-P.trs",21565),rep("P.dav-P.trs",21565))
fst_bal_fill=c(rep("bal_non_outlier",nrow(bal_non_outlier)),
               rep("bal_outlier",nrow(bal_outlier)),
               rep("bal_non_outlier",nrow(bal_non_outlier)),
               rep("bal_outlier",nrow(bal_outlier)),
               rep("bal_non_outlier",nrow(bal_non_outlier)),
               rep("bal_outlier",nrow(bal_outlier)))

fst_bal_table=as.data.frame(cbind(fst_bal,fst_bal_group,fst_bal_fill))

fst_bal_table$fst_bal=as.numeric(as.character(fst_bal_table$fst_bal))
fst_bal_table$fst_bal_group=factor(fst_bal_table$fst_bal_group,levels=c("P.tra-P.dav","P.tra-P.trs","P.dav-P.trs"))
fst_bal_table$fst_bal_fill=factor(fst_bal_table$fst_bal_fill,levels=c("bal_non_outlier","bal_outlier"))

ggplot(fst_bal_table, aes(x=fst_bal_group, y=fst_bal,fill=fst_bal_fill))+ 
  geom_boxplot(outlier.shape=NA)+scale_fill_manual(values=c("grey60","red"))+
  scale_y_continuous(limits=c(0,1))+
  #xlab("")+ylab(expression(paste(bal," (P.tra)",sep="")))+
  xlab("")+ylab("")+
  annotate(geom="text",x=c(1,2,3),y=c(0.98,0.98,0.98),label=c("***","***","***"),size=8)+
  theme(legend.position="none",axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),axis.text.x = element_blank(),axis.text.y = element_text(size=16))

ggsave("bal.fst.boxplot.png",width=3.5,height=4.5,units='in',dpi=500)

##############################################################################
###4. Divergence, dxy
dxy_bal=c(bal_non_outlier$tra_dav_dxy,bal_outlier$tra_dav_dxy,
          bal_non_outlier$tra_trs_dxy,bal_outlier$tra_trs_dxy,
          bal_non_outlier$dav_trs_dxy,bal_outlier$dav_trs_dxy)

wilcox.test(bal_non_outlier$tra_dav_dxy,bal_outlier$tra_dav_dxy)
wilcox.test(bal_non_outlier$tra_trs_dxy,bal_outlier$tra_trs_dxy)
wilcox.test(bal_non_outlier$dav_trs_dxy,bal_outlier$dav_trs_dxy)


dxy_bal_group=c(rep("P.tra-P.dav",21565),rep("P.tra-P.trs",21565),rep("P.dav-P.trs",21565))
dxy_bal_fill=c(rep("bal_non_outlier",nrow(bal_non_outlier)),
               rep("bal_outlier",nrow(bal_outlier)),
               rep("bal_non_outlier",nrow(bal_non_outlier)),
               rep("bal_outlier",nrow(bal_outlier)),
               rep("bal_non_outlier",nrow(bal_non_outlier)),
               rep("bal_outlier",nrow(bal_outlier)))

dxy_bal_table=as.data.frame(cbind(dxy_bal,dxy_bal_group,dxy_bal_fill))

dxy_bal_table$dxy_bal=as.numeric(as.character(dxy_bal_table$dxy_bal))
dxy_bal_table$dxy_bal_group=factor(dxy_bal_table$dxy_bal_group,levels=c("P.tra-P.dav","P.tra-P.trs","P.dav-P.trs"))
dxy_bal_table$dxy_bal_fill=factor(dxy_bal_table$dxy_bal_fill,levels=c("bal_non_outlier","bal_outlier"))

ggplot(dxy_bal_table, aes(x=dxy_bal_group, y=dxy_bal,fill=dxy_bal_fill))+ 
  geom_boxplot(outlier.shape=NA)+scale_fill_manual(values=c("grey60","red"))+
  scale_y_continuous(limits=c(0,0.05))+
  #xlab("")+ylab(expression(paste(bal," (P.tra)",sep="")))+
  xlab("")+ylab("")+
  annotate(geom="text",x=c(1,2,3),y=c(0.049,0.049,0.049),label=c("n.s.","*","*"),size=8)+
  theme(legend.position="none",axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),axis.text.x = element_blank(),axis.text.y = element_text(size=16))

ggsave("bal.dxy.boxplot.png",width=3.5,height=4.5,units='in',dpi=500)


##############################################################################
###5. topology
topo_bal=c(bal_non_outlier$top1,bal_outlier$top1,
           bal_non_outlier$top2,bal_outlier$top2,
           bal_non_outlier$top3,bal_outlier$top3)

wilcox.test(bal_non_outlier$top1,bal_outlier$top1)
wilcox.test(bal_non_outlier$top2,bal_outlier$top2)
wilcox.test(bal_non_outlier$top3,bal_outlier$top3)

# * P<0.05
# ** P<0.01
# *** P<0.001


topo_bal_group=c(rep("topo1",21565),rep("topo2",21565),rep("topo3",21565))
topo_bal_fill=c(rep("bal_non_outlier",nrow(bal_non_outlier)),
                rep("bal_outlier",nrow(bal_outlier)),
                rep("bal_non_outlier",nrow(bal_non_outlier)),
                rep("bal_outlier",nrow(bal_outlier)),
                rep("bal_non_outlier",nrow(bal_non_outlier)),
                rep("bal_outlier",nrow(bal_outlier)))

topo_bal_table=as.data.frame(cbind(topo_bal,topo_bal_group,topo_bal_fill))

topo_bal_table$topo_bal=as.numeric(as.character(topo_bal_table$topo_bal))
topo_bal_table$topo_bal_group=factor(topo_bal_table$topo_bal_group,levels=c("topo1","topo2","topo3"))
topo_bal_table$topo_bal_fill=factor(topo_bal_table$topo_bal_fill,levels=c("bal_non_outlier","bal_outlier"))

ggplot(topo_bal_table, aes(x=topo_bal_group, y=topo_bal,fill=topo_bal_fill))+ 
  geom_boxplot(outlier.shape=NA)+scale_fill_manual(values=c("grey60","red"))+
  scale_y_continuous(limits=c(0,1.1))+
  #xlab("")+ylab(expression(paste(bal," (P.tra)",sep="")))+
  xlab("")+ylab("")+
  annotate(geom="text",x=c(1,2,3),y=c(1.08,1.08,1.08),label=c("n.s.","n.s.","*"),size=8)+
  theme(legend.position="none",axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),axis.text.x = element_blank(),axis.text.y = element_text(size=16))

ggsave("bal.topo.boxplot.png",width=3.5,height=4.5,units='in',dpi=500)

##############################################################################
###6.ILS

ILS_bal=as.data.frame(cbind(bal_non_outlier$nls_n,bal_outlier$nls_n))
names(ILS_bal)=c("non_outlier","outlier")
ILS_melt_bal=melt(ILS_bal)
wilcox.test(bal_non_outlier$nls_n,bal_outlier$nls_n)

ggplot(ILS_melt_bal, aes(x=variable, y=value,fill=variable))+ 
  geom_boxplot(outlier.shape=NA)+scale_fill_manual(values=c("grey60","red"))+
  scale_y_continuous(limits=c(0,0.3))+
  #xlab("")+ylab(expression(paste(bal," (P.tra)",sep="")))+
  xlab("")+ylab("")+
  annotate(geom="text",x=1.5,y=0.29,label=c("**"),size=8)+
  theme(legend.position="none",axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),axis.text.x = element_blank(),axis.text.y = element_text(size=16))

ggsave("bal.ILS.boxplot.png",width=2,height=4,units='in',dpi=500)


##############################################################################
###7.fd

fd_bal=as.data.frame(cbind(bal_non_outlier$fd,bal_outlier$fd))
names(fd_bal)=c("non_outlier","outlier")
fd_melt_bal=melt(fd_bal)
wilcox.test(bal_non_outlier$fd,bal_outlier$fd)

ggplot(fd_melt_bal, aes(x=variable, y=value,fill=variable))+ 
  geom_boxplot(outlier.shape=NA)+scale_fill_manual(values=c("grey60","red"))+
  scale_y_continuous(limits=c(-0.3,0.5))+
  #xlab("")+ylab(expression(paste(bal," (P.tra)",sep="")))+
  xlab("")+ylab("")+
  annotate(geom="text",x=1.5,y=0.49,label=c("*"),size=8)+
  theme(legend.position="none",axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),axis.text.x = element_blank(),axis.text.y = element_text(size=16))

ggsave("bal.fd.boxplot.png",width=2,height=4,units='in',dpi=500)



#########################################################################
########Output the SNPs in the data
snp_data=data.frame(do.call(rbind, strsplit(common_two, ":", fixed=TRUE)))
names(snp_data)=c("chromo","SNP")

gene_data=fread("~/Dropbox/davidiana_paper/data/SnIPRE/scripts/out/trichocarpa.gene.filter.bed.txt",header=F)
names(gene_data)=c("Chr","Start","End","Gene")

gene_data$Chr=gsub("Chr0","",gene_data$Chr)
gene_data$Chr=gsub("Chr","",gene_data$Chr)

snp_data$Gene="NA"

###Define a function to summarize how many common SNPs (balancing selected SNPs detected in each gene
bal_win=function(data,snp){
  snp_data=data.frame(do.call(rbind, strsplit(snp, ":", fixed=TRUE)))
  names(snp_data)=c("chromo","SNP")
  n=nrow(snp_data)
  n_data=nrow(data)
  for (i in 1:n) {
    gene=data[which(data$Start <= as.numeric(as.character(snp_data[i,]$SNP)) & data$End >= as.numeric(as.character(snp_data[i,]$SNP)) & data$Chr == as.numeric(as.character(snp_data[i,]$chromo))),]$Gene
    if (identical(gene,character(0))) {
      next
    } else {
      snp_data[i,]$Gene=gene
    }
  }
  return(snp_data)
}

### only find three genes  "Potri.001G261600" "Potri.001G350800" "Potri.005G181600"


##############################################################################
##############################################################################
##Extract the genes within the selected regions and do GO enrichment analysis

n_gene=nrow(bal_outlier)
genes_bal_outlier=c()
for (i in 1:n_gene){
  genes_sf2=(strsplit(bal_outlier$geneIDs[i],","))[[1]]
  if (genes_sf2==".") {
    next
  }
  genes_bal_outlier=c(genes_bal_outlier,genes_sf2)
}

outlier_gene=unique(genes_bal_outlier) ##166 genes

##GO enrichment
source("https://bioconductor.org/biocLite.R")
#biocLite("topGO")
library(topGO)
#library(dplyr)
d=read.table("~/Dropbox/davidiana_paper/data/SnIPRE/scripts/out/trichocarpa.gene_go.txt",header=T)
d$gene_id=as.character(d$gene_id)
d$go_id=as.character(d$go_id)
all_go <- lapply(split(d, sub("\\.\\d+$", "", d[, 1])), function(x) unique(x[, 2]))

geneNames=names(all_go)
geneList=factor(as.integer(geneNames %in% outlier_gene))
names(geneList)=geneNames

GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList,annot=annFUN.gene2GO, gene2GO = all_go)
restRes=runTest(GOdata,algorithm="classic",statistic="fisher")
GenTable(GOdata, p.value = restRes, orderBy = "p.value")

###choose GO with pvalue <0.01
gene_table=GenTable(GOdata,Fisher.p=restRes,topNodes=6)

###because there is 1395 GO terms and we can calculated the adjust P-values
gene_table_total=GenTable(GOdata,Fisher.p=restRes,topNodes=1395)
allGO=usedGO(GOdata)
pvalues=gene_table_total$Fisher.p
gene_table$adjust.p=round(head(p.adjust(pvalues,method="BH"),6),10)

write.table(gene_table,file="betascan.outlier.GO.enrich.txt",sep="\t", quote=F, row.names=F, col.names=T)



