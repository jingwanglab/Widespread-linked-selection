library(data.table)
library(RColorBrewer)
library(dplyr)
library(reshape2)
library(ggplot2)
library("VennDiagram")


##set working directory
setwd("~/Dropbox/davidiana_paper/data/selection/XPCLR/results/")
###############
#read in xpclr results
xpclr_trs_dav_trs=fread('cat trs_dav.trs.xpclr.txt.gz |gunzip',sep=" ",header=F)
names(xpclr_trs_dav_trs)=c("chr","grid","snp_n","physical_pos","genetic_pos","xpclr","max_s")
xpclr_trs_dav_dav=fread('cat trs_dav.dav.xpclr.txt.gz |gunzip',sep=" ",header=F)
names(xpclr_trs_dav_dav)=c("chr","grid","snp_n","physical_pos","genetic_pos","xpclr","max_s")

xpclr_tra_dav_tra=fread('cat tra_dav.tra.xpclr.txt.gz |gunzip',sep=" ",header=F)
names(xpclr_tra_dav_tra)=c("chr","grid","snp_n","physical_pos","genetic_pos","xpclr","max_s")
xpclr_tra_dav_dav=fread('cat tra_dav.dav.xpclr.txt.gz |gunzip',sep=" ",header=F)
names(xpclr_tra_dav_dav)=c("chr","grid","snp_n","physical_pos","genetic_pos","xpclr","max_s")

xpclr_tra_trs_tra=fread('cat tra_trs.tra.xpclr.txt.gz |gunzip',sep=" ",header=F)
names(xpclr_tra_trs_tra)=c("chr","grid","snp_n","physical_pos","genetic_pos","xpclr","max_s")
xpclr_tra_trs_trs=fread('cat tra_trs.trs.xpclr.txt.gz |gunzip',sep=" ",header=F)
names(xpclr_tra_trs_trs)=c("chr","grid","snp_n","physical_pos","genetic_pos","xpclr","max_s")

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

summary$chr=gsub("Chr0","",summary$Chr)
summary$chr=gsub("Chr","",summary$chr)

##############################################################################
##############################################################################
#define the function the merge the XPCLR results into the summary filesh

merge_xpclr=function(xpclr,summary){

xpclr_vector=c()
for (i in 1:19) {
  chr_v=which(as.numeric(as.vector(xpclr$chr))==i)  ###choose the vector of the specific chromosome
  xpclr_chr=xpclr[chr_v,]  ###the data from of xpclr for the chromosome
  chr_s=which(as.numeric(as.vector(summary$chr))==i)
  summary_chr=summary[chr_s,]
  
  summary_pos=summary_chr$Pos
  summary_length=length(summary_pos)
  
  for (j in 1:summary_length) {
    xpclr_chr_win=xpclr_chr[which(xpclr_chr$physical_pos>=summary_pos[j]-5000 & xpclr_chr$physical_pos<=summary_pos[j]+5000),]
    xpclr_vector=c(xpclr_vector,mean(xpclr_chr_win[!is.infinite(xpclr_chr_win$xpclr)]$xpclr))
  }

}
return(xpclr_vector)
}

summary$xpclr_trs_dav_trs=merge_xpclr(xpclr_trs_dav_trs,summary)
summary$xpclr_trs_dav_dav=merge_xpclr(xpclr_trs_dav_dav,summary)
summary$xpclr_tra_dav_tra=merge_xpclr(xpclr_tra_dav_tra,summary)
summary$xpclr_tra_dav_dav=merge_xpclr(xpclr_tra_dav_dav,summary)
summary$xpclr_tra_trs_tra=merge_xpclr(xpclr_tra_trs_tra,summary)
summary$xpclr_tra_trs_trs=merge_xpclr(xpclr_tra_trs_trs,summary)


##############################################################################
##############################################################################
#The script has three aims:
#1. For selective regions detected by XPCLR, choose the top 1% windows from all three aspen species(those detected by the other two reference), and compare the 

###for the above three groups, compare: diveristy(pi), divergence (fst,dxy),ILS,topology(%),fd


###making the venndiagram plot
colors <- brewer.pal(12,"Paired")
cols=colors[c(1,2,3,4,5,6,7,8)]


##############################################################################
source("../script/xpclr_manhanttan_plot.R")
##############################################################################


############################################################
###tra
##tra-dav(ref)
png("tra.tra_dav.manhattan.png",width=16,height=2.5,units='in',res=400)
par(mar=c(4,5,1,1))
tra_tra_dav=data.frame(cbind(CHR=as.numeric(summary$chr),BP=as.numeric(as.character(summary$Pos)),P=as.numeric(as.character(summary$xpclr_tra_dav_tra))))
tra_tra_dav_new=tra_tra_dav[!is.nan(tra_tra_dav$P),]
tra_tra_dav_new$SNP=paste(tra_tra_dav_new$CHR,tra_tra_dav_new$BP,sep=":")
manhattan(tra_tra_dav_new,chr="CHR",bp="BP",p="P",highligt=NA,xlabel="Chromosome",ylabel="XPCLR",col=c(cols[4],cols[3]),mar_value=c(4,5,1,1)+0.1,cex_value=1.5)
tra_tra_dav_cutoff=quantile(tra_tra_dav_new$P,0.99)
abline(h=tra_tra_dav_cutoff,lty=1,col="red",lwd=2)
dev.off()

##tra-trs(ref)
png("tra.tra_trs.manhattan.png",width=16,height=2.5,units='in',res=400)
par(mar=c(4,5,1,1))
tra_tra_trs=data.frame(cbind(CHR=as.numeric(summary$chr),BP=as.numeric(as.character(summary$Pos)),P=as.numeric(as.character(summary$xpclr_tra_trs_tra))))
tra_tra_trs_new=tra_tra_trs[!is.nan(tra_tra_trs$P),]
tra_tra_trs_new$SNP=paste(tra_tra_trs_new$CHR,tra_tra_trs_new$BP,sep=":")
manhattan(tra_tra_trs_new,chr="CHR",bp="BP",p="P",highligt=NA,xlabel="Chromosome",ylabel="XPCLR",col=c(cols[4],cols[3]),mar_value=c(4,5,1,1)+0.1,cex_value=1.5)
tra_tra_trs_cutoff=quantile(tra_tra_trs_new$P,0.99)
abline(h=tra_tra_trs_cutoff,lty=1,col="red",lwd=2)
dev.off()


####################
###check intersect of the two dataset
tra_tra_dav_outlier=tra_tra_dav_new[which(tra_tra_dav_new$P>tra_tra_dav_cutoff),]
tra_tra_trs_outlier=tra_tra_trs_new[which(tra_tra_trs_new$P>tra_tra_trs_cutoff),]
tra=unique(c(tra_tra_dav_outlier$SNP,tra_tra_trs_outlier$SNP))

############################################################
###dav
##dav-tra(ref)
png("dav.tra_dav.manhattan.png",width=16,height=2.5,units='in',res=400)
par(mar=c(4,5,1,1))
dav_tra_dav=data.frame(cbind(CHR=as.numeric(summary$chr),BP=as.numeric(as.character(summary$Pos)),P=as.numeric(as.character(summary$xpclr_tra_dav_dav))))
dav_tra_dav_new=dav_tra_dav[!is.nan(dav_tra_dav$P),]
dav_tra_dav_new$SNP=paste(dav_tra_dav_new$CHR,dav_tra_dav_new$BP,sep=":")
manhattan(dav_tra_dav_new,chr="CHR",bp="BP",p="P",highligt=NA,xlabel="Chromosome",ylabel="XPCLR",col=c(cols[2],cols[1]),mar_value=c(4,5,1,1)+0.1,cex_value=1.5)
dav_tra_dav_cutoff=quantile(dav_tra_dav_new$P,0.99)
abline(h=dav_tra_dav_cutoff,lty=1,col="red",lwd=2)
dev.off()

##dav-trs(ref)
png("dav.dav_trs.manhattan.png",width=16,height=2.5,units='in',res=400)
par(mar=c(4,5,1,1))
dav_dav_trs=data.frame(cbind(CHR=as.numeric(summary$chr),BP=as.numeric(as.character(summary$Pos)),P=as.numeric(as.character(summary$xpclr_trs_dav_dav))))
dav_dav_trs_new=dav_dav_trs[!is.nan(dav_dav_trs$P),]
dav_dav_trs_new$SNP=paste(dav_dav_trs_new$CHR,dav_dav_trs_new$BP,sep=":")
manhattan(dav_dav_trs_new,chr="CHR",bp="BP",p="P",highligt=NA,xlabel="Chromosome",ylabel="XPCLR",col=c(cols[2],cols[1]),mar_value=c(4,5,1,1)+0.1,cex_value=1.5)
dav_dav_trs_cutoff=quantile(dav_dav_trs_new$P,0.99)
abline(h=dav_dav_trs_cutoff,lty=1,col="red",lwd=2)
dev.off()

############################################
####################
dav_tra_dav_outlier=dav_tra_dav_new[which(dav_tra_dav_new$P>dav_tra_dav_cutoff),]
dav_dav_trs_outlier=dav_dav_trs_new[which(dav_dav_trs_new$P>dav_dav_trs_cutoff),]
dav=unique(c(dav_tra_dav_outlier$SNP,dav_dav_trs_outlier$SNP))

############################################################
###trs
##trs-tra(ref)
png("trs.tra_trs.manhattan.png",width=16,height=2.5,units='in',res=400)
par(mar=c(4,5,1,1))
trs_tra_trs=data.frame(cbind(CHR=as.numeric(summary$chr),BP=as.numeric(as.character(summary$Pos)),P=as.numeric(as.character(summary$xpclr_tra_trs_trs))))
trs_tra_trs_new=trs_tra_trs[!is.nan(trs_tra_trs$P),]
trs_tra_trs_new$SNP=paste(trs_tra_trs_new$CHR,trs_tra_trs_new$BP,sep=":")
manhattan(trs_tra_trs_new,chr="CHR",bp="BP",p="P",highligt=NA,xlabel="Chromosome",ylabel="XPCLR",col=c(cols[6],cols[5]),mar_value=c(4,5,1,1)+0.1,cex_value=1.5)
trs_tra_trs_cutoff=quantile(trs_tra_trs_new$P,0.99)
abline(h=trs_tra_trs_cutoff,lty=1,col="red",lwd=2)
dev.off()

##trs-dav(ref)
png("trs.dav_trs.manhattan.png",width=16,height=2.5,units='in',res=400)
par(mar=c(4,5,1,1))
trs_dav_trs=data.frame(cbind(CHR=as.numeric(summary$chr),BP=as.numeric(as.character(summary$Pos)),P=as.numeric(as.character(summary$xpclr_trs_dav_trs))))
trs_dav_trs_new=trs_dav_trs[!is.nan(trs_dav_trs$P),]
trs_dav_trs_new$SNP=paste(trs_dav_trs_new$CHR,trs_dav_trs_new$BP,sep=":")
manhattan(trs_dav_trs_new,chr="CHR",bp="BP",p="P",highligt=NA,xlabel="Chromosome",ylabel="XPCLR",col=c(cols[6],cols[5]),mar_value=c(4,5,1,1)+0.1,cex_value=1.5)
trs_dav_trs_cutoff=quantile(trs_dav_trs_new$P,0.99)
abline(h=trs_dav_trs_cutoff,lty=1,col="red",lwd=2)
dev.off()

############################################
####################
trs_tra_trs_outlier=trs_tra_trs_new[which(trs_tra_trs_new$P>trs_tra_trs_cutoff),]
trs_dav_trs_outlier=trs_dav_trs_new[which(trs_dav_trs_new$P>trs_dav_trs_cutoff),]
trs=unique(c(trs_tra_trs_outlier$SNP,trs_dav_trs_outlier$SNP))


#outlier SNPs with extremely high beta for each species
three_outlier=venn.diagram(list(a=tra,
                              b=dav,
                              c=trs),
                         fill=c(cols[3],cols[1],cols[5]),
                         category.names=c("P.tra","P.dav","P.trs"),
                         filename=NULL,cat.cex = 2.5,cex = 2.5,margin = 0.05,
                         cat.fontface = "bold",fontface = "bold")

png("VennDiagram.3XPCLR.png",height=6,width=6,units='in',res=600)
grid.draw(three_outlier)
dev.off()


##############################################################################
setwd("/Users/jiwa0011/Dropbox/davidiana_paper/data/selection/XPCLR/results/significant_snps/")

###994windows
xpclr_outlier=summary[which(summary$xpclr_tra_dav_tra>tra_tra_dav_cutoff | 
          summary$xpclr_tra_dav_dav>dav_tra_dav_cutoff |
          summary$xpclr_tra_trs_tra>tra_tra_trs_cutoff |
          summary$xpclr_tra_trs_trs>trs_tra_trs_cutoff |
          summary$xpclr_trs_dav_dav>dav_dav_trs_cutoff |
          summary$xpclr_trs_dav_trs>trs_dav_trs_cutoff
          ),]

###20,571windows
xpclr_non_outlier=summary[!which(summary$xpclr_tra_dav_tra>tra_tra_dav_cutoff | 
                                 summary$xpclr_tra_dav_dav>dav_tra_dav_cutoff |
                                 summary$xpclr_tra_trs_tra>tra_tra_trs_cutoff |
                                 summary$xpclr_tra_trs_trs>trs_tra_trs_cutoff |
                                 summary$xpclr_trs_dav_dav>dav_dav_trs_cutoff |
                                 summary$xpclr_trs_dav_trs>trs_dav_trs_cutoff
),]

########################
#write the windows and genes into txt file
write.table(xpclr_outlier,file="xpclr.outlier.txt",sep="\t", quote=F, row.names=F, col.names=T)


##############################################################################
##############################################################################
# * P<0.05
# ** P<0.01
# *** P<0.001

###1. Diversity
pi_xpclr=c(xpclr_non_outlier$tra_tP,xpclr_outlier$tra_tP,
         xpclr_non_outlier$dav_tP,xpclr_outlier$dav_tP,
         xpclr_non_outlier$trs_tP,xpclr_outlier$trs_tP)

wilcox.test(xpclr_non_outlier$tra_tP,xpclr_outlier$tra_tP)
wilcox.test(xpclr_non_outlier$dav_tP,xpclr_outlier$dav_tP)
wilcox.test(xpclr_non_outlier$trs_tP,xpclr_outlier$trs_tP)


pi_xpclr_group=c(rep("P.tra",21565),rep("P.dav",21565),rep("P.trs",21565))
pi_xpclr_fill=c(rep("xpclr_non_outlier",nrow(xpclr_non_outlier)),
              rep("xpclr_outlier",nrow(xpclr_outlier)),
              rep("xpclr_non_outlier",nrow(xpclr_non_outlier)),
              rep("xpclr_outlier",nrow(xpclr_outlier)),
              rep("xpclr_non_outlier",nrow(xpclr_non_outlier)),
              rep("xpclr_outlier",nrow(xpclr_outlier)))

pi_xpclr_table=as.data.frame(cbind(pi_xpclr,pi_xpclr_group,pi_xpclr_fill))

pi_xpclr_table$pi_xpclr=as.numeric(as.character(pi_xpclr_table$pi_xpclr))
pi_xpclr_table$pi_xpclr_group=factor(pi_xpclr_table$pi_xpclr_group,levels=c("P.tra","P.dav","P.trs"))
pi_xpclr_table$pi_xpclr_fill=factor(pi_xpclr_table$pi_xpclr_fill,levels=c("xpclr_non_outlier","xpclr_outlier"))

ggplot(pi_xpclr_table, aes(x=pi_xpclr_group, y=pi_xpclr,fill=pi_xpclr_fill))+ 
  geom_boxplot(outlier.shape=NA)+scale_fill_manual(values=c("grey60","red"))+
  scale_y_continuous(limits=c(0,0.03))+
  #xlab("")+ylab(expression(paste(xpclr," (P.tra)",sep="")))+
  xlab("")+ylab("")+
  annotate(geom="text",x=c(1,2,3),y=c(0.029,0.029,0.029),label=c(rep("***",3)),size=8)+
  theme(legend.position="none",axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),axis.text.x = element_blank(),axis.text.y = element_text(size=16))

ggsave("xpclr.pi.boxplot.png",width=3.5,height=4.5,units='in',dpi=500)

##############################################################################
###2. recombination
rho_xpclr=c(xpclr_non_outlier$tra_rho_mean,xpclr_outlier$tra_rho_mean,
          xpclr_non_outlier$dav_rho_mean,xpclr_outlier$dav_rho_mean,
          xpclr_non_outlier$trs_rho_mean,xpclr_outlier$trs_rho_mean)

wilcox.test(xpclr_outlier$tra_rho_mean,xpclr_non_outlier$tra_rho_mean)
wilcox.test(xpclr_outlier$dav_rho_mean,xpclr_non_outlier$dav_rho_mean)
wilcox.test(xpclr_outlier$trs_rho_mean,xpclr_non_outlier$trs_rho_mean)


rho_xpclr_group=c(rep("P.tra",21565),rep("P.dav",21565),rep("P.trs",21565))
rho_xpclr_fill=c(rep("xpclr_non_outlier",nrow(xpclr_non_outlier)),
               rep("xpclr_outlier",nrow(xpclr_outlier)),
               rep("xpclr_non_outlier",nrow(xpclr_non_outlier)),
               rep("xpclr_outlier",nrow(xpclr_outlier)),
               rep("xpclr_non_outlier",nrow(xpclr_non_outlier)),
               rep("xpclr_outlier",nrow(xpclr_outlier)))

rho_xpclr_table=as.data.frame(cbind(rho_xpclr,rho_xpclr_group,rho_xpclr_fill))

rho_xpclr_table$rho_xpclr=as.numeric(as.character(rho_xpclr_table$rho_xpclr))
rho_xpclr_table$rho_xpclr_group=factor(rho_xpclr_table$rho_xpclr_group,levels=c("P.tra","P.dav","P.trs"))
rho_xpclr_table$rho_xpclr_fill=factor(rho_xpclr_table$rho_xpclr_fill,levels=c("xpclr_non_outlier","xpclr_outlier"))

ggplot(rho_xpclr_table, aes(x=rho_xpclr_group, y=rho_xpclr,fill=rho_xpclr_fill))+ 
  geom_boxplot(outlier.shape=NA)+scale_fill_manual(values=c("grey60","red"))+
  scale_y_continuous(limits=c(0,0.07))+
  #xlab("")+ylab(expression(paste(xpclr," (P.tra)",sep="")))+
  xlab("")+ylab("")+
  annotate(geom="text",x=c(1,2,3),y=c(0.068,0.068,0.068),label=c(rep("***",3)),size=8)+
  theme(legend.position="none",axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),axis.text.x = element_blank(),axis.text.y = element_text(size=16))

ggsave("xpclr.rho.boxplot.png",width=3.5,height=4.5,units='in',dpi=500)


##############################################################################
###3. Divergence, Fst
fst_xpclr=c(xpclr_non_outlier$tra_dav_fst,xpclr_outlier$tra_dav_fst,
          xpclr_non_outlier$tra_trs_fst,xpclr_outlier$tra_trs_fst,
          xpclr_non_outlier$dav_trs_fst,xpclr_outlier$dav_trs_fst)

wilcox.test(xpclr_non_outlier$tra_dav_fst,xpclr_outlier$tra_dav_fst)
wilcox.test(xpclr_non_outlier$tra_trs_fst,xpclr_outlier$tra_trs_fst)
wilcox.test(xpclr_non_outlier$dav_trs_fst,xpclr_outlier$dav_trs_fst)


fst_xpclr_group=c(rep("P.tra-P.dav",21565),rep("P.tra-P.trs",21565),rep("P.dav-P.trs",21565))
fst_xpclr_fill=c(rep("xpclr_non_outlier",nrow(xpclr_non_outlier)),
               rep("xpclr_outlier",nrow(xpclr_outlier)),
               rep("xpclr_non_outlier",nrow(xpclr_non_outlier)),
               rep("xpclr_outlier",nrow(xpclr_outlier)),
               rep("xpclr_non_outlier",nrow(xpclr_non_outlier)),
               rep("xpclr_outlier",nrow(xpclr_outlier)))

fst_xpclr_table=as.data.frame(cbind(fst_xpclr,fst_xpclr_group,fst_xpclr_fill))

fst_xpclr_table$fst_xpclr=as.numeric(as.character(fst_xpclr_table$fst_xpclr))
fst_xpclr_table$fst_xpclr_group=factor(fst_xpclr_table$fst_xpclr_group,levels=c("P.tra-P.dav","P.tra-P.trs","P.dav-P.trs"))
fst_xpclr_table$fst_xpclr_fill=factor(fst_xpclr_table$fst_xpclr_fill,levels=c("xpclr_non_outlier","xpclr_outlier"))

ggplot(fst_xpclr_table, aes(x=fst_xpclr_group, y=fst_xpclr,fill=fst_xpclr_fill))+ 
  geom_boxplot(outlier.shape=NA)+scale_fill_manual(values=c("grey60","red"))+
  scale_y_continuous(limits=c(0,1))+
  #xlab("")+ylab(expression(paste(xpclr," (P.tra)",sep="")))+
  xlab("")+ylab("")+
  annotate(geom="text",x=c(1,2,3),y=c(0.98,0.98,0.98),label=c("***","***","***"),size=8)+
  theme(legend.position="none",axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),axis.text.x = element_blank(),axis.text.y = element_text(size=16))

ggsave("xpclr.fst.boxplot.png",width=3.5,height=4.5,units='in',dpi=500)

##############################################################################
###4. Divergence, dxy
dxy_xpclr=c(xpclr_non_outlier$tra_dav_dxy,xpclr_outlier$tra_dav_dxy,
          xpclr_non_outlier$tra_trs_dxy,xpclr_outlier$tra_trs_dxy,
          xpclr_non_outlier$dav_trs_dxy,xpclr_outlier$dav_trs_dxy)

wilcox.test(xpclr_non_outlier$tra_dav_dxy,xpclr_outlier$tra_dav_dxy)
wilcox.test(xpclr_non_outlier$tra_trs_dxy,xpclr_outlier$tra_trs_dxy)
wilcox.test(xpclr_non_outlier$dav_trs_dxy,xpclr_outlier$dav_trs_dxy)


dxy_xpclr_group=c(rep("P.tra-P.dav",21565),rep("P.tra-P.trs",21565),rep("P.dav-P.trs",21565))
dxy_xpclr_fill=c(rep("xpclr_non_outlier",nrow(xpclr_non_outlier)),
               rep("xpclr_outlier",nrow(xpclr_outlier)),
               rep("xpclr_non_outlier",nrow(xpclr_non_outlier)),
               rep("xpclr_outlier",nrow(xpclr_outlier)),
               rep("xpclr_non_outlier",nrow(xpclr_non_outlier)),
               rep("xpclr_outlier",nrow(xpclr_outlier)))

dxy_xpclr_table=as.data.frame(cbind(dxy_xpclr,dxy_xpclr_group,dxy_xpclr_fill))

dxy_xpclr_table$dxy_xpclr=as.numeric(as.character(dxy_xpclr_table$dxy_xpclr))
dxy_xpclr_table$dxy_xpclr_group=factor(dxy_xpclr_table$dxy_xpclr_group,levels=c("P.tra-P.dav","P.tra-P.trs","P.dav-P.trs"))
dxy_xpclr_table$dxy_xpclr_fill=factor(dxy_xpclr_table$dxy_xpclr_fill,levels=c("xpclr_non_outlier","xpclr_outlier"))

ggplot(dxy_xpclr_table, aes(x=dxy_xpclr_group, y=dxy_xpclr,fill=dxy_xpclr_fill))+ 
  geom_boxplot(outlier.shape=NA)+scale_fill_manual(values=c("grey60","red"))+
  scale_y_continuous(limits=c(0,0.05))+
  #xlab("")+ylab(expression(paste(xpclr," (P.tra)",sep="")))+
  xlab("")+ylab("")+
  annotate(geom="text",x=c(1,2,3),y=c(0.049,0.049,0.049),label=c("n.s.","n.s.","n.s."),size=8)+
  theme(legend.position="none",axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),axis.text.x = element_blank(),axis.text.y = element_text(size=16))

ggsave("xpclr.dxy.boxplot.png",width=3.5,height=4.5,units='in',dpi=500)


##############################################################################
###5. topology
topo_xpclr=c(xpclr_non_outlier$top1,xpclr_outlier$top1,
           xpclr_non_outlier$top2,xpclr_outlier$top2,
           xpclr_non_outlier$top3,xpclr_outlier$top3)

wilcox.test(xpclr_non_outlier$top1,xpclr_outlier$top1)
wilcox.test(xpclr_non_outlier$top2,xpclr_outlier$top2)
wilcox.test(xpclr_non_outlier$top3,xpclr_outlier$top3)

# * P<0.05
# ** P<0.01
# *** P<0.001


topo_xpclr_group=c(rep("topo1",21565),rep("topo2",21565),rep("topo3",21565))
topo_xpclr_fill=c(rep("xpclr_non_outlier",nrow(xpclr_non_outlier)),
                rep("xpclr_outlier",nrow(xpclr_outlier)),
                rep("xpclr_non_outlier",nrow(xpclr_non_outlier)),
                rep("xpclr_outlier",nrow(xpclr_outlier)),
                rep("xpclr_non_outlier",nrow(xpclr_non_outlier)),
                rep("xpclr_outlier",nrow(xpclr_outlier)))

topo_xpclr_table=as.data.frame(cbind(topo_xpclr,topo_xpclr_group,topo_xpclr_fill))

topo_xpclr_table$topo_xpclr=as.numeric(as.character(topo_xpclr_table$topo_xpclr))
topo_xpclr_table$topo_xpclr_group=factor(topo_xpclr_table$topo_xpclr_group,levels=c("topo1","topo2","topo3"))
topo_xpclr_table$topo_xpclr_fill=factor(topo_xpclr_table$topo_xpclr_fill,levels=c("xpclr_non_outlier","xpclr_outlier"))

ggplot(topo_xpclr_table, aes(x=topo_xpclr_group, y=topo_xpclr,fill=topo_xpclr_fill))+ 
  geom_boxplot(outlier.shape=NA)+scale_fill_manual(values=c("grey60","red"))+
  scale_y_continuous(limits=c(0,1.1))+
  #xlab("")+ylab(expression(paste(xpclr," (P.tra)",sep="")))+
  xlab("")+ylab("")+
  annotate(geom="text",x=c(1,2,3),y=c(1.08,1.08,1.08),label=c("***","n.s.","***"),size=8)+
  theme(legend.position="none",axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),axis.text.x = element_blank(),axis.text.y = element_text(size=16))

ggsave("xpclr.topo.boxplot.png",width=3.5,height=4.5,units='in',dpi=500)

##############################################################################
###6.ILS
wilcox.test(xpclr_non_outlier$nls_n,xpclr_outlier$nls_n)
ILS_melt_xpclr=as.data.frame(rbind(cbind(xpclr_non_outlier$nls_n,rep("non_outlier",20571)),
                                    cbind(xpclr_outlier$nls_n,rep("outlier",994))))
names(ILS_melt_xpclr)=c("value","variable")
ILS_melt_xpclr$value=as.numeric(as.character(ILS_melt_xpclr$value))

ggplot(ILS_melt_xpclr, aes(x=variable, y=value,fill=variable))+ 
  geom_boxplot(outlier.shape=NA)+scale_fill_manual(values=c("grey60","red"))+
  scale_y_continuous(limits=c(0,0.3))+
  #xlab("")+ylab(expression(paste(xpclr," (P.tra)",sep="")))+
  xlab("")+ylab("")+
  annotate(geom="text",x=1.5,y=0.29,label=c("***"),size=8)+
  theme(legend.position="none",axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),axis.text.x = element_blank(),axis.text.y = element_text(size=16))

ggsave("xpclr.ILS.boxplot.png",width=2,height=4,units='in',dpi=500)


##############################################################################
###7.fd

xpclr_outlier_fd=xpclr_outlier[which(xpclr_outlier$fd>=0),]
xpclr_non_outlier_fd=xpclr_non_outlier[which(xpclr_non_outlier$fd>=0),]
wilcox.test(xpclr_non_outlier_fd$fd,xpclr_outlier_fd$fd)

fd_melt_xpclr=as.data.frame(rbind(cbind(xpclr_non_outlier_fd$nls_n,rep("non_outlier",12894)),
                                   cbind(xpclr_outlier_fd$nls_n,rep("outlier",609))))
names(fd_melt_xpclr)=c("value","variable")
fd_melt_xpclr$value=as.numeric(as.character(fd_melt_xpclr$value))

wilcox.test(xpclr_non_outlier$fd,xpclr_outlier$fd)
fd_melt_xpclr=as.data.frame(rbind(cbind(xpclr_non_outlier$fd,rep("non_outlier",20571)),
                                   cbind(xpclr_outlier$fd,rep("outlier",994))))
names(fd_melt_xpclr)=c("value","variable")
fd_melt_xpclr$value=as.numeric(as.character(fd_melt_xpclr$value))

ggplot(fd_melt_xpclr, aes(x=variable, y=value,fill=variable))+ 
  geom_boxplot(outlier.shape=NA)+scale_fill_manual(values=c("grey60","red"))+
  scale_y_continuous(limits=c(-0.3,0.5))+
  #xlab("")+ylab(expression(paste(xpclr," (P.tra)",sep="")))+
  xlab("")+ylab("")+
  annotate(geom="text",x=1.5,y=0.49,label=c("n.s."),size=8)+
  theme(legend.position="none",axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),axis.text.x = element_blank(),axis.text.y = element_text(size=16))

ggsave("xpclr.fd.boxplot.png",width=2,height=4,units='in',dpi=500)



#########################################################################
########Output the windows in the data
gene_data=fread("~/Dropbox/davidiana_paper/data/SnIPRE/scripts/out/trichocarpa.gene.filter.bed.txt",header=F)
names(gene_data)=c("Chr","Start","End","Gene")

gene_data$Chr=gsub("Chr0","",gene_data$Chr)
gene_data$Chr=gsub("Chr","",gene_data$Chr)



##############################################################################
##############################################################################
##Extract the genes within the selected regions and do GO enrichment analysis

n_gene=nrow(xpclr_outlier)
genes_xpclr_outlier=c()
for (i in 1:n_gene){
  genes_sf2=(strsplit(xpclr_outlier$geneIDs[i],","))[[1]]
  if (genes_sf2==".") {
    next
  }
  genes_xpclr_outlier=c(genes_xpclr_outlier,genes_sf2)
}

outlier_gene=unique(genes_xpclr_outlier) ##976 genes

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
gene_table=GenTable(GOdata,Fisher.p=restRes,topNodes=69)
###because there is 1395 GO terms and we can calculated the adjust P-values
gene_table_total=GenTable(GOdata,Fisher.p=restRes,topNodes=1395)
allGO=usedGO(GOdata)
pvalues=gene_table_total$Fisher.p
gene_table$adjust.p=round(head(p.adjust(pvalues,method="BH"),69),10)

write.table(gene_table,file="xpclr.outlier.GO.enrich.txt",sep="\t", quote=F, row.names=F, col.names=T)




