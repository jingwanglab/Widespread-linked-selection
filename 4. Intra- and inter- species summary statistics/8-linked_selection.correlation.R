library(data.table)
library(RColorBrewer)
library(dplyr)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(ppcor)

colors <- brewer.pal(12,"Paired")

setwd("~/Dropbox/davidiana_paper/data/linked_selection")

##########################################################################
###read 10kb
#####Read in the results from many analysis and to estimate the correlations among them

###1.Summary of windows containing the selection and constraint effects of genes located within windows
###2.Introgression detection: D,fd
###3.Phylogenetic tree(twisst):weightings, top1,top2,top3
###4.Incomplete lineage sorting estimated from allele freqeuncy data


###read1
summary_10kb=fread("data/summary.24235win.div.fd.gene.txt",header=T)

###read2
#gatk
D_gatk_10kb=read.csv("../abba_baba/gatk_data/4species.gatk.hap.snp.rm_indel.bed.biallelic.DP5.GQ10.rm_missing.recode.gt.beagle.10kb.abbababa.csv",header=T)
#angsd
abbababa_angsd_10kb=fread("../abba_baba/angsd/10kb/4species.all.outDstat.Observed.txt",header=T)
abbababa_angsd_10kb_sites=fread(input = 'zcat < ../abba_baba/angsd/10kb/4species.abbababa.all.10000.out.abbababa2.gz',header=T)
abbababa_angsd_10kb$numSites=abbababa_angsd_10kb_sites$numSites
D_angsd_10kb=abbababa_angsd_10kb %>% filter(numSites>1000)

###read3
twisst_10kb=read.table("../twisst/data/4species.gatk.hap.snp.rm_indel.bed.biallelic.DP5.GQ10.rm_missing.recode.gt.beagle.twisst.weights.10kb.csv")
twisst_10kb_new=twisst_10kb[-1,]
names(twisst_10kb_new)=c("top1","top2","top3")

indx <- sapply(twisst_10kb_new, is.factor)
twisst_10kb_new[indx] <- lapply(twisst_10kb_new[indx], function(x) as.numeric(as.character(x)))
colSums(twisst_10kb_new)/sum(colSums(twisst_10kb_new))

twisst_10kb_prop=prop.table(as.matrix(twisst_10kb_new),1)
###twisst sites
twisst_10kb_sites=read.table("../twisst/data/4species.gatk.hap.snp.rm_indel.bed.biallelic.DP5.GQ10.rm_missing.recode.gt.beagle.raxml.w10kb.data.tsv",header=T)
twisst_10kb_sites_new=twisst_10kb_sites[which(twisst_10kb_sites$sites>50),]
twisst_10kb_prop_sites=cbind(twisst_10kb_sites_new,twisst_10kb_prop)

###read4
nls_10kb=fread("../snp_allele_freq/allele_frq_aspen_10kb.nls.txt",header=T)
nls_10kb_new=nls_10kb %>% filter(n_sites>50)


##########################################################################
###choose common windows
summary_10kb$win=paste(summary_10kb$Chr,summary_10kb$Pos+5000,sep=":")
D_gatk_10kb$win=paste(D_gatk_10kb$scaffold,D_gatk_10kb$end,sep=":")
D_angsd_10kb$win=paste(D_angsd_10kb$CHR,D_angsd_10kb$BLOCKend,sep=":")
twisst_10kb_prop_sites$win=paste(twisst_10kb_prop_sites$scaffold,twisst_10kb_prop_sites$end,sep=":")
nls_10kb_new$win=paste(nls_10kb_new$chro,nls_10kb_new$end,sep=":")

common_win=Reduce(intersect,list(summary_10kb$win,
                                 D_gatk_10kb$win,
                                 D_angsd_10kb$win,
                                 twisst_10kb_prop_sites$win,
                                 nls_10kb_new$win
))

###new dataset results with common windows
summary_10kb_w=summary_10kb[which(summary_10kb$win %in% common_win),]
D_gatk_10kb_w=D_gatk_10kb[which(D_gatk_10kb$win %in% common_win),]
D_angsd_10kb_w=D_angsd_10kb[which(D_angsd_10kb$win %in% common_win),]
twisst_10kb_prop_sites_w=twisst_10kb_prop_sites[which(twisst_10kb_prop_sites$win %in% common_win),]
nls_10kb_new_w=nls_10kb_new[which(nls_10kb_new$win %in% common_win),]


total_10kb=cbind(summary_10kb_w,D_gatk_10kb_w$D,D_gatk_10kb_w$fd,
                 twisst_10kb_prop_sites_w$top1,twisst_10kb_prop_sites_w$top2,twisst_10kb_prop_sites_w$top3,
                 nls_10kb_new_w$nls_n)
names(total_10kb)=c(names(summary_10kb_w),"D","fd","top1","top2","top3","nls")

##########################################################################
##########################################################################

###read100kb
summary_100kb=fread("data/summary.div.diverg.rho.100kb.txt",header=T)

###read2
#gatk
D_gatk_100kb=read.csv("../abba_baba/gatk_data/4species.gatk.hap.snp.rm_indel.bed.biallelic.DP5.GQ10.rm_missing.recode.gt.beagle.100kb.abbababa.csv",header=T)
#angsd
abbababa_angsd_100kb=fread("../abba_baba/angsd/100kb/4species.all.outDstat.Observed.txt",header=T)
abbababa_angsd_100kb_sites=fread(input = 'zcat < ../abba_baba/angsd/100kb/4species.abbababa.all.100000.out.abbababa2.gz',header=T)
abbababa_angsd_100kb$numSites=abbababa_angsd_100kb_sites$numSites
D_angsd_100kb=abbababa_angsd_100kb %>% filter(numSites>10000)

###read3
twisst_100kb=read.table("../twisst/data/4species.gatk.hap.snp.rm_indel.bed.biallelic.DP5.GQ10.rm_missing.recode.gt.beagle.twisst.weights.100kb.csv")
twisst_100kb_new=twisst_100kb[-1,]
names(twisst_100kb_new)=c("top1","top2","top3")

indx <- sapply(twisst_100kb_new, is.factor)
twisst_100kb_new[indx] <- lapply(twisst_100kb_new[indx], function(x) as.numeric(as.character(x)))
colSums(twisst_100kb_new)/sum(colSums(twisst_100kb_new))

twisst_100kb_prop=prop.table(as.matrix(twisst_100kb_new),1)
###twisst sites
twisst_100kb_sites=read.table("../twisst/data/4species.gatk.hap.snp.rm_indel.bed.biallelic.DP5.GQ10.rm_missing.recode.gt.beagle.raxml.w100kb.data.tsv",header=T)
twisst_100kb_sites_new=twisst_100kb_sites[which(twisst_100kb_sites$sites>200),]
twisst_100kb_prop_sites=cbind(twisst_100kb_sites_new,twisst_100kb_prop)

###read4
nls_100kb=fread("../snp_allele_freq/allele_frq_aspen_100kb.nls.txt",header=T)
nls_100kb_new=nls_100kb %>% filter(n_sites>200)


##########################################################################
###choose common windows
summary_100kb$win=paste(summary_100kb$Chr,summary_100kb$Pos+50000,sep=":")
D_gatk_100kb$win=paste(D_gatk_100kb$scaffold,D_gatk_100kb$end,sep=":")
D_angsd_100kb$win=paste(D_angsd_100kb$CHR,D_angsd_100kb$BLOCKend,sep=":")
twisst_100kb_prop_sites$win=paste(twisst_100kb_prop_sites$scaffold,twisst_100kb_prop_sites$end,sep=":")
nls_100kb_new$win=paste(nls_100kb_new$chro,nls_100kb_new$end,sep=":")

common_win=Reduce(intersect,list(summary_100kb$win,
                                 D_gatk_100kb$win,
                                 D_angsd_100kb$win,
                                 twisst_100kb_prop_sites$win,
                                 nls_100kb_new$win
))

###new dataset results with common windows
summary_100kb_w=summary_100kb[which(summary_100kb$win %in% common_win),]
D_gatk_100kb_w=D_gatk_100kb[which(D_gatk_100kb$win %in% common_win),]
D_angsd_100kb_w=D_angsd_100kb[which(D_angsd_100kb$win %in% common_win),]
twisst_100kb_prop_sites_w=twisst_100kb_prop_sites[which(twisst_100kb_prop_sites$win %in% common_win),]
nls_100kb_new_w=nls_100kb_new[which(nls_100kb_new$win %in% common_win),]

total_100kb=cbind(summary_100kb_w,D_gatk_100kb_w$D,D_gatk_100kb_w$fd,
                 twisst_100kb_prop_sites_w$top1,twisst_100kb_prop_sites_w$top2,twisst_100kb_prop_sites_w$top3,
                 nls_100kb_new_w$nls_n)
names(total_100kb)=c(names(summary_100kb_w),"D","fd","top1","top2","top3","nls")


##########################################################################
##########################################################################
total_10kb <- na.omit(total_10kb)
total_100kb <- na.omit(total_100kb)

##########################################################################
##############Correlation with coding density##############
#1.1cd vs. tra-tP
cor.test(total_10kb$Coding_prop,total_10kb$tra_tP,method="spearman")
pcor.test(total_10kb$Coding_prop,total_10kb$tra_tP,total_10kb[,c("tra_rho_mean","GC")],method="spearman")

cor.test(total_100kb$Coding_prop,total_100kb$tra_tP,method="spearman")
pcor.test(total_100kb$Coding_prop,total_100kb$tra_tP,total_100kb[,c("tra_rho_mean","GC")],method="spearman")

#1.1cd vs. dav-tP
cor.test(total_10kb$Coding_prop,total_10kb$dav_tP,method="spearman")
pcor.test(total_10kb$Coding_prop,total_10kb$dav_tP,total_10kb[,c("dav_rho_mean","GC")],method="spearman")

cor.test(total_100kb$Coding_prop,total_100kb$dav_tP,method="spearman")
pcor.test(total_100kb$Coding_prop,total_100kb$dav_tP,total_100kb[,c("dav_rho_mean","GC")],method="spearman")


#1.1cd vs. trs-tP
cor.test(total_10kb$Coding_prop,total_10kb$trs_tP,method="spearman")
pcor.test(total_10kb$Coding_prop,total_10kb$trs_tP,total_10kb[,c("trs_rho_mean","GC")],method="spearman")

cor.test(total_100kb$Coding_prop,total_100kb$trs_tP,method="spearman")
pcor.test(total_100kb$Coding_prop,total_100kb$trs_tP,total_100kb[,c("trs_rho_mean","GC")],method="spearman")

#1.1cd vs. tri-tP
cor.test(total_10kb$Coding_prop,total_10kb$tri_tP,method="spearman")
pcor.test(total_10kb$Coding_prop,total_10kb$tri_tP,total_10kb[,c("tri_rho_mean","GC")],method="spearman")

cor.test(total_100kb$Coding_prop,total_100kb$tri_tP,method="spearman")
pcor.test(total_100kb$Coding_prop,total_100kb$tri_tP,total_100kb[,c("tri_rho_mean","GC")],method="spearman")


##########################################################################
##1.2 cd vs. fst
#tra-dav
cor.test(total_10kb$Coding_prop,total_10kb$tra_dav_fst,method="spearman")
pcor.test(total_10kb$Coding_prop,total_10kb$tra_dav_fst,total_10kb[,c("tra_rho_mean","GC")],method="spearman")

cor.test(total_100kb$Coding_prop,total_100kb$tra_dav_fst,method="spearman")
pcor.test(total_100kb$Coding_prop,total_100kb$tra_dav_fst,total_100kb[,c("tra_rho_mean","GC")],method="spearman")

#tra-trs
cor.test(total_10kb$Coding_prop,total_10kb$tra_trs_fst,method="spearman")
pcor.test(total_10kb$Coding_prop,total_10kb$tra_trs_fst,total_10kb[,c("tra_rho_mean","GC")],method="spearman")

cor.test(total_100kb$Coding_prop,total_100kb$tra_trs_fst,method="spearman")
pcor.test(total_100kb$Coding_prop,total_100kb$tra_trs_fst,total_100kb[,c("tra_rho_mean","GC")],method="spearman")

#dav-trs
cor.test(total_10kb$Coding_prop,total_10kb$dav_trs_fst,method="spearman")
pcor.test(total_10kb$Coding_prop,total_10kb$dav_trs_fst,total_10kb[,c("tra_rho_mean","GC")],method="spearman")

cor.test(total_100kb$Coding_prop,total_100kb$dav_trs_fst,method="spearman")
pcor.test(total_100kb$Coding_prop,total_100kb$dav_trs_fst,total_100kb[,c("tra_rho_mean","GC")],method="spearman")

#tra-tri
cor.test(total_10kb$Coding_prop,total_10kb$tra_tri_fst,method="spearman")
pcor.test(total_10kb$Coding_prop,total_10kb$tra_tri_fst,total_10kb[,c("tra_rho_mean","GC")],method="spearman")

cor.test(total_100kb$Coding_prop,total_100kb$tra_tri_fst,method="spearman")
pcor.test(total_100kb$Coding_prop,total_100kb$tra_tri_fst,total_100kb[,c("tra_rho_mean","GC")],method="spearman")

#dav-tri
cor.test(total_10kb$Coding_prop,total_10kb$dav_tri_fst,method="spearman")
pcor.test(total_10kb$Coding_prop,total_10kb$dav_tri_fst,total_10kb[,c("tra_rho_mean","GC")],method="spearman")

cor.test(total_100kb$Coding_prop,total_100kb$dav_tri_fst,method="spearman")
pcor.test(total_100kb$Coding_prop,total_100kb$dav_tri_fst,total_100kb[,c("tra_rho_mean","GC")],method="spearman")

#trs-tri
cor.test(total_10kb$Coding_prop,total_10kb$trs_tri_fst,method="spearman")
pcor.test(total_10kb$Coding_prop,total_10kb$trs_tri_fst,total_10kb[,c("tra_rho_mean","GC")],method="spearman")

cor.test(total_100kb$Coding_prop,total_100kb$trs_tri_fst,method="spearman")
pcor.test(total_100kb$Coding_prop,total_100kb$trs_tri_fst,total_100kb[,c("tra_rho_mean","GC")],method="spearman")


##########################################################################
##1.3 cd vs. dxy
#tra-dav
cor.test(total_10kb$Coding_prop,total_10kb$tra_dav_dxy,method="spearman")
pcor.test(total_10kb$Coding_prop,total_10kb$tra_dav_dxy,total_10kb[,c("tra_rho_mean","GC")],method="spearman")

cor.test(total_100kb$Coding_prop,total_100kb$tra_dav_dxy,method="spearman")
pcor.test(total_100kb$Coding_prop,total_100kb$tra_dav_dxy,total_100kb[,c("tra_rho_mean","GC")],method="spearman")

#tra-trs
cor.test(total_10kb$Coding_prop,total_10kb$tra_trs_dxy,method="spearman")
pcor.test(total_10kb$Coding_prop,total_10kb$tra_trs_dxy,total_10kb[,c("tra_rho_mean","GC")],method="spearman")

cor.test(total_100kb$Coding_prop,total_100kb$tra_trs_dxy,method="spearman")
pcor.test(total_100kb$Coding_prop,total_100kb$tra_trs_dxy,total_100kb[,c("tra_rho_mean","GC")],method="spearman")

#dav-trs
cor.test(total_10kb$Coding_prop,total_10kb$dav_trs_dxy,method="spearman")
pcor.test(total_10kb$Coding_prop,total_10kb$dav_trs_dxy,total_10kb[,c("tra_rho_mean","GC")],method="spearman")

cor.test(total_100kb$Coding_prop,total_100kb$dav_trs_dxy,method="spearman")
pcor.test(total_100kb$Coding_prop,total_100kb$dav_trs_dxy,total_100kb[,c("tra_rho_mean","GC")],method="spearman")

#tra-tri
cor.test(total_10kb$Coding_prop,total_10kb$tra_tri_dxy,method="spearman")
pcor.test(total_10kb$Coding_prop,total_10kb$tra_tri_dxy,total_10kb[,c("tra_rho_mean","GC")],method="spearman")

cor.test(total_100kb$Coding_prop,total_100kb$tra_tri_dxy,method="spearman")
pcor.test(total_100kb$Coding_prop,total_100kb$tra_tri_dxy,total_100kb[,c("tra_rho_mean","GC")],method="spearman")

#dav-tri
cor.test(total_10kb$Coding_prop,total_10kb$dav_tri_dxy,method="spearman")
pcor.test(total_10kb$Coding_prop,total_10kb$dav_tri_dxy,total_10kb[,c("tra_rho_mean","GC")],method="spearman")

cor.test(total_100kb$Coding_prop,total_100kb$dav_tri_dxy,method="spearman")
pcor.test(total_100kb$Coding_prop,total_100kb$dav_tri_dxy,total_100kb[,c("tra_rho_mean","GC")],method="spearman")

#trs-tri
cor.test(total_10kb$Coding_prop,total_10kb$trs_tri_dxy,method="spearman")
pcor.test(total_10kb$Coding_prop,total_10kb$trs_tri_dxy,total_10kb[,c("tra_rho_mean","GC")],method="spearman")

cor.test(total_100kb$Coding_prop,total_100kb$trs_tri_dxy,method="spearman")
pcor.test(total_100kb$Coding_prop,total_100kb$trs_tri_dxy,total_100kb[,c("tra_rho_mean","GC")],method="spearman")

##########################################################################
##1.4 cd vs. nls

cor.test(total_10kb$Coding_prop,total_10kb$nls,method="spearman")
pcor.test(total_10kb$Coding_prop,total_10kb$nls,total_10kb[,c("tra_rho_mean","GC")],method="spearman")

cor.test(total_100kb$Coding_prop,total_100kb$nls,method="spearman")
pcor.test(total_100kb$Coding_prop,total_100kb$nls,total_100kb[,c("tra_rho_mean","GC")],method="spearman")

##########################################################################
##1.5 cd vs. top2

cor.test(total_10kb$Coding_prop,total_10kb$top2,method="spearman")
pcor.test(total_10kb$Coding_prop,total_10kb$top2,total_10kb[,c("tra_rho_mean","GC")],method="spearman")

cor.test(total_100kb$Coding_prop,total_100kb$top2,method="spearman")
pcor.test(total_100kb$Coding_prop,total_100kb$top2,total_100kb[,c("tra_rho_mean","GC")],method="spearman")

##########################################################################
##1.6 cd vs. fd
total_10kb_fd=total_10kb[which(total_10kb$fd>0),]
total_100kb_fd=total_100kb[which(total_100kb$fd>0),]

cor.test(total_10kb_fd$Coding_prop,total_10kb_fd$fd,method="spearman")
pcor.test(total_10kb_fd$Coding_prop,total_10kb_fd$fd,total_10kb_fd[,c("tra_rho_mean","GC")],method="spearman")

cor.test(total_100kb_fd$Coding_prop,total_100kb_fd$fd,method="spearman")
pcor.test(total_100kb_fd$Coding_prop,total_100kb_fd$fd,total_100kb_fd[,c("tra_rho_mean","GC")],method="spearman")



##########################################################################
##########################################################################

##############Correlation with recombination rate##############

#1.1rho vs. tra-tP
cor.test(total_10kb$tra_rho_mean,total_10kb$tra_tP,method="spearman")
pcor.test(total_10kb$tra_rho_mean,total_10kb$tra_tP,total_10kb[,c("Coding_prop","GC")],method="spearman")

cor.test(total_100kb$tra_rho_mean,total_100kb$tra_tP,method="spearman")
pcor.test(total_100kb$tra_rho_mean,total_100kb$tra_tP,total_100kb[,c("Coding_prop","GC")],method="spearman")

#1.1rho vs. dav-tP
cor.test(total_10kb$dav_rho_mean,total_10kb$dav_tP,method="spearman")
pcor.test(total_10kb$dav_rho_mean,total_10kb$dav_tP,total_10kb[,c("Coding_prop","GC")],method="spearman")

cor.test(total_100kb$dav_rho_mean,total_100kb$dav_tP,method="spearman")
pcor.test(total_100kb$dav_rho_mean,total_100kb$dav_tP,total_100kb[,c("Coding_prop","GC")],method="spearman")


#1.1 rho vs. trs-tP
cor.test(total_10kb$trs_rho_mean,total_10kb$trs_tP,method="spearman")
pcor.test(total_10kb$trs_rho_mean,total_10kb$trs_tP,total_10kb[,c("Coding_prop","GC")],method="spearman")

cor.test(total_100kb$trs_rho_mean,total_100kb$trs_tP,method="spearman")
pcor.test(total_100kb$trs_rho_mean,total_100kb$trs_tP,total_100kb[,c("Coding_prop","GC")],method="spearman")

#1.1 rho vs. tri-tP
cor.test(total_10kb$tri_rho_mean,total_10kb$tri_tP,method="spearman")
pcor.test(total_10kb$tri_rho_mean,total_10kb$tri_tP,total_10kb[,c("Coding_prop","GC")],method="spearman")

cor.test(total_100kb$tri_rho_mean,total_100kb$tri_tP,method="spearman")
pcor.test(total_100kb$tri_rho_mean,total_100kb$tri_tP,total_100kb[,c("Coding_prop","GC")],method="spearman")


##########################################################################
##1.2 rho vs. fst
#tra-dav
cor.test(total_10kb$tra_rho_mean,total_10kb$tra_dav_fst,method="spearman")
pcor.test(total_10kb$tra_rho_mean,total_10kb$tra_dav_fst,total_10kb[,c("Coding_prop","GC")],method="spearman")

cor.test(total_100kb$tra_rho_mean,total_100kb$tra_dav_fst,method="spearman")
pcor.test(total_100kb$tra_rho_mean,total_100kb$tra_dav_fst,total_100kb[,c("Coding_prop","GC")],method="spearman")

#tra-trs
cor.test(total_10kb$tra_rho_mean,total_10kb$tra_trs_fst,method="spearman")
pcor.test(total_10kb$tra_rho_mean,total_10kb$tra_trs_fst,total_10kb[,c("Coding_prop","GC")],method="spearman")

cor.test(total_100kb$tra_rho_mean,total_100kb$tra_trs_fst,method="spearman")
pcor.test(total_100kb$tra_rho_mean,total_100kb$tra_trs_fst,total_100kb[,c("Coding_prop","GC")],method="spearman")

#dav-trs
cor.test(total_10kb$tra_rho_mean,total_10kb$dav_trs_fst,method="spearman")
pcor.test(total_10kb$tra_rho_mean,total_10kb$dav_trs_fst,total_10kb[,c("Coding_prop","GC")],method="spearman")

cor.test(total_100kb$tra_rho_mean,total_100kb$dav_trs_fst,method="spearman")
pcor.test(total_100kb$tra_rho_mean,total_100kb$dav_trs_fst,total_100kb[,c("Coding_prop","GC")],method="spearman")

#tra-tri
cor.test(total_10kb$tra_rho_mean,total_10kb$tra_tri_fst,method="spearman")
pcor.test(total_10kb$tra_rho_mean,total_10kb$tra_tri_fst,total_10kb[,c("Coding_prop","GC")],method="spearman")

cor.test(total_100kb$tra_rho_mean,total_100kb$tra_tri_fst,method="spearman")
pcor.test(total_100kb$tra_rho_mean,total_100kb$tra_tri_fst,total_100kb[,c("Coding_prop","GC")],method="spearman")

#dav-tri
cor.test(total_10kb$tra_rho_mean,total_10kb$dav_tri_fst,method="spearman")
pcor.test(total_10kb$tra_rho_mean,total_10kb$dav_tri_fst,total_10kb[,c("Coding_prop","GC")],method="spearman")

cor.test(total_100kb$tra_rho_mean,total_100kb$dav_tri_fst,method="spearman")
pcor.test(total_100kb$tra_rho_mean,total_100kb$dav_tri_fst,total_100kb[,c("Coding_prop","GC")],method="spearman")

#trs-tri
cor.test(total_10kb$tra_rho_mean,total_10kb$trs_tri_fst,method="spearman")
pcor.test(total_10kb$tra_rho_mean,total_10kb$trs_tri_fst,total_10kb[,c("Coding_prop","GC")],method="spearman")

cor.test(total_100kb$tra_rho_mean,total_100kb$trs_tri_fst,method="spearman")
pcor.test(total_100kb$tra_rho_mean,total_100kb$trs_tri_fst,total_100kb[,c("Coding_prop","GC")],method="spearman")


##########################################################################
##1.3 rho vs. dxy
#tra-dav
cor.test(total_10kb$tra_rho_mean,total_10kb$tra_dav_dxy,method="spearman")
pcor.test(total_10kb$tra_rho_mean,total_10kb$tra_dav_dxy,total_10kb[,c("Coding_prop","GC")],method="spearman")

cor.test(total_100kb$tra_rho_mean,total_100kb$tra_dav_dxy,method="spearman")
pcor.test(total_100kb$tra_rho_mean,total_100kb$tra_dav_dxy,total_100kb[,c("Coding_prop","GC")],method="spearman")

#tra-trs
cor.test(total_10kb$tra_rho_mean,total_10kb$tra_trs_dxy,method="spearman")
pcor.test(total_10kb$tra_rho_mean,total_10kb$tra_trs_dxy,total_10kb[,c("Coding_prop","GC")],method="spearman")

cor.test(total_100kb$tra_rho_mean,total_100kb$tra_trs_dxy,method="spearman")
pcor.test(total_100kb$tra_rho_mean,total_100kb$tra_trs_dxy,total_100kb[,c("Coding_prop","GC")],method="spearman")

#dav-trs
cor.test(total_10kb$tra_rho_mean,total_10kb$dav_trs_dxy,method="spearman")
pcor.test(total_10kb$tra_rho_mean,total_10kb$dav_trs_dxy,total_10kb[,c("Coding_prop","GC")],method="spearman")

cor.test(total_100kb$tra_rho_mean,total_100kb$dav_trs_dxy,method="spearman")
pcor.test(total_100kb$tra_rho_mean,total_100kb$dav_trs_dxy,total_100kb[,c("Coding_prop","GC")],method="spearman")

#tra-tri
cor.test(total_10kb$tra_rho_mean,total_10kb$tra_tri_dxy,method="spearman")
pcor.test(total_10kb$tra_rho_mean,total_10kb$tra_tri_dxy,total_10kb[,c("Coding_prop","GC")],method="spearman")

cor.test(total_100kb$tra_rho_mean,total_100kb$tra_tri_dxy,method="spearman")
pcor.test(total_100kb$tra_rho_mean,total_100kb$tra_tri_dxy,total_100kb[,c("Coding_prop","GC")],method="spearman")

#dav-tri
cor.test(total_10kb$tra_rho_mean,total_10kb$dav_tri_dxy,method="spearman")
pcor.test(total_10kb$tra_rho_mean,total_10kb$dav_tri_dxy,total_10kb[,c("Coding_prop","GC")],method="spearman")

cor.test(total_100kb$tra_rho_mean,total_100kb$dav_tri_dxy,method="spearman")
pcor.test(total_100kb$tra_rho_mean,total_100kb$dav_tri_dxy,total_100kb[,c("Coding_prop","GC")],method="spearman")

#trs-tri
cor.test(total_10kb$tra_rho_mean,total_10kb$trs_tri_dxy,method="spearman")
pcor.test(total_10kb$tra_rho_mean,total_10kb$trs_tri_dxy,total_10kb[,c("Coding_prop","GC")],method="spearman")

cor.test(total_100kb$tra_rho_mean,total_100kb$trs_tri_dxy,method="spearman")
pcor.test(total_100kb$tra_rho_mean,total_100kb$trs_tri_dxy,total_100kb[,c("Coding_prop","GC")],method="spearman")

##########################################################################
##1.4 rho vs. nls

cor.test(total_10kb$tra_rho_mean,total_10kb$nls,method="spearman")
pcor.test(total_10kb$tra_rho_mean,total_10kb$nls,total_10kb[,c("Coding_prop","GC")],method="spearman")

cor.test(total_100kb$tra_rho_mean,total_100kb$nls,method="spearman")
pcor.test(total_100kb$tra_rho_mean,total_100kb$nls,total_100kb[,c("Coding_prop","GC")],method="spearman")

##########################################################################
##1.5 rho vs. top2

cor.test(total_10kb$tra_rho_mean,total_10kb$top2,method="spearman")
pcor.test(total_10kb$tra_rho_mean,total_10kb$top2,total_10kb[,c("Coding_prop","GC")],method="spearman")

cor.test(total_100kb$tra_rho_mean,total_100kb$top2,method="spearman")
pcor.test(total_100kb$tra_rho_mean,total_100kb$top2,total_100kb[,c("Coding_prop","GC")],method="spearman")

##########################################################################
##1.6 rho vs. fd
total_10kb_fd=total_10kb[which(total_10kb$fd>0),]
total_100kb_fd=total_100kb[which(total_100kb$fd>0),]

cor.test(total_10kb_fd$tra_rho_mean,total_10kb_fd$fd,method="spearman")
pcor.test(total_10kb_fd$tra_rho_mean,total_10kb_fd$fd,total_10kb_fd[,c("Coding_prop","GC")],method="spearman")

cor.test(total_100kb_fd$tra_rho_mean,total_100kb_fd$fd,method="spearman")
pcor.test(total_100kb_fd$tra_rho_mean,total_100kb_fd$fd,total_100kb_fd[,c("Coding_prop","GC")],method="spearman")




