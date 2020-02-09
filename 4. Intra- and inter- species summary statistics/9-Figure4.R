library(data.table)
library(RColorBrewer)
library(dplyr)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(tidyverse)

colors <- brewer.pal(12,"Paired")

setwd("~/Dropbox/davidiana_paper/data/linked_selection")

##########################################################################
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


##########################################################################
#####The purpose of this script is to make a matrix plot, with relationship between
#1. Coding prop with pi, fst, dxy, ILS, top2, fd
#2. Recombination rate with pi, fst, dxy, ILS, top2, fd
#3. CLR from SF2 with pi, fst, dxy, ILS, top2, fd
#4. Selection effects of genes with pi, fst, dxy, ILS, top2, fd

###notes: here recombination rate, CLR, Selection effects, I only use the results from P.tra as representative
###notes: for pi, fst, dxy, I use Ptra, Ptra-Pdav as representative 


colors <- brewer.pal(12,"Paired")
cols=colors[c(2,4,8,10)]

setwd("./plot/")

####################################################
###1. Coding prop

###Using coding density
summary_10kb_w_cd=summary_10kb_w %>% mutate(cd_quantile=ntile(Coding_prop,10))

nls_cd=data.frame(cbind(summary_10kb_w_cd$cd_quantile,summary_10kb_w_cd$Coding_prop,nls_10kb_new_w$nls_n,twisst_10kb_prop_sites_w$top1,twisst_10kb_prop_sites_w$top2,twisst_10kb_prop_sites_w$top3))
names(nls_cd)=c("cd_quantile","coding_prop","nls","top1","top2","top3")

nls_cd_mean=as.vector(by(nls_cd$nls, nls_cd$cd_quantile, function(x) mean(x,na.rm=T)))
nls_cd_sd=as.vector(by(nls_cd$nls, nls_cd$cd_quantile, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

top2_cd_mean=as.vector(by(nls_cd$top2, nls_cd$cd_quantile, function(x) mean(x,na.rm=T)))
top2_cd_sd=as.vector(by(nls_cd$top2, nls_cd$cd_quantile, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

non_top2_cd_mean=as.vector(by(1-nls_cd$top2, nls_cd$cd_quantile, function(x) mean(x,na.rm=T)))
non_top2_cd_sd=as.vector(by(1-nls_cd$top2, nls_cd$cd_quantile, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

cd_mean=as.vector(by(nls_cd$coding_prop, nls_cd$cd_quantile, function(x) mean(x,na.rm=T)))
cd_sd=as.vector(by(nls_cd$coding_prop, nls_cd$cd_quantile, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

nls_cd.df=data.frame(cd_quantile=1:10,cd_mean=cd_mean,cd_sd=cd_sd,nls_cd_mean=nls_cd_mean,nls_cd_sd=nls_cd_sd)
top2_cd.df=data.frame(cd_quantile=1:10,cd_mean=cd_mean,cd_sd=cd_sd,top2_cd_mean=top2_cd_mean,top2_cd_sd=top2_cd_sd)
non_top2_cd.df=data.frame(cd_quantile=1:10,cd_mean=cd_mean,cd_sd=cd_sd,non_top2_cd_mean=non_top2_cd_mean,non_top2_cd_sd=non_top2_cd_sd)


#####1.1 cd vs. pi
###pi of P.tra
pi_tra_cd=data.frame(cbind(summary_10kb_w_cd$cd_quantile,summary_10kb_w_cd$Coding_prop,summary_10kb_w_cd$tra_tP))
names(pi_tra_cd)=c("cd_quantile","coding_prop","tra_tP")

###pi
pi_tra_cd_mean=as.vector(by(pi_tra_cd$tra_tP, pi_tra_cd$cd_quantile, function(x) mean(x,na.rm=T)))
pi_tra_cd_sd=as.vector(by(pi_tra_cd$tra_tP, pi_tra_cd$cd_quantile, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

tra_cd_mean=as.vector(by(pi_tra_cd$coding_prop, pi_tra_cd$cd_quantile, function(x) mean(x,na.rm=T)))
tra_cd_sd=as.vector(by(pi_tra_cd$coding_prop, pi_tra_cd$cd_quantile, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

pi_tra_cd.df=data.frame(cd_quantile=1:10,tra_cd_mean=tra_cd_mean,tra_cd_sd=tra_cd_sd,pi_tra_cd_mean=pi_tra_cd_mean,pi_tra_cd_sd=pi_tra_cd_sd)

#pi
ggplot(pi_tra_cd.df,aes(as.factor(cd_quantile),y=pi_tra_cd_mean,ymin=pi_tra_cd_mean-pi_tra_cd_sd,ymax=pi_tra_cd_mean+pi_tra_cd_sd))+
  geom_pointrange(color=cols[1],size=0.5)+
  scale_x_discrete(labels= c("0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%"))+
  xlab("")+ylab("")+ylim(c(0.009,0.019))+
 # xlab("Percentile of coding density")+ylab("Scaled rate of ILS sites")+
  theme_bw()+theme(plot.margin = margin(1, 1, 0.5, 0.5, "cm"),axis.text=element_text(size=18),axis.text.x = element_blank())
ggsave("cd.tra_pi.pdf",width=4,height=3)


#####1.2 cd vs.fst
###fst of P.tra-Pdav
fst_tra_dav_cd=data.frame(cbind(summary_10kb_w_cd$cd_quantile,summary_10kb_w_cd$Coding_prop,summary_10kb_w_cd$tra_dav_fst))
names(fst_tra_dav_cd)=c("cd_quantile","coding_prop","tra_dav_fst")


fst_tra_dav_cd_mean=as.vector(by(fst_tra_dav_cd$tra_dav_fst, fst_tra_dav_cd$cd_quantile, function(x) mean(x,na.rm=T)))
fst_tra_dav_cd_sd=as.vector(by(fst_tra_dav_cd$tra_dav_fst, fst_tra_dav_cd$cd_quantile, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

tra_cd_mean=as.vector(by(fst_tra_dav_cd$coding_prop, fst_tra_dav_cd$cd_quantile, function(x) mean(x,na.rm=T)))
tra_cd_sd=as.vector(by(fst_tra_dav_cd$coding_prop, fst_tra_dav_cd$cd_quantile, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

fst_tra_dav_cd.df=data.frame(cd_quantile=1:10,tra_cd_mean=tra_cd_mean,tra_cd_sd=tra_cd_sd,fst_tra_dav_cd_mean=fst_tra_dav_cd_mean,fst_tra_dav_cd_sd=fst_tra_dav_cd_sd)

#fst
ggplot(fst_tra_dav_cd.df,aes(as.factor(cd_quantile),y=fst_tra_dav_cd_mean,ymin=fst_tra_dav_cd_mean-fst_tra_dav_cd_sd,ymax=fst_tra_dav_cd_mean+fst_tra_dav_cd_sd))+
  geom_pointrange(color=cols[1],size=0.5)+
  scale_x_discrete(labels= c("0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%"))+
  xlab("")+ylab("")+ylim(c(0.27,0.41))+
  # xlab("Percentile of coding density")+ylab("Scaled rate of ILS sites")+
  theme_bw()+theme(plot.margin = margin(1, 1, 0.5, 0.5, "cm"),axis.text=element_text(size=18),axis.text.x = element_blank())
ggsave("cd.tra_dav.fst.pdf",width=4,height=3)


#####1.3 cd vs.dxy
###dxy of P.tra-Pdav
dxy_tra_dav_cd=data.frame(cbind(summary_10kb_w_cd$cd_quantile,summary_10kb_w_cd$Coding_prop,summary_10kb_w_cd$tra_dav_dxy))
names(dxy_tra_dav_cd)=c("cd_quantile","coding_prop","tra_dav_dxy")


dxy_tra_dav_cd_mean=as.vector(by(dxy_tra_dav_cd$tra_dav_dxy, dxy_tra_dav_cd$cd_quantile, function(x) mean(x,na.rm=T)))
dxy_tra_dav_cd_sd=as.vector(by(dxy_tra_dav_cd$tra_dav_dxy, dxy_tra_dav_cd$cd_quantile, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

tra_cd_mean=as.vector(by(dxy_tra_dav_cd$coding_prop, dxy_tra_dav_cd$cd_quantile, function(x) mean(x,na.rm=T)))
tra_cd_sd=as.vector(by(dxy_tra_dav_cd$coding_prop, dxy_tra_dav_cd$cd_quantile, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

dxy_tra_dav_cd.df=data.frame(cd_quantile=1:10,tra_cd_mean=tra_cd_mean,tra_cd_sd=tra_cd_sd,dxy_tra_dav_cd_mean=dxy_tra_dav_cd_mean,dxy_tra_dav_cd_sd=dxy_tra_dav_cd_sd)

#dxy
ggplot(dxy_tra_dav_cd.df,aes(as.factor(cd_quantile),y=dxy_tra_dav_cd_mean,ymin=dxy_tra_dav_cd_mean-dxy_tra_dav_cd_sd,ymax=dxy_tra_dav_cd_mean+dxy_tra_dav_cd_sd))+
  geom_pointrange(color=cols[1],size=0.5)+
  scale_x_discrete(labels= c("0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%"))+
  xlab("")+ylab("")+ylim(c(0.0185,0.0215))+
  # xlab("Percentile of coding density")+ylab("Scaled rate of ILS sites")+
  theme_bw()+theme(plot.margin = margin(1, 1, 0.5, 0.5, "cm"),axis.text.x=element_blank(),axis.text=element_text(size=18))
ggsave("cd.tra_dav.dxy.pdf",width=4,height=3)

#####1.4 cd vs.nls
#nls
ggplot(nls_cd.df,aes(as.factor(cd_quantile),y=nls_cd_mean,ymin=nls_cd_mean-nls_cd_sd,ymax=nls_cd_mean+nls_cd_sd))+
  geom_pointrange(color=cols[1],size=0.5)+
  scale_x_discrete(labels= c("0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%"))+
  xlab("")+ylab("")+ylim(c(0.095,0.128))+
  #xlab("Percentile of coding density")+ylab("Scaled rate of ILS sites")+
  theme_bw()+theme(plot.margin = margin(1, 1, 0.5, 0.5, "cm"),axis.text=element_text(size=18),axis.text.x = element_blank())
ggsave("cd.ILS.pdf",width=4,height=3)

#####1.5 cd vs.top2
#the percentile of topology2 
ggplot(top2_cd.df,aes(as.factor(cd_quantile),y=top2_cd_mean,ymin=top2_cd_mean-top2_cd_sd,ymax=top2_cd_mean+top2_cd_sd))+
  geom_pointrange(color=cols[1],size=0.5)+
  scale_x_discrete(labels= c("0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%"))+
  xlab("")+ylab("")+ylim(c(0.48,0.63))+
  #xlab("Percentile of coding density")+ylab("Species topology (%)")+
  theme_bw()+theme(plot.margin = margin(1, 1, 0.5, 0.5, "cm"),axis.text=element_text(size=18),axis.text.x = element_blank())
ggsave("cd.top2.pdf",width=4,height=3)

#####1.6 cd vs.fd
fd_cd_old=data.frame(cbind(summary_10kb_w_cd$cd_quantile,summary_10kb_w_cd$Coding_prop,D_gatk_10kb_w$fd))
names(fd_cd_old)=c("cd_quantile","coding_prop","fd")
fd_cd=fd_cd_old[which(fd_cd_old$fd>=0),]
cor.test(fd_cd$coding_prop,fd_cd$fd,method="spearman")

fd_cd_mean=as.vector(by(fd_cd$fd, fd_cd$cd_quantile, function(x) mean(x,na.rm=T)))
fd_cd_sd=as.vector(by(fd_cd$fd, fd_cd$cd_quantile, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

cd_mean=as.vector(by(fd_cd$coding_prop, fd_cd$cd_quantile, function(x) mean(x,na.rm=T)))
cd_sd=as.vector(by(fd_cd$coding_prop, fd_cd$cd_quantile, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

fd_cd.df=data.frame(cd_quantile=1:10,cd_mean=cd_mean,cd_sd=cd_sd,fd_cd_mean=fd_cd_mean,fd_cd_sd=fd_cd_sd)

#fd
ggplot(fd_cd.df,aes(as.factor(cd_quantile),y=fd_cd_mean,ymin=fd_cd_mean-fd_cd_sd,ymax=fd_cd_mean+fd_cd_sd))+
  geom_pointrange(color=cols[1],size=0.5)+
  scale_x_discrete(labels= c("0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%"))+
  xlab("")+ylab("")+ylim(c(0.12,0.19))+
  #xlab("Percentile of coding density")+ylab(expression(f[d]))+
  theme_bw()+theme(plot.margin = margin(1, 1, 0.5, 0.5, "cm"),axis.text=element_text(size=18),axis.text.x = element_blank())
ggsave("cd.fd.pdf",width=4,height=3)




####################################################
###2 recombination rate from P.tra

###Using coding density
summary_10kb_w_rho=summary_10kb_w %>% mutate(rho_quantile=ntile(tra_rho_mean,10))

nls_rho=data.frame(cbind(summary_10kb_w_rho$rho_quantile,summary_10kb_w_rho$tra_rho_mean,nls_10kb_new_w$nls_n,twisst_10kb_prop_sites_w$top1,twisst_10kb_prop_sites_w$top2,twisst_10kb_prop_sites_w$top3))
names(nls_rho)=c("rho_quantile","tra_rho","nls","top1","top2","top3")

nls_rho_mean=as.vector(by(nls_rho$nls, nls_rho$rho_quantile, function(x) mean(x,na.rm=T)))
nls_rho_sd=as.vector(by(nls_rho$nls, nls_rho$rho_quantile, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

top2_rho_mean=as.vector(by(nls_rho$top2, nls_rho$rho_quantile, function(x) mean(x,na.rm=T)))
top2_rho_sd=as.vector(by(nls_rho$top2, nls_rho$rho_quantile, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

non_top2_rho_mean=as.vector(by(1-nls_rho$top2, nls_rho$rho_quantile, function(x) mean(x,na.rm=T)))
non_top2_rho_sd=as.vector(by(1-nls_rho$top2, nls_rho$rho_quantile, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

rho_mean=as.vector(by(nls_rho$tra_rho, nls_rho$rho_quantile, function(x) mean(x,na.rm=T)))
rho_sd=as.vector(by(nls_rho$tra_rho, nls_rho$rho_quantile, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

nls_rho.df=data.frame(rho_quantile=1:10,rho_mean=rho_mean,rho_sd=rho_sd,nls_rho_mean=nls_rho_mean,nls_rho_sd=nls_rho_sd)
top2_rho.df=data.frame(rho_quantile=1:10,rho_mean=rho_mean,rho_sd=rho_sd,top2_rho_mean=top2_rho_mean,top2_rho_sd=top2_rho_sd)
non_top2_rho.df=data.frame(rho_quantile=1:10,rho_mean=rho_mean,rho_sd=rho_sd,non_top2_rho_mean=non_top2_rho_mean,non_top2_rho_sd=non_top2_rho_sd)


#####1.1 rho vs. pi
###pi of P.tra
pi_tra_rho=data.frame(cbind(summary_10kb_w_rho$rho_quantile,summary_10kb_w_rho$Coding_prop,summary_10kb_w_rho$tra_tP))
names(pi_tra_rho)=c("rho_quantile","tra_rho","tra_tP")

###pi
pi_tra_rho_mean=as.vector(by(pi_tra_rho$tra_tP, pi_tra_rho$rho_quantile, function(x) mean(x,na.rm=T)))
pi_tra_rho_sd=as.vector(by(pi_tra_rho$tra_tP, pi_tra_rho$rho_quantile, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

tra_rho_mean=as.vector(by(pi_tra_rho$tra_rho, pi_tra_rho$rho_quantile, function(x) mean(x,na.rm=T)))
tra_rho_sd=as.vector(by(pi_tra_rho$tra_rho, pi_tra_rho$rho_quantile, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

pi_tra_rho.df=data.frame(rho_quantile=1:10,tra_rho_mean=tra_rho_mean,tra_rho_sd=tra_rho_sd,pi_tra_rho_mean=pi_tra_rho_mean,pi_tra_rho_sd=pi_tra_rho_sd)

#pi
ggplot(pi_tra_rho.df,aes(as.factor(rho_quantile),y=pi_tra_rho_mean,ymin=pi_tra_rho_mean-pi_tra_rho_sd,ymax=pi_tra_rho_mean+pi_tra_rho_sd))+
  geom_pointrange(color=cols[2],size=0.5)+
  scale_x_discrete(labels= c("0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%"))+
  xlab("")+ylab("")+ylim(c(0.009,0.019))+
  # xlab("Percentile of coding density")+ylab("Scaled rate of ILS sites")+
  theme_bw()+theme(plot.margin = margin(1, 1, 0.5, 0.5, "cm"),axis.text=element_text(size=18),axis.text.x = element_blank())
ggsave("rho.tra_pi.pdf",width=4,height=3)


#####1.2 rho vs.fst
###fst of P.tra-Pdav
fst_tra_dav_rho=data.frame(cbind(summary_10kb_w_rho$rho_quantile,summary_10kb_w_rho$tra_rho_mean,summary_10kb_w_rho$tra_dav_fst))
names(fst_tra_dav_rho)=c("rho_quantile","tra_rho","tra_dav_fst")


fst_tra_dav_rho_mean=as.vector(by(fst_tra_dav_rho$tra_dav_fst, fst_tra_dav_rho$rho_quantile, function(x) mean(x,na.rm=T)))
fst_tra_dav_rho_sd=as.vector(by(fst_tra_dav_rho$tra_dav_fst, fst_tra_dav_rho$rho_quantile, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

tra_rho_mean=as.vector(by(fst_tra_dav_rho$tra_rho, fst_tra_dav_rho$rho_quantile, function(x) mean(x,na.rm=T)))
tra_rho_sd=as.vector(by(fst_tra_dav_rho$tra_rho, fst_tra_dav_rho$rho_quantile, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

fst_tra_dav_rho.df=data.frame(rho_quantile=1:10,tra_rho_mean=tra_rho_mean,tra_rho_sd=tra_rho_sd,fst_tra_dav_rho_mean=fst_tra_dav_rho_mean,fst_tra_dav_rho_sd=fst_tra_dav_rho_sd)

#fst
ggplot(fst_tra_dav_rho.df,aes(as.factor(rho_quantile),y=fst_tra_dav_rho_mean,ymin=fst_tra_dav_rho_mean-fst_tra_dav_rho_sd,ymax=fst_tra_dav_rho_mean+fst_tra_dav_rho_sd))+
  geom_pointrange(color=cols[2],size=0.5)+
  scale_x_discrete(labels= c("0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%"))+
  xlab("")+ylab("")+ylim(c(0.27,0.41))+
  # xlab("Percentile of coding density")+ylab("Scaled rate of ILS sites")+
  theme_bw()+theme(plot.margin = margin(1, 1, 0.5, 0.5, "cm"),axis.text=element_text(size=18),axis.text.x = element_blank())
ggsave("rho.tra_dav.fst.pdf",width=4,height=3)


#####1.3 rho vs.dxy
###dxy of P.tra-Pdav
dxy_tra_dav_rho=data.frame(cbind(summary_10kb_w_rho$rho_quantile,summary_10kb_w_rho$tra_rho_mean,summary_10kb_w_rho$tra_dav_dxy))
names(dxy_tra_dav_rho)=c("rho_quantile","tra_rho","tra_dav_dxy")


dxy_tra_dav_rho_mean=as.vector(by(dxy_tra_dav_rho$tra_dav_dxy, dxy_tra_dav_rho$rho_quantile, function(x) mean(x,na.rm=T)))
dxy_tra_dav_rho_sd=as.vector(by(dxy_tra_dav_rho$tra_dav_dxy, dxy_tra_dav_rho$rho_quantile, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

tra_rho_mean=as.vector(by(dxy_tra_dav_rho$tra_rho, dxy_tra_dav_rho$rho_quantile, function(x) mean(x,na.rm=T)))
tra_rho_sd=as.vector(by(dxy_tra_dav_rho$tra_rho, dxy_tra_dav_rho$rho_quantile, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

dxy_tra_dav_rho.df=data.frame(rho_quantile=1:10,tra_rho_mean=tra_rho_mean,tra_rho_sd=tra_rho_sd,dxy_tra_dav_rho_mean=dxy_tra_dav_rho_mean,dxy_tra_dav_rho_sd=dxy_tra_dav_rho_sd)

#dxy
ggplot(dxy_tra_dav_rho.df,aes(as.factor(rho_quantile),y=dxy_tra_dav_rho_mean,ymin=dxy_tra_dav_rho_mean-dxy_tra_dav_rho_sd,ymax=dxy_tra_dav_rho_mean+dxy_tra_dav_rho_sd))+
  geom_pointrange(color=cols[2],size=0.5)+
  scale_x_discrete(labels= c("0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%"))+
  xlab("")+ylab("")+ylim(c(0.0185,0.0215))+
  # xlab("Percentile of coding density")+ylab("Scaled rate of ILS sites")+
  theme_bw()+theme(plot.margin = margin(1, 1, 0.5, 0.5, "cm"),axis.text.x=element_blank(),axis.text=element_text(size=18))
ggsave("rho.tra_dav.dxy.pdf",width=4,height=3)

#####1.4 rho vs.nls
#nls
ggplot(nls_rho.df,aes(as.factor(rho_quantile),y=nls_rho_mean,ymin=nls_rho_mean-nls_rho_sd,ymax=nls_rho_mean+nls_rho_sd))+
  geom_pointrange(color=cols[2],size=0.5)+
  scale_x_discrete(labels= c("0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%"))+
  xlab("")+ylab("")+ylim(c(0.095,0.128))+
  #xlab("Percentile of coding density")+ylab("Scaled rate of ILS sites")+
  theme_bw()+theme(plot.margin = margin(1, 1, 0.5, 0.5, "cm"),axis.text=element_text(size=18),axis.text.x = element_blank())
ggsave("rho.ILS.pdf",width=4,height=3)

#####1.5 rho vs.top2
#the percentile of topology2 
ggplot(top2_rho.df,aes(as.factor(rho_quantile),y=top2_rho_mean,ymin=top2_rho_mean-top2_rho_sd,ymax=top2_rho_mean+top2_rho_sd))+
  geom_pointrange(color=cols[2],size=0.5)+
  scale_x_discrete(labels= c("0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%"))+
  xlab("")+ylab("")+ylim(c(0.48,0.63))+
  #xlab("Percentile of coding density")+ylab("Species topology (%)")+
  theme_bw()+theme(plot.margin = margin(1, 1, 0.5, 0.5, "cm"),axis.text=element_text(size=18),axis.text.x = element_blank())
ggsave("rho.top2.pdf",width=4,height=3)

#####1.6 rho vs.fd
fd_rho_old=data.frame(cbind(summary_10kb_w_rho$rho_quantile,summary_10kb_w_rho$tra_rho_mean,D_gatk_10kb_w$fd))
names(fd_rho_old)=c("rho_quantile","tra_rho","fd")
fd_rho=fd_rho_old[which(fd_rho_old$fd>=0),]
cor.test(fd_rho$tra_rho,fd_rho$fd,method="spearman")

fd_rho_mean=as.vector(by(fd_rho$fd, fd_rho$rho_quantile, function(x) mean(x,na.rm=T)))
fd_rho_sd=as.vector(by(fd_rho$fd, fd_rho$rho_quantile, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

rho_mean=as.vector(by(fd_rho$tra_rho, fd_rho$rho_quantile, function(x) mean(x,na.rm=T)))
rho_sd=as.vector(by(fd_rho$tra_rho, fd_rho$rho_quantile, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

fd_rho.df=data.frame(rho_quantile=1:10,rho_mean=rho_mean,rho_sd=rho_sd,fd_rho_mean=fd_rho_mean,fd_rho_sd=fd_rho_sd)

#fd
ggplot(fd_rho.df,aes(as.factor(rho_quantile),y=fd_rho_mean,ymin=fd_rho_mean-fd_rho_sd,ymax=fd_rho_mean+fd_rho_sd))+
  geom_pointrange(color=cols[2],size=0.5)+
  scale_x_discrete(labels= c("0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%"))+
  xlab("")+ylab("")+ylim(c(0.12,0.19))+
  #xlab("Percentile of coding density")+ylab(expression(f[d]))+
  theme_bw()+theme(plot.margin = margin(1, 1, 0.5, 0.5, "cm"),axis.text=element_text(size=18),axis.text.x = element_blank())
ggsave("rho.fd.pdf",width=4,height=3)


####################################################
###3 CLR of SweepFinder2 from P.tra

###Using CLR of SweepFinder2
summary_10kb_w_clr=summary_10kb_w %>% mutate(clr_quantile=ntile(sf2_tra,10))

nls_clr=data.frame(cbind(summary_10kb_w_clr$clr_quantile,summary_10kb_w_clr$sf2_tra,nls_10kb_new_w$nls_n,twisst_10kb_prop_sites_w$top1,twisst_10kb_prop_sites_w$top2,twisst_10kb_prop_sites_w$top3))
names(nls_clr)=c("clr_quantile","tra_clr","nls","top1","top2","top3")

nls_clr_mean=as.vector(by(nls_clr$nls, nls_clr$clr_quantile, function(x) mean(x,na.rm=T)))
nls_clr_sd=as.vector(by(nls_clr$nls, nls_clr$clr_quantile, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

top2_clr_mean=as.vector(by(nls_clr$top2, nls_clr$clr_quantile, function(x) mean(x,na.rm=T)))
top2_clr_sd=as.vector(by(nls_clr$top2, nls_clr$clr_quantile, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

non_top2_clr_mean=as.vector(by(1-nls_clr$top2, nls_clr$clr_quantile, function(x) mean(x,na.rm=T)))
non_top2_clr_sd=as.vector(by(1-nls_clr$top2, nls_clr$clr_quantile, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

clr_mean=as.vector(by(nls_clr$tra_clr, nls_clr$clr_quantile, function(x) mean(x,na.rm=T)))
clr_sd=as.vector(by(nls_clr$tra_clr, nls_clr$clr_quantile, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

nls_clr.df=data.frame(clr_quantile=1:10,clr_mean=clr_mean,clr_sd=clr_sd,nls_clr_mean=nls_clr_mean,nls_clr_sd=nls_clr_sd)
top2_clr.df=data.frame(clr_quantile=1:10,clr_mean=clr_mean,clr_sd=clr_sd,top2_clr_mean=top2_clr_mean,top2_clr_sd=top2_clr_sd)
non_top2_clr.df=data.frame(clr_quantile=1:10,clr_mean=clr_mean,clr_sd=clr_sd,non_top2_clr_mean=non_top2_clr_mean,non_top2_clr_sd=non_top2_clr_sd)


#####3.1 clr vs. pi
###pi of P.tra
pi_tra_clr=data.frame(cbind(summary_10kb_w_clr$clr_quantile,summary_10kb_w_clr$sf2_tra,summary_10kb_w_clr$tra_tP))
names(pi_tra_clr)=c("clr_quantile","tra_clr","tra_tP")

###pi
pi_tra_clr_mean=as.vector(by(pi_tra_clr$tra_tP, pi_tra_clr$clr_quantile, function(x) mean(x,na.rm=T)))
pi_tra_clr_sd=as.vector(by(pi_tra_clr$tra_tP, pi_tra_clr$clr_quantile, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

tra_clr_mean=as.vector(by(pi_tra_clr$tra_clr, pi_tra_clr$clr_quantile, function(x) mean(x,na.rm=T)))
tra_clr_sd=as.vector(by(pi_tra_clr$tra_clr, pi_tra_clr$clr_quantile, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

pi_tra_clr.df=data.frame(clr_quantile=1:10,tra_clr_mean=tra_clr_mean,tra_clr_sd=tra_clr_sd,pi_tra_clr_mean=pi_tra_clr_mean,pi_tra_clr_sd=pi_tra_clr_sd)

#pi
ggplot(pi_tra_clr.df,aes(as.factor(clr_quantile),y=pi_tra_clr_mean,ymin=pi_tra_clr_mean-pi_tra_clr_sd,ymax=pi_tra_clr_mean+pi_tra_clr_sd))+
  geom_pointrange(color=cols[3],size=0.5)+
  scale_x_discrete(labels= c("0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%"))+
  xlab("")+ylab("")+ylim(c(0.009,0.019))+
  # xlab("Percentile of coding density")+ylab("Scaled rate of ILS sites")+
  theme_bw()+theme(plot.margin = margin(1, 1, 0.5, 0.5, "cm"),axis.text=element_text(size=18),axis.text.x = element_blank())
ggsave("clr.tra_pi.pdf",width=4,height=3)


#####3.2 clr vs.fst
###fst of P.tra-Pdav
fst_tra_dav_clr=data.frame(cbind(summary_10kb_w_clr$clr_quantile,summary_10kb_w_clr$sf2_tra,summary_10kb_w_clr$tra_dav_fst))
names(fst_tra_dav_clr)=c("clr_quantile","tra_clr","tra_dav_fst")


fst_tra_dav_clr_mean=as.vector(by(fst_tra_dav_clr$tra_dav_fst, fst_tra_dav_clr$clr_quantile, function(x) mean(x,na.rm=T)))
fst_tra_dav_clr_sd=as.vector(by(fst_tra_dav_clr$tra_dav_fst, fst_tra_dav_clr$clr_quantile, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

tra_clr_mean=as.vector(by(fst_tra_dav_clr$tra_clr, fst_tra_dav_clr$clr_quantile, function(x) mean(x,na.rm=T)))
tra_clr_sd=as.vector(by(fst_tra_dav_clr$tra_clr, fst_tra_dav_clr$clr_quantile, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

fst_tra_dav_clr.df=data.frame(clr_quantile=1:10,tra_clr_mean=tra_clr_mean,tra_clr_sd=tra_clr_sd,fst_tra_dav_clr_mean=fst_tra_dav_clr_mean,fst_tra_dav_clr_sd=fst_tra_dav_clr_sd)

#fst
ggplot(fst_tra_dav_clr.df,aes(as.factor(clr_quantile),y=fst_tra_dav_clr_mean,ymin=fst_tra_dav_clr_mean-fst_tra_dav_clr_sd,ymax=fst_tra_dav_clr_mean+fst_tra_dav_clr_sd))+
  geom_pointrange(color=cols[3],size=0.5)+
  scale_x_discrete(labels= c("0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%"))+
  xlab("")+ylab("")+ylim(c(0.27,0.41))+
  # xlab("Percentile of coding density")+ylab("Scaled rate of ILS sites")+
  theme_bw()+theme(plot.margin = margin(1, 1, 0.5, 0.5, "cm"),axis.text=element_text(size=18),axis.text.x = element_blank())
ggsave("clr.tra_dav.fst.pdf",width=4,height=3)


#####3.3 rho vs.dxy
###dxy of P.tra-Pdav
dxy_tra_dav_clr=data.frame(cbind(summary_10kb_w_clr$clr_quantile,summary_10kb_w_clr$sf2_tra,summary_10kb_w_clr$tra_dav_dxy))
names(dxy_tra_dav_clr)=c("clr_quantile","tra_clr","tra_dav_dxy")


dxy_tra_dav_clr_mean=as.vector(by(dxy_tra_dav_clr$tra_dav_dxy, dxy_tra_dav_clr$clr_quantile, function(x) mean(x,na.rm=T)))
dxy_tra_dav_clr_sd=as.vector(by(dxy_tra_dav_clr$tra_dav_dxy, dxy_tra_dav_clr$clr_quantile, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

tra_clr_mean=as.vector(by(dxy_tra_dav_clr$tra_clr, dxy_tra_dav_clr$clr_quantile, function(x) mean(x,na.rm=T)))
tra_clr_sd=as.vector(by(dxy_tra_dav_clr$tra_clr, dxy_tra_dav_clr$clr_quantile, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

dxy_tra_dav_clr.df=data.frame(clr_quantile=1:10,tra_clr_mean=tra_clr_mean,tra_clr_sd=tra_clr_sd,dxy_tra_dav_clr_mean=dxy_tra_dav_clr_mean,dxy_tra_dav_clr_sd=dxy_tra_dav_clr_sd)

#dxy
ggplot(dxy_tra_dav_clr.df,aes(as.factor(clr_quantile),y=dxy_tra_dav_clr_mean,ymin=dxy_tra_dav_clr_mean-dxy_tra_dav_clr_sd,ymax=dxy_tra_dav_clr_mean+dxy_tra_dav_clr_sd))+
  geom_pointrange(color=cols[3],size=0.5)+
  scale_x_discrete(labels= c("0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%"))+
  xlab("")+ylab("")+ylim(c(0.0185,0.0215))+
  # xlab("Percentile of coding density")+ylab("Scaled rate of ILS sites")+
  theme_bw()+theme(plot.margin = margin(1, 1, 0.5, 0.5, "cm"),axis.text.x=element_blank(),axis.text=element_text(size=18))
ggsave("clr.tra_dav.dxy.pdf",width=4,height=3)

#####3.4 clr vs.nls
#nls
ggplot(nls_clr.df,aes(as.factor(clr_quantile),y=nls_clr_mean,ymin=nls_clr_mean-nls_clr_sd,ymax=nls_clr_mean+nls_clr_sd))+
  geom_pointrange(color=cols[3],size=0.5)+
  scale_x_discrete(labels= c("0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%"))+
  xlab("")+ylab("")+ylim(c(0.095,0.128))+
  #xlab("Percentile of coding density")+ylab("Scaled rate of ILS sites")+
  theme_bw()+theme(plot.margin = margin(1, 1, 0.5, 0.5, "cm"),axis.text=element_text(size=18),axis.text.x = element_blank())
ggsave("clr.ILS.pdf",width=4,height=3)

#####3.5 clr vs.top2
#the percentile of topology2 
ggplot(top2_clr.df,aes(as.factor(clr_quantile),y=top2_clr_mean,ymin=top2_clr_mean-top2_clr_sd,ymax=top2_clr_mean+top2_clr_sd))+
  geom_pointrange(color=cols[3],size=0.5)+
  scale_x_discrete(labels= c("0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%"))+
  xlab("")+ylab("")+ylim(c(0.48,0.63))+
  #xlab("Percentile of coding density")+ylab("Species topology (%)")+
  theme_bw()+theme(plot.margin = margin(1, 1, 0.5, 0.5, "cm"),axis.text=element_text(size=18),axis.text.x = element_blank())
ggsave("clr.top2.pdf",width=4,height=3)

#####3.6 clr vs.fd
fd_clr_old=data.frame(cbind(summary_10kb_w_clr$clr_quantile,summary_10kb_w_clr$sf2_tra,D_gatk_10kb_w$fd))
names(fd_clr_old)=c("clr_quantile","tra_clr","fd")
fd_clr=fd_clr_old[which(fd_clr_old$fd>=0),]
cor.test(fd_clr$tra_clr,fd_clr$fd,method="spearman")

fd_clr_mean=as.vector(by(fd_clr$fd, fd_clr$clr_quantile, function(x) mean(x,na.rm=T)))
fd_clr_sd=as.vector(by(fd_clr$fd, fd_clr$clr_quantile, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

clr_mean=as.vector(by(fd_clr$tra_clr, fd_clr$clr_quantile, function(x) mean(x,na.rm=T)))
clr_sd=as.vector(by(fd_clr$tra_clr, fd_clr$clr_quantile, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

fd_clr.df=data.frame(clr_quantile=1:10,clr_mean=clr_mean,clr_sd=clr_sd,fd_clr_mean=fd_clr_mean,fd_clr_sd=fd_clr_sd)

#fd
ggplot(fd_clr.df,aes(as.factor(clr_quantile),y=fd_clr_mean,ymin=fd_clr_mean-fd_clr_sd,ymax=fd_clr_mean+fd_clr_sd))+
  geom_pointrange(color=cols[3],size=0.5)+
  scale_x_discrete(labels= c("0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%"))+
  xlab("")+ylab("")+ylim(c(0.12,0.19))+
  #xlab("Percentile of coding density")+ylab(expression(f[d]))+
  theme_bw()+theme(plot.margin = margin(1, 1, 0.5, 0.5, "cm"),axis.text=element_text(size=18),axis.text.x = element_blank())
ggsave("clr.fd.pdf",width=4,height=3)


####################################################
###4 selection effects og genes from P.tra

###Using selection effects og genes from P.tra
summary_10kb_w_selec_eff=summary_10kb_w %>% mutate(selec_eff_quantile=ntile(selection_eff_tra,10))

nls_selec_eff=data.frame(cbind(summary_10kb_w_selec_eff$selec_eff_quantile,summary_10kb_w_selec_eff$selection_eff_tra,nls_10kb_new_w$nls_n,twisst_10kb_prop_sites_w$top1,twisst_10kb_prop_sites_w$top2,twisst_10kb_prop_sites_w$top3))
names(nls_selec_eff)=c("selec_eff_quantile","tra_selec_eff","nls","top1","top2","top3")

nls_selec_eff_mean=as.vector(by(nls_selec_eff$nls, nls_selec_eff$selec_eff_quantile, function(x) mean(x,na.rm=T)))
nls_selec_eff_sd=as.vector(by(nls_selec_eff$nls, nls_selec_eff$selec_eff_quantile, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

top2_selec_eff_mean=as.vector(by(nls_selec_eff$top2, nls_selec_eff$selec_eff_quantile, function(x) mean(x,na.rm=T)))
top2_selec_eff_sd=as.vector(by(nls_selec_eff$top2, nls_selec_eff$selec_eff_quantile, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

non_top2_selec_eff_mean=as.vector(by(1-nls_selec_eff$top2, nls_selec_eff$selec_eff_quantile, function(x) mean(x,na.rm=T)))
non_top2_selec_eff_sd=as.vector(by(1-nls_selec_eff$top2, nls_selec_eff$selec_eff_quantile, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

selec_eff_mean=as.vector(by(nls_selec_eff$tra_selec_eff, nls_selec_eff$selec_eff_quantile, function(x) mean(x,na.rm=T)))
selec_eff_sd=as.vector(by(nls_selec_eff$tra_selec_eff, nls_selec_eff$selec_eff_quantile, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

nls_selec_eff.df=data.frame(selec_eff_quantile=1:10,selec_eff_mean=selec_eff_mean,selec_eff_sd=selec_eff_sd,nls_selec_eff_mean=nls_selec_eff_mean,nls_selec_eff_sd=nls_selec_eff_sd)
top2_selec_eff.df=data.frame(selec_eff_quantile=1:10,selec_eff_mean=selec_eff_mean,selec_eff_sd=selec_eff_sd,top2_selec_eff_mean=top2_selec_eff_mean,top2_selec_eff_sd=top2_selec_eff_sd)
non_top2_selec_eff.df=data.frame(selec_eff_quantile=1:10,selec_eff_mean=selec_eff_mean,selec_eff_sd=selec_eff_sd,non_top2_selec_eff_mean=non_top2_selec_eff_mean,non_top2_selec_eff_sd=non_top2_selec_eff_sd)


#####3.1 selec_eff vs. pi
###pi of P.tra
pi_tra_selec_eff=data.frame(cbind(summary_10kb_w_selec_eff$selec_eff_quantile,summary_10kb_w_selec_eff$selection_eff_tra,summary_10kb_w_selec_eff$tra_tP))
names(pi_tra_selec_eff)=c("selec_eff_quantile","tra_selec_eff","tra_tP")

###pi
pi_tra_selec_eff_mean=as.vector(by(pi_tra_selec_eff$tra_tP, pi_tra_selec_eff$selec_eff_quantile, function(x) mean(x,na.rm=T)))
pi_tra_selec_eff_sd=as.vector(by(pi_tra_selec_eff$tra_tP, pi_tra_selec_eff$selec_eff_quantile, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

tra_selec_eff_mean=as.vector(by(pi_tra_selec_eff$tra_selec_eff, pi_tra_selec_eff$selec_eff_quantile, function(x) mean(x,na.rm=T)))
tra_selec_eff_sd=as.vector(by(pi_tra_selec_eff$tra_selec_eff, pi_tra_selec_eff$selec_eff_quantile, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

pi_tra_selec_eff.df=data.frame(selec_eff_quantile=1:10,tra_selec_eff_mean=tra_selec_eff_mean,tra_selec_eff_sd=tra_selec_eff_sd,pi_tra_selec_eff_mean=pi_tra_selec_eff_mean,pi_tra_selec_eff_sd=pi_tra_selec_eff_sd)

#pi
ggplot(pi_tra_selec_eff.df,aes(as.factor(selec_eff_quantile),y=pi_tra_selec_eff_mean,ymin=pi_tra_selec_eff_mean-pi_tra_selec_eff_sd,ymax=pi_tra_selec_eff_mean+pi_tra_selec_eff_sd))+
  geom_pointrange(color=cols[4],size=0.5)+
  scale_x_discrete(labels= c("0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%"))+
  xlab("")+ylab("")+ylim(c(0.009,0.019))+
  # xlab("Percentile of coding density")+ylab("Scaled rate of ILS sites")+
  theme_bw()+theme(plot.margin = margin(1, 1, 0.5, 0.5, "cm"),axis.text=element_text(size=18),axis.text.x = element_blank())
ggsave("selec_eff.tra_pi.pdf",width=4,height=3)


#####3.2 selec_eff vs.fst
###fst of P.tra-Pdav
fst_tra_dav_selec_eff=data.frame(cbind(summary_10kb_w_selec_eff$selec_eff_quantile,summary_10kb_w_selec_eff$selection_eff_tra,summary_10kb_w_selec_eff$tra_dav_fst))
names(fst_tra_dav_selec_eff)=c("selec_eff_quantile","tra_selec_eff","tra_dav_fst")


fst_tra_dav_selec_eff_mean=as.vector(by(fst_tra_dav_selec_eff$tra_dav_fst, fst_tra_dav_selec_eff$selec_eff_quantile, function(x) mean(x,na.rm=T)))
fst_tra_dav_selec_eff_sd=as.vector(by(fst_tra_dav_selec_eff$tra_dav_fst, fst_tra_dav_selec_eff$selec_eff_quantile, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

tra_selec_eff_mean=as.vector(by(fst_tra_dav_selec_eff$tra_selec_eff, fst_tra_dav_selec_eff$selec_eff_quantile, function(x) mean(x,na.rm=T)))
tra_selec_eff_sd=as.vector(by(fst_tra_dav_selec_eff$tra_selec_eff, fst_tra_dav_selec_eff$selec_eff_quantile, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

fst_tra_dav_selec_eff.df=data.frame(selec_eff_quantile=1:10,tra_selec_eff_mean=tra_selec_eff_mean,tra_selec_eff_sd=tra_selec_eff_sd,fst_tra_dav_selec_eff_mean=fst_tra_dav_selec_eff_mean,fst_tra_dav_selec_eff_sd=fst_tra_dav_selec_eff_sd)

#fst
ggplot(fst_tra_dav_selec_eff.df,aes(as.factor(selec_eff_quantile),y=fst_tra_dav_selec_eff_mean,ymin=fst_tra_dav_selec_eff_mean-fst_tra_dav_selec_eff_sd,ymax=fst_tra_dav_selec_eff_mean+fst_tra_dav_selec_eff_sd))+
  geom_pointrange(color=cols[4],size=0.5)+
  scale_x_discrete(labels= c("0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%"))+
  xlab("")+ylab("")+ylim(c(0.27,0.41))+
  # xlab("Percentile of coding density")+ylab("Scaled rate of ILS sites")+
  theme_bw()+theme(plot.margin = margin(1, 1, 0.5, 0.5, "cm"),axis.text=element_text(size=18),axis.text.x = element_blank())
ggsave("selec_eff.tra_dav.fst.pdf",width=4,height=3)


#####4.3 select_eff vs.dxy
###dxy of P.tra-Pdav
dxy_tra_dav_selec_eff=data.frame(cbind(summary_10kb_w_selec_eff$selec_eff_quantile,summary_10kb_w_selec_eff$selection_eff_tra,summary_10kb_w_selec_eff$tra_dav_dxy))
names(dxy_tra_dav_selec_eff)=c("selec_eff_quantile","tra_selec_eff","tra_dav_dxy")


dxy_tra_dav_selec_eff_mean=as.vector(by(dxy_tra_dav_selec_eff$tra_dav_dxy, dxy_tra_dav_selec_eff$selec_eff_quantile, function(x) mean(x,na.rm=T)))
dxy_tra_dav_selec_eff_sd=as.vector(by(dxy_tra_dav_selec_eff$tra_dav_dxy, dxy_tra_dav_selec_eff$selec_eff_quantile, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

tra_selec_eff_mean=as.vector(by(dxy_tra_dav_selec_eff$tra_selec_eff, dxy_tra_dav_selec_eff$selec_eff_quantile, function(x) mean(x,na.rm=T)))
tra_selec_eff_sd=as.vector(by(dxy_tra_dav_selec_eff$tra_selec_eff, dxy_tra_dav_selec_eff$selec_eff_quantile, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

dxy_tra_dav_selec_eff.df=data.frame(selec_eff_quantile=1:10,tra_selec_eff_mean=tra_selec_eff_mean,tra_selec_eff_sd=tra_selec_eff_sd,dxy_tra_dav_selec_eff_mean=dxy_tra_dav_selec_eff_mean,dxy_tra_dav_selec_eff_sd=dxy_tra_dav_selec_eff_sd)

#dxy
ggplot(dxy_tra_dav_selec_eff.df,aes(as.factor(selec_eff_quantile),y=dxy_tra_dav_selec_eff_mean,ymin=dxy_tra_dav_selec_eff_mean-dxy_tra_dav_selec_eff_sd,ymax=dxy_tra_dav_selec_eff_mean+dxy_tra_dav_selec_eff_sd))+
  geom_pointrange(color=cols[4],size=0.5)+
  scale_x_discrete(labels= c("0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%"))+
  xlab("")+ylab("")+ylim(c(0.0185,0.0215))+
  # xlab("Percentile of coding density")+ylab("Scaled rate of ILS sites")+
  theme_bw()+theme(plot.margin = margin(1, 1, 0.5, 0.5, "cm"),axis.text.x=element_blank(),axis.text=element_text(size=18))
ggsave("selec_eff.tra_dav.dxy.pdf",width=4,height=3)

#####4.4 selec_eff vs.nls
#nls
ggplot(nls_selec_eff.df,aes(as.factor(selec_eff_quantile),y=nls_selec_eff_mean,ymin=nls_selec_eff_mean-nls_selec_eff_sd,ymax=nls_selec_eff_mean+nls_selec_eff_sd))+
  geom_pointrange(color=cols[4],size=0.5)+
  scale_x_discrete(labels= c("0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%"))+
  xlab("")+ylab("")+ylim(c(0.095,0.128))+
  #xlab("Percentile of coding density")+ylab("Scaled rate of ILS sites")+
  theme_bw()+theme(plot.margin = margin(1, 1, 0.5, 0.5, "cm"),axis.text=element_text(size=18),axis.text.x = element_blank())
ggsave("selec_eff.ILS.pdf",width=4,height=3)

#####4.5 selec_eff vs.top2
#the percentile of topology2 
ggplot(top2_selec_eff.df,aes(as.factor(selec_eff_quantile),y=top2_selec_eff_mean,ymin=top2_selec_eff_mean-top2_selec_eff_sd,ymax=top2_selec_eff_mean+top2_selec_eff_sd))+
  geom_pointrange(color=cols[4],size=0.5)+
  scale_x_discrete(labels= c("0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%"))+
  xlab("")+ylab("")+ylim(c(0.48,0.63))+
  #xlab("Percentile of coding density")+ylab("Species topology (%)")+
  theme_bw()+theme(plot.margin = margin(1, 1, 0.5, 0.5, "cm"),axis.text=element_text(size=18),axis.text.x = element_blank())
ggsave("selec_eff.top2.pdf",width=4,height=3)

#####4.6 selec_eff vs.fd
fd_selec_eff_old=data.frame(cbind(summary_10kb_w_selec_eff$selec_eff_quantile,summary_10kb_w_selec_eff$selection_eff_tra,D_gatk_10kb_w$fd))
names(fd_selec_eff_old)=c("selec_eff_quantile","tra_selec_eff","fd")
fd_selec_eff=fd_selec_eff_old[which(fd_selec_eff_old$fd>=0),]
cor.test(fd_selec_eff$tra_selec_eff,fd_selec_eff$fd,method="spearman")

fd_selec_eff_mean=as.vector(by(fd_selec_eff$fd, fd_selec_eff$selec_eff_quantile, function(x) mean(x,na.rm=T)))
fd_selec_eff_sd=as.vector(by(fd_selec_eff$fd, fd_selec_eff$selec_eff_quantile, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

selec_eff_mean=as.vector(by(fd_selec_eff$tra_selec_eff, fd_selec_eff$selec_eff_quantile, function(x) mean(x,na.rm=T)))
selec_eff_sd=as.vector(by(fd_selec_eff$tra_selec_eff, fd_selec_eff$selec_eff_quantile, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

fd_selec_eff.df=data.frame(selec_eff_quantile=1:10,selec_eff_mean=selec_eff_mean,selec_eff_sd=selec_eff_sd,fd_selec_eff_mean=fd_selec_eff_mean,fd_selec_eff_sd=fd_selec_eff_sd)

#fd
ggplot(fd_selec_eff.df,aes(as.factor(selec_eff_quantile),y=fd_selec_eff_mean,ymin=fd_selec_eff_mean-fd_selec_eff_sd,ymax=fd_selec_eff_mean+fd_selec_eff_sd))+
  geom_pointrange(color=cols[4],size=0.5)+
  scale_x_discrete(labels= c("0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%"))+
  xlab("")+ylab("")+ylim(c(0.12,0.19))+
  #xlab("Percentile of coding density")+ylab(expression(f[d]))+
  theme_bw()+theme(plot.margin = margin(1, 1, 0.5, 0.5, "cm"),axis.text=element_text(size=18),axis.text.x = element_blank())
ggsave("selec_eff.fd.pdf",width=4,height=3)


#############################################################################################
#############################################################################################
###making the plot and combining the coding density and recombination rate together to make the plot

setwd("./cd_rho/")
####################################################
###1. Coding prop

###Using coding density
summary_10kb_w_cd=summary_10kb_w %>% mutate(cd_quantile=ntile(Coding_prop,10))

nls_cd=data.frame(cbind(summary_10kb_w_cd$cd_quantile,summary_10kb_w_cd$Coding_prop,nls_10kb_new_w$nls_n,twisst_10kb_prop_sites_w$top1,twisst_10kb_prop_sites_w$top2,twisst_10kb_prop_sites_w$top3))
names(nls_cd)=c("cd_quantile","coding_prop","nls","top1","top2","top3")

nls_cd_mean=as.vector(by(nls_cd$nls, nls_cd$cd_quantile, function(x) mean(x,na.rm=T)))
nls_cd_sd=as.vector(by(nls_cd$nls, nls_cd$cd_quantile, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

top2_cd_mean=as.vector(by(nls_cd$top2, nls_cd$cd_quantile, function(x) mean(x,na.rm=T)))
top2_cd_sd=as.vector(by(nls_cd$top2, nls_cd$cd_quantile, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

non_top2_cd_mean=as.vector(by(1-nls_cd$top2, nls_cd$cd_quantile, function(x) mean(x,na.rm=T)))
non_top2_cd_sd=as.vector(by(1-nls_cd$top2, nls_cd$cd_quantile, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

cd_mean=as.vector(by(nls_cd$coding_prop, nls_cd$cd_quantile, function(x) mean(x,na.rm=T)))
cd_sd=as.vector(by(nls_cd$coding_prop, nls_cd$cd_quantile, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

nls_cd.df=data.frame(cd_quantile=1:10,cd_mean=cd_mean,cd_sd=cd_sd,nls_cd_mean=nls_cd_mean,nls_cd_sd=nls_cd_sd)
top2_cd.df=data.frame(cd_quantile=1:10,cd_mean=cd_mean,cd_sd=cd_sd,top2_cd_mean=top2_cd_mean,top2_cd_sd=top2_cd_sd)
non_top2_cd.df=data.frame(cd_quantile=1:10,cd_mean=cd_mean,cd_sd=cd_sd,non_top2_cd_mean=non_top2_cd_mean,non_top2_cd_sd=non_top2_cd_sd)


####################################################
###2 recombination rate from P.tra

###Using coding density
summary_10kb_w_rho=summary_10kb_w %>% mutate(rho_quantile=ntile(tra_rho_mean,10))

nls_rho=data.frame(cbind(summary_10kb_w_rho$rho_quantile,summary_10kb_w_rho$tra_rho_mean,nls_10kb_new_w$nls_n,twisst_10kb_prop_sites_w$top1,twisst_10kb_prop_sites_w$top2,twisst_10kb_prop_sites_w$top3))
names(nls_rho)=c("rho_quantile","tra_rho","nls","top1","top2","top3")

nls_rho_mean=as.vector(by(nls_rho$nls, nls_rho$rho_quantile, function(x) mean(x,na.rm=T)))
nls_rho_sd=as.vector(by(nls_rho$nls, nls_rho$rho_quantile, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

top2_rho_mean=as.vector(by(nls_rho$top2, nls_rho$rho_quantile, function(x) mean(x,na.rm=T)))
top2_rho_sd=as.vector(by(nls_rho$top2, nls_rho$rho_quantile, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

non_top2_rho_mean=as.vector(by(1-nls_rho$top2, nls_rho$rho_quantile, function(x) mean(x,na.rm=T)))
non_top2_rho_sd=as.vector(by(1-nls_rho$top2, nls_rho$rho_quantile, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

rho_mean=as.vector(by(nls_rho$tra_rho, nls_rho$rho_quantile, function(x) mean(x,na.rm=T)))
rho_sd=as.vector(by(nls_rho$tra_rho, nls_rho$rho_quantile, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

nls_rho.df=data.frame(rho_quantile=1:10,rho_mean=rho_mean,rho_sd=rho_sd,nls_rho_mean=nls_rho_mean,nls_rho_sd=nls_rho_sd)
top2_rho.df=data.frame(rho_quantile=1:10,rho_mean=rho_mean,rho_sd=rho_sd,top2_rho_mean=top2_rho_mean,top2_rho_sd=top2_rho_sd)
non_top2_rho.df=data.frame(rho_quantile=1:10,rho_mean=rho_mean,rho_sd=rho_sd,non_top2_rho_mean=non_top2_rho_mean,non_top2_rho_sd=non_top2_rho_sd)

############################
###pi of P.tra
#####1.1 cd vs. pi
pi_tra_cd=data.frame(cbind(summary_10kb_w_cd$cd_quantile,summary_10kb_w_cd$Coding_prop,summary_10kb_w_cd$tra_tP))
names(pi_tra_cd)=c("cd_quantile","coding_prop","tra_tP")

pi_tra_cd_mean=as.vector(by(pi_tra_cd$tra_tP, pi_tra_cd$cd_quantile, function(x) mean(x,na.rm=T)))
pi_tra_cd_sd=as.vector(by(pi_tra_cd$tra_tP, pi_tra_cd$cd_quantile, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

tra_cd_mean=as.vector(by(pi_tra_cd$coding_prop, pi_tra_cd$cd_quantile, function(x) mean(x,na.rm=T)))
tra_cd_sd=as.vector(by(pi_tra_cd$coding_prop, pi_tra_cd$cd_quantile, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

pi_tra_cd.df=data.frame(cd_quantile=1:10,tra_cd_mean=tra_cd_mean,tra_cd_sd=tra_cd_sd,pi_tra_cd_mean=pi_tra_cd_mean,pi_tra_cd_sd=pi_tra_cd_sd)


#####1.2 rho vs. pi

pi_tra_rho=data.frame(cbind(summary_10kb_w_rho$rho_quantile,summary_10kb_w_rho$Coding_prop,summary_10kb_w_rho$tra_tP))
names(pi_tra_rho)=c("rho_quantile","tra_rho","tra_tP")

###pi
pi_tra_rho_mean=as.vector(by(pi_tra_rho$tra_tP, pi_tra_rho$rho_quantile, function(x) mean(x,na.rm=T)))
pi_tra_rho_sd=as.vector(by(pi_tra_rho$tra_tP, pi_tra_rho$rho_quantile, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

tra_rho_mean=as.vector(by(pi_tra_rho$tra_rho, pi_tra_rho$rho_quantile, function(x) mean(x,na.rm=T)))
tra_rho_sd=as.vector(by(pi_tra_rho$tra_rho, pi_tra_rho$rho_quantile, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

pi_tra_rho.df=data.frame(rho_quantile=1:10,tra_rho_mean=tra_rho_mean,tra_rho_sd=tra_rho_sd,pi_tra_rho_mean=pi_tra_rho_mean,pi_tra_rho_sd=pi_tra_rho_sd)
pi_tra_rho.df$class="Recombination"

######create new file for cd and rho and merge them together
pi_tra_cd.df.new=pi_tra_cd.df[,c("cd_quantile","pi_tra_cd_mean","pi_tra_cd_sd","class")]
names(pi_tra_cd.df.new)=c("quantile","pi_mean","pi_sd","class")
pi_tra_rho.df.new=pi_tra_rho.df[,c("rho_quantile","pi_tra_rho_mean","pi_tra_rho_sd","class")]
names(pi_tra_rho.df.new)=c("quantile","pi_mean","pi_sd","class")
#merge
pi_tra_cd_rho.df=data.frame(rbind(pi_tra_cd.df.new,pi_tra_rho.df.new))
pi_tra_cd_rho.df$quantile=factor(pi_tra_cd_rho.df$quantile,levels=c("1","2","3","4","5","6","7","8","9","10"))


#plot
ggplot(pi_tra_cd_rho.df,aes(as.factor(quantile),y=pi_mean,ymin=pi_mean-pi_sd,ymax=pi_mean+pi_sd))+
 # geom_pointrange(color=cols[2],size=0.5)+
  geom_pointrange(aes(color=class))+
  scale_x_discrete(labels= c("0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%"))+
  xlab("")+ylab("")+ylim(c(0.009,0.019))+
  # xlab("Percentile of coding density")+ylab("Scaled rate of ILS sites")+
 theme_bw()+theme(legend.position="none",plot.margin = margin(1, 1, 0.5, 0.5, "cm"),axis.text=element_text(size=16),axis.text.x = element_blank())
 # theme_bw()+theme(legend.position="none",plot.margin = margin(1, 1, 0.5, 0.5, "cm"),axis.text=element_text(size=16),axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave("cd.rho.tra_pi.pdf",width=5,height=4)

############################################################
###2.fst of P.tra-Pdav
#####2.1 cd vs.fst
fst_tra_dav_cd=data.frame(cbind(summary_10kb_w_cd$cd_quantile,summary_10kb_w_cd$Coding_prop,summary_10kb_w_cd$tra_dav_fst))
names(fst_tra_dav_cd)=c("cd_quantile","coding_prop","tra_dav_fst")


fst_tra_dav_cd_mean=as.vector(by(fst_tra_dav_cd$tra_dav_fst, fst_tra_dav_cd$cd_quantile, function(x) mean(x,na.rm=T)))
fst_tra_dav_cd_sd=as.vector(by(fst_tra_dav_cd$tra_dav_fst, fst_tra_dav_cd$cd_quantile, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

tra_cd_mean=as.vector(by(fst_tra_dav_cd$coding_prop, fst_tra_dav_cd$cd_quantile, function(x) mean(x,na.rm=T)))
tra_cd_sd=as.vector(by(fst_tra_dav_cd$coding_prop, fst_tra_dav_cd$cd_quantile, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

fst_tra_dav_cd.df=data.frame(cd_quantile=1:10,tra_cd_mean=tra_cd_mean,tra_cd_sd=tra_cd_sd,fst_tra_dav_cd_mean=fst_tra_dav_cd_mean,fst_tra_dav_cd_sd=fst_tra_dav_cd_sd)
fst_tra_dav_cd.df$class="Coding"

#####2.2 rho vs.fst

fst_tra_dav_rho=data.frame(cbind(summary_10kb_w_rho$rho_quantile,summary_10kb_w_rho$tra_rho_mean,summary_10kb_w_rho$tra_dav_fst))
names(fst_tra_dav_rho)=c("rho_quantile","tra_rho","tra_dav_fst")

fst_tra_dav_rho_mean=as.vector(by(fst_tra_dav_rho$tra_dav_fst, fst_tra_dav_rho$rho_quantile, function(x) mean(x,na.rm=T)))
fst_tra_dav_rho_sd=as.vector(by(fst_tra_dav_rho$tra_dav_fst, fst_tra_dav_rho$rho_quantile, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

tra_rho_mean=as.vector(by(fst_tra_dav_rho$tra_rho, fst_tra_dav_rho$rho_quantile, function(x) mean(x,na.rm=T)))
tra_rho_sd=as.vector(by(fst_tra_dav_rho$tra_rho, fst_tra_dav_rho$rho_quantile, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

fst_tra_dav_rho.df=data.frame(rho_quantile=1:10,tra_rho_mean=tra_rho_mean,tra_rho_sd=tra_rho_sd,fst_tra_dav_rho_mean=fst_tra_dav_rho_mean,fst_tra_dav_rho_sd=fst_tra_dav_rho_sd)
fst_tra_dav_rho.df$class="Recombination"

######create new file for cd and rho and merge them together
fst_tra_dav_cd.df.new=fst_tra_dav_cd.df[,c("cd_quantile","fst_tra_dav_cd_mean","fst_tra_dav_cd_sd","class")]
names(fst_tra_dav_cd.df.new)=c("quantile","fst_mean","fst_sd","class")
fst_tra_dav_rho.df.new=fst_tra_dav_rho.df[,c("rho_quantile","fst_tra_dav_rho_mean","fst_tra_dav_rho_sd","class")]
names(fst_tra_dav_rho.df.new)=c("quantile","fst_mean","fst_sd","class")
#merge
fst_cd_rho.df=data.frame(rbind(fst_tra_dav_cd.df.new,fst_tra_dav_rho.df.new))
fst_cd_rho.df$quantile=factor(fst_cd_rho.df$quantile,levels=c("1","2","3","4","5","6","7","8","9","10"))


#fst
ggplot(fst_cd_rho.df,aes(as.factor(quantile),y=fst_mean,ymin=fst_mean-fst_sd,ymax=fst_mean+fst_sd))+
  geom_pointrange(aes(color=class))+
  scale_x_discrete(labels= c("0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%"))+
  xlab("")+ylab("")+ylim(c(0.28,0.37))+
  # xlab("Percentile of coding density")+ylab("Scaled rate of ILS sites")+
  theme_bw()+theme(legend.position="none",plot.margin = margin(1, 1, 0.5, 0.5, "cm"),axis.text=element_text(size=16),axis.text.x = element_blank())
ggsave("cd.rho.fst.pdf",width=5,height=4)

#################################################
###3dxy of P.tra-Pdav
#####3.1 cd vs.dxy
dxy_tra_dav_cd=data.frame(cbind(summary_10kb_w_cd$cd_quantile,summary_10kb_w_cd$Coding_prop,summary_10kb_w_cd$tra_dav_dxy))
names(dxy_tra_dav_cd)=c("cd_quantile","coding_prop","tra_dav_dxy")

dxy_tra_dav_cd_mean=as.vector(by(dxy_tra_dav_cd$tra_dav_dxy, dxy_tra_dav_cd$cd_quantile, function(x) mean(x,na.rm=T)))
dxy_tra_dav_cd_sd=as.vector(by(dxy_tra_dav_cd$tra_dav_dxy, dxy_tra_dav_cd$cd_quantile, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

tra_cd_mean=as.vector(by(dxy_tra_dav_cd$coding_prop, dxy_tra_dav_cd$cd_quantile, function(x) mean(x,na.rm=T)))
tra_cd_sd=as.vector(by(dxy_tra_dav_cd$coding_prop, dxy_tra_dav_cd$cd_quantile, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

dxy_tra_dav_cd.df=data.frame(cd_quantile=1:10,tra_cd_mean=tra_cd_mean,tra_cd_sd=tra_cd_sd,dxy_tra_dav_cd_mean=dxy_tra_dav_cd_mean,dxy_tra_dav_cd_sd=dxy_tra_dav_cd_sd)
dxy_tra_dav_cd.df$class="Coding"

#####3.2 rho vs.dxy
dxy_tra_dav_rho=data.frame(cbind(summary_10kb_w_rho$rho_quantile,summary_10kb_w_rho$tra_rho_mean,summary_10kb_w_rho$tra_dav_dxy))
names(dxy_tra_dav_rho)=c("rho_quantile","tra_rho","tra_dav_dxy")

dxy_tra_dav_rho_mean=as.vector(by(dxy_tra_dav_rho$tra_dav_dxy, dxy_tra_dav_rho$rho_quantile, function(x) mean(x,na.rm=T)))
dxy_tra_dav_rho_sd=as.vector(by(dxy_tra_dav_rho$tra_dav_dxy, dxy_tra_dav_rho$rho_quantile, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

tra_rho_mean=as.vector(by(dxy_tra_dav_rho$tra_rho, dxy_tra_dav_rho$rho_quantile, function(x) mean(x,na.rm=T)))
tra_rho_sd=as.vector(by(dxy_tra_dav_rho$tra_rho, dxy_tra_dav_rho$rho_quantile, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

dxy_tra_dav_rho.df=data.frame(rho_quantile=1:10,tra_rho_mean=tra_rho_mean,tra_rho_sd=tra_rho_sd,dxy_tra_dav_rho_mean=dxy_tra_dav_rho_mean,dxy_tra_dav_rho_sd=dxy_tra_dav_rho_sd)
dxy_tra_dav_rho.df$class="Recombination"

######create new file for cd and rho and merge them together
dxy_tra_dav_cd.df.new=dxy_tra_dav_cd.df[,c("cd_quantile","dxy_tra_dav_cd_mean","dxy_tra_dav_cd_sd","class")]
names(dxy_tra_dav_cd.df.new)=c("quantile","dxy_mean","dxy_sd","class")
dxy_tra_dav_rho.df.new=dxy_tra_dav_rho.df[,c("rho_quantile","dxy_tra_dav_rho_mean","dxy_tra_dav_rho_sd","class")]
names(dxy_tra_dav_rho.df.new)=c("quantile","dxy_mean","dxy_sd","class")
#merge
dxy_cd_rho.df=data.frame(rbind(dxy_tra_dav_cd.df.new,dxy_tra_dav_rho.df.new))
dxy_cd_rho.df$quantile=factor(dxy_cd_rho.df$quantile,levels=c("1","2","3","4","5","6","7","8","9","10"))

###Plot
ggplot(dxy_cd_rho.df,aes(as.factor(quantile),y=dxy_mean,ymin=dxy_mean-dxy_sd,ymax=dxy_mean+dxy_sd))+
  geom_pointrange(aes(color=class))+
  scale_x_discrete(labels= c("0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%"))+
  xlab("")+ylab("")+ylim(c(0.0185,0.0215))+
  # xlab("Percentile of coding density")+ylab("Scaled rate of ILS sites")+
  theme_bw()+theme(legend.position="none",plot.margin = margin(1, 1, 0.5, 0.5, "cm"),axis.text=element_text(size=16),axis.text.x = element_blank())
ggsave("cd.rho.dxy.pdf",width=5,height=4)


#########################################################
#4.nls
#####4.1 cd vs.nls
nls_cd.df$class="Coding"
nls_rho.df$class="Recombination"

######create new file for cd and rho and merge them together
nls_cd.df.new=nls_cd.df[,c("cd_quantile","nls_cd_mean","nls_cd_sd","class")]
names(nls_cd.df.new)=c("quantile","nls_mean","nls_sd","class")
nls_rho.df.new=nls_rho.df[,c("rho_quantile","nls_rho_mean","nls_rho_sd","class")]
names(nls_rho.df.new)=c("quantile","nls_mean","nls_sd","class")
#merge
nls_cd_rho.df=data.frame(rbind(nls_cd.df.new,nls_rho.df.new))
nls_cd_rho.df$quantile=factor(nls_cd_rho.df$quantile,levels=c("1","2","3","4","5","6","7","8","9","10"))

###Plot
ggplot(nls_cd_rho.df,aes(as.factor(quantile),y=nls_mean,ymin=nls_mean-nls_sd,ymax=nls_mean+nls_sd))+
  geom_pointrange(aes(color=class))+
  scale_x_discrete(labels= c("0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%"))+
  xlab("")+ylab("")+ylim(c(0.095,0.128))+
  # xlab("Percentile of coding density")+ylab("Scaled rate of ILS sites")+
  theme_bw()+theme(legend.position="none",plot.margin = margin(1, 1, 0.5, 0.5, "cm"),axis.text=element_text(size=16),axis.text.x = element_blank())
ggsave("cd.rho.nls.pdf",width=5,height=4)

#######################################
#5.the percentile of topology2
#####5.1 cd vs.top2
top2_cd.df$class="Coding" 
top2_rho.df$class="Recombination"

######create new file for cd and rho and merge them together
top2_cd.df.new=top2_cd.df[,c("cd_quantile","top2_cd_mean","top2_cd_sd","class")]
names(top2_cd.df.new)=c("quantile","top2_mean","top2_sd","class")
top2_rho.df.new=top2_rho.df[,c("rho_quantile","top2_rho_mean","top2_rho_sd","class")]
names(top2_rho.df.new)=c("quantile","top2_mean","top2_sd","class")
#merge
top2_cd_rho.df=data.frame(rbind(top2_cd.df.new,top2_rho.df.new))
top2_cd_rho.df$quantile=factor(top2_cd_rho.df$quantile,levels=c("1","2","3","4","5","6","7","8","9","10"))

###Plot
ggplot(top2_cd_rho.df,aes(as.factor(quantile),y=top2_mean,ymin=top2_mean-top2_sd,ymax=top2_mean+top2_sd))+
  geom_pointrange(aes(color=class))+
  scale_x_discrete(labels= c("0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%"))+
  xlab("")+ylab("")+ylim(c(0.49,0.61))+
  # xlab("Percentile of coding density")+ylab("Scaled rate of ILS sites")+
  theme_bw()+theme(legend.position="none",plot.margin = margin(1, 1, 0.5, 0.5, "cm"),axis.text=element_text(size=16),axis.text.x = element_blank())
ggsave("cd.rho.top2.pdf",width=5,height=4)


#############################################
###6.fd
#####6.1 cd vs.fd
fd_cd_old=data.frame(cbind(summary_10kb_w_cd$cd_quantile,summary_10kb_w_cd$Coding_prop,D_gatk_10kb_w$fd))
names(fd_cd_old)=c("cd_quantile","coding_prop","fd")
fd_cd=fd_cd_old[which(fd_cd_old$fd>=0),]
cor.test(fd_cd$coding_prop,fd_cd$fd,method="spearman")

fd_cd_mean=as.vector(by(fd_cd$fd, fd_cd$cd_quantile, function(x) mean(x,na.rm=T)))
fd_cd_sd=as.vector(by(fd_cd$fd, fd_cd$cd_quantile, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

cd_mean=as.vector(by(fd_cd$coding_prop, fd_cd$cd_quantile, function(x) mean(x,na.rm=T)))
cd_sd=as.vector(by(fd_cd$coding_prop, fd_cd$cd_quantile, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

fd_cd.df=data.frame(cd_quantile=1:10,cd_mean=cd_mean,cd_sd=cd_sd,fd_cd_mean=fd_cd_mean,fd_cd_sd=fd_cd_sd)
fd_cd.df$class="Coding"

#####6.2 rho vs.fd
fd_rho_old=data.frame(cbind(summary_10kb_w_rho$rho_quantile,summary_10kb_w_rho$tra_rho_mean,D_gatk_10kb_w$fd))
names(fd_rho_old)=c("rho_quantile","tra_rho","fd")
fd_rho=fd_rho_old[which(fd_rho_old$fd>=0),]
cor.test(fd_rho$tra_rho,fd_rho$fd,method="spearman")

fd_rho_mean=as.vector(by(fd_rho$fd, fd_rho$rho_quantile, function(x) mean(x,na.rm=T)))
fd_rho_sd=as.vector(by(fd_rho$fd, fd_rho$rho_quantile, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

rho_mean=as.vector(by(fd_rho$tra_rho, fd_rho$rho_quantile, function(x) mean(x,na.rm=T)))
rho_sd=as.vector(by(fd_rho$tra_rho, fd_rho$rho_quantile, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

fd_rho.df=data.frame(rho_quantile=1:10,rho_mean=rho_mean,rho_sd=rho_sd,fd_rho_mean=fd_rho_mean,fd_rho_sd=fd_rho_sd)
fd_rho.df$class="Recombination"


######create new file for cd and rho and merge them together
fd_cd.df.new=fd_cd.df[,c("cd_quantile","fd_cd_mean","fd_cd_sd","class")]
names(fd_cd.df.new)=c("quantile","fd_mean","fd_sd","class")
fd_rho.df.new=fd_rho.df[,c("rho_quantile","fd_rho_mean","fd_rho_sd","class")]
names(fd_rho.df.new)=c("quantile","fd_mean","fd_sd","class")
#merge
fd_cd_rho.df=data.frame(rbind(fd_cd.df.new,fd_rho.df.new))
fd_cd_rho.df$quantile=factor(fd_cd_rho.df$quantile,levels=c("1","2","3","4","5","6","7","8","9","10"))

###Plot
ggplot(fd_cd_rho.df,aes(as.factor(quantile),y=fd_mean,ymin=fd_mean-fd_sd,ymax=fd_mean+fd_sd))+
  geom_pointrange(aes(color=class))+
  scale_x_discrete(labels= c("0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%"))+
  xlab("")+ylab("")+ylim(c(0.12,0.19))+
  # xlab("Percentile of coding density")+ylab("Scaled rate of ILS sites")+
  theme_bw()+theme(legend.position="none",plot.margin = margin(1, 1, 0.5, 0.5, "cm"),axis.text=element_text(size=16),axis.text.x = element_blank())
ggsave("cd.rho.fd.pdf",width=5,height=4)


