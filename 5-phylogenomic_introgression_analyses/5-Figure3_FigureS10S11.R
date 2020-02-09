library(circlize)
library(data.table)
library(RColorBrewer)
library(dplyr)
library(reshape2)
library(ggplot2)
#install.packages("gridExtra")
library(gridExtra)

colors <- brewer.pal(12,"Paired")

setwd("~/Dropbox/davidiana_paper/data/abba_baba/")

#######################ABBA-BABA##########################
###Main aim: Determine the relative frequency of phylogenetic tree: D-statistics (ABBA-BABA),Twisst weighting tree for the 4 species
##To determine the effects of window size on the results, I used two window size: (1)10kb (2)100kb

####For results from ANGSD
win<-function(data){
  data$Win=paste(data$CHR,data$BLOCKend,sep=":")
  return(data)
}

###For results from GATK
win2<-function(data){ ####for divergence estimates, the Chr is replaced by Pop
  data$Win=paste(data$scaffold,data$end,sep=":")
  return(data)
}



##1. 10kb
##1.1 abba-baba
###1.1.1 gatk results
abbababa_gatk_10kb=read.csv("gatk_data/4species.gatk.hap.snp.rm_indel.bed.biallelic.DP5.GQ10.rm_missing.recode.gt.beagle.10kb.abbababa.csv")
abbababa_gatk_10kb_new=abbababa_gatk_10kb %>% filter(sites>50)
abbababa_gatk_10kb_new_w=win2(abbababa_gatk_10kb_new)

###distribution plot
#mean D=-0.1076366
ggplot(abbababa_gatk_10kb_new,aes(-D))+geom_histogram(binwidth=0.01)+xlab("D-statistic value")+
geom_vline(aes(xintercept=mean(-D,na.rm=T)),color="dark red",linetype="dashed",size=0.8)
 
ggsave(filename="4species.gatk.abbababa.10kb.png",width=5,height=5,dpi=300)

###1.1.2 ANGSD results
abbababa_angsd_10kb=fread("angsd/10kb/4species.all.outDstat.Observed.txt",header=T)
abbababa_angsd_10kb_sites=fread(input = 'zcat < angsd/10kb/4species.abbababa.all.10000.out.abbababa2.gz',header=T)
abbababa_angsd_10kb$numSites=abbababa_angsd_10kb_sites$numSites
abbababa_angsd_10kb_new=abbababa_angsd_10kb %>% filter(numSites>1000)
abbababa_angsd_10kb_new_w=win(abbababa_angsd_10kb_new)

####common windows of ABBA-BABA results from ANGSD and GATK
common_win=Reduce(intersect,list(abbababa_gatk_10kb_new_w$Win,
                                 abbababa_angsd_10kb_new_w$Win))

abbababa_gatk_10kb_new_w_c=abbababa_gatk_10kb_new_w[which(abbababa_gatk_10kb_new_w$Win %in% common_win),]
abbababa_angsd_10kb_new_w_c=abbababa_angsd_10kb_new_w[which(abbababa_angsd_10kb_new_w$Win %in% common_win),]

#####making the correlation plots between D-statistics estimated from ANGSD and GATK
png("abbababa.gatk_angsd.10kb.png",width = 6, height = 6, units = 'in',res=300)
par(mar=c(4,4,1,1))
pal <- colorRampPalette(c("light blue", "yellow", "red"))
colors_abbababa=densCols(abbababa_gatk_10kb_new_w_c$D,abbababa_angsd_10kb_new_w_c$`mean(D)`,colramp=pal)
plot(-abbababa_gatk_10kb_new_w_c$D,abbababa_angsd_10kb_new_w_c$`mean(D)`,pch=19,col=colors_abbababa,cex=.5,cex.lab=1,cex.axis=0.7,xlab="D-statistic (SNPs)",ylab="D-statistic (ANGSD)")
cor.test(-abbababa_gatk_10kb_new_w_c$D,abbababa_angsd_10kb_new_w_c$`mean(D)`,method="spearman")
#text(0.0175,0.06,expression(paste(rho,"=0.829")^"***"))
legend("topleft",lty=NA,lwd=2,pch=NA, bty="n",text.font=1,legend=expression(paste(rho,"=0.758")^"***"))
dev.off()

#####making the correlation plots between D adn fdM estimated by GATK
png("D_fdM.10kb.png",width = 6, height = 6, units = 'in',res=300)
par(mar=c(4,4,1,1))
pal <- colorRampPalette(c("light blue", "yellow", "red"))
colors_abbababa=densCols(-abbababa_gatk_10kb_new_w_c$D,-abbababa_gatk_10kb_new_w_c$fdM,colramp=pal)
plot(-abbababa_gatk_10kb_new_w_c$D,-abbababa_gatk_10kb_new_w_c$fdM,pch=19,col=colors_abbababa,cex=.5,cex.lab=1,cex.axis=0.7,xlab=expression(italic(D)),ylab=expression(italic(f[d])))
cor.test(-abbababa_gatk_10kb_new_w_c$D,-abbababa_gatk_10kb_new_w_c$fdM,method="spearman")
#text(0.0175,0.06,expression(paste(rho,"=0.829")^"***"))
legend("topleft",lty=NA,lwd=2,pch=NA, bty="n",text.font=1,legend=expression(paste(rho,"=0.941")^"***"))
dev.off()

###Summary of the statistis according to chromsomes
mean_d_gatk=aggregate(-D~scaffold,data=abbababa_gatk_10kb_new_w_c,function(x){mean(x)})
se_d_gatk=aggregate(-D~scaffold,data=abbababa_gatk_10kb_new_w_c,function(x){sd(x)/sqrt(length(x))})

mean_fd_gatk=aggregate(-fdM~scaffold,data=abbababa_gatk_10kb_new_w_c,function(x){mean(x)})
se_fd_gatk=aggregate(-fdM~scaffold,data=abbababa_gatk_10kb_new_w_c,function(x){sd(x)/sqrt(length(x))})

mean_D_angsd=aggregate(`mean(D)`~CHR,data=abbababa_angsd_10kb_new_w_c,function(x){mean(x)})
names(mean_D_angsd)=c("scaffold","D")
se_D_angsd=aggregate(`mean(D)`~CHR,data=abbababa_angsd_10kb_new_w_c,function(x){sd(x)/sqrt(length(x))})
names(se_D_angsd)=c("scaffold","D_se")

abbababa_summary=as.data.frame(cbind(as.character(mean_d_gatk$scaffold),as.numeric(as.character(mean_d_gatk$`-D`)),as.numeric(se_d_gatk$`-D`),as.numeric(mean_fd_gatk$`-fdM`),as.numeric(se_fd_gatk$`-fdM`),as.numeric(mean_D_angsd$D),as.numeric(se_D_angsd$D_se)))
abbababa_summary$V2=as.numeric(as.character(abbababa_summary$V2))
abbababa_summary$V3=as.numeric(as.character(abbababa_summary$V3))
abbababa_summary$V4=as.numeric(as.character(abbababa_summary$V4))
abbababa_summary$V5=as.numeric(as.character(abbababa_summary$V5))
abbababa_summary$V6=as.numeric(as.character(abbababa_summary$V6))
abbababa_summary$V7=as.numeric(as.character(abbababa_summary$V7))

is.num <- sapply(abbababa_summary, is.numeric)
abbababa_summary[is.num] <- lapply(abbababa_summary[is.num], round, 4)
names(abbababa_summary)=c("CHR","D_gatk","D(se)_gatk","fd_gatk","fd(se)_gatk","D_angsd","D(se)_angsd")


###1.2 Twisst weightings
##topologies
#topo1 ((P.tremula,P.tremuloides),P.davidiana,P.trichocarpa);
#topo2 ((P.tremula,P.davidiana),P.tremuloides,P.trichocarpa);
#topo3 ((P.tremula,P.trichocarpa),P.tremuloides,P.davidiana);

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

####histgram plot of the proportion of incomplete lineage sorting (top1 and top3) across windows
ggplot(twisst_10kb_prop_sites,aes(1-top2))+geom_histogram(binwidth=0.01)+xlab("ILS")+
  geom_vline(aes(xintercept=mean(1-top2,na.rm=T)),color="dark red",linetype="dashed",size=0.8)
ggsave(filename="4species.twisst.ILS.10kb.png",width=5,height=5,dpi=300)

##top1
ggplot(twisst_10kb_prop_sites,aes(top1))+geom_histogram(binwidth=0.01)+xlab("Top1")+
  geom_vline(aes(xintercept=mean(top1,na.rm=T)),color="dark red",linetype="dashed",size=0.8)
ggsave(filename="4species.twisst.top1.10kb.png",width=5,height=5,dpi=300)
##top2
ggplot(twisst_10kb_prop_sites,aes(top2))+geom_histogram(binwidth=0.01)+xlab("Top2")+
  geom_vline(aes(xintercept=mean(top2,na.rm=T)),color="dark red",linetype="dashed",size=0.8)
ggsave(filename="4species.twisst.top2.10kb.png",width=5,height=5,dpi=300)
##top3
ggplot(twisst_10kb_prop_sites,aes(top3))+geom_histogram(binwidth=0.01)+xlab("Top3")+
  geom_vline(aes(xintercept=mean(top3,na.rm=T)),color="dark red",linetype="dashed",size=0.8)
ggsave(filename="4species.twisst.top3.10kb.png",width=5,height=5,dpi=300)

####plot top1,2,3 in a single plot

topology_10kb=twisst_10kb_prop_sites[,c("top1","top2","top3")]
names(topology_10kb)=c("Top1","Top2","Top3")
topology_10kb_melt=melt(topology_10kb)
names(topology_10kb_melt)=c("Topology","value")

ggplot(topology_10kb_melt,aes(x=value,fill=Topology))+geom_histogram(binwidth=.03, position="dodge")+
  scale_fill_manual(values=c(colors[1], colors[7], colors[9]))+xlab("Proportion")+ylab("Count")
ggsave(filename="4species.twisst.topology.prop.10kb.pdf",width=6,height=5,dpi=300)



###Summary across chromosomes
mean_top1_10kb_twisst=aggregate(top1~scaffold,data=twisst_10kb_prop_sites,function(x){mean(x)})
mean_top2_10kb_twisst=aggregate(top2~scaffold,data=twisst_10kb_prop_sites,function(x){mean(x)})
mean_top3_10kb_twisst=aggregate(top3~scaffold,data=twisst_10kb_prop_sites,function(x){mean(x)})

mean_top_10kb_twisst=cbind(mean_top1_10kb_twisst,mean_top2_10kb_twisst[,"top2"],mean_top3_10kb_twisst[,"top3"])
names(mean_top_10kb_twisst)=c("scaffold","top1","top2","top3")

#########################
##compare the candidate "sex" chrom (chr19) with other chroms
chr19_10kb=twisst_10kb_prop_sites[which(twisst_10kb_prop_sites$scaffold=="Chr19"),]
other_chr_10kb=twisst_10kb_prop_sites[which(twisst_10kb_prop_sites$scaffold!="Chr19"),]


###2. 100kb
##2.1 abbababa
###P1: P.tremula; P2:P.davidiana; P3: P. tremuloides; O: P.trichocarpa
###In order to make the gatk results and angsd resutls consistent, I will transfer the D and f value to -D value for GATK

abbababa_gatk_100kb=read.csv("gatk_data/4species.gatk.hap.snp.rm_indel.bed.biallelic.DP5.GQ10.rm_missing.recode.gt.beagle.100kb.abbababa.csv")
abbababa_gatk_100kb_new=abbababa_gatk_100kb %>% filter(sites>=200)

##mean D=-0.1103929
ggplot(abbababa_gatk_100kb_new,aes(-D))+geom_histogram(binwidth=0.01)+xlab("D-statistic value")+
  geom_vline(aes(xintercept=mean(-D,na.rm=T)),color="dark red",linetype="dashed",size=0.8)
ggsave(filename="4species.gatk.abbababa.100kb.png",width=5,height=5,dpi=300)


##2.2 Twisst
twisst_100kb=read.table("../twisst/data/4species.gatk.hap.snp.rm_indel.bed.biallelic.DP5.GQ10.rm_missing.recode.gt.beagle.twisst.weights.100kb.csv")
twisst_100kb_new=twisst_100kb[-1,]
names(twisst_100kb_new)=c("top1","top2","top3")

indx <- sapply(twisst_100kb_new, is.factor)
twisst_100kb_new[indx] <- lapply(twisst_100kb_new[indx], function(x) as.numeric(as.character(x)))
colSums(twisst_100kb_new)/sum(colSums(twisst_100kb_new))
head(rowSums(twisst_100kb_new))

twisst_100kb_prop=prop.table(as.matrix(twisst_100kb_new),1)
twisst_100kb_prop

###twisst sites
twisst_100kb_sites=read.table("../twisst/data/4species.gatk.hap.snp.rm_indel.bed.biallelic.DP5.GQ10.rm_missing.recode.gt.beagle.raxml.w100kb.data.tsv",header=T)
twisst_100kb_sites_new=twisst_100kb_sites[which(twisst_100kb_sites$sites>200),]
twisst_100kb_prop_sites=cbind(twisst_100kb_sites_new,twisst_100kb_prop)


#########################
##compare the candidate "sex" chrom (chr19) with other chroms
chr19_100kb=twisst_100kb_prop_sites[which(twisst_100kb_prop_sites$scaffold=="Chr19"),]
other_chr_100kb=twisst_100kb_prop_sites[which(twisst_100kb_prop_sites$scaffold!="Chr19"),]

length(which(chr19_100kb$top2==1))/nrow(chr19_100kb)
length(which(other_chr_100kb$top2==1))/nrow(other_chr_100kb)

####histgram plot of the proportion of incomplete lineage sorting (top1 and top3) across windows
ggplot(twisst_100kb_prop_sites,aes(1-top2))+geom_histogram(binwidth=0.01)+xlab("ILS")+
  geom_vline(aes(xintercept=mean(1-top2,na.rm=T)),color="dark red",linetype="dashed",size=0.8)
ggsave(filename="4species.twisst.ILS.100kb.png",width=5,height=5,dpi=300)

##top1
ggplot(twisst_100kb_prop_sites,aes(top1))+geom_histogram(binwidth=0.01)+xlab("Top1")+
  geom_vline(aes(xintercept=mean(top1,na.rm=T)),color="dark red",linetype="dashed",size=0.8)
ggsave(filename="4species.twisst.top1.100kb.png",width=5,height=5,dpi=300)
##top2
ggplot(twisst_100kb_prop_sites,aes(top2))+geom_histogram(binwidth=0.01)+xlab("Top2")+
  geom_vline(aes(xintercept=mean(top2,na.rm=T)),color="dark red",linetype="dashed",size=0.8)
ggsave(filename="4species.twisst.top2.100kb.png",width=5,height=5,dpi=300)
##top3
ggplot(twisst_100kb_prop_sites,aes(top3))+geom_histogram(binwidth=0.01)+xlab("Top3")+
  geom_vline(aes(xintercept=mean(top3,na.rm=T)),color="dark red",linetype="dashed",size=0.8)
ggsave(filename="4species.twisst.top3.100kb.png",width=5,height=5,dpi=300)

####plot top1,2,3 in a single plot

topology_100kb=twisst_100kb_prop_sites[,c("top1","top2","top3")]
names(topology_100kb)=c("Top1","Top2","Top3")
topology_100kb_melt=melt(topology_100kb)
names(topology_100kb_melt)=c("Topology","value")


###make the overall summary
colSums(twisst_10kb_new)/sum(colSums(twisst_10kb_new))
colSums(twisst_100kb_new)/sum(colSums(twisst_100kb_new))






####calculating the numbers of genomic windows supporting each of four topologies (Top1,Top2,Top3,unresolved trees)
length(which(topology_100kb_melt$Topology=="Top1" & topology_100kb_melt$value=="1"))
length(which(topology_100kb_melt$Topology=="Top2" & topology_100kb_melt$value=="1"))
length(which(topology_100kb_melt$Topology=="Top3" & topology_100kb_melt$value=="1"))

ggplot(topology_100kb_melt,aes(x=value,fill=Topology))+geom_histogram(binwidth=.03, position="dodge")+
  scale_fill_manual(values=c(colors[1], colors[7], colors[9]))+xlab("Proportion")+ylab("Count")
ggsave(filename="4species.twisst.topology.prop.100kb.pdf",width=6,height=5,dpi=300)


###Summary across chromosomes
mean_top1_100kb_twisst=aggregate(top1~scaffold,data=twisst_100kb_prop_sites,function(x){mean(x)})
mean_top2_100kb_twisst=aggregate(top2~scaffold,data=twisst_100kb_prop_sites,function(x){mean(x)})
mean_top3_100kb_twisst=aggregate(top3~scaffold,data=twisst_100kb_prop_sites,function(x){mean(x)})
mean_top_100kb_twisst=cbind(mean_top1_100kb_twisst,mean_top2_100kb_twisst[,"top2"],mean_top3_100kb_twisst[,"top3"])
names(mean_top_100kb_twisst)=c("scaffold","top1","top2","top3")
  
  

###add the sites information here
twisst_100kb_sites=fread("../twisst/data/4species.gatk.hap.snp.rm_indel.bed.biallelic.DP5.GQ10.rm_missing.recode.gt.beagle.raxml.w100kb.data.tsv")
twisst_100kb_sites_new=twisst_100kb_sites[which(twisst_100kb_sites$sites>200),]

hist(as.data.frame(twisst_100kb_prop)$top1,100)
hist(as.data.frame(twisst_10kb_prop)$top1,100)

twisst_100kb_prop_sites=cbind(twisst_100kb_sites_new,twisst_100kb_prop)



#######################################
###common windows between abbababa and twisst results
win3_100kb<-function(data){
  data$center=data$start+50000
  data$Win=paste(data$scaffold,data$center,sep=":")
  return(data)
}

twisst_100kb_prop_sites_win=win3_100kb(twisst_100kb_prop_sites)
abbababa_gatk_100kb_win=win3_100kb(abbababa_gatk_100kb_new)

common_win=Reduce(intersect,list(twisst_100kb_prop_sites_win$Win,
                                 abbababa_gatk_100kb_win$Win))

twisst_100kb_prop_sites_win=twisst_100kb_prop_sites_win[which(twisst_100kb_prop_sites_win$Win %in% common_win),]
abbababa_gatk_100kb_win=abbababa_gatk_100kb_win[which(abbababa_gatk_100kb_win$Win %in% common_win),]


###############Creating the sliding window plot for both D-statistic and Topology weighing test from Twisst########
abbababa_twisst_100kb=cbind(abbababa_gatk_100kb_win,twisst_100kb_prop_sites_win[,c("top1","top2","top3")])

####plot
colors <- brewer.pal(12,"Paired")[c(1,2,3,4,5,6,7,8,9,10,11,12)]

chr=as.list(matrix(,19))
chr_pos=as.list(matrix(,19))
chr_pos_new=as.list(matrix(,19))
D=as.list(matrix(,19))
fdM=as.list(matrix(,19))
top1=as.list(matrix(,19))
top2=as.list(matrix(,19))
top3=as.list(matrix(,19))

D_new=as.list(matrix(,19))
fdM_new=as.list(matrix(,19))
top1_new=as.list(matrix(,19))
top2_new=as.list(matrix(,19))
top3_new=as.list(matrix(,19))


for (i in 1:19) { 
  if (i<10) {
    chr[[i]]=paste("Chromosome ","0",i,sep="")
  }
  else {
    chr[[i]]=paste("Chromosome ",i,sep="")
  }
  
  if(i<10){
  chrom=paste("Chr0",i,sep="")
  }else{chrom=paste("Chr",i,sep="")}
  
  chr_pos[[i]]=abbababa_twisst_100kb$center[abbababa_twisst_100kb$scaffold==chrom]
  D[[i]]=abbababa_twisst_100kb$D[abbababa_twisst_100kb$scaffold==chrom]
  fdM[[i]]=abbababa_twisst_100kb$fdM[abbababa_twisst_100kb$scaffold==chrom]
  top1[[i]]=abbababa_twisst_100kb$top1[abbababa_twisst_100kb$scaffold==chrom]
  top2[[i]]=abbababa_twisst_100kb$top2[abbababa_twisst_100kb$scaffold==chrom]
  top3[[i]]=abbababa_twisst_100kb$top3[abbababa_twisst_100kb$scaffold==chrom]
  
  ###create new dataset to include all windows (not only just those windows with sites> N), and replace all missing dataset by NA
  start_pos=chr_pos[[i]][1]
  end_pos=chr_pos[[i]][length(chr_pos[[i]])]
  chr_pos_new[[i]]=seq(start_pos,end_pos,by=100000)
  new_dataset=data.frame(matrix(nrow=length(chr_pos_new[[i]]),ncol=6))
  names(new_dataset)=c("center","D","fdM","top1","top2","top3")
  new_dataset$center=chr_pos_new[[i]]
  
  new_dataset[which(new_dataset$center %in% chr_pos[[i]]),]$D=D[[i]]
  new_dataset[which(new_dataset$center %in% chr_pos[[i]]),]$fdM=fdM[[i]]
  new_dataset[which(new_dataset$center %in% chr_pos[[i]]),]$top1=top1[[i]]
  new_dataset[which(new_dataset$center %in% chr_pos[[i]]),]$top2=top2[[i]]
  new_dataset[which(new_dataset$center %in% chr_pos[[i]]),]$top3=top3[[i]]
  
  D_new[[i]]=-new_dataset$D
  fdM_new[[i]]=-new_dataset$fdM
  top1_new[[i]]=new_dataset$top1
  top2_new[[i]]=new_dataset$top2
  top3_new[[i]]=new_dataset$top3
}
  


##########plot the distribution of D-statistics and topology weightings 

pdf("twisst.D_statistic.4species.pdf",width=15,height=10)
  mat=matrix(c(rep(1,49),rep(2,25),rep(0,1),rep(3,21),rep(4,24),rep(5,26),rep(0,4),rep(6,28),rep(7,16),rep(8,20),rep(0,11),rep(9,13),rep(10,22),rep(11,14),rep(12,16),rep(0,10),rep(13,16),rep(14,13),rep(15,15),rep(16,15),rep(0,16),rep(17,16),rep(18,17),rep(19,16),rep(0,26)),byrow=T,nrow=6)
  layout(mat)

  
  for (i in 1:19) { 
    if(i<10){
      chrom=paste("Chr0",i,sep="")
    }else{chrom=paste("Chr",i,sep="")}
    
    if (i %in% c(1,3,6,9,13,17)){
    par(mar=c(1.5,4,0.5,1))
    plot(chr_pos_new[[i]]/1e6,D_new[[i]],axes=F,col="black",xlab="",ylab="",xaxt="n",type="n",ylim=c(-1,1))
    }
    else{ par(mar=c(1.5,1,0.5,1))
    plot(chr_pos_new[[i]]/1e6,D_new[[i]],axes=F,col="black",xlab="",ylab="",xaxt="n",type="n",ylim=c(-1,1))
    }
    
    ##add segments of topologies
    #top1
    segments(x0=c(chr_pos_new[[i]][which(top1_new[[i]]==1)]/1e6),x1=c(chr_pos_new[[i]][which(top1_new[[i]]==1)]/1e6),y0=c(rep(0.4,length(which(top1_new[[i]]==1)))),y1=c(rep(1,length(which(top1_new[[i]]==1)))),col=colors[1],lwd=1.5)
    #top2
    segments(x0=c(chr_pos_new[[i]][which(top2_new[[i]]==1)]/1e6),x1=c(chr_pos_new[[i]][which(top2_new[[i]]==1)]/1e6),y0=c(rep(-0.3,length(which(top2_new[[i]]==1)))),y1=c(rep(0.3,length(which(top2_new[[i]]==1)))),col=colors[7],lwd=1.5)
    #top3
    segments(x0=c(chr_pos_new[[i]][which(top3_new[[i]]==1)]/1e6),x1=c(chr_pos_new[[i]][which(top3_new[[i]]==1)]/1e6),y0=c(rep(-1,length(which(top3_new[[i]]==1)))),y1=c(rep(-0.4,length(which(top3_new[[i]]==1)))),col=colors[9],lwd=1.5)
    ##no clear topologies
    #segments(x0=c(chr_pos_new[[i]][which(top1_new[[i]]!=1 & top2_new[[i]]!=1 & top3_new[[i]]!=1 )]/1e6),
    #         x1=c(chr_pos_new[[i]][which(top1_new[[i]]!=1 & top2_new[[i]]!=1 & top3_new[[i]]!=1 )]/1e6),
    #         y0=c(rep(-1,length(which(top1_new[[i]]!=1 & top2_new[[i]]!=1 & top3_new[[i]]!=1 )))),
    #         y1=c(rep(1,length(which(top1_new[[i]]!=1 & top2_new[[i]]!=1 & top3_new[[i]]!=1 )))),
    #         col=colors[11])
    #add axis
    axis(labels=NA,side=1,tck=-0.015,pos=0)
    axis(side=2,at=seq(-1,1,0.5),cex.axis=0.7)
    mtext(chrom,side=1,cex=1,line=0.2)
    ##D-statistics
    lines(chr_pos_new[[i]]/1e6,D_new[[i]],axes=F,col="grey10",xlab="",ylab="",xaxt="n",ylim=c(-1,1))

   # if (i == 6){mtext("D",side=3,cex=1.1,line=-0.3,padj=6.3,adj=-0.1)}
}
    
dev.off()  

#########plot the distribution of fdM-statistics and topology weightings 

pdf("twisst.fdM_statistic.4species.pdf",width=15,height=10)
mat=matrix(c(rep(1,49),rep(2,25),rep(0,1),rep(3,21),rep(4,24),rep(5,26),rep(0,4),rep(6,28),rep(7,16),rep(8,20),rep(0,11),rep(9,13),rep(10,22),rep(11,14),rep(12,16),rep(0,10),rep(13,16),rep(14,13),rep(15,15),rep(16,15),rep(0,16),rep(17,16),rep(18,17),rep(19,16),rep(0,26)),byrow=T,nrow=6)
layout(mat)


for (i in 1:19) { 
  if(i<10){
    chrom=paste("Chr0",i,sep="")
  }else{chrom=paste("Chr",i,sep="")}
  
  if (i %in% c(1,3,6,9,13,17)){
    par(mar=c(1.5,4,0.5,1))
    plot(chr_pos_new[[i]]/1e6,fdM_new[[i]],axes=F,col="black",xlab="",ylab="",xaxt="n",type="n",ylim=c(-0.5,0.5))
  }
  else{ par(mar=c(1.5,1,0.5,1))
    plot(chr_pos_new[[i]]/1e6,fdM_new[[i]],axes=F,col="black",xlab="",ylab="",xaxt="n",type="n",ylim=c(-0.5,0.5))
  }
  
  ##add segments of topologies
  #top1
  segments(x0=c(chr_pos_new[[i]][which(top1_new[[i]]==1)]/1e6),x1=c(chr_pos_new[[i]][which(top1_new[[i]]==1)]/1e6),y0=c(rep(0.17,length(which(top1_new[[i]]==1)))),y1=c(rep(0.5,length(which(top1_new[[i]]==1)))),col=colors[1],lwd=1.5)
  #top2
  segments(x0=c(chr_pos_new[[i]][which(top2_new[[i]]==1)]/1e6),x1=c(chr_pos_new[[i]][which(top2_new[[i]]==1)]/1e6),y0=c(rep(-0.18,length(which(top2_new[[i]]==1)))),y1=c(rep(0.15,length(which(top2_new[[i]]==1)))),col=colors[7],lwd=1.5)
  #top3
  segments(x0=c(chr_pos_new[[i]][which(top3_new[[i]]==1)]/1e6),x1=c(chr_pos_new[[i]][which(top3_new[[i]]==1)]/1e6),y0=c(rep(-0.5,length(which(top3_new[[i]]==1)))),y1=c(rep(-0.17,length(which(top3_new[[i]]==1)))),col=colors[9],lwd=1.5)
  ##no clear topologies
  #segments(x0=c(chr_pos_new[[i]][which(top1_new[[i]]!=1 & top2_new[[i]]!=1 & top3_new[[i]]!=1 )]/1e6),
  #         x1=c(chr_pos_new[[i]][which(top1_new[[i]]!=1 & top2_new[[i]]!=1 & top3_new[[i]]!=1 )]/1e6),
  #         y0=c(rep(-0.5,length(which(top1_new[[i]]!=1 & top2_new[[i]]!=1 & top3_new[[i]]!=1 )))),
  #         y1=c(rep(0.5,length(which(top1_new[[i]]!=1 & top2_new[[i]]!=1 & top3_new[[i]]!=1 )))),
  #         col=colors[11])
  #add axis
  axis(labels=NA,side=1,tck=-0.015,pos=0)
  axis(side=2,at=seq(-0.5,0.5,0.25),cex.axis=0.7)
  mtext(chrom,side=1,cex=1,line=0.2)
  ##fdM-statistics
  lines(chr_pos_new[[i]]/1e6,fdM_new[[i]],axes=F,col="grey10",xlab="",ylab="",xaxt="n",ylim=c(-1,1))
  
  # if (i == 6){mtext("D",side=3,cex=1.1,line=-0.3,padj=6.3,adj=-0.1)}
}

dev.off()  



######branch length
setwd("../twisst/data/")
brach_length=read.table("/Users/Jing/Dropbox/davidiana_paper/data/twisst/data/4species.gatk.hap.snp.rm_indel.bed.biallelic.DP5.GQ10.rm_missing.recode.gt.beagle.twisst.distance.100kb.csv",header=T)

###brach length for P.tra-Ptrs
pdf("twisst.brach_length.100kb.ptra_ptrs.pdf",width=5,height=5)
par(mar=c(3,4,1,1))
boxplot(brach_length$topo1_P.tremula_P.tremuloides,brach_length$topo2_P.tremula_P.tremuloides,brach_length$topo3_P.tremula_P.tremuloides,
        col=c(colors[1],colors[7],colors[9]),ylab="Distance",xlab="",ylim=c(0.05,0.4),
        outline=FALSE,frame=F,axes=FALSE)
axis(2,at=c(0.1,0.2,0.3,0.4))
abline(h=0.1,lty=2,col="grey")
abline(h=0.15,lty=2,col="grey")
abline(h=0.2,lty=2,col="grey")
abline(h=0.25,lty=2,col="grey")
abline(h=0.3,lty=2,col="grey")
abline(h=0.35,lty=2,col="grey")
abline(h=0.4,lty=2,col="grey")
dev.off()


###brach length for Ptrs-Pdav
pdf("twisst.brach_length.100kb.ptrs_pdav.pdf",width=5,height=5)
par(mar=c(3,4,1,1))
boxplot(brach_length$topo1_P.tremuloides_P.davidiana,brach_length$topo2_P.tremuloides_P.davidiana,brach_length$topo3_P.tremuloides_P.davidiana,
        col=c(colors[1],colors[7],colors[9]),ylab="Distance",xlab="",ylim=c(0.05,0.4),
        outline=FALSE,frame=F,axes=FALSE)
axis(2,at=c(0.1,0.2,0.3,0.4))
abline(h=0.1,lty=2,col="grey")
abline(h=0.15,lty=2,col="grey")
abline(h=0.2,lty=2,col="grey")
abline(h=0.25,lty=2,col="grey")
abline(h=0.3,lty=2,col="grey")
abline(h=0.35,lty=2,col="grey")
abline(h=0.4,lty=2,col="grey")
dev.off()




###brach length for P.tra-Pdav
pdf("twisst.brach_length.100kb.ptra_pdav.pdf",width=5,height=5)
par(mar=c(3,4,1,1))
boxplot(brach_length$topo1_P.tremula_P.davidiana,brach_length$topo2_P.tremula_P.davidiana,brach_length$topo3_P.tremula_P.davidiana,
        col=c(colors[1],colors[7],colors[9]),ylab="Distance",xlab="",ylim=c(0.05,0.4),
        outline=FALSE,frame=F,axes=FALSE)
axis(2,at=c(0.1,0.2,0.3,0.4))
abline(h=0.1,lty=2,col="grey")
abline(h=0.15,lty=2,col="grey")
abline(h=0.2,lty=2,col="grey")
abline(h=0.25,lty=2,col="grey")
abline(h=0.3,lty=2,col="grey")
abline(h=0.35,lty=2,col="grey")
abline(h=0.4,lty=2,col="grey")
dev.off()


###2. density plot of brach length for P.tra-Ptrs
pdf("topo.density.pdf",width=7,height=5)
tra_trs=brach_length[,c("topo1_P.tremula_P.tremuloides","topo2_P.tremula_P.davidiana","topo3_P.tremuloides_P.davidiana")]
tra_trs_melt=melt(tra_trs)
ggplot(tra_trs_melt, aes(x=value, group=variable,fill=variable)) + 
  geom_density(alpha=.3) + 
  xlab("Topology") +
  ylab("Density")
dev.off()

###brach length 10kb
brach_length_10kb=read.table("/Users/Jing/Dropbox/davidiana_paper/data/twisst/data/4species.gatk.hap.snp.rm_indel.bed.biallelic.DP5.GQ10.rm_missing.recode.gt.beagle.twisst.distance.10kb.csv",header=T)

###brach length for P.tra-Ptrs
pdf("twisst.brach_length.10kb.ptra_ptrs.pdf",width=5,height=5)
par(mar=c(3,4,1,1))
boxplot(brach_length_10kb$topo1_P.tremula_P.tremuloides,brach_length_10kb$topo2_P.tremula_P.tremuloides,brach_length_10kb$topo3_P.tremula_P.tremuloides,
        col=c(colors[1],colors[7],colors[9]),ylab="Distance",xlab="",ylim=c(0.00,0.4),
        outline=FALSE,frame=F,axes=FALSE)
axis(2,at=c(0.0,0.1,0.2,0.3,0.4))
abline(h=0.05,lty=2,col="grey")
abline(h=0.1,lty=2,col="grey")
abline(h=0.15,lty=2,col="grey")
abline(h=0.2,lty=2,col="grey")
abline(h=0.25,lty=2,col="grey")
abline(h=0.3,lty=2,col="grey")
abline(h=0.35,lty=2,col="grey")
abline(h=0.4,lty=2,col="grey")
dev.off()

###brach length for Ptrs-Pdav
pdf("twisst.brach_length.10kb.ptrs_pdav.pdf",width=5,height=5)
par(mar=c(3,4,1,1))
boxplot(brach_length_10kb$topo1_P.tremuloides_P.davidiana,brach_length_10kb$topo2_P.tremuloides_P.davidiana,brach_length_10kb$topo3_P.tremuloides_P.davidiana,
        col=c(colors[1],colors[7],colors[9]),ylab="Distance",xlab="",ylim=c(0.0,0.4),
        outline=FALSE,frame=F,axes=FALSE)
axis(2,at=c(0,0.1,0.2,0.3,0.4))
abline(h=0.05,lty=2,col="grey")
abline(h=0.1,lty=2,col="grey")
abline(h=0.15,lty=2,col="grey")
abline(h=0.2,lty=2,col="grey")
abline(h=0.25,lty=2,col="grey")
abline(h=0.3,lty=2,col="grey")
abline(h=0.35,lty=2,col="grey")
abline(h=0.4,lty=2,col="grey")
dev.off()

###brach length for P.tra-Pdav
pdf("twisst.brach_length.10kb.ptra_pdav.pdf",width=5,height=5)
par(mar=c(3,4,1,1))
boxplot(brach_length_10kb$topo1_P.tremula_P.davidiana,brach_length_10kb$topo2_P.tremula_P.davidiana,brach_length_10kb$topo3_P.tremula_P.davidiana,
        col=c(colors[1],colors[7],colors[9]),ylab="Distance",xlab="",ylim=c(0.0,0.4),
        outline=FALSE,frame=F,axes=FALSE)
axis(2,at=c(0,0.1,0.2,0.3,0.4))
abline(h=0.05,lty=2,col="grey")
abline(h=0.1,lty=2,col="grey")
abline(h=0.15,lty=2,col="grey")
abline(h=0.2,lty=2,col="grey")
abline(h=0.25,lty=2,col="grey")
abline(h=0.3,lty=2,col="grey")
abline(h=0.35,lty=2,col="grey")
abline(h=0.4,lty=2,col="grey")
dev.off()




