library(data.table)
library(RColorBrewer)
#install.packages("gplots")
library(gplots)
library(lattice)
library(ggplot2)


setwd("~/Dropbox/davidiana_paper/data/beagle_ibd")
sample=read.table("4species.samples.txt",header=F)

ibd=fread('cat 4species.gatk.hap.snp.rm_indel.bed.biallelic.DP5.GQ10.rm_missing.recode.gt.beagle.ibd.all.ibd.gz | gunzip')
#ibd=ibd_old[which(ibd_old$V8>0),]
hbd=fread('cat 4species.gatk.hap.snp.rm_indel.bed.biallelic.DP5.GQ10.rm_missing.recode.gt.beagle.ibd.all.hbd.gz | gunzip')
#hbd=hbd_old[which(hbd_old$V8>0),]

ibd$length=ibd$V7-ibd$V6

###etimate the length of shared haplotypes
###potra
potra=as.character(sample[which(sample$V1=="P.tremula"),]$V2)
potra_ibd=ibd[which(ibd$V1 %in% potra & ibd$V3 %in% potra),]
potra_ibd$species="Ptra"
summary(potra_ibd$length)
nrow(potra_ibd)

###podav
podav=as.character(sample[which(sample$V1=="P.davidiana"),]$V2)
podav_ibd=ibd[which(ibd$V1 %in% podav & ibd$V3 %in% podav),]
podav_ibd$species="Pdav"
summary(podav_ibd$length)
nrow(podav_ibd)

###potrs
potrs=as.character(sample[which(sample$V1=="P.tremuloides"),]$V2)
potrs_ibd=ibd[which(ibd$V1 %in% potrs & ibd$V3 %in% potrs),]
potrs_ibd$species="Ptrs"
summary(potrs_ibd$length)
nrow(potrs_ibd)

###potri
potri=as.character(sample[which(sample$V1=="P.trichocarpa"),]$V2)
potri_ibd=ibd[which(ibd$V1 %in% potri & ibd$V3 %in% potri),]
potri_ibd$species="Ptri"
summary(potri_ibd$length)
nrow(potri_ibd)

#####Summarize the length of shared IBD haplotypes within species
ibd_species=rbind(potra_ibd,podav_ibd,potrs_ibd,potri_ibd)
ibd_species$species=factor(ibd_species$species,levels=c("Ptra","Pdav","Ptrs","Ptri"))


pdf("ibd.length.4species.pdf",width=6,height=6)
ggplot(ibd_species,aes(x=length,fill=species,color=species))+
  geom_histogram(position="identity",alpha=0.2)+
  scale_x_continuous(limits = c(0, 100000))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  xlab("IBD block length, bp")+ylab("Count")
dev.off()

##potra-podav
potra_podav_ibd=ibd[which((ibd$V1 %in% potra & ibd$V3 %in% podav) | (ibd$V1 %in% podav & ibd$V3 %in% potra)),]
summary(potra_podav_ibd$length)
nrow(potra_podav_ibd)
potra_podav_ibd$species="Ptra-Pdav"
  
##potra-potrs
potra_potrs_ibd=ibd[which((ibd$V1 %in% potra & ibd$V3 %in% potrs) | (ibd$V1 %in% potrs & ibd$V3 %in% potra)),]
summary(potra_potrs_ibd$length)
nrow(potra_potrs_ibd)
potra_potrs_ibd$species="Ptra-Ptrs"

##podav-potrs
podav_potrs_ibd=ibd[which((ibd$V1 %in% podav & ibd$V3 %in% potrs) | (ibd$V1 %in% potrs & ibd$V3 %in% podav)),]
summary(podav_potrs_ibd$length)
nrow(podav_potrs_ibd)
podav_potrs_ibd$species="Pdav-Ptrs"

#####
#####Summarize the length of shared IBD haplotypes between species
ibd_between_species=rbind(potra_podav_ibd,potra_potrs_ibd,podav_potrs_ibd)
ibd_between_species$species=factor(ibd_between_species$species,levels=c("Ptra-Pdav","Ptra-Ptrs","Pdav-Ptrs"))


pdf("ibd.length.between.species.pdf",width=6,height=6)
ggplot(ibd_between_species,aes(x=length,fill=species,color=species))+
  geom_histogram(position="identity",alpha=0.2)+
  scale_x_continuous(limits = c(0, 100000))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  xlab("IBD block length, bp")+ylab("Count")
dev.off()



##potra-potri
potra_potri_ibd=ibd[which((ibd$V1 %in% potra & ibd$V3 %in% potri) | (ibd$V1 %in% potri & ibd$V3 %in% potra)),]
summary(potra_potri_ibd$length)

##potra-potri
potrs_potri_ibd=ibd[which((ibd$V1 %in% potrs & ibd$V3 %in% potri) | (ibd$V1 %in% potri & ibd$V3 %in% potrs)),]
summary(potrs_potri_ibd$length)

##podav-potri
podav_potri_ibd=ibd[which((ibd$V1 %in% podav & ibd$V3 %in% potri) | (ibd$V1 %in% potri & ibd$V3 %in% podav)),]
summary(podav_potri_ibd$length)


