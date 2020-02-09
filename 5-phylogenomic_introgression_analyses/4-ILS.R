library(data.table)
library(RColorBrewer)
library(dplyr)
library(reshape2)
library(ggplot2)
#install.packages("Hmisc")
library(Hmisc)

setwd("~/Dropbox/davidiana_paper/data/snp_allele_freq")


#########################################################################
###the allele frequency are non-reference allele frequency, when treating P.trichocarpa as ancestral species, the allele on the other three aspen species can be used as "derived" allele.
allele_frq=fread('cat species.frq.gz |gunzip',header=T)

###The possible phylogenetic relationship assuming bi-allelic SNPs
#(((P.tra,Pdav),P.trs),P.tri)
#1.BBBA
#2.AAAB
#3.BBAA
#4.BABA
#5.ABBA
#6.BAAA
#7.ABAA
#8.AABA
###the tree 3-8 are the tree under our focus


###remove cases of AAAB and BBBA
allele_frq_aspen=allele_frq %>% filter((tremula_frq+davidiana_frq+tremuloides_frq)>0) %>% filter((tremula_frq+davidiana_frq+tremuloides_frq)<3)
###Then, I also found that for cases that at least one of the three aspen species shared ancestral allele with P.trichocarpa, there is always something wrong with the computation of incomplete lineage sorting, so I removed these sites


allele_frq_aspen$bbaa=round(allele_frq_aspen$tremula_frq*allele_frq_aspen$davidiana_frq*(1-allele_frq_aspen$tremuloides_frq)*(1-allele_frq_aspen$trichocarpa_frq),digits=6)
allele_frq_aspen$baba=round(allele_frq_aspen$tremula_frq*(1-allele_frq_aspen$davidiana_frq)*allele_frq_aspen$tremuloides_frq*(1-allele_frq_aspen$trichocarpa_frq),digits=6)
allele_frq_aspen$abba=round((1-allele_frq_aspen$tremula_frq)*allele_frq_aspen$davidiana_frq*allele_frq_aspen$tremuloides_frq*(1-allele_frq_aspen$trichocarpa_frq),digits=6)
allele_frq_aspen$baaa=round(allele_frq_aspen$tremula_frq*(1-allele_frq_aspen$davidiana_frq)*(1-allele_frq_aspen$tremuloides_frq)*(1-allele_frq_aspen$trichocarpa_frq),digits=6)
allele_frq_aspen$abaa=round((1-allele_frq_aspen$tremula_frq)*allele_frq_aspen$davidiana_frq*(1-allele_frq_aspen$tremuloides_frq)*(1-allele_frq_aspen$trichocarpa_frq),digits=6)
allele_frq_aspen$aaba=round((1-allele_frq_aspen$tremula_frq)*(1-allele_frq_aspen$davidiana_frq)*allele_frq_aspen$tremuloides_frq*(1-allele_frq_aspen$trichocarpa_frq),digits=6)

##calculate the incomplete lineage sorting
allele_frq_aspen$h=round((allele_frq_aspen$baaa+allele_frq_aspen$aaba+allele_frq_aspen$abaa+(2*(allele_frq_aspen$bbaa+allele_frq_aspen$baba+allele_frq_aspen$abba)))/3,digits=6)
allele_frq_aspen$nls=round((allele_frq_aspen$baba+allele_frq_aspen$abba)/allele_frq_aspen$h,digits=6)


#########################################################################
###summarize across windows
sum_win=function(data,window,step){
  data_out=c()
  for (chr in 1:19){
    if(chr<10){chrom=paste("Chr0",chr,sep="")}else{chrom=paste("Chr",chr,sep="")}
    print(chrom)
    data_chr=data %>% filter(CHROM==chrom)
    
    end_pos=data_chr[nrow(data_chr),]$POS
    idx1<-seq(0,end_pos,by=step);
    idx1+window->idx2;
   
###define new columns for the table    
    
    chro=c()
    start=c()
    end=c()
    n_sites=c()
    bbaa_n=c()
    baba_n=c()
    abba_n=c()
    baaa_n=c()
    aaba_n=c()
    abaa_n=c()
    nls_n=c()
    
    
    for (i in 1:length(idx1)){
      chro[i]=chrom
      start[i]=idx1[i]
      end[i]=idx2[i]
      n_sites[i]=length(which(data_chr$POS<=idx2[i] & data_chr$POS>idx1[i]))
      
      bbaa_n[i]=mean(data_chr[which(data_chr$POS<=idx2[i] & data_chr$POS>idx1[i]),]$bbaa,na.rm=T)
      baba_n[i]=mean(data_chr[which(data_chr$POS<=idx2[i] & data_chr$POS>idx1[i]),]$baba,na.rm=T)
      abba_n[i]=mean(data_chr[which(data_chr$POS<=idx2[i] & data_chr$POS>idx1[i]),]$abba,na.rm=T)
      baaa_n[i]=mean(data_chr[which(data_chr$POS<=idx2[i] & data_chr$POS>idx1[i]),]$baaa,na.rm=T)
      aaba_n[i]=mean(data_chr[which(data_chr$POS<=idx2[i] & data_chr$POS>idx1[i]),]$aaba,na.rm=T)
      abaa_n[i]=mean(data_chr[which(data_chr$POS<=idx2[i] & data_chr$POS>idx1[i]),]$abaa,na.rm=T)
      nls_n[i]=mean(data_chr[which(data_chr$POS<=idx2[i] & data_chr$POS>idx1[i]),]$nls,na.rm=T)
    }
    data_chr_new=data.frame(cbind(chro,start,end,n_sites,bbaa_n,baba_n,abba_n,baaa_n,aaba_n,abaa_n,nls_n))
    data_out=rbind(data_out,data_chr_new)
  }
  return(data_out)
}

allele_frq_aspen_100kb=sum_win(allele_frq_aspen,100000,100000)
allele_frq_aspen_10kb=sum_win(allele_frq_aspen,10000,10000)

###write out the table
write.table(allele_frq_aspen_100kb,file="allele_frq_aspen_100kb.nls.txt",sep="\t", quote=F, row.names=F, col.names=T)
write.table(allele_frq_aspen_10kb,file="allele_frq_aspen_10kb.nls.txt",sep="\t", quote=F, row.names=F, col.names=T)


colors <- brewer.pal(12,"Paired")

#########################################################################
###plot the pattern of incomplete lineage sorting along with the distance of genes
allele_frq_aspen$SNP=paste(allele_frq_aspen$CHROM,allele_frq_aspen$POS,sep=":")


snp_gene_distance=fread('cat 4species.gatk.hap.snp.rm_indel.bed.biallelic.DP5.GQ10.rm_missing.recode.gene.distance.txt.gz |gunzip',header=F)
names(snp_gene_distance)=c("Chr","Pos","Gene","Distance")
snp_gene_distance$SNP=paste(snp_gene_distance$Chr,snp_gene_distance$Pos,sep=":")
snp_gene_distance_new=snp_gene_distance[which(snp_gene_distance$SNP %in% allele_frq_aspen$SNP),]




dis_0=allele_frq_aspen[which(abs(snp_gene_distance_new$Distance)=="0"),]
gene_larger_0=allele_frq_aspen[which(abs(snp_gene_distance_new$Distance)>0),]
gene_larger_0$dis=snp_gene_distance_new[which(abs(snp_gene_distance_new$Distance)>0),]$Distance

gene_larger_0$group <- as.numeric(cut2(abs(gene_larger_0$dis), g=20))

by(gene_larger_0$nls, gene_larger_0$group, function(x) mean(x,na.rm=T))

gene_larger_0_nls_g=as.vector(by(gene_larger_0$nls, gene_larger_0$group, function(x) mean(x,na.rm=T)))
gene_larger_0_dis_g=as.vector(by(abs(gene_larger_0$dis), gene_larger_0$group, function(x) mean(x,na.rm=T)))
gene_larger_0_sd_g=as.vector(by(gene_larger_0$nls, gene_larger_0$group, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

###summary into a table
gene_larger_0.df=data.frame(mean=gene_larger_0_nls_g,std=gene_larger_0_sd_g,dis=gene_larger_0_dis_g)


####making the plot
ggplot(gene_larger_0.df,aes(dis,y=mean,ymin=mean-std,ymax=mean+std))+
  geom_pointrange(color=colors[4])+xlab("Physical distance from gene")+ylab("Scaled rate of ILS sites")+
  theme_bw()+theme(plot.margin = margin(1, 1, 0.5, 0.5, "cm"),axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))
ggsave("Dis_gene.ILS.pdf",width=6,height=5)



#########################################################################
####exon

snp_exon_distance=fread('cat 4species.gatk.hap.snp.rm_indel.bed.biallelic.DP5.GQ10.rm_missing.recode.exon.distance.txt.gz |gunzip',header=F)
names(snp_exon_distance)=c("Chr","Pos","Distance")
snp_exon_distance$SNP=paste(snp_exon_distance$Chr,snp_exon_distance$Pos,sep=":")
snp_exon_distance_new=snp_exon_distance[which(snp_exon_distance$SNP %in% allele_frq_aspen$SNP),]

exon_dis_0=allele_frq_aspen[which(abs(snp_exon_distance_new$Distance)=="0"),]
exon_larger_0=allele_frq_aspen[which(abs(snp_exon_distance_new$Distance)>0),]
exon_larger_0$dis=snp_exon_distance_new[which(abs(snp_exon_distance_new$Distance)>0),]$Distance

exon_larger_0$group <- as.numeric(cut2(abs(exon_larger_0$dis), g=20))

exon_larger_0_nls_g=as.vector(by(exon_larger_0$nls, exon_larger_0$group, function(x) mean(x,na.rm=T)))
exon_larger_0_dis_g=as.vector(by(abs(exon_larger_0$dis), exon_larger_0$group, function(x) mean(x,na.rm=T)))
exon_larger_0_sd_g=as.vector(by(exon_larger_0$nls, exon_larger_0$group, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

###summary into a table
exon_larger_0.df=data.frame(mean=exon_larger_0_nls_g,std=exon_larger_0_sd_g,dis=exon_larger_0_dis_g)


####making the plot
ggplot(exon_larger_0.df,aes(dis,y=mean,ymin=mean-std,ymax=mean+std))+
  geom_pointrange(color=colors[2])+xlab("Physical distance from exon")+ylab("Scaled rate of ILS sites")+
  theme_bw()+theme(plot.margin = margin(1, 1, 0.5, 0.5, "cm"),axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))
ggsave("Dis_exon.ILS.pdf",width=6,height=5)


#########################################################################
####cds

snp_cds_distance=fread('cat 4species.gatk.hap.snp.rm_indel.bed.biallelic.DP5.GQ10.rm_missing.recode.cds.distance.txt.gz |gunzip',header=F)
names(snp_cds_distance)=c("Chr","Pos","Distance")
snp_cds_distance$SNP=paste(snp_cds_distance$Chr,snp_cds_distance$Pos,sep=":")
snp_cds_distance_new=snp_cds_distance[which(snp_cds_distance$SNP %in% allele_frq_aspen$SNP),]

cds_dis_0=allele_frq_aspen[which(abs(snp_cds_distance_new$Distance)=="0"),]
cds_larger_0=allele_frq_aspen[which(abs(snp_cds_distance_new$Distance)>0),]
cds_larger_0$dis=snp_cds_distance_new[which(abs(snp_cds_distance_new$Distance)>0),]$Distance

cds_larger_0$group <- as.numeric(cut2(abs(cds_larger_0$dis), g=20))

cds_larger_0_nls_g=as.vector(by(cds_larger_0$nls, cds_larger_0$group, function(x) mean(x,na.rm=T)))
cds_larger_0_dis_g=as.vector(by(abs(cds_larger_0$dis), cds_larger_0$group, function(x) mean(x,na.rm=T)))
cds_larger_0_sd_g=as.vector(by(cds_larger_0$nls, cds_larger_0$group, function(x) 1.96*(sd(x,na.rm=T)/sqrt(length(x)))))

###summary into a table
cds_larger_0.df=data.frame(mean=cds_larger_0_nls_g,std=cds_larger_0_sd_g,dis=cds_larger_0_dis_g)


colors <- brewer.pal(12,"Paired")

####making the plot
ggplot(cds_larger_0.df,aes(dis,y=mean,ymin=mean-std,ymax=mean+std))+
  geom_pointrange(color=colors[6])+xlab("Physical distance from coding sites")+ylab("Scaled rate of ILS sites")+
  theme_bw()+theme(plot.margin = margin(1, 1, 0.5, 0.5, "cm"),axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))
ggsave("Dis_cds.ILS.pdf",width=6,height=5)





