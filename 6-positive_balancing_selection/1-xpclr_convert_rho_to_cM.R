#! /usr/bin/Rscript --no-save --no-restore

library(data.table)

setwd("/proj/sllstore2017050/nobackup/milou_files/PaperIV-phylogenomics/GATK/HC/total3/XPCLR/recombination_cM")

#Size of sliding windows used in ld_helmet
win_size<-1e4

##Total chromosome lengths in map units (cM), taken from our recent tremula map
chrMapUnits<-c(498.77,266.10,235.24,219.79,271.90,296.41,154.65,214.53,172.22,273.71,182.01,143.34,178.67,194.34,159.22,145.19,149.44,165.47,135.32)

##Read ld_helmet results
rec<-fread("gunzip -c /proj/sllstore2017050/nobackup/milou_files/PaperIV-phylogenomics/recombination_rates/ldhelmet/tremula/out/tremula.window10000.ldhelmet.summary.txt.gz")
rec$chr<-as.numeric(substr(rec$Chr,4,5))

##Assign rho=0 to all windows with no SNPs
rec[is.na(rec$rho_mean),]$rho_mean<-0

##Calculate cumulative sum of rho/window
rec$rho_cum=0
for(i in 1:19){
rec[which(rec$chr==i),]$rho_cum<-cumsum(rec[which(rec$chr==i),]$rho_mean*win_size)
}
#Get chromosome lenths in rho
chrRhoUnits<-tapply(rec$rho_cum,rec$chr,max)

##Calculate conversion factor for the two chromsome scales (rho and cM)
scale<-chrMapUnits/chrRhoUnits

##Calculate cM value for all windows across the 19 chromosomes
rec$cM_cum<-0
for(i in 1:19){
  rec[rec$chr==i,]$cM_cum<-rec[rec$chr==i,]$rho_cum*scale[i]
}

##create the genetic position file used for XPCLR
rec$gene_pos=rec$cM_cum/100

#write.table
write.table(rec,file="tremula.ldhelmet_10kb_cM.txt",sep="\t", quote=F, row.names=F, col.names=T)

