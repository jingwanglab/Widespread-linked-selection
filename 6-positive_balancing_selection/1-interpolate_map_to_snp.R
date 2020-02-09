#! /usr/bin/Rscript --no-save --no-restore

.libPaths("/home/jingwang/R/x86_64-redhat-linux-gnu-library/3.6")

###This script is used to interpolate the genetic map position to snp position

library("data.table")

###read the SNP file that needs to transfer into genetic map position
args=(commandArgs(TRUE))
snp=args[1]   ##snp position Chr01.tra_dav.dav.snp

mis.sub<-fread(snp,header=F)
names(mis.sub)<-c("snp","chr","gen_pos","pos","reference","alternative")

map <- fread("/proj/sllstore2017050/nobackup/milou_files/PaperIV-phylogenomics/GATK/HC/total3/XPCLR/recombination_cM/tremula.ldhelmet_10kb_cM.txt",header=T)
names(map) <- c("Chr","pos","n_snp","rho_mean","rho_p0.025","rho_p0.5","rho_p0.975","chr","rho_cum","cM","gene_pos")

tri_bed=read.table("/proj/sllstore2017050/nobackup/milou_files/PaperIV-phylogenomics/bed/Potri_bed/trichocarpa.bed",header=F)
chrLengths <- tri_bed$V3[1:19]

rows<-dim(mis.sub)[1]

for(i in 1:nrow(mis.sub)){
  cat( 100*i/rows, "% done", "\r")
  lowerIndex <- which(map$chr == mis.sub$chr[i] & map$pos <= mis.sub$pos[i]) #find index of map anchors smaller than observed position
  belowPhys <- max(map$pos[lowerIndex[length(lowerIndex)]],1,na.rm=T) #take largest position of anchor that is smaller than observed position
  belowGen <- max(map$cM[lowerIndex[length(lowerIndex)]],map$cM[which(map$chr==mis.sub$chr[i])][1],na.rm=T) #take corresponding genetic position

  higherIndex <- which(map$chr == mis.sub$chr[i] & map$pos >= mis.sub$pos[i]) #find index of map anchors larger than observed position
  abovePhys <- min(map$pos[higherIndex[1]],chrLengths[mis.sub$chr[i]],na.rm=T) #take smallest position of anchor that is larger than observed position
  aboveGen <- min(map$cM[higherIndex[1]],map$cM[which(map$chr==mis.sub$chr[i])][length(which(map$chr==mis.sub$chr[i]))]+1,na.rm=T) #take corresponding genetic position

  scale <- {mis.sub$pos[i]-belowPhys}/{abovePhys-belowPhys} #compute linear scale for position of observed relative to anchors
  newGen <- {aboveGen-belowGen}*scale + belowGen # compute genetic position for observed position

  mis.sub$cM[i] <- newGen
  mis.sub$gen_pos[i]=newGen/100 ###assign to the new genetic position to the mis.sub position
#write.table
write.table(mis.sub[i,1:6],file=paste(snp,".new",sep=""),sep="\t",quote=F,row.names=F,col.names=F,append=T)
}




