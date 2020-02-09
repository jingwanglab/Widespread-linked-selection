#! /usr/bin/Rscript --no-save --no-restore


setwd("/proj/sllstore2017050/nobackup/milou_files/PaperIV-phylogenomics/demographic_inference/fastsimcoal2/input/angsd/chr")

###1.tra_dav
tra_dav=read.table("tremula.davidiana.Chr01.2dsfs.sfs",header=F)
###read in all other sfs datasets and add into the chr01
for (chr in 2:19){
	if (chr <10){
		chrom=paste("0",chr,sep="")
	}else{  chrom=chr}
	tra_dav_chr=paste("tremula.davidiana.Chr",chrom,".2dsfs.sfs",sep="")
	tra_dav_chr_data=read.table(tra_dav_chr,header=F)
	tra_dav=tra_dav+tra_dav_chr_data
} 
tra_dav_m=matrix(tra_dav,ncol=17)
##write the matrix into the table
write.table(tra_dav_m, file="2d.ptra_pdav.matrix.txt", sep="\t", quote=F, row.names=F, col.names=F)
#write the fastsimcoal input into the file
colnames(tra_dav_m)=paste("d0_",0:16,sep="") ###colnames corresponds to d0:ptra
rownames(tra_dav_m)=paste("d1_",0:16,sep="") ###rownames corresponds to d1:pdav
cat("1 observation\n", file="2d.ptra_pdav.pop0_1.obs")
write.table(tra_dav_m, file="2d.ptra_pdav.pop0_1.obs", sep="\t", quote=F, row.names=T, col.names=T,append=T)

#2.tra_trs
tra_trs=read.table("tremula.tremuloides.Chr01.2dsfs.sfs",header=F)
###read in all other sfs datasets and add into the chr01
for (chr in 2:19){
        if (chr<10){
                chrom=paste("0",chr,sep="")
        }else{  chrom=chr}
        tra_trs_chr=paste("tremula.tremuloides.Chr",chrom,".2dsfs.sfs",sep="")
        tra_trs_chr_data=read.table(tra_trs_chr,header=F)
        tra_trs=tra_trs+tra_trs_chr_data
}

tra_trs_m=matrix(tra_trs,ncol=17)
##write the matrix into the table
write.table(tra_trs_m, file="2d.ptra_ptrs.matrix.txt", sep="\t", quote=F, row.names=F, col.names=F)
#write the fastsimcoal input into the file
colnames(tra_trs_m)=paste("d0_",0:16,sep="") ###colnames corresponds to d0:ptra
rownames(tra_trs_m)=paste("d2_",0:16,sep="") ###rownames corresponds to d1:pdav
cat("1 observation\n", file="2d.ptra_ptrs.pop0_2.obs")
write.table(tra_trs_m, file="2d.ptra_ptrs.pop0_2.obs", sep="\t", quote=F, row.names=T, col.names=T,append=T)


dav_trs=read.table("davidiana.tremuloides.Chr01.2dsfs.sfs",header=F)
for (chr in 2:19){
        if (chr<10){
                chrom=paste("0",chr,sep="")
        }else{  chrom=chr}
        dav_trs_chr=paste("davidiana.tremuloides.Chr",chrom,".2dsfs.sfs",sep="")
        dav_trs_chr_data=read.table(dav_trs_chr,header=F)
        dav_trs=dav_trs+dav_trs_chr_data
}
dav_trs_m=matrix(dav_trs,ncol=17)
##write the matrix into the table
write.table(dav_trs_m, file="2d.pdav_ptrs.matrix.txt", sep="\t", quote=F, row.names=F, col.names=F)
#write the fastsimcoal input into the file
colnames(dav_trs_m)=paste("d1_",0:16,sep="") ###colnames corresponds to d0:ptra
rownames(dav_trs_m)=paste("d2_",0:16,sep="") ###rownames corresponds to d1:pdav
cat("1 observation\n", file="2d.pdav_ptrs.pop1_2.obs")
write.table(dav_trs_m, file="2d.pdav_ptrs.pop1_2.obs", sep="\t", quote=F, row.names=T, col.names=T,append=T)




