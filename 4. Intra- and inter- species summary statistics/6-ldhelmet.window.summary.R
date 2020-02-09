#! /usr/bin/Rscript --no-save --no-restore
.libPaths("/pica/h1/jingwang/R/x86_64-redhat-linux-gnu-library/3.2")
library(data.table)
library(dplyr)

args=(commandArgs(TRUE))
setwd(args[1])

species=args[2]
ldhelmet_file=paste(species,".ldhelmet.txt",sep="")
window=as.numeric(args[3])
step=as.numeric(args[4])

outfile=paste(species,".window",window,".ldhelmet.summary.txt",sep="")


ldhelmet <- fread(ldhelmet_file,header=F)
names(ldhelmet)=c("chrom","left_snp","right_snp","mean","p0.025","p0.500","p0.098")

for (chr in 1:19){
	if (chr<10) {chromosome=paste("Chr0",chr,sep="")}else{chromosome=paste("Chr",chr,sep="")}
	ldhelmet_chr=filter(ldhelmet,chrom==chromosome)
	start=as.integer(min(ldhelmet_chr$left_snp)/window,0)
	if(start>=1){start_real=start*window}else{start_real=start}  
	end=max(ldhelmet_chr$right_snp)
	start_chr=seq(start_real,end,step)  ##a set of start position
	end_chr=start_chr+window-1   ##a set of end position
	#if(start<1){wpos=start_chr+(window/2)}else{wpos=start_chr+(window/2)-1}   ##a set of middle position
	wpos=start_chr+(window/2)	

	snp_n=c()  ##counting the number of SNPs in each window, which could be used for later filtering
	rho_mean=c()  ##the average rho value across the SNPs in the window
	rho_p0_025=c()  ##the 2.5 percentile of the rho value
	rho_p0_5=c()   ##the 50th percentile of the rho value
	rho_p0_975=c()   ##the 97.5th percentile of the rho value
	
	for (i in 1:length(start_chr)){
		ldhelmet_chr_win=filter(ldhelmet_chr,left_snp>=start_chr[i],right_snp<=end_chr[i])
		snp_n[i]=nrow(ldhelmet_chr_win)
		rho_mean[i]=mean(ldhelmet_chr_win$mean)
		rho_p0_025[i]=mean(ldhelmet_chr_win$p0.025)
		rho_p0_5[i]=mean(ldhelmet_chr_win$p0.500)
		rho_p0_975[i]=mean(ldhelmet_chr_win$p0.098)
	}
	rho_summary=data.frame(cbind(Chr=rep(chromosome,length(start_chr)),Pos=wpos,N_SNP=snp_n,rho_mean=rho_mean,rho_p0.025=rho_p0_025,rho_p0.5=rho_p0_5,rho_p0.975=rho_p0_975))
	if (chr == 1){write.table(rho_summary,file=outfile,append=T,sep="\t", quote=F, row.names=F, col.names=T)}else{
	write.table(rho_summary,file=outfile,append=T,sep="\t", quote=F, row.names=F, col.names=F)}
}

