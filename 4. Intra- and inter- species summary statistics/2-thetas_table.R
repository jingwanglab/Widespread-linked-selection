#! /usr/bin/Rscript --no-save --no-restore

# File: plotThetas.R, the following is one example of how to use this script
#for file in {01..19}; do Rscript thetas_table.R trichocarpa $file; done


args=(commandArgs(TRUE))
species <- args[1]
chr <- args[2]
window <- as.character(args[3])            ###100000
step <- as.character(args[4])             ####100000


#set working directory
wd=paste("/proj/b2011141/nobackup/PaperIV-phylogenomics/ANGSD/",species,"/",species,"_",chr,sep="")
setwd(wd)

args <- commandArgs(TRUE)
if (length(args) != 4) {
    message("Usage: plotThetas.R tremuloides 01 100000 100000!")
    quit("yes")
}

#set window size and step size
#window <- gsub("^[a-z]+_[0-9]+.[a-z]+([0-9]+)[a-z]+.*","\\1",tremula_input)
#step <- gsub("^[a-z]+_[0-9]+.[a-z]+([0-9]+)[a-z]+.*","\\1",tremula_input)

#read table
tremula_thetas <- read.table(paste(species,"_",chr,".thetas",window,".window.gz.pestPG",sep=""))

#set names of columns
names(tremula_thetas) <- c("range", "chr", "pos", "tW", "tP", "tF", "tH", "tL", "tajD", "fulif", "fuliD", "fayH", "zengsE", "numSites")

#normarize theta by numSites
tremula_thetas <- cbind(tremula_thetas, tremula_thetas[,c("tW", "tP", "tF", "tH", "tL")] / tremula_thetas$numSites)
names(tremula_thetas) <- c("range", "chr", "pos", "tW", "tP", "tF", "tH", "tL", "tajD", "fulif", "fuliD", "fayH", "zengsE", "numSites", "tW.norm", "tP.norm", "tF.norm", "tH.norm", "tL.norm")


#only selected the sites that have to been 20% left for analysis
window=as.numeric(window)
minimum_sites=0.1*window
tremula_thetas[which(tremula_thetas$numSites < minimum_sites),c(9,10,11,12,13,15,16,17,18,19)]=rep(NA,10)

tremula_theta=data.frame(cbind(Chr=as.character(tremula_thetas$chr), Pos=tremula_thetas$pos, numSites=tremula_thetas$numSites,tW.norm=tremula_thetas$tW.norm,tP.norm=tremula_thetas$tP.norm,tajD=tremula_thetas$tajD,fulif=tremula_thetas$fulif,fuliD=tremula_thetas$fuliD,fayH=tremula_thetas$fayH,zengsE=tremula_thetas$zengsE))

window=as.character(args[3])
write.table(tremula_theta, file=paste(species,"_",chr,".",window,"bp",step,"bp",".thetas.txt",sep="",collapse=""), sep="\t", quote=F, row.names=F, col.names=T)



