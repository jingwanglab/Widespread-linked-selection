#! /usr/bin/Rscript --no-save --no-restore

library(data.table)

args=(commandArgs(TRUE))
setwd(args[1])

compare <- args[2]


fst_file=paste(compare,".all",".fst_stat",sep="")
fst_threshd="0.99"

outfile=paste(compare,".win.fst.stat",sep="")
window=args[3]
step=args[4]


values=fread(fst_file,header=F)
names(values)=c("CHR","POS","A","AB","f","FST","Pvar","fixed","dxy")
chr=as.list(matrix(,19))

for (n in 1:19) { 
  if (n<10) {
    chr[[n]]=paste("Chr","0",n,sep="")
  }
  else {
    chr[[n]]=paste("Chr",n,sep="")
  }


ind <- which(values$Pvar>=as.numeric(fst_threshd) & values$CHR==chr[[n]])
pos <- as.numeric(values$POS)
pos_chr <- as.numeric(values$POS[values$CHR==chr[[n]]])

fst_values <- values[ind,]
fst_pos <- pos[ind]
#cat("Removed",length(ind),"sites; now there are",nrow(values),"sites going from",min(fst_pos),"to",max(fst_pos),"\n")
#cat("Overall FST:",sum(fst_values$A)/sum(fst_values$AB),"\n")

#for Fst and fixed sites
fst_len=max(pos_chr)
window=as.numeric(window)
step=as.numeric(step)
#fst_start=seq(min(pos), fst_len, step)

fst_start=seq(step,fst_len, step)
fst_end=fst_start+window-1
fst_wpos=round(fst_start+(window/2)); # position of the window in the plot (center)
fst=c() 
fix=c()
dxy=c()

for (i in 1:length(fst_start)) {
	fst_ipos=which(fst_pos>=fst_start[i] & fst_pos<=fst_end[i])
	ipos=which(pos_chr>=fst_start[i] & pos_chr<=fst_end[i])
	if (length(ipos)<0.1*window){
	fst[i]=NA
	}
	else{
	fst[i]=sum(fst_values$A[fst_ipos])/sum(fst_values$AB[fst_ipos])
}
}

for (i in 1:length(fst_start)) {
        ipos=which(pos_chr>=fst_start[i] & pos_chr<=fst_end[i])
        if (length(ipos)<0.1*window){
        fix[i]=NA
	dxy[i]=NA
        }
        else{
        fix[i]=sum(values$fixed[ipos])/length(ipos)
        dxy[i]=sum(values$dxy[ipos])/length(ipos)
	
}
}

#fst[is.na(fst)]<-NA
#fix[is.na(fix)]<-NA

# Data
fst_df=data.frame(cbind(Pop=rep(chr[[n]],length(fst_wpos)), Pos=fst_wpos, fst=fst,fixed=fix,dxy=dxy));
fst_df[,2:5]=sapply(fst_df[,2:5], as.character)
fst_df[,2:5]=sapply(fst_df[,2:5], as.numeric)
write.table(fst_df, file=paste(outfile,".",chr[[n]],".window",window,"_step",step,".fst.dxy.txt",sep="",collapse=""), sep="\t", quote=F, row.names=F, col.names=T)
}

#For derived fixed sites in each species (P.tremula and P.tremuloides)

#tremula_fix_ind=which(values$Pvar>=as.numeric(fst_threshd) & values$FST>=as.numeric(fst_threshd) & values$MAF_tremula>=as.numeric(fst_threshd))
#tremuloides_fix_ind=which(values$Pvar>=as.numeric(fst_threshd) & values$FST>=as.numeric(fst_threshd) & values$MAF_tremuloides>=as.numeric(fst_threshd))

#tremula_fix_fst=data.frame(cbind(POS=values$POS[tremula_fix_ind],FST=values$FST[tremula_fix_ind],MAF_tremula=values$MAF_tremula[tremula_fix_ind]))
#tremuloides_fix_fst=data.frame(cbind(POS=values$POS[tremuloides_fix_ind],FST=values$FST[tremuloides_fix_ind],MAF_tremuloides=values$MAF_tremuloides[tremuloides_fix_ind]))

#tremula_fix_fst[,1:3]=sapply(tremula_fix_fst[,1:3],as.numeric)
#tremuloides_fix_fst[,1:3]=sapply(tremuloides_fix_fst[,1:3],as.numeric)

#write.table(tremula_fix_fst,file=paste("Chr",chr_n,".tremula_fixed.fst.txt",sep="",collapse=""),sep="\t",quote=F,row.names=F,col.names=T)
#write.table(tremuloides_fix_fst,file=paste("Chr",chr_n,".tremuloides_fixed.fst.txt",sep="",collapse=""),sep="\t",quote=F,row.names=F,col.names=T)


