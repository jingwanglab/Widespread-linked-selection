library(dplyr)
library(ape)
#install.packages("phangorn")
library(phangorn)
library(ggtree)
#source("http://bioconductor.org/biocLite.R")
#biocLite("ggtree")
library(RColorBrewer)

setwd("~/Dropbox/davidiana_paper/data/population_structure/chloroplast")
species=read.table("4species.annotation.txt",header=T)


######cls_group->data.frame

cls=list(a1=as.character(species[which(species$CLUSTER=="P.tremula"),]$FID),
         a2=as.character(species[which(species$CLUSTER=="P.davidiana"),]$FID),
         a3=as.character(species[which(species$CLUSTER=="P.tremuloides"),]$FID),
         a4=as.character(species[which(species$CLUSTER=="P.trichocarpa"),]$FID)
         )


tr=read.tree("outtree.txt")
tr_group=groupOTU(tr,cls)

##ggplot default color values
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
n = 4
cols = gg_color_hue(n)

pdf("ggtree.ml.chloroplast.4species.pdf",width=7,height=6)
    #units='in',res=300)
par(oma=c(2,2,2,2))
par(mar=c(1,1,1,1))
ggtree(midpoint(tr_group))+geom_tippoint(aes(color=group),alpha=1,size=2)+geom_tiplab(align=F,size=2,linesize=2, hjust = 0)+ggplot2::xlim(0, 0.35)+
geom_treescale(offset=-2)+
  theme(legend.position=c(0.1,0.9),legend.key.size=unit(0.1,"cm"),legend.text=element_text(size=12))+
  scale_colour_manual(name="Species",values=cols[c(2,1,3,4)],labels=c(expression(italic("P. davidiana")),expression(italic("P. tremula")),expression(italic("P. tremuloides")),expression(italic("P. trichocarpa"))))
dev.off()

pdf("ggtree.ml.chloroplast.4species.2.pdf",width=4,height=6)
#units='in',res=300)
par(oma=c(2,2,2,2))
par(mar=c(1,1,1,1))
ggtree(midpoint(tr_group))+geom_tippoint(aes(color=group),alpha=1,size=2)+ggplot2::xlim(0, 0.35)+
  geom_treescale(offset=-2)+scale_colour_manual(values=cols[c(2,1,3,4)])
 # theme(legend.position=c(0.2,0.9),legend.key.size=unit(0.1,"cm"),legend.text=element_text(size=12))+
 # scale_colour_manual(name="Species",values=cols[c(2,1,3,4)],labels=c(expression(italic("P. davidiana")),expression(italic("P. tremula")),expression(italic("P. tremuloides")),expression(italic("P. trichocarpa"))))
dev.off()



pdf("ggtree.ml.chloroplast.4species.unrooted.pdf",width=7,height=6)
ggtree(tr_group,layout="unrooted")+geom_tippoint(aes(color=group),alpha=1,size=2)+geom_tiplab(align=F,size=2,linesize=2, hjust = 0)+
  theme(legend.position=c(0.1,0.9),legend.key.size=unit(0.1,"cm"),legend.text=element_text(size=12))+
  scale_colour_manual(name="Species",values=cols[c(2,1,3,4)],labels=c(expression(italic("P. davidiana")),expression(italic("P. tremula")),expression(italic("P. tremuloides")),expression(italic("P. trichocarpa"))))
dev.off()
  
