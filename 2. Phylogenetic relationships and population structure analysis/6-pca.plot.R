library(ggplot2)

setwd("~/Dropbox/davidiana_paper/data/population_structure/pca")


##ggplot default color values
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
n = 4
cols = gg_color_hue(n)

pca=read.table("angsd_pca.PC.txt")
eval=read.table("angsd_pca.eval.txt",header=T)
eval$x[1]/sum(eval$x)

#########making the plot
#pc1-pc2
pdf("pca.4species.pc1_2.pdf",width=5,height=4)
ggplot(pca,aes(x=PC1,y=PC2,colour=Species))+geom_point(size=2)+
  xlab("PC1 (28.79%)")+ylab("PC2 (7.52%)")+
  theme_bw()+
  scale_colour_manual(values=cols,labels=c(expression(italic("P. tremula")),expression(italic("P. davidiana")),expression(italic("P. tremuloides")),expression(italic("P. trichocarpa"))))
dev.off()

#pc2-pc3
pdf("pca.4species.pc2_3.pdf",width=5,height=4)
ggplot(pca,aes(x=PC2,y=PC3,colour=Species))+geom_point(size=)+
  xlab("PC2 (7.52%)")+ylab("PC3 (5.33%)")+
  theme_bw()+
  scale_colour_manual(values=cols,labels=c(expression(italic("P. tremula")),expression(italic("P. davidiana")),expression(italic("P. tremuloides")),expression(italic("P. trichocarpa"))))
dev.off()

