#library(data.table)
library(RColorBrewer)
library(dplyr)
library(reshape2)
library(ggplot2)
library(survival)

colors <- brewer.pal(10,"Paired")[c(1,2,3,4,5,6,7,8,9,10,11,12)]
##ggplot default color values
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
n = 4
cols = gg_color_hue(n)


####Define function to summarize the .final.txt dataset into the simple table for making plot

scale_func<-function(dataset){
  ###summarize the dataset across different combinations
  dataset_sum=as.data.frame(dataset %>% group_by(time_index) %>% 
                              summarise(left_time_mean=mean(left_time_boundary),
                                        right_time_mean=mean(right_time_boundary),
                                        lambda_median=median(lambda),
                                        lambda_mean=mean(lambda),
                                        lambda_sd=sd(lambda),
                                        lambda_se=sd(lambda)/sqrt(length(lambda))))
  
  ##scaling to real time and population size
  mu=3.75e-8 
  g=15
  
  dataset_sum_new=mutate(dataset_sum,left_time_real=(left_time_mean/mu)*g, ###scale to real time (left)
                         right_time_real=(right_time_mean/mu)*g,  ###scale to real time (right)
                         ne_median=(1/lambda_median)/(2*mu),
                         ne_mean=(1/lambda_mean)/(2*mu),
                         ne_sd=(lambda_sd/lambda_mean)*(1/lambda_mean)/(2*mu),
                         ne_se=(lambda_sd/lambda_mean)*(1/lambda_mean)/(2*mu))

  
  ##remove the top and bottom two rows, since they may not very accurate 
  #dataset_sum_new2=dataset_sum_new[-c(1,2,nrow(dataset_sum_new)-2,nrow(dataset_sum_new)-1,nrow(dataset_sum_new)),]
  ###do not remove the top and bottome rows
  dataset_sum_new2=dataset_sum_new
  return(dataset_sum_new2)
}


###read in dataset
setwd("~/Dropbox/davidiana_paper/data/msmc2/")

######################################################################

#######hap8
#tremula
tremula_hap8=read.table("~/Dropbox/davidiana_paper/data/msmc2/tremula/hap8/tremula.hap8.msmc2.summary.txt",header=T)
tremula_hap8_sum=scale_func(tremula_hap8)
tremula_hap8_sum$Species="P.tremula"
#tremuloides
tremuloides_hap8=read.table("~/Dropbox/davidiana_paper/data/msmc2/tremuloides/hap8/tremuloides.hap8.msmc2.summary.txt",header=T)
tremuloides_hap8_sum=scale_func(tremuloides_hap8)
tremuloides_hap8_sum$Species="P.tremuloides"
#davidiana
davidiana_hap8=read.table("~/Dropbox/davidiana_paper/data/msmc2/davidiana/hap8/davidiana.hap8.msmc2.summary.txt",header=T)
davidiana_hap8_sum=scale_func(davidiana_hap8)
davidiana_hap8_sum$Species="P.davidiana"
#trichocapra
trichocarpa_hap8=read.table("~/Dropbox/davidiana_paper/data/msmc2/trichocarpa/hap8/trichocarpa.hap8.msmc2.summary.txt",header=T)
trichocarpa_hap8_sum=scale_func(trichocarpa_hap8)
trichocarpa_hap8_sum$Species="P.trichocarpa"

all_species_hap8=rbind(tremula_hap8_sum,tremuloides_hap8_sum,davidiana_hap8_sum,trichocarpa_hap8_sum)
all_species_hap8$Species=factor(all_species_hap8$Species,levels=c("P.tremula","P.davidiana","P.tremuloides","P.trichocarpa"))

all_species_hap8$ne_left=all_species_hap8$ne_mean-all_species_hap8$ne_sd
all_species_hap8$ne_right=all_species_hap8$ne_mean+all_species_hap8$ne_sd
all_species_hap8[which(all_species_hap8$ne_left<=0),]$ne_left=0
#all_species_hap8_new=all_species_hap8[which(all_species_hap8$time_index>1),]
##ggplot
ggplot(all_species_hap8)+
  geom_line(aes(x=log10(left_time_real),y=log10(ne_median),colour=Species))+
  geom_ribbon(aes(x=log10(left_time_real),ymin=log10((ne_median-ne_sd)),ymax=log10((ne_median+ne_sd)),fill=Species),alpha=0.2)+
  theme_bw()+scale_x_continuous(expand = c(0, 0), limits=c(2.4,6.5),breaks=c(3,4,5,6),labels=c(expression(10^3),expression(10^4),expression(10^5),expression(10^6)))+
  xlab("Time (years ago)")+ylab("Effective population size")+
  scale_y_continuous(expand = c(0, 0), limits=c(3.5,7),breaks=c(3,4,5,6),labels=c(expression(10^3),expression(10^4),expression(10^5),expression(10^6)))+
  annotation_logticks()+
  scale_colour_manual(values=cols,labels=c(expression(italic("P. tremula")),expression(italic("P. davidiana")),expression(italic("P. tremuloides")),expression(italic("P. trichocarpa"))))+
  scale_fill_manual(values=cols,labels=c(expression(italic("P. tremula")),expression(italic("P. davidiana")),expression(italic("P. tremuloides")),expression(italic("P. trichocarpa"))))

ggsave("msms2.4species.hap8.pdf",width=7,height=6)



######################################################################

#######hap2
#tremula
tremula_hap2=read.table("~/Dropbox/davidiana_paper/data/msmc2/tremula/hap2/tremula.hap2.msmc2.summary.txt",header=T)
tremula_hap2_sum=scale_func(tremula_hap2)
tremula_hap2_sum$Species="P.tremula"
#tremuloides
tremuloides_hap2=read.table("~/Dropbox/davidiana_paper/data/msmc2/tremuloides/hap2/tremuloides.hap2.msmc2.summary.txt",header=T)
tremuloides_hap2_sum=scale_func(tremuloides_hap2)
tremuloides_hap2_sum$Species="P.tremuloides"
#davidiana
davidiana_hap2=read.table("~/Dropbox/davidiana_paper/data/msmc2/davidiana/hap2/davidiana.hap2.msmc2.summary.txt",header=T)
davidiana_hap2_sum=scale_func(davidiana_hap2)
davidiana_hap2_sum$Species="P.davidiana"
#trichocapra
trichocarpa_hap2=read.table("~/Dropbox/davidiana_paper/data/msmc2/trichocarpa/hap2/trichocarpa.hap2.msmc2.summary.txt",header=T)
trichocarpa_hap2_sum=scale_func(trichocarpa_hap2)
trichocarpa_hap2_sum$Species="P.trichocarpa"

all_species_hap2=rbind(tremula_hap2_sum,tremuloides_hap2_sum,davidiana_hap2_sum,trichocarpa_hap2_sum)
all_species_hap2$Species=factor(all_species_hap2$Species,levels=c("P.tremula","P.davidiana","P.tremuloides","P.trichocarpa"))

##ggplot
ggplot(all_species_hap2)+
  geom_line(aes(x=log10(left_time_real),y=log10(ne_median),colour=Species))+
  geom_ribbon(aes(x=log10(left_time_real),ymin=log10((ne_median-ne_sd)),ymax=log10((ne_median+ne_sd)),fill=Species),alpha=0.2)+
  theme_bw()+scale_x_continuous(expand = c(0, 0), limits=c(4,6.7),breaks=c(4,5,6,7),labels=c(expression(10^4),expression(10^5),expression(10^6),expression(10^7)))+
  xlab("Time (years ago)")+ylab("Effective population size")+
  scale_y_continuous(expand = c(0, 0), limits=c(3.5,6),breaks=c(4,5,6),labels=c(expression(10^4),expression(10^5),expression(10^6)))+
  annotation_logticks()+
  scale_colour_manual(values=cols,labels=c(expression(italic("P. tremula")),expression(italic("P. davidiana")),expression(italic("P. tremuloides")),expression(italic("P. trichocarpa"))))+
  scale_fill_manual(values=cols,labels=c(expression(italic("P. tremula")),expression(italic("P. davidiana")),expression(italic("P. tremuloides")),expression(italic("P. trichocarpa"))))

ggsave("msms2.4species.hap2.pdf",width=7,height=6)

######################################################################

#######hap4
#tremula
tremula_hap4=read.table("~/Dropbox/davidiana_paper/data/msmc2/tremula/hap4/tremula.hap4.msmc2.summary.txt",header=T)
tremula_hap4_sum=scale_func(tremula_hap4)
tremula_hap4_sum$Species="P.tremula"
#tremuloides
tremuloides_hap4=read.table("~/Dropbox/davidiana_paper/data/msmc2/tremuloides/hap4/tremuloides.hap4.msmc2.summary.txt",header=T)
tremuloides_hap4_sum=scale_func(tremuloides_hap4)
tremuloides_hap4_sum$Species="P.tremuloides"
#davidiana
davidiana_hap4=read.table("~/Dropbox/davidiana_paper/data/msmc2/davidiana/hap4/davidiana.hap4.msmc2.summary.txt",header=T)
davidiana_hap4_sum=scale_func(davidiana_hap4)
davidiana_hap4_sum$Species="P.davidiana"
#trichocapra
trichocarpa_hap4=read.table("~/Dropbox/davidiana_paper/data/msmc2/trichocarpa/hap4/trichocarpa.hap4.msmc2.summary.txt",header=T)
trichocarpa_hap4_sum=scale_func(trichocarpa_hap4)
trichocarpa_hap4_sum$Species="P.trichocarpa"

all_species_hap4=rbind(tremula_hap4_sum,tremuloides_hap4_sum,davidiana_hap4_sum,trichocarpa_hap4_sum)
all_species_hap4$Species=factor(all_species_hap4$Species,levels=c("P.tremula","P.davidiana","P.tremuloides","P.trichocarpa"))

##ggplot
ggplot(all_species_hap4)+
  geom_line(aes(x=log10(left_time_real),y=log10(ne_median),colour=Species))+
  geom_ribbon(aes(x=log10(left_time_real),ymin=log10((ne_median-ne_sd)),ymax=log10((ne_median+ne_sd)),fill=Species),alpha=0.2)+
  theme_bw()+scale_x_continuous(expand = c(0, 0), limits=c(3.5,6.5),breaks=c(4,5,6),labels=c(expression(10^4),expression(10^5),expression(10^6)))+
  xlab("Time (years ago)")+ylab("Effective population size")+
  scale_y_continuous(expand = c(0, 0), limits=c(3.5,5.5),breaks=c(4,5),labels=c(expression(10^4),expression(10^5)))+
  annotation_logticks()+
  scale_colour_manual(values=cols,labels=c(expression(italic("P. tremula")),expression(italic("P. davidiana")),expression(italic("P. tremuloides")),expression(italic("P. trichocarpa"))))+
  scale_fill_manual(values=cols,labels=c(expression(italic("P. tremula")),expression(italic("P. davidiana")),expression(italic("P. tremuloides")),expression(italic("P. trichocarpa"))))

ggsave("msms2.4species.hap4.pdf",width=7,height=6)




