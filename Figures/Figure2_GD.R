

################################################
#  File Name:Heatmap_Figure2_GD.R
#  Author: Baoyan Bai


#################################################



##load required libraries
library("plyr")
library(dplyr)
library(reshape)
library(distances)

library(ggplot2)
library(xlsx)
library("gridExtra")


setwd("~/Desktop/FL/fl_latest")
#####Patients with mutiple biopsies


#####using clonal mutation from pyclone to calcuate genomic distance
#####Input file: output file from pyclone
setwd("~/Desktop/FL/11_pyclone/BED_bam/counts_new/pyclone_output_301018/ouput_pyclone/output_200X_i2")

empty<- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("Pair", "Pat",  "distance"))

for (i in pat) {
  
  cluster<- read.table(file=paste(i, "/tables/cluster.tsv",sep=""), sep="\t", header=T, stringsAsFactors = F)
  cluster.cast<- data.frame(cast(cluster, cluster_id ~ sample_id, value= "mean"))
  cluster.cast$max<- do.call(pmax, cluster.cast[c(-1)])
  
  cluster.filter<- cluster.cast%>% dplyr::filter(max>0.4)
  
  loci<- read.table(file=paste(i, "/tables/loci.tsv",sep=""), header=T, sep="\t", stringsAsFactors = F)
  
  loci.fil<- subset(loci, cluster_id %in% cluster.filter$cluster_id)
  
  loci.fil.tmp<- loci.fil%>% select(mutation_id, sample_id, cellular_prevalence)
  loci.cast<- cast(loci.fil.tmp, mutation_id ~ sample_id, value="cellular_prevalence")
  
  loci.cast<- as.matrix(loci.cast)
  dis<- melt(matrix(dist(t(loci.cast))))
  names(dis)<- c("Pair", "Pat", "distance")
  dis$Pat<- i
  empty<- rbind(empty, dis)
}
write.table(empty, file="~/Desktop/FL/fl_latest/34_GDist/FL_genomicdist_clonal.txt", sep="\t", row.name=F,quote=F)


####Plot genomic distance, with clinical information
library(cowplot)
library(forcats)
library(plyr)
fldist<- read.table(file="FL_gDistance_update.txt",sep="\t",header=T, stringsAsFactors = F)



fldist$Group[fldist$Group=="clinial tFL"]<- "tFL"


pat.order<- ddply(fldist, .(Pat), summarise, largest_dist=max(Distance_CF))

pat.order$Pat<- reorder(pat.order$Pat, -pat.order$largest_dist)
fldist$Pat<- factor(fldist$Pat, levels=levels(pat.order$Pat))
names(fldist)[names(fldist)=="classify"]<- "Comparison"
fldist$Comparison[fldist$Comparison=="FL3B-tFL"]<- "FL-tFL"
fldist$Comparison[fldist$Comparison=="tFL-tFL"]<- "FL-tFL"

library(tidyverse)
dis_plot<-ggplot(fldist, aes(x=Pat, y=Distance_CF, color=Comparison,size=Interval))+ 
  geom_point(position=position_jitter(h=0.1, w=0.1),
             shape = 19,alpha=0.7)+
  scale_color_manual(values=c("dark blue", "dark red"))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=12),
        axis.ticks.x=element_blank())+ ylab("Genomic Distance")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


library(reshape)
clinic<- xlsx2dfs("../2_clinic_infor/FL_Baoyan_clinical_file_FINAL_020719 - survival.xlsx")
clinic<- clinic$FL_Baoyan_clinical_file_FINAL_0
clinic$Patient_ID<- row.names(clinic)
clinic$Transf<- ifelse(is.na(clinic$Date_transf), "N", "Y")
clinic.info<- clinic%>%select(Patient_ID, Transf, POD24)
clinic.info.melt<- melt(clinic.info, id.vars = "Patient_ID")
names(clinic.info.melt)<- c("sample", "variable", "value")

clinic.info.melt$sample<- as.character(clinic.info.melt$sample)

pats.info<- unique(fldist%>%dplyr::select(Pat, Group))
pats.info$variable<- "Group"


clinic.info.30 <- subset(clinic.info.melt, sample %in% pats.info$Pat)
clinic.info.30$sample<- factor(clinic.info.30$sample, levels=levels(fldist$Pat))

clinic.info.30$value[clinic.info.30$value==0]<- "No"
clinic.info.30$value[clinic.info.30$value==1]<- "Yes"
clinic.info.30$value[clinic.info.30$value=="N"]<- "nFL"
clinic.info.30$value[clinic.info.30$value=="Y"]<- "tFL"
clinic.info.30$variable<- as.character(clinic.info.30$variable)

clinic.info.30$variable[clinic.info.30$variable=="Transf"]<- "Group"
clinic.info.30$variable<- factor(clinic.info.30$variable, levels=c( "Group", "POD24"))
clinic.info.30$value<- factor(clinic.info.30$value, levels=c("No", "Yes", "nFL", "tFL"))


clinic.anno<- ggplot(clinic.info.30, aes(x=sample,y=variable, fill=value))+ geom_tile(color="white")+
  scale_fill_manual(values = c("grey","black", "#998EC3","#F1A340"))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()
  )+theme_void()+
  guides(fill=guide_legend(ncol=2))+ylab("")

pdf(file="fl_dist_sorted_20211213.pdf", width =7, height=5)
plot_grid(dis_plot, clinic.anno, ncol = 1, align = "v", rel_heights=c(6,2))
dev.off()

#########Violin plot
####violin plot to show the genomic distance between nFL and tFL groups

p1<- ggplot(fldist, aes(x=Group, y=Distance_CF)) + 
  geom_violin(trim=FALSE)

pdf(file="fldist_group_violin_20210609.pdf", width=5, height=5)
p1+geom_dotplot( binaxis='y', stackdir='center', dotsize=1)+ theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size=12),
        axis.title.x =element_text(size=14),
        axis.title.y =element_text(size=14),
        axis.text.y = element_text(size=12))+ylab("Genomic Distance")+stat_compare_means()
# print: p-value 0.6

dev.off()
####violin plot to show the genomic distance between different comparisons
fldist$classify[fldist$classify=="tFL-tFL"]<- "FL-tFL"
fldist$classify[fldist$classify=="FL3B-tFL"]<- "tFL-tFL"

p2<- ggplot(fldist, aes(x=classify, y=Distance_CF)) + 
  geom_violin(trim=FALSE)

######Type of comparison
pdf(file="fldist_type_violin_202100609.pdf", width=5, height=5)
p2+geom_dotplot( binaxis='y', stackdir='center', dotsize=1)+ theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        axis.text.x = element_text(size=12),
        axis.text.y =element_text(size=12),
        axis.title.x = element_text(size=14),
        axis.title.y =element_text(size=14),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ 
  xlab("Type of Comparison")+ylab("Genomic Distance")+stat_compare_means()
#print: pvalue =0.37
dev.off()

###Based on POD status


p3<-ggplot(fldist, aes(x=POD24, y=Distance_CF)) + 
  geom_violin(trim=FALSE)
pdf(file="fldist_POD24_violin_20230913.pdf", width=5, height=5)
p3+geom_dotplot( binaxis='y', stackdir='center', dotsize=1)+ theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        axis.text.x = element_text(size=12),
        axis.text.y =element_text(size=12),
        axis.title.x = element_text(size=14),
        axis.title.y =element_text(size=14),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ 
  xlab("POD24 status")+ylab("Genomic Distance")+stat_compare_means()
dev.off()

pdf(file="fldist_POD24_violin_20230913_wPval.pdf", width=5, height=5)
p3+geom_dotplot( binaxis='y', stackdir='center', dotsize=1)+ theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        axis.text.x = element_text(size=12),
        axis.text.y =element_text(size=12),
        axis.title.x = element_text(size=14),
        axis.title.y =element_text(size=14),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ 
  xlab("POD24 status")+ylab("Genomic Distance")+stat_compare_means()
dev.off()

p1<- ggplot(fl.dist, aes(x=Defi_2020, y=Distance_CF))+geom_violin(trim=FALSE)
pdf(file="evolution_gDIST_group_wPvalue.pdf", width=5, height=5)
p1+geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+ theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        axis.text.x = element_text(size=12),
        axis.text.y =element_text(size=12),
        axis.title.x = element_text(size=14),
        axis.title.y =element_text(size=14),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ xlab("Evolution patten")+ylab("Genomic Distance")+
  facet_wrap(~Group, scale="free")+stat_compare_means()
dev.off()

pdf(file="evolution_gDIST_group_woPvalue.pdf", width=8, height=5)
p1+geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+ theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        axis.text.x = element_text(size=12),
        axis.text.y =element_text(size=12),
        axis.title.x = element_text(size=14),
        axis.title.y =element_text(size=14),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ xlab("Evolution patten")+ylab("Genomic Distance")+
  facet_wrap(~Group, scale="free")
dev.off()

pdf(file="evolution_gDIST_all_wPvalue.pdf", width=5, height=5)
p1+geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+ theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        axis.text.x = element_text(size=12),
        axis.text.y =element_text(size=12),
        axis.title.x = element_text(size=14),
        axis.title.y =element_text(size=14),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ xlab("Evolution patten")+ylab("Genomic Distance")+
  stat_compare_means()
dev.off()




pdf(file="evolution_gDIST_all_woPvalue.pdf", width=5, height=5)
p1+geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+ theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        axis.text.x = element_text(size=12),
        axis.text.y =element_text(size=12),
        axis.title.x = element_text(size=14),
        axis.title.y =element_text(size=14),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ xlab("Evolution patten")+ylab("Genomic Distance")

dev.off()

#pval : 6.9e-05
