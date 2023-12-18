
################################################
#  File Name:Figure3.R
#  Author: Baoyan Bai


#################################################

####Integrate three different mutation files:
#1. original mutation calling file from Sigve
#   this files has been filtered by Ole-Christian'script and by Jillian's TCGA PON 
#2. raw vcf files from TSD. Many mutations in serial biopsies were in fact called, but rejcted for different reasons. 
#   I started to work this. Noticed some of mutations were still not called. But if looked at coverage, variant allele can be detected

#3. Pyclone result. Some mutations are missing due to CN=0 by ascat. I tried to recovery those mutations
#
# two mutations targeting on TBL1XR1 in P46 have wrong AF by IGV inspection. Those mutations were also manually added

#######Aim of 

#######1. calcuate which mutation is clonal, subclonal
#######2. Check whether's possible to detect the variants in control DNA and which group of mutations can be detected in blood DNAs

####data.frame empty.raw: raw mutation file.
####this file: from raw VCF files. add additional columns to indicate the mutation is shared by all biopsy for a given patient.

####data.frame empty pyclone result
###flcombined imported from 16_combined/fl42_combined_updated010719.txt the latest mutation file, original file from Sigve

library(plyr)
library(dplyr)
library(ggplot2)


setwd("G:/FL_resequncing/FL_exome_final/fl_latest/11_pyclone/BED_bam/counts_new/CPC_codingPlusUTR_pyclone/All_codingPlusUTR/All_codingPlusUTR_i1")
flcombined<- read.csv(file="../fl_latest/16_combined/fl44_reseq_final_060619.txt", sep="\t",header=T, stringsAsFactors = F)

mutation<- unique(flcombined%>%select(Patient, Gene, SYMBOL,VAR_ID, CDS_CHANGE,HGVSp_short,Region,
                                      VARIANT,VARIANT_CLASS))
mutation$chr<- as.character(lapply(strsplit(as.character(mutation$VAR_ID), "\\_"),"[",1))
mutation$star<- as.character(lapply(strsplit(as.character(mutation$VAR_ID), "\\_"),"[",2))
mutation$mutation_id<- paste(mutation$chr, mutation$star,sep="_")


###process JP6 seperately because JP6 does not have raw VCF file
JP6_anno<- mutation%>%filter(Patient=="P6B")
#########annotate pyclone result
pyclone.anno<- merge(empty, mutation, by=c("Patient", "mutation_id"))

###raw vcf df
###mutation calling
empty.raw$mutation_id<- paste(empty.raw$V1, empty.raw$V2,sep="_")

raw.vcf<- empty.raw%>%select(Patient,SAMPLE,mutation_id, VAR_ID, V7,V8, num_biopsy,sharedbyAll )

raw.vcf.anno<- merge(raw.vcf, mutation, by=c("Patient", "VAR_ID"))

##extract only coding
raw.vcf.anno.coding<- raw.vcf.anno%>%filter((Region=="Coding")| (VARIANT=="3'UTR") | (VARIANT=="5'UTR"))

#########
fl32<- subset(raw.vcf.anno.coding, Patient %in% unique(pyclone.anno$ Patient))
names(fl32)[names(fl32)=="mutation_id.x"]<- "mutation_id"

pyclone.df<- pyclone.anno%>%select(SAMPLE, mutation_id, cellular_prevalence,variant_allele_frequency, mean, class)


fl32.pyclone<- merge(fl32, pyclone.df, by=c("SAMPLE", "mutation_id"), all.x = T)

fl32.pyclone.removed<- fl32.pyclone%>%filter(is.na(cellular_prevalence))
missed_gene<- c("TAS1R3","ATAD3C","MMP23B","MORN1","HES2","TRAV8-6","SNRPA1","DLL4","VPS39","PAQR5","CEMIP","AMELY","ANHX","ATAD3C","B3GALT6","BRCA2","CBLB","CEMIP","DLL4","ESPN","FAM46A","HES2","IL5RA","IL9R","LCE1E","LEPROTL1","MAP1A","MMP23B","MORN1","PAQR5",
                "PPIAL4G","RAB7A","RORC","SNRPA1","TAS1R3","TNFRSF14","TRAV8-6","TRMT10A","TSPY1","USF3","USP9Y","VPREB1","VPS39","ZBTB20","ZNF268")

fl32.missed<- subset(fl32.pyclone.removed, SYMBOL %in%missed_gene)
fl32.ig<- setdiff(fl32.pyclone.removed, fl32.missed)

#######have to manually add clonal informations for missed genes
write.table(fl32.missed,file="clonal_gene/fl32.missed_genes.txt",sep="\t",row.names=F, quote=F )
########read manually revised genes
fl32.missed<- read.table(file="clonal_gene/fl32.missed_genes_manulrevised.txt",sep="\t", quote="",header=T,stringsAsFactors = F)


fl32.pyclone.kept<- fl32.pyclone%>%filter(!is.na(cellular_prevalence))

fl32.pyclone.new<- rbind(fl32.pyclone.kept, fl32.missed)

########defining mutation order: early dominant, dominant, subclonal


###
biopsy_empty_new<- fl32.pyclone.new[-FALSE,]
for (i in unique(fl32.pyclone.new$SAMPLE)){
  
  biopsy<- subset(fl32.pyclone.new, SAMPLE==i)
  biopsy_specific<- biopsy%>%filter(variant_allele_frequency>0)
  biopsy_specific$clonal<- ifelse(biopsy_specific$cellular_prevalence>=0.5 | biopsy_specific$mean>=0.5, "dominant", "subclonal")
  biopsy_empty_new<- rbind(biopsy_empty_new, biopsy_specific)
  
}


removed<- fl32.pyclone.new%>%filter(variant_allele_frequency==0)
##remove"TBL1XR1", discard_rescue_TBL1XR1

recover<- subset(removed, class=="discard_rescue_TBL1XR1")
recover$clonal<- "dominant"

fl32.missed$clonal<- fl32.missed$class

######
fl32.all<- rbind( biopsy_empty_new, recover,fl32.missed)

mutation_af<- ddply(fl32.all, .(mutation_id, Patient), summarise, min_CF=min(cellular_prevalence))

fl32.all.af<- merge(fl32.all, mutation_af,by=c("Patient","mutation_id"))


fl32.all.af$clonal[fl32.all.af$min_CF>0.5 & fl32.all.af$sharedbyAll=="Y" ]<-"early_dominant"

###corrected on April 5th 2019
fl32.all.af$clonal[fl32.all.af$class=="CPC_dominant_lost_LOH" & fl32.all.af$clonal=="dominant"]<-"early_dominant"
fl32.all.af$clonal[fl32.all.af$class=="CPC_dominant"]<- "early_dominant"


fl32.final<- fl32.all.af[,-c(16,17,18)]
##########process JP6 ,see previous

JP6_final<- JP6_df[,-c(4,5,7,9,11,13,14,15,16,17,18,19)]

JP6_final<- ddply(JP6_final, .(mutation_id), mutate, min_CF=min(cellular_prevalence))
fl32.all.final<-rbind(fl32.final,JP6_final)

##########
write.table(fl32.all.final, file="clonal_gene/fl32.all.coding_040319.txt",sep="\t", row.names=F,quote=F)
write.table(fl32.all.final, file="clonal_gene/fl32.all.coding_040519.txt",sep="\t", row.names=F,quote=F)


####plotting
fl32<- read.table(file="clonal_gene/fl32.all.coding_040519.txt",sep="\t",header=T,quote="",stringsAsFactors = F)
fl32.dominant<- fl32%>%filter(clonal !="subclonal")
fl32.dominant.ns<- fl32.dominant%>%filter(Region=="Coding" & VARIANT!="Silent")%>%filter(class!="discard")

fl32.dominant.ns.df<- unique(fl32.dominant.ns%>%select(Patient, SYMBOL,clonal))

fl32.dominant.ns.df.count<- ddply(fl32.dominant.ns.df, .(SYMBOL),mutate, num_pat=length(unique(Patient)))

clinic<- read.xlsx(file="../../../../../../2_clinic_infor/FL_Baoyan_clinical_file_FINAL_020719 - survival.xlsx",sheetIndex = 1)
clinic$group<- ifelse(is.na(clinic$Date_transf), "nFL", "tFL")
clinic.reduced<- clinic%>%select(Patient_ID, group, POD24)
names(clinic.reduced)<- c("Patient", "group", "POD24")

fl32.dominant.clinic<- merge(fl32.dominant.ns.df.count,clinic.reduced,by="Patient", all.x=T )
#########
discard_genes<- c("TTN","MUC16","MUC12")

fl32.dominant.clinic.clean<- subset(fl32.dominant.clinic, !SYMBOL %in% discard_genes)
#######taking genes mutated in at least 4 cases


dominant_4cases<- fl32.dominant.clinic.clean%>%filter(num_pat>=4)
dominant_4cases$SYMBOL<- reorder(dominant_4cases$SYMBOL, dominant_4cases$num_pat)


dominant_4cases$SYMBOL_ordered <- with(dominant_4cases, reorder(SYMBOL, SYMBOL, function(x) -length(x)))
####sepearte nFL and tFL

dominant_4cases.nFL<- dominant_4cases%>%filter(group=="nFL")

tmp.nFL<-ddply(dominant_4cases.nFL, .(Patient,SYMBOL), mutate, status=length(Patient))
tmp.nFL$clonal_correct<-tmp.nFL$clonal
tmp.nFL$clonal_correct[tmp.nFL$status=="2"]<-"Both"

nFL_plot<- unique(tmp.nFL[, -c(3)])
nFL_plot$clonal_correct[nFL_plot$clonal_correct=="dominant"]<- "Late"
nFL_plot$clonal_correct[nFL_plot$clonal_correct=="early_dominant"]<- "CPC"
nFL_plot$clonal_correct<-factor(nFL_plot$clonal_correct, levels=c("Both","Late","CPC"))


nFL<-ggplot(nFL_plot, aes(x=SYMBOL,fill=clonal_correct,drop=FALSE))+geom_bar(width=0.5) +
  scale_fill_manual(values=c("dark blue","pink","purple" ))+theme_bw()+coord_flip()+
  theme(legend.position = "none",
        axis.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 18,face="bold"),
        legend.text=element_text(size=8,face="bold"),
        legend.title = element_text(size=8, face="bold"),
        plot.margin = unit(c(0.2,0,0.2,0.2), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border= element_rect(color="blue",size=2))+
  scale_x_discrete(drop=FALSE)+
  scale_y_continuous(breaks = c(0,5,10,15))+xlab("")+ylab("")



############separate tFL group

dominant_4cases.tFL<- dominant_4cases%>%filter(group=="tFL")
tmp.tFL<-ddply(dominant_4cases.tFL, .(Patient,SYMBOL), mutate, status=length(Patient))
tmp.tFL$clonal_correct<-tmp.tFL$clonal
tmp.tFL$clonal_correct[tmp.tFL$status=="2"]<-"Both"

tFL_plot<- unique(tmp.tFL[, -c(3)])
tFL_plot$clonal_correct[tFL_plot$clonal_correct=="dominant"]<- "Late"
tFL_plot$clonal_correct[tFL_plot$clonal_correct=="early_dominant"]<- "CPC"
tFL_plot$clonal_correct<-factor(tFL_plot$clonal_correct, levels=c("Both","Late","CPC"))

tFL<-ggplot(tFL_plot, aes(x=SYMBOL,fill=clonal_correct))+geom_bar(width=0.5) +
  scale_fill_manual(values=c("dark blue","pink","purple" ))+coord_flip()+theme_bw()+
  theme(legend.position = "none",
        axis.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 18,face="bold"),
        legend.text=element_text(size=8,face="bold"),
        legend.title = element_text(size=8, face="bold"),
        plot.margin = unit(c(0.2,0.2,0.2,-0.1), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border=element_rect(color="red", size=2))+xlab("")+ylab("")+
  scale_x_discrete(drop=FALSE)+scale_y_continuous(breaks = c(0,5,10,15))

pdf(file="CPC_plots/24genes_freq_groups_20191203.pdf",width=4, height=6)
plot_grid(nFL, 
          tFL, align="h"
          
          , nrow = 1,rel_widths=c(1,1))
dev.off()
#######try to plot them together, waterfall plot, frequency plot, clinic data
all<- rbind(nFL_plot, tFL_plot)

patient_order_nFL<- unique(fl32.dominant.clinic.clean$Patient[fl32.dominant.clinic.clean$group=="nFL"])

patient_order_tFL<- unique(fl32.dominant.clinic.clean$Patient[fl32.dominant.clinic.clean$group=="tFL"])

all$Patient<- factor(all$Patient, levels=c(patient_order_nFL,patient_order_tFL))
all$clonal_correct<- factor(all$clonal_correct,levels=c("Both", "CPC", "Late", "NA"))
all.tmp<- all%>%select(Patient, SYMBOL, clonal_correct)
all.tmp.cast<- cast(all.tmp, Patient~SYMBOL)

all.melt.again<- melt(all.tmp.cast, id.vars="Patient")
all.melt.again$value<- factor(all.melt.again$value, levels=c("Both", "CPC", "Late"))
all.plot<- ggplot(all.melt.again, aes(x=Patient, y=SYMBOL, fill=value))+
  geom_tile(color="gray", size=0.1)+ scale_fill_manual(values=c("dark blue","purple","pink"),drop=FALSE)+scale_x_discrete(drop=FALSE)+
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y =element_text(size=15,face="bold"),
        legend.position = "NULL",
        plot.margin = unit(c(0.5,-0.5,-0.5,0.5), "cm"))+ xlab("")+ylab("")

clinic.32<-unique(fl32.dominant.clinic.clean%>%select(Patient, group, POD24))

names(clinic.32)<- c("Patient", "Group", "POD24")
clinic.32$POD24<- factor(clinic.32$POD24, levels=c(0,1))
clinic.32$Patient<- factor(clinic.32$Patient, levels=c(patient_order_nFL,patient_order_tFL))
clinic.32.melt<- melt(clinic.32, id.vars = "Patient")

clinic_plot<- ggplot(clinic.32.melt, aes(x=Patient, y=variable, fill=value))+
  geom_tile()+
  scale_fill_manual(values=c("blue", "red", "grey","black"))+
  
  theme(panel.background = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y =element_text(size=16,face="bold"),
        axis.ticks.x=element_blank(),
        legend.position = "none",
        plot.margin = unit(c(-0.1,0.5,-0.5,0.5), "cm"))+
  xlab("")+ylab("")

library(egg)



####legend
tFL_plot$clonal_correct<- factor(tFL_plot$clonal_correct, levels=c("Both", "CPC", "Late"))
pdf(file="CPC_plots/Figure4B_20191112_legend.pdf",width=12,height=13)
ggplot(tFL_plot, aes(x=SYMBOL,fill=clonal_correct))+geom_bar(width=0.5) +
  scale_fill_manual(values=c("dark blue","purple","pink"))+coord_flip()+theme_bw()+
  theme(
    axis.text.y=element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 18,face="bold"),
    legend.text=element_text(size=14,face="bold"),
    legend.title = element_text(size=16, face="bold"),
    plot.margin = unit(c(0.2,0.2,0.2,-0.2), "cm"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank())+xlab("")+ylab("")+labs(fill = "Timing")+
  scale_x_discrete(drop=FALSE)+panel_border(color="red", size=2)+scale_y_continuous(breaks = c(0,5,10,15))

dev.off()
save.image("G:/FL_resequncing/FL_exome_final/fl_latest/11_pyclone/BED_bam/counts_new/CPC_codingPlusUTR_pyclone/All_codingPlusUTR/All_codingPlusUTR_i1/Figure4B_20191112.RData")

###updated on 2019.11.15




pdf(file="CPC_plots/24genes_groups_20191203.pdf",width=6, height=5.8)
plot_grid(all.plot, 
          clinic_plot, align="v"
          
          , nrow = 2,rel_heights =c(14,2))
dev.off()
