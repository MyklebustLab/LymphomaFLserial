
################################################
#  File Name:Figure1_drivers.R
#  Author: Baoyan Bai


#################################################


########Manually plot driver genes (DO NOT USE GENVISR)
########date: 20211012, plot 51 genes with VAF adjusted for tumor% >=0.15
########add clinic features: group, POD24 infor
########genes ordered based on mutation frequency

library(maftools)
library(plyr)
library(dplyr)
library(pheatmap)
library(xlsx)
library(reshape)
library(readxl)
library(ggplot2)
library(RColorBrewer)
library(egg)
setwd("~/Desktop/FL/fl_latest")
flreseq<- read.table(file="16_combined/fl44_reseq_final_060619.txt",sep="\t",header = T, 
                     stringsAsFactors = F,quote="")

drivers<- read.table(file="28_drivergenes/51drivers.txt",header=T,stringsAsFactors = F)


###clinical information

clinic_biopsy<- read.table(file="2_clinic_infor/FL94biopsies_Sigve_20211125.txt",sep="\t",header=T,stringsAsFactors = F)
clinic_biopsy$type[clinic_biopsy$type=="pre&transform"]<- "transform"




clinic_biopsy_redced<- clinic_biopsy%>%select(Sample_name, type, group, POD24.1)
names(clinic_biopsy_redced)<- c("sample", "type", "group", "POD24")
clinic_biopsy_redced.melt<- melt(clinic_biopsy_redced, id.vars = "sample")
clinic_biopsy_redced.melt$variable<- factor(clinic_biopsy_redced.melt$variable, levels=c("group","POD24","type"))

####clinic information




###making data.frame for waterfall plot

plot_df_all<- unique(flreseq%>%filter(Region=="Coding"& VARIANT!="Silent")%>%select(SAMPLE, SYMBOL,VARIANT))
names(plot_df_all)<- c("sample", "gene", "variant_class")
plot_df_all$variant_class[plot_df_all$variant_class=="Frame_Shift_Ins" |plot_df_all$variant_class=="Frame_Shift_Del"]<- "Frame_Shift"
plot_df_all$variant_class[plot_df_all$variant_class=="In_Frame_Ins" |plot_df_all$variant_class=="In_Frame_Del"]<- "In_Frame"

sample_order<- c(as.character(clinic_biopsy_redced$sample[clinic_biopsy_redced$group=="nFL"]),
                 as.character(clinic_biopsy_redced$sample[clinic_biopsy_redced$group=="tFL"& clinic_biopsy_redced$type!="transform"]) ,
                 as.character(clinic_biopsy_redced$sample[clinic_biopsy_redced$group=="tFL" & clinic_biopsy_redced$type=="transform"]))

remove_sample<- c("P34_1", "P35_1", "P58_2")
sample_order_94<- setdiff(sample_order, remove_sample )







plot_df_all$variant_class[plot_df_all$variant_class=="Translation_Start_Site" |plot_df_all$variant_class=="Nonstop_Mutation"]<- "Translation_Site"

plot_df_all$variant_class[plot_df_all$variant_class=="Missense_Mutation"]<- "Missense"

plot_df_all$variant_class[plot_df_all$variant_class=="Nonsense_Mutation"]<- "Nonsense"

plot_df_all$variant_class<- factor(plot_df_all$variant_class, 
                                   levels=c("Frame_Shift","Nonsense", "Splice_Site", "In_Frame", "Translation_Site", "Missense"))

##########
plot_df<- ddply(plot_df_all, .(sample, gene), mutate, variant_class_correct=min(as.integer( variant_class )))



plot_df$variant_1<-levels(plot_df$variant_class)[plot_df$variant_class_correct]


plot_df.reduced<- unique(plot_df%>%select(sample, gene, variant_1))
plot_df.cast<- cast(plot_df.reduced, gene~ sample)

plot_df.melt<- melt(plot_df.cast, id.vars="gene")

plot_df.drivers<- subset(plot_df.melt, gene %in%drivers$SYMBOL)


plot_df.drivers$sample<-factor(plot_df.drivers$sample, levels=sample_order_94)
plot_df.drivers$value<- as.character(plot_df.drivers$value)


plot_df.drivers$value<- factor(plot_df.drivers$value, 
                               levels=c("Frame_Shift", "Nonsense","Splice_Site", "In_Frame", "Translation_Site","Missense"))


#####sort genes based on patient-level frequency

plot_df_pat<- unique(flreseq%>%filter(Region=="Coding"& VARIANT!="Silent")%>%select(Patient, SYMBOL)%>% filter(SYMBOL %in% drivers$SYMBOL))

###########

drivers.sum.pat<- ddply(plot_df_pat, .(SYMBOL), summarise,, Patient=length(Patient))

drivers.sum.pat$SYMBOL<- reorder(drivers.sum.pat$SYMBOL, drivers.sum.pat$Patient)



###
plot_df.drivers$gene<- factor(plot_df.drivers$gene, levels=levels(drivers.sum.pat$SYMBOL))


mutation<-ggplot(plot_df.drivers, aes(x=sample, y=gene,fill=value))+geom_tile(color="gray", size=0.1)+
  scale_fill_manual(values=c("#e7298a","firebrick","navyblue","#1b9e77", "#7570b3","#66a61e","white"),na.translate=FALSE)+scale_x_discrete(drop=FALSE)+
  theme(panel.background = element_blank(),
        legend.text = element_text(size=14,face="bold"),
        legend.title = element_text(size=16,face="bold"),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "none",
        axis.text.x =element_blank(),
        axis.text.y =element_text(size=10,hjust =0),
        plot.margin = unit(c(0.5,0.5,-0.2,0.5), "cm"))+ xlab("")+ylab("")
guides(fill=guide_legend(title="Mutation Type"))



#####frequency plot
pat<- unique(clinic_biopsy%>% select(Patient, group))

drivers.pat.clinic<-merge(plot_df_pat, pat,by="Patient")
drivers.pat.clinic$SYMBOL<- factor(drivers.pat.clinic$SYMBOL, levels=levels(plot_df.drivers$gene))
drivers.pat.clinic$group<- factor(drivers.pat.clinic$group, levels=c("tFL", "nFL"))

plot_freq<- ggplot(drivers.pat.clinic, aes(x=SYMBOL, fill=group))+ geom_bar()+  scale_fill_manual(values=c("red","blue"))+
  coord_flip()+theme_bw()+ theme(axis.text.y = element_blank() ,
                                 axis.text.x=element_text(size=12),
                                 legend.position = c(0.80, 0.1),
                                 axis.ticks.y=element_blank(),
                                 plot.margin = unit(c(0.3,0.3,-0.8,-0.5), "cm"))+
  ylab("Number of patients")+xlab("")

freq<-plot_freq+ guides(fill = guide_legend(reverse=TRUE))





#######clinic information
clinic_biopsy_redced.melt$sample<- factor(clinic_biopsy_redced.melt$sample, levels=sample_order_94)
clinic.final<- subset(clinic_biopsy_redced.melt,sample %in%sample_order_94)
clinic.final$variable<- as.character(clinic.final$variable)

clinic.final$variable[clinic.final$variable=="group"]<- "Group"
clinic.final$variable[clinic.final$variable=="type"]<- "Biopsy"
clinic.final$variable<- factor(clinic.final$variable,levels=c("Group", "POD24", "Biopsy"))

clinic.final$value<- factor(clinic.final$value, levels=c("pre", "relapse", "transform", "nFL", "tFL","0", "1"))
clinic.plot<- ggplot(clinic.final, aes(x=sample, y=variable, fill=value))+geom_tile()+
  scale_fill_manual(values=c("dark cyan","#fbb4ae", "purple","#998EC3","#F1A340", "grey", "black"))+
  theme(panel.background = element_blank(),
        axis.ticks.x=element_blank(),
        legend.text = element_text(size=12,face="bold"),
        legend.title = element_text(size=14,face="bold"),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        legend.position="right",
        axis.text.y =element_text(size=14,face="bold"),
        plot.margin = unit(c(0.1,0.2,0,2,0.5), "cm"))+ xlab("")+ylab("")+
  guides(fill=guide_legend(ncol=2, title="Biopsy      Group"))

###select colors

"#998EC3","#F1A340",
     
clinic.plot.nolegend<- ggplot(clinic.final, aes(x=sample, y=variable, fill=value))+geom_tile()+
  scale_fill_manual(values=c("dark cyan","#fbb4ae", "purple","#998EC3","#F1A340", "grey", "black"))+
  theme(panel.background = element_blank(),
        axis.ticks.x=element_blank(),
        legend.text = element_text(size=12,face="bold"),
        legend.title = element_text(size=14,face="bold"),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        legend.position="null",
        axis.text.y =element_text(size=12),
        plot.margin = unit(c(-0.2,0.5,0.5,0.5), "cm"))+ xlab("")+ylab("")

##########

#####arrange the plots


pdf(file="figures/Final_figures/51drivers_complete_clinic_20211130.pdf", width=13,height=8)
ggarrange(mutation,  clinic.plot.nolegend, nrow=2, ncol=1,  heights=c(10,1))

dev.off()

###########stop

#####making forest plot
clinic<- read_excel("./2_clinic_infor/FL_Baoyan_clinical_file_FINAL_020719_use.xlsx", sheet = 1)

clinic$group<- ifelse(is.na(clinic$Date_transf), "nFL", "tFL")

clinic.info<- clinic%>%select(Patient_ID, group)

names(clinic.info)<- c("Patient", "Group")

flreseq$chrom<- as.character(lapply(strsplit(as.character(flreseq$VAR_ID), "\\_"),"[",1))
flreseq$pos<-as.character(lapply(strsplit(as.character(flreseq$VAR_ID), "\\_"),"[",2))
flreseq$ref<-as.character(lapply(strsplit(as.character(flreseq$VAR_ID), "\\_"),"[",3))
flreseq$alt<-as.character(lapply(strsplit(as.character(flreseq$VAR_ID), "\\_"),"[",4))



flreseq_sample.maf<- unique(flreseq%>%select(Patient, SYMBOL, chrom, pos, ref, alt, VARIANT,VARIANT_CLASS, HGVSp_short))



names(flreseq_sample.maf)<- c("Tumor_Sample_Barcode",
                              "Hugo_Symbol", "Chromosome", "Start_Position", "Reference_Allele", 
                              "Tumor_Seq_Allele2", "Variant_Classification", "Variant_Type", "HGVSp_short")

flreseq_sample.maf$End_Position<- flreseq_sample.maf$Start_Position  
flreseq_sample.maf$Variant_Type[flreseq_sample.maf$Variant_Type=="snv"]<- "SNP"
flreseq_sample.maf$Variant_Classification[flreseq_sample.maf$Variant_Classification=="3' UTR"]<- "3'UTR"
flreseq_sample.maf$Variant_Classification[flreseq_sample.maf$Variant_Classification=="5' UTR"]<- "5'UTR"


flreseq_sample.maf$Variant_Type[flreseq_sample.maf$Variant_Type=="deletion"]<- "DEL"
flreseq_sample.maf$Variant_Type[flreseq_sample.maf$Variant_Type=="inDEL"]<- "DEL"
flreseq_sample.maf$Variant_Type[flreseq_sample.maf$Variant_Type=="insertion"]<- "INS"
flreseq_sample.maf$Variant_Type[flreseq_sample.maf$Variant_Type=="INdel"]<- "INS"
flreseq_sample.maf$Variant_Type[flreseq_sample.maf$Variant_Type=="substitution"]<- "SUB"
flreseq_sample.maf$Variant_Type[flreseq_sample.maf$Variant_Type=="SNV"]<- "SNP"
#######seperate nFL and tFL



flreseq.nFL<- subset(flreseq_sample.maf, Tumor_Sample_Barcode%in% clinic.info$Patient[clinic.info$Group=="nFL"])
flreseq.tFL<- subset(flreseq_sample.maf, Tumor_Sample_Barcode%in% clinic.info$Patient[clinic.info$Group=="tFL"])

###taking only drivers
flreseq.nFL.drivers<- flreseq.nFL%>%filter (Hugo_Symbol%in% drivers$SYMBOL)
flreseq.tFL.drivers<- flreseq.tFL%>%filter (Hugo_Symbol%in% drivers$SYMBOL)

#####


nFL.maf<- read.maf(flreseq.nFL.drivers)
tFL.maf<- read.maf(flreseq.tFL.drivers)
nFL.vs.tFL<- mafCompare(m1=nFL.maf, m2=tFL.maf, m1Name = "nFL group", m2Name = "tFL group", minMut=0)
print(nFL.vs.tFL)



drivers.sum.pat$SYMBOL<- reorder(drivers.sum.pat$SYMBOL, drivers.sum.pat$Patient)

trace(forestPlot, edit=TRUE)
###change
#m.sigs$Hugo_Symbol <- factor(m.sigs$Hugo_Symbol, levels = levels(drivers.sum.pat$SYMBOL))
#m.sigs = m.sigs[order(m.sigs$Hugo_Symbol)]

#m.sigs$or_new = ifelse(test = m.sigs$or > 3, yes = 3, no = m.sigs$or)
#m.sigs$upper = ifelse(test = m.sigs$ci.up > 3, yes = 3, 
#                      no = m.sigs$ci.up)
#m.sigs$lower = ifelse(test = m.sigs$ci.low > 3, yes = 3, 
#                      no = m.sigs$ci.low)

#####arrage

#if (is.null(fdr)) {
#  m.sigs = mafCompareRes$results[pval < pVal]
#}
#else {
#  m.sigs = mafCompareRes$results[adjPval < fdr]
#}

pdf(file="2021_codes/figures/51drivers_complete_clinic_20211125.pdf", width=13, height=8)
ggarrange(mutation,  clinic.plot.nolegend, nrow=2,   heights=c(8.5,1))

dev.off()

pdf(file="figures/51drivers_forest_20211130.pdf",width=6, height=10)
forestPlot(mafCompareRes = nFL.vs.tFL, color=c("#F1A340","#998EC3"),geneFontSize = 1,pVal = 1)
dev.off()



##################################
flreseq.ns<- unique(flreseq%>%filter(Region=="Coding" & VARIANT !="Silent")%>%select(Patient, SYMBOL))
clinic<- read.xlsx(file="2_clinic_infor/FL_Baoyan_clinical_file_FINAL_020719 - survival.xlsx",sheetIndex = 1,header=T)
fl.drivers<- subset(flreseq.ns, SYMBOL %in% drivers$Genes)

clinic$trans<- ifelse(is.na(clinic$Date_transf),"0","1")
patient<- clinic%>%select(Patient_ID, trans)
names(patient)<- c("Patient", "Trans_status")

fl.drivers.clinic<- merge(fl.drivers, patient, by="Patient")
fl.drivers.clinic$SYMBOL_ordered<- factor(fl.drivers.clinic$SYMBOL, levels=drivers$Genes)

fl.drivers.summary<-ddply(fl.drivers.clinic,.(SYMBOL, Trans_status), summarise, num_pat=length(Patient))

fl.drivers.sum.cast<- cast(fl.drivers.summary, SYMBOL~ Trans_status)
fl.drivers.sum.cast[is.na(fl.drivers.sum.cast)]<-0

fl.drivers.melt<- melt(fl.drivers.sum.cast, id.vars = 'SYMBOL')

fl.drivers.melt$value[fl.drivers.melt$Trans_status=="0"]<- fl.drivers.melt$value[fl.drivers.melt$Trans_status=="0"]*-1

fl.drivers.melt$SYMBOL<- factor(fl.drivers.melt$SYMBOL, levels=rev(drivers$Genes))





mutation_nolegend<-mutation<-ggplot(plot_df.drivers, aes(x=sample, y=gene,fill=value))+geom_tile(color="gray", size=0.1)+
  scale_fill_manual(values=c("#e7298a","firebrick","navyblue","#1b9e77", "#7570b3","#66a61e","white"),na.translate=FALSE)+scale_x_discrete(drop=FALSE)+
  theme(panel.background = element_blank(),
        legend.text = element_text(size=14,face="bold"),
        legend.title = element_text(size=16,face="bold"),
        axis.ticks.x=element_blank(),
        legend.position = "",
        axis.text.x =element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        
        plot.margin = unit(c(0.2,-0.2,-0.5,0.1), "cm"))+ xlab("")+ylab("")+
  
  guides(fill=guide_legend(title="Mutation Type"))

clinic.nolegend<- ggplot(clinic.final, aes(x=sample, y=variable, fill=value))+geom_tile()+
  scale_fill_manual(values=c("dark cyan","#fbb4ae", "purple","blue", "red"))+
  theme(panel.background = element_blank(),
        axis.ticks.x=element_blank(),
        legend.text = element_text(size=14,face="bold"),
        legend.title = element_text(size=16,face="bold"),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        legend.position="",
        axis.text.y =element_text(size=16,face="bold"),
        plot.margin = unit(c(-0.5,0.2,0,2,0.1), "cm"))+ xlab("")+ylab("")+
  guides(fill=guide_legend(ncol=2, title="Biopsy      Group"))










mutation_legend<-mutation<-ggplot(plot_df.drivers, aes(x=sample, y=gene,fill=value))+geom_tile(color="gray", size=0.1)+
  scale_fill_manual(values=c("#e7298a","firebrick","navyblue","#1b9e77", "#7570b3","#66a61e","white"),na.translate=FALSE)+scale_x_discrete(drop=FALSE)+
  theme(panel.background = element_blank(),
        legend.text = element_text(size=14,face="bold"),
        legend.title = element_text(size=16,face="bold"),
        axis.ticks.x=element_blank(),
        legend.position = "left",
        axis.text.x =element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        
        plot.margin = unit(c(0.2,-0.2,-0.5,0.1), "cm"))+ xlab("")+ylab("")+
  
  guides(fill=guide_legend(title="Mutation Type"))


clinic.legend<- ggplot(clinic.final, aes(x=sample, y=variable, fill=value))+geom_tile()+
  scale_fill_manual(values=c("dark cyan","#fbb4ae", "purple","blue", "red"))+
  theme(panel.background = element_blank(),
        axis.ticks.x=element_blank(),
        legend.text = element_text(size=14,face="bold"),
        legend.title = element_text(size=16,face="bold"),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        legend.position="left",
        axis.text.y =element_text(size=16,face="bold"),
        plot.margin = unit(c(-0.5,0.2,0,2,0.1), "cm"))+ xlab("")+ylab("")+
  guides(fill=guide_legend(ncol=2, title="Biopsy      Group"))


