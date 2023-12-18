

################################################
#  File Name:Heatmap_Figure5C.R
#  Author: Baoyan Bai


#################################################

###plot heatmap for de-histone genes identified with the new filter (Figure 5C)
###load required libraries
###expression data: output from CIBERSORTx


library(pheatmap)
library(dplyr)
library(plyr)
library(RColorBrewer)
setwd("~/Desktop/FL/Ankush_CIBERSORT")

###import expression data and remove samples not relative
match_names<- read.csv(file="hgnc-symbol-check.csv",sep=";", header=T)
flciber<- read.csv(file="CIBERSORTxHiRes_LM4_protein_coding_fl_Normal_Bcells_Window20.txt",header=T, sep="\t", stringsAsFactors = F)
samples.exp<- flciber[!grepl("BC", names(flciber))]
samples.exp.76<- select(samples.exp,-c("JP6_2", "P34_1", "P35_1"))

samples.hist.76<- subset(samples.exp.76, Gene %in% match_names$Input)


##update gene names


samples.hist.76$Gene[match(match_names$Input, samples.hist.76$Gene)]<- match_names$Approved.symbol



hist.new.matrix<- samples.hist.76[,2:77]
row.names(hist.new.matrix)<- samples.hist.76$Gene

###clinical information
clinic.anno.clean<- read.table(file="clinic_infor.text",sep="\t")


###correct P23 POD24 status
clinic.anno.clean$POD24[rownames(clinic.anno.clean)=="P23_1"]<- "No" 
clinic.anno.clean$POD24[rownames(clinic.anno.clean)=="P23_2"]<- "No" 
clinic.anno.clean$POD24[rownames(clinic.anno.clean)=="P23_3"]<- "No" 

aka3 = list(POD24 = c(No = "grey", Yes="black"),
            Type= c(pre = "dark cyan", relapse="#fbb4ae", transform="purple", Naive= "#7570B3",GC="#E7298A", Memory="#66A61E"),
            
            Group =c(nFL="#998EC3",tFL= "#F1A340", Tonsil="white"),
            EZH2=c(Wt = "grey", Mut="black"))

brewer.pal(n = 8, name = "Dark2")

#####new colors
paletteLength <- 100
mycolor<- colorRampPalette(c("blue","white",  "yellow"))(7)
pdf(file="hist_heatmap_20230912_scale_color.pdf", width=10, height=6.5)
pheatmap(log(hist.new.matrix), annotation_col = clinic.anno.clean, annotation_colors   = aka3, 
         fontsize = 6, show_colnames = F, scale ="row",
         color = mycolor)

dev.off()

###Add mutation status of histone genes

flreseq<- read.csv(file="../fl_latest/16_combined/fl44_reseq_final_060619.txt", sep="\t",header=T, stringsAsFactors = F)
flreseq.76<-subset(flreseq, SAMPLE %in% names(hist.new.matrix))
flreseq.76.hist<- subset(flreseq.76, SYMBOL %in% match_names$Input)
flreseq.76.hist.coding<- flreseq.76.hist %>% filter(Region=="Coding")%>% filter(VARIANT!="Silent")
summ.1<-ddply(flreseq.64.4gene, .(SAMPLE, SYMBOL), summarise, n=length(SYMBOL))
summ.cast<-cast(summ.1, SAMPLE~SYMBOL, value="n")



