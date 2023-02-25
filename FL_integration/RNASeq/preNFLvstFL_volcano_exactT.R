################################################
#  File Name:preNFL_tFLvolcano_exact.R
#  Author: Ankush Sharma

#################################################
#' Script using EdgeR exact test to perform differential expression and plots enhanced volcano plots for preNFL vs preTFL RNA Seq cohort for CibersortX imputed B cell transcriptome
#'
#

library(edgeR)
library(dplyr)
library(tidyverse)
library(factoextra)
library(tibble)
library(scales)
currentjob="prenFLvstFL_Exacttest_"

rawcounts<-read.csv("LM4_IMPUTED_BCELLS_FL_NORMAL_CIBERSORTxHiRes_NA_Ribosomal_variance_0_omited_protein_coding_Window20.txt", sep= "\t",header=T,row.names = 1)
counts=  subset(rawcounts, select = c(P1_1,P1_2,P2_1,P2_2,P4_1,P6_3,P8_2,P9_2,P10_2,P11_2,P12_1,P14_1,P14_2,P23_1,P27_1,P27_2,P28_1,P29_1,P29_2,P30_2,P31_1,P31_4,P33_1,P34_2,P35_2,P36_2,P40_1,P41_1,P42_1,P42_2,P43_1,P44_1,P44_2,P46_1,P46_3,P48_1,P49_1,P49_2,P50_1,P56_1,P56_2,P75_1,P11_1,P9_1,P3_1,P3_2,P46_2,P6_2,P57_1,P58_1,P10_1,P21_1,P21_2,P23_2,P23_3,P25_2,P27_3,P31_3,P32_1,P43_2,P59_1,P61_1,P75_2,P75_3))

#samples showing fraction less than 0.4 in cibersortx imputation
#counts1<- select(counts,-c("P30_2","P2_1","P58_1","P61_1","P11_2","P27_3","P44_2","P3_2","P25_2","P42_2","P41_1"))

experiment_design <-read.table("1_metadata_flcohort_mut_wt_driver_1_51.txt", sep="\t",header=TRUE, check.names=FALSE)


# Set the row names of the experiment design to be the sample names
rownames(experiment_design) <- experiment_design$Sample_name

# Subset the experiment design to only include the columns that are in the counts matrix
experiment_design.ord <- experiment_design[colnames(counts), ]

# Convert the sample names column in the experiment design to a character vector
samples <- as.character(experiment_design.ord$Sample_name)

# Create a factor variable based on the given gene
group<-factor(experiment_design.ord$Group)


# Filter the counts matrix to only include genes with at least 1 count in at least one sample and a minimum total count of 300
keep <- filterByExpr(counts, min.count = 1, min.total.count = 300, large.n = 52)

# Subset the counts matrix to only include the filtered genes
mycounts <- counts[keep, ]

# Load the limma library
library(limma)

# Create a design matrix for the differential expression analysis
design <- model.matrix(~0 + group)

# Remove the "group" prefix from the column names of the design matrix
colnames(design) <- gsub("group", "", colnames(design))

# Create a DGEList object for the differential expression analysis
y <- DGEList(counts = mycounts, group = group, lib.size = colSums(mycounts))

# Normalize the count data in the DGEList object
y_et <- calcNormFactors(y)

# Estimate the dispersion in the normalized count data
y_et <- estimateDisp(y_et)

# Perform the exact test for differential expression
et <- exactTest(y_et, pair = c(2, 1))

# Convert the exact test results to a data frame
et <- as.data.frame(et)


filename = paste("filtrbyexpr_stringentremove1Classicmode_qcml_exacttest_tFlvsnFL_wt",".txt",sep="")
write.table(et,file=filename,row.names=T,quote=T,sep="\t")

###Volcano Plot######
  # Convert the rownames of the exact test results to a column in the data frame
  resdata <- rownames_to_column(et, var = "GeneSymbol")
  
  # Subset the exact test results to only include genes with "HIST" in their name
  Hist <- resdata %>% filter(grepl("HIST", resdata$GeneSymbol))

# Filter for genes with significant log fold changes
  Hist.SIG <- Hist %>% filter(PValue < 0.05) %>% filter(logFC < -0.3 | logFC > 0.3)
  
  # Filter for genes with significant log fold changes greater than 0.5
  Hist.label <- Hist %>% filter(PValue < 0.05) %>% filter(logFC < -0.5 | logFC > 0.5)
  
  # Read in a list of approved gene symbols
  match_name <- read.csv(file = "../../drivergenes_de/LM4_Signature/FL_DifferentialExp/Baoyan_volcano/hgnc-symbol-check.csv", sep = ";", header = TRUE)
  
  # Combine the list of approved gene symbols with a list of predefined genes
  genes <- c("HRK","PCDH9","FCGR2B","LMO2","LIMS1","PTPN22","SKAP2","PRKCD", as.character(match_name$Approved.symbol[match_name$Input == Hist.SIG$GeneSymbol]))
  
  # Replace the gene symbols in the result data with the approved symbols
  resdata$GeneSymbol[match(match_name$Input, resdata$GeneSymbol)] <- match_name$Approved.symbol
  
resdata <- rownames_to_column(et, var = "GeneSymbol")
library(EnhancedVolcano)
keyvals <- ifelse(resdata$PValue> 0.005, "grey", 
                  ifelse( resdata$logFC < -0.3 , 'royalblue',
                          ifelse(resdata$logFC > 0.3, 'red', 'grey')))
keyvals[is.na(keyvals)] <-  'grey'
  names(keyvals)[keyvals == 'red'] <- 'upregulated'
  names(keyvals)[keyvals == 'grey'] <- 'non significant'
  names(keyvals)[keyvals == 'royalblue'] <- 'downregulated'
  
  EnhancedVolcano(resdata,
                  lab = as.character(resdata$GeneSymbol),
                  x = 'logFC',
                  y = 'PValue',
                  #selectLab = genes,
                  selectLab = as.character(resdata$GeneSymbol)[which(names(keyvals) %in% c('upregulated', 'downregulated'))],
                  xlim = c(-0.8,1),
                  ylim=c(1,6),
                  xlab = bquote(~Log[2]~ 'fold change'),
                  title = "prenFLvstFL: Wt",
                  pCutoff = 0.005,
                  FCcutoff = 0.3,
                  pointSize = 2,
                  labSize = 5,
                  colCustom = keyvals,
                  colAlpha = 1,
                  legendPosition = 'right',
                  legendLabels=c('Not significant','log2-fold change','PValue','PValue & log2-fold change'),
                  legendLabSize = 12,
                  legendIconSize = 5.0,
                  drawConnectors = TRUE,
                  widthConnectors = 0.20,
                  #colConnectors = 'grey30',
                  gridlines.major = FALSE,
                  gridlines.minor = FALSE,
                  border = 'full',
                  borderWidth = 1.0,
                  borderColour = 'black')
  dev.copy(pdf,paste0(currentjob, "_20230122_stringent_ExactTest_prenFLvstFL_EnhancedVolcano.pdf"), width=10, height=12,paper='special')
  dev.off()
  saveRDS(et, file="prenFLvstFL_QCML.RDS")