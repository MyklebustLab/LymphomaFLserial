################################################
#  File Name:DeEx_FLCohort_CibersortXimputed.R
#  Author: Ankush Sharma


#################################################
#' Function "analyze_gene" to perform differential expression using EdgeR exact test and functionmakeVolcanoPlot plot enhanced volcano plots for FL RNA Seq cohort for CibersortX imputed B cell transcriptome
#'
#


 Load required libraries
library(edgeR)
library(dplyr)
library(tidyverse)
library(factoextra)
library(tibble)
library(scales)
library(limma)
library(EnhancedVolcano)

# Define function to analyze a gene
analyze_gene <- function(gene) {
  

  rawcounts<-read.csv("/LymphomaBiology/FL1/FL1_RNASEQ/0_INPUTFILES/2.DOWNSTREAM_INPUT/Cibersortx/LM4_IMPUTED_BCELLS_FL_NORMAL_CIBERSORTxHiRes_NA_Ribosomal_variance_0_omited_protein_coding_Window20.txt", sep= "\t",header=T,row.names = 1)
# Subset the count data to include only specific columns
  counts=  subset(rawcounts, select = c(P1_1,P1_2,P2_1,P2_2,P4_1,P6_3,P8_2,P9_2,P10_2,P11_2,P12_1,P14_1,P14_2,P23_1,P27_1,P27_2,P28_1,P29_1,P29_2,P30_2,P31_1,P31_4,P33_1,P34_2,P35_2,P36_2,P40_1,P41_1,P42_1,P42_2,P43_1,P44_1,P44_2,P46_1,P46_3,P48_1,P49_1,P49_2,P50_1,P56_1,P56_2,P75_1,P11_1,P9_1,P3_1,P3_2,P46_2,P6_2,P57_1,P58_1,P10_1,P21_1,P21_2,P23_2,P23_3,P25_2,P27_3,P31_3,P32_1,P43_2,P59_1,P61_1,P75_2,P75_3))  
# Read in experiment design information 
  experiment_design <- read.table("/LymphomaBiology/FL1/FL1_RNASEQ/0_INPUTFILES/2.DOWNSTREAM_INPUT/Metadata/1_metadata_flcohort_mut_wt_driver_1_51.txt", sep="\t",header=TRUE, check.names=FALSE)
  
  # Set the row names of the experiment design data to the sample names
  rownames(experiment_design) <- experiment_design$Sample_name
  
  # Subset the experiment design data to include only relevant columns
  experiment_design.ord <- experiment_design[colnames(counts), ]
  
  samples <- as.character(experiment_design.ord$Sample_name)
  group <- factor(experiment_design.ord[, gene])
    
#Filter the count data based on certain criteria
  keep <- filterByExpr(counts, min.count = 1, min.total.count = 300, large.n = 52)
  mycounts <- counts[keep, ]

#Load the limma library
library(limma)

#Create a design matrix
  design <- model.matrix(~0 + group)

#Remove the "group" prefix from the column names of the design matrix
  colnames(design) <- gsub("group", "", colnames(design))

#Create a DGEList object
  y <- DGEList(counts = mycounts, group = group, lib.size = colSums(mycounts))

#Normalize the count data
  y_et <- calcNormFactors(y)
  y_et <- estimateDisp(y_et)

#Perform the exact test
  et <- exactTest(y_et, pair = c(2, 1))

#Convert the results to a data frame
  et <- as.data.frame(et)

#Create a file name for the results
  filename <- paste(gene, "_exact_test_results.txt", sep = "")

#Write the results to a tab-delimited file
  write.table(et, file = filename, row.names = TRUE, quote = TRUE, sep = "\t")
      }
                 

##############
##############
# Function to create a volcano plot
makeVolcanoPlot <- function(gene) {
  
  # Read in the exact test results from file
  et <- read.table(paste0(gene, "_exact_test_results.txt"))
  
  # Load the EnhancedVolcano library
  library(EnhancedVolcano)
  
  # Convert row names to a column in the data frame
  resdata <- rownames_to_column(et, var = "GeneSymbol")
  
  # Filter for genes containing "HIST" in their name
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
  
  # Assign colors based on significance and log fold change
  keyvals <- ifelse(resdata$PValue > 0.005, "grey",
                    ifelse(resdata$logFC < -0.3, 'royalblue',
                           ifelse(resdata$logFC > 0.3, 'red', 'grey')))
  keyvals[is.na(keyvals)] <-  'grey'
  names(keyvals)[keyvals == 'red'] <- 'higher in Mutated'
  names(keyvals)[keyvals == 'grey'] <- 'Mid'
  names(keyvals)[keyvals == 'royalblue'] <- 'higher in WT'
  
  # Save the volcano plot as a PNG file
  ggsave(file = paste0(getwd(), "/", gene, "_20230122_stringent_ExactTest_mutvswT_Enhanced_Volcano.png"), 
         width = 10, height = 12, dpi = 300)
}
###################
###################
# Function to create a volcano plot
makeVolcanoPlot <- function(gene) {
  
  # Read in the exact test results from file
  et <- read.table(paste0(gene, "_exact_test_results.txt"))
  
  # Load the EnhancedVolcano library
  library(EnhancedVolcano)
  
  # Convert row names to a column in the data frame
  resdata <- rownames_to_column(et, var = "GeneSymbol")
  
  # Filter for genes containing "HIST" in their name
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
  
  # Assign colors based on significance and log fold change
  keyvals <- ifelse(resdata$PValue > 0.005, "grey",
                    ifelse(resdata$logFC < -0.3, 'royalblue',
                           ifelse(resdata$logFC > 0.3, 'red', 'grey')))
  keyvals[is.na(keyvals)] <-  'grey'
  names(keyvals)[keyvals == 'red'] <- 'higher in Mutated'
  names(keyvals)[keyvals == 'grey'] <- 'Mid'
  names(keyvals)[keyvals == 'royalblue'] <- 'higher in WT'
  
  # Save the volcano plot as a PNG file
  ggsave(file = paste0(getwd(), "/", gene, "_20230122_stringent_ExactTest_mutvswT_Enhanced_Volcano.png"), 
         width = 10, height = 12, dpi = 300)


  EnhancedVolcano(resdata,
                    lab = as.character(resdata$GeneSymbol),
                    x = 'logFC',
                    y = 'PValue',
                    selectLab = as.character(resdata$GeneSymbol)[which(names(keyvals) %in% c('higher in Mutated', 'higher in WT'))],
                    xlim = c(-0.8,1),
                    ylim = c(1,6),
                    xlab = bquote(~Log[2]~ 'fold change'),
                    title = paste("Mutated vs Wt"),
                    pCutoff = 0.05,
                    FCcutoff = 0.3,
                    pointSize = 2,
                    labSize = 5,
                    colCustom = keyvals,
                    colAlpha = 1,
                    legendPosition = 'right',
                    legendLabels = c('Not significant','log2-fold change','PValue','PValue & log2-fold change'),
                    legendLabSize = 12,
                    legendIconSize = 5.0,
                    drawConnectors = TRUE,
                    widthConnectors = 0.20,
                    gridlines.major = FALSE,
                    gridlines.minor = FALSE,
                    border = 'full',
                    borderWidth = 1.0,
                    borderColour = 'black')
     dev.off() 
     dev.set(dev.next())
} 
    
  
  analyze_gene("KMT2D")
  makeVolcanoPlot("KMT2D") 

  analyze_gene("EZH2")
  makeVolcanoPlot("EZH2") 

  analyze_gene("CREBBP")
  makeVolcanoPlot("CREBBP") 
  analyze_gene("EP300")
  makeVolcanoPlot("EP300") 
  