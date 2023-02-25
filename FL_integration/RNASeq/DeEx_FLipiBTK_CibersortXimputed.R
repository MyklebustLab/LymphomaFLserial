#################################################
#' Function "analyze_gene" to perform differential expression using EdgeR exact test and functionmakeVolcanoPlot plot enhanced volcano plots for Validation RNA Seq cohort for CibersortX imputed B cell transcriptome
#'
#



library(edgeR)
library(dplyr)
library(tidyverse)
library(factoextra)
library(tibble)
library(scales)
currentjob="_Exacttest_"

analyze_gene <- function(gene) {
  rawcounts<-read.csv("~/flipi_btk_validation/Annotation/LM4_flipi_btk_gencode_annotated_CIBERSORTxHiRes_NA_Ribosomal_variance_0_omited_Window20.txt", sep= "\t",header=T,row.names = 1)
  write.csv(colnames(rawcounts), "samplesnamesflipibtk.csv", sep="\t")
  
  counts=  subset(rawcounts, select = c(II_01.CEL,II_16.CEL,II_2.CEL,II_21.CEL,II_24_A.CEL,II_46.CEL,II_48.CEL,II_55.CEL,II_57.CEL,II_68.CEL,II_74.CEL,II_79.CEL,II_8.CEL,II_86_A.CEL,II_97.CEL,III_1003.CEL,III_1012.CEL,III_1019.CEL,III_1028.CEL,III_2001.CEL,III_2002.CEL,III_2005.CEL,III_2019.CEL,III_2030.CEL,III_2033.CEL,III_2034.CEL,III_2039.CEL,III_2048.CEL,III_2049.CEL,III_2052.CEL,III_2053.CEL,III_2060.CEL,III_2061.CEL,III_2064_A.CEL,III_2069_A.CEL,III_2070.CEL,III_2076.CEL,III_2079.CEL,III_2082.CEL,III_2085.CEL,III_2093.CEL,III_2097.CEL,III_2098.CEL,III_2103.CEL,III_2110.CEL,III_2111.CEL,III_2113.CEL,III_2114.CEL,III_2121.CEL,III_2128.CEL,III_2129.CEL,III_2130.CEL,III_2135.CEL,III_2138.CEL,III_2140.CEL,III_2142.CEL,III_2144.CEL,III_2145.CEL,III_2146.CEL,III_2149.CEL,III_2152.CEL,III_5006_A.CEL,III_5006_B.CEL,III_5007.CEL,III_6002.CEL,III_6034_A.CEL,III_6034_B.CEL,III_6043.CEL,III_6062.CEL,III_6074.CEL,III_6082.CEL,III_6087.CEL))
  
  filename = paste("sample",".txt",sep="")
  write.table(out,file=filename,row.names=T,quote=T,sep="\t")
  
  #samples showing fraction less than 0.4 in cibersortx imputation
  #counts1<- select(counts,-c("P30_2","P2_1","P58_1","P61_1","P11_2","P27_3","P44_2","P3_2","P25_2","P42_2","P41_1"))
  #experiment_design <-read.table("/Users/ankushs/Dropbox (UiO)/LymphomaBiology/FL1/FL1_RNASEQ/0_INPUTFILES/2.DOWNSTREAM_INPUT/Metadata/1_metadata_flcohort_mut_unmut_driver_1_51.txt", sep="\t",header=TRUE, check.names=FALSE)
  experiment_design <-read.table("~/LymphomaBiology/FL1/figures_baietal/ChloeBTS_FLIPI_Validation/METADATA/EZH2_MUT_UNMUT_METADATA.TXT", sep="\t",header=TRUE, check.names=FALSE)
  

# Set the row names of the experiment design to be the sample names
rownames(experiment_design) <- experiment_design$Sample_name

# Subset the experiment design to only include the columns that are in the counts matrix
experiment_design.ord <- experiment_design[colnames(counts), ]

# Convert the sample names column in the experiment design to a character vector
samples <- as.character(experiment_design.ord$Sample_name)

# Create a factor variable based on the given gene
group <- factor(experiment_design.ord[, gene])

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

# Create a filename for the exact test results based on the given gene
filename <- paste(gene, "_exact_test_results.txt", sep = "")

# Write the exact test results to a tab-delimited file
write.table(et, file = filename, row.names = TRUE, quote = TRUE, sep = "\t")  


############
############


# Function for creating a volcano plot
makeVolcanoPlot <- function(gene) {
  # Read in the exact test results for the given gene
  et <- read.table(paste0(gene, "_exact_test_results.txt"))
  
  # Load the EnhancedVolcano library
  library(EnhancedVolcano)
  
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



analyze_gene("EZH2")
makeVolcanoPlot("EZH2") 

