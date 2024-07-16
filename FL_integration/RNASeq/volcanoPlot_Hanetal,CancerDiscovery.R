
setwd("~/UiO Dropbox/Ankush Sharma/LymphomaBiology/2024/MyklebustLab/FL_Multiomic_Revision/MichaelGreenPaper")
library(Seurat)
library(dplyr)
library(ggplot2)
library(EnhancedVolcano)
currentjob="Green_FL"
#Single Cell RNASeq Data from Han et al, Cancer Discovery article, GEO ID GSE203610_FL-gex-meta  
GreenFL1 <- readRDS("../Data/FL_MultiomicRevision/Data_MGreenPaper/Green_FL_CZIbioscience.rds")

#Extract gene names
gene_names <- GreenFL1@assays[["RNA"]]@meta.features[["feature_name"]]
head(feature_names)

# Ensure that the length of gene_names matches the number of rows in your counts matrix
if(length(gene_names) == nrow(GreenFL1@assays[["RNA"]]@counts)) {
  # Step 2: Replace Ensembl IDs with gene names in the counts matrix
  rownames(GreenFL1@assays[["RNA"]]@counts) <- gene_names
} else {
  stop("The number of gene names does not match the number of Ensembl IDs in counts.")
}
rownames(GreenFL1@assays[["RNA"]]@data) <- gene_names


#Define the group based on sample_id
#Corresponding Sample name "FL-013", "FL-044", "FL-030", "FL-039") mapped for samples (FL03,FL11,FL06,FL10) that are EZH2 Mutated in Metadata for HAN et al  manuscript GSE203610_FL-gex-meta 


# Subset to include only 'Malignant' cells based on the 'cell_type' metadata
GreenFL1_Malignant <- subset(GreenFL1, subset = cell_type == "malignant cell")
fl_group_ids <- c("FL-013", "FL-044", "FL-030", "FL-039")
GreenFL1_Malignant$Group <- ifelse(GreenFL1_Malignant$sample_id %in% fl_group_ids, "EZH2Mutated", "EZH2WT")
table(GreenFL1_Malignant$Group)

GreenFL1_Malignant <- SetIdent(GreenFL1_Malignant, value = "Group")

# Perform FindMarkers between activated Tregs and naive Tregs
de_results <- FindMarkers(
  GreenFL1_Malignant,
  ident.1 = "EZH2Mutated", ident.2 = "EZH2WT",
  min.pct = 0.1
)

write.table(de_results,"GreenPaper_MalignantCells_FL_EZH2Mutated_WT_updated20240508.csv")


Histone <- c("H3C8","H2AC17", "H2AC7", "H2AC8", "H2BC10", "H2BC5", "H3C4", "H2AC19", "H2BC4", "H3C1", "H2BC14", "H4C4", "H2BC9", "H4C13", "H4C16", "H1-2", "H1-3", "H2AC4", "H2BC13", "H2AC16", "H2BC11", "H2AC6", "H4C12", "H4C8", "H4C2", "H2AC18", "H3C13", "H2AC20", "H2AC15", "H1-5", "H2AC12", "H2BC8", "H4C1", "H2BC3", "H2BC17", "H2AC13", "H2AC11", "H4C3", "H3C12", "H2BC7", "H2BC12", "H4C6", "H2BC18", "H4C14", "H2AC21", "H2BC15", "HLA-DQB1", "HLA-DQB2", "HLA-DRA", "TCL1A", "HLA-DRB1", "HLA-DOB", "HLA-B", "HLA-F", "LGALS1", "FOS", "JUN", "HLA_DQA", "SOX5", "CLDN5", "H2AC25", "SYT11", "CLEC4A", "IL2RA", "SLC25A23", "HLA-DQA2", "ALKAL2", "LMO7", "HLA-DMB", "HLA-E","HLA-DMB")



  # Assigning basic colors based on significance and fold change
  keyvals <- ifelse(de_results$p_val_adj > 0.000005, "#434c54", 
                    ifelse(de_results$avg_log2FC < -0.2 , '#5870B3',
                           ifelse(de_results$avg_log2FC > 0.2, '#C86065', 'lightgray')))
  keyvals[is.na(keyvals)] <- '#292524' # Assigning 'darkgrey' to NA values
    
  # Updated gene names based on color coding
  names(keyvals)[keyvals == '#5870B3' & de_results$avg_log2FC < -0.2] <- 'downregulated'
  names(keyvals)[keyvals == '#C86065' & de_results$avg_log2FC > 0.2] <- 'upregulated'
  names(keyvals)[keyvals == '#434c54'] <- 'non significant'
  names(keyvals)[keyvals == 'lightgray'] <- 'significant'
  Green_Volc <- EnhancedVolcano(de_results,
                                lab = rownames(de_results),
                                x = 'avg_log2FC',
                                y = 'p_val_adj',
                                selectLab = Histone,
                                
                                xlim = c(-4,2),
                                xlab = bquote(~Log[2]~ 'fold change'),
                                title = 'EZH2 Mutated vs WT',
                                subtitle = 'Han et al.(2022) Cancer Discovery',
                                pCutoff = 0.0005,
                                FCcutoff = 0.2,
                                pointSize = 0.5,
                                labSize = 1.5,
                                axisLabSize=6,
                                
                                titleLabSize = 8,  # Smaller font size for the title
                                subtitleLabSize = 6,
                                colCustom = keyvals,
                                colAlpha = 1,
                                legendLabSize =6,
                                legendIconSize = 3.0,
                                legendPosition = 'right',
                                legendLabels=c('Not significant','Avg log2FC','Adj.Pval','Adj.Pval & Avg. log2FC'),
                                drawConnectors = TRUE,
                                widthConnectors = 0.1,
                                colConnectors = 'grey10',
                                gridlines.major = FALSE,
                                gridlines.minor = FALSE,
                                border = 'full',
                                borderWidth = 0.3,
                                borderColour = 'black')
  
Green_Volc
# Close the PDF device to finalize and save the file
dev.copy(pdf,paste0(currentjob, "_Volcano_EZH2Mutated_vsWT_MalignantCells_default_20240508_v1.pdf"), width=6, height=7,paper='special')
dev.off()
dev.set(dev.next())
saveRDS(GreenFL1,"GREENFL1_group_AS.rds") 

#########################
