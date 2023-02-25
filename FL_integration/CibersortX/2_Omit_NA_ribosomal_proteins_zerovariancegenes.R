
################################################
#  File Name:2_Omit_NA_proteins_zerovariancegenes.R
#  Author: Ankush Sharma

#################################################
#' Remove NA and genes Not imputed in Cibersortx Imputed BCell Trancriptome gene expression based om LM4 Merged class signature matrix
#'

fl_exp<-read.csv("../CIbersortX_local_Hires/LM4_FL_NORMAL_PROTEIN_CODING/CIBERSORTxHiRes_LM4_protein_coding_fl_Normal_Bcells_Window20.txt", "\t",header=T,row.names=1)

data<-fl_exp
vf = read.table("../0_INPUTFILES/2.DOWNSTREAM_INPUT/BCELL_SPECIFIC_GENE/listofribosomalproteins_toremove.txt", header=T, sep="\t")

remove = vf$ribosomalproteins # list of rownames I would like to remove from file "data"

head(remove)
datawithoutVF = data[-which(rownames(data) %in% remove), ]

dim(datawithoutVF)
datawithoutVF<-na.omit(datawithoutVF)
fl_exp<-datawithoutVF[apply(datawithoutVF, 1, var) != 0, ]
dim (fl_exp)
write.table(fl_exp,"../0_INPUTFILES/2.DOWNSTREAM_INPUT/Cibersortx/LM4_IMPUTED_BCELLS_FL_NORMAL_CIBERSORTxHiRes_NA_Ribosomal_variance_0_omited_protein_coding_Window20.txt",sep="\t",quote = FALSE)
