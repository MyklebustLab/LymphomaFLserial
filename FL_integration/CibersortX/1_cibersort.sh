################################################
#  File Name:1_cibersort.sh
#  Author: Ankush Sharma

#################################################
#' Script to run Docker CibersortX to impute fraction and Gene expression
#'

#!/bin/sh
#Computes Fraction of celltypes on LM22 signature matrix
docker run -v  ~/Desktop/CIBERSORTX/inputfiles/salmon:/src/data -v  ~/Desktop/CIBERSORTX/output/salmon_lm4_hires:/src/outdir cibersortx/fractions --username email@email.com --token SecretKey --mixture inputfiles/salmon/FL_cibersort_fraction_input.txt --sigmatrix LM22_signature_matrix.txt --label all_transcripts_fl --fraction 0 --rmbatchBmode TRUE
# imputed gene expression on merged LM4 Classes 
docker run -v  ~/Desktop/CIBERSORTX/inputfiles/salmon:/src/data -v  ~/Desktop/CIBERSORTX/output/salmon_lm4_hires:/src/outdir cibersortx/hires --username email@email.com --token SecretKey --mixture ~/Desktop/CIBERSORTX/inputfiles/FL_cibersort_hires_input_protein_coding.txt --sigmatrix LM22_signature_matrix.txt --label protein_coding_fl --classes merged_classes_LM8.txt --window 20 --rmbatchBmode TRUE

