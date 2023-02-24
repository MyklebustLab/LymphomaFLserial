################################################
#  File Name:determine_allCNbreaks.r
#  Author: jillianwise
#  Mail: jwise7@mgh.harvard.edu
#  Created Time: Mon 22 August 2022 9:25 AM EST
#################################################
#' Given multiple samples determine the CN state at all site breaks present in the dataset
#'
#' Given a data frame of CNV segments of DNA with rows representing each unique segment and columns representing the CN states and sample information break each sample into all available segments
#' @chr column of a dataframe that lists the numeric (please reformat X and Y to 23 and 24) chromosome number of the alignment
#' @start column of a dataframe representing the start site of a segment
#' @end column of a dataframe representing the end site of a segment
#' @SegCN The CN state in Numeric CN, either integer or summed raw copy states of both alleles
#' @sample the sample name used in segmentation and expression files
#' @case the case name used for identifying multiple related samples
#' @ploidy the estimated ploidy from any number of CN tools (ASCAT, facets, etc)
#' Given an a directory to write the dataframe to


determine_allCNbreaks <- function(df,outdir) {
  print("converting CN to log2-1 ploidy adjusted states")
  ###check numeric SegCN
df$SegCN <- as.numeric(df$SegCN)
  ###Now create SegCNploidyadjusted
df$tumor_ploidyc <- round(df$ploidy-2)
df$SegCNa <- (log2(df$SegCN-df$tumor_ploidyc)-1)
  ###Lets replace -inf or other negative errors with -10 a
df$SegCNa <- ifelse(df$SegCNa=="NaN", "-10", paste(df$SegCNa)) 
df$SegCNa <- ifelse(df$SegCNa=="-Inf", "-10", paste(df$SegCNa))
df$SegCNa <- as.numeric(df$SegCNa)

#Make list of every unique start and end for every chr, regardless of case, so one set of segments for all cases to be compared across#
df2 <- NULL
print("Finding All Break Points")
for (chr in unique(df$chr)) {
  print(chr)
  x=unique(df$start[df$chr == chr])
  y=unique(df$end[df$chr == chr])
  pos=sort(c(x,y))
  df2 = rbind(df2, data.frame(chr=rep(chr,length(pos)), pos=pos))
}

df2$end <- dplyr::lead(df2$pos)
for(i in 1:nrow(df2)){
  df2$end[i] <- ifelse(df2$chr[i] != df2$chr[i+1], 0, df2$end[i])}
df2 <- filter(df2, df2$end !=0)

###For each patient make a column of their CN state in that segment
require(GenomicRanges)

all <- makeGRangesFromDataFrame(df2,
                                keep.extra.columns = TRUE,
                                ignore.strand = FALSE,
                                start.field="pos",
                                end.field="end", 
                                seqnames = "chr")
###################
df3 <- list()
print("Making GRanges Object")
for (i in unique(df$sample)){
  print(i)
  tmp <- makeGRangesFromDataFrame(df[df$sample == i,],
                                  keep.extra.columns=TRUE,
                                  ignore.strand=FALSE,
                                  start.field="start",
                                  end.field="end", 
                                  seqnames = "chr") 
  df3[[i]]=tmp
}

colnames <- names(df3)
df2$idx <- rownames(df2)
df2$idx <- as.numeric(df2$idx)

#####Run granges overlap between the total segment site list and each patient segmentation granges object
tmp2 <- NULL
tmp <- NULL
print("Finding Copy State For Each Break")
for (i in (1:length(df3))){
  print(i)
  over <- findOverlaps(all, df3[[i]], minoverlap=2)
  tmp <- all[queryHits(over)]
  tmp$SegCNa <- df3[[i]]$SegCNa[subjectHits(over)]
  tmp$sample <- df3[[i]]$sample[subjectHits(over)]
  tmp <- as.data.frame(tmp, row.names = NULL)
  #####Next is to index which results segments match the total segmentation list
  for (j in (1:length(tmp$seqnames))){
    tmp2[j] <- which(tmp$seqnames[j]==df2$chr & tmp$start[j]==df2$pos & tmp$end[j]==df2$end)
  }
  tmp2 <- as.data.frame(tmp2)
  tmp2$SegCNa <- as.numeric(tmp$SegCNa)
  tmp2 <- tmp2 %>% arrange(tmp2)
  #####For each sample make a column and for the index matches add in the correct CNstate
  for (z in (1:length(tmp2$tmp2))){
    df2[tmp2$tmp2[z],paste(colnames[i],sep="")] <- tmp2$SegCNa[z]}
  tmp <- NULL
  tmp2 <- NULL
}
write.table(df2, file = paste(outdir,"allCNbreaks.txt",sep=""), append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", col.names = TRUE, row.names = FALSE)
return(df2)
}



