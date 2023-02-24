################################################
#  File Name:find_earlystable.r
#  Author: jillianwise
#  Mail: jwise7@mgh.harvard.edu
#  Created Time: Mon 22 August 2022 9:25 AM EST
#################################################
#' Given a list of CN states at all avaliable site breaks determine which CN States are stable throughout sampling/time
#'
#' Given a data frame of CNV segments in which all breaks are reported for all cases (see determine_allCNbreaks.R) and a sample information file- output the early stable CN regions for each case
#' allCNbreaks file
#' @chr column of a dataframe that lists the numeric (please reformat X and Y to 23 and 24) chromosome number of the alignment
#' @start column of a dataframe representing the start site of a segment
#' @end column of a dataframe representing the end site of a segment
#' @idx determined from determine_allCNbreaks.R this gives a unique id to each break region for each sample
#' The following columns should be the SegCN adjusted for ploidy values for every sample as a column (log2(df$SegCN-round(df$ploidy-2))-1)
#' The second file should represent a sample manifest sheet
#' @sample each samples name, matching those used to find all CN breaks
#' @case the case name used for identifying multiple related samples
#' @ploidy the estimated ploidy from any number of CN tools (ASCAT, facets, etc)
#' Given an output directory
#' Given the choice of modeling for Stable regions: 
#' plusminus allows for CN states to range between +- a threshold typically set at 0.3 while allowing for a copy state abberation (+-0.3) at initial biopsy to further increase or decrease.
#' mahalanobis uses a calculated distribution for each sample and determines the distance of each point from that distribution
#' linearpc1 calculates distances from the first principal component
#' Given a threshold (thresh) for plusminus modeling (default=0.3)
#' Given an option for returning a merged/unionized case segmentation file representing the early stable copy states


find_earlystable <- function(df,names,outdir,change="plusminus",thresh=0.3,return_segmentfile=FALSE) {
  require(dplyr)
  df2 <- subset(df, select = c("chr", "pos", "end", "idx"))
  df2[,unique(names$case)] <- NA
  if(change != "plusminus"){
    ####I report ploidy adjust CN as (log2(df$SegCN-round(df$ploidy-2))-1)
    ####To run Mahalanobis and PC1 We need total CN
    tmp <- df[,-c(1:4)]
    tmp <- (2^(tmp+1))
    tmp3 <- tmp 
    tmp3[,unique(names$sample)] <- NA
    for (i in colnames(tmp)){
      tmp2 <- tmp[,i]
      tmp2 <- tmp2+(names[names$sample==i,3]-2)
      tmp3[,i] <- tmp2
    }
    tmp3 <- cbind(df[,c(1:4)],tmp3)
    df <- tmp3
  }
  
  ######Check overlap between names files and CN file
  print("Checking Sample Intersection")
  rem <- setdiff(names$sample,colnames(df))
  if (length(rem)!=0){
    print(paste("removing ",rem, " from sample file",sep=""))
    names <- dplyr::filter(names, names$sample != rem)
  }
  
  ########################
  ##Check for any columns with all NAs
  len <- dim(df)[1]
  df <- df[,colSums(is.na(df)) != len]
  ########################
  
  #########################
  ###Remove any unpaired samples
  print("Checking for Unpaired Samples")
  tmpnames <- dplyr::filter(names, names$sample %in% colnames(df))
  tmp <- tmpnames %>% dplyr::group_by(case) %>% dplyr::count()
  rem2 <- as.data.frame(tmp[tmp$n < 2,1])
  print(paste("Sample", rem2$case," has no Pairs, Removing"))
  rem3 <- as.data.frame(names[names$case==rem2$case,1])
  colnames(rem3) <- c("sample")
  drop <- c(rem3$sample)
  df <- df[,!(names(df) %in% drop)]
  
  #####################
  ##Update names after sample removal
  names <- dplyr::filter(names, names$sample %in% colnames(df))
  
  ###For each case (based on the first biopsy) make a loss dataframe and gain (If one uses absolute value a switch to gain or loss would be counted)
  if(change == "plusminus"){
    if (is.na(plot)) {
      plot = F
    } else {
      fname = plot
      plot = T
    }
  for (i in unique(names$case)){
    print(i)
    namestmp <- dplyr::filter(names, names$case == i)  
    tmp <- unique(namestmp$sample)
    tmp <- c(tmp)
    tmpdf2 <- df %>% dplyr::select(tmp)
    tmpdf3 <- tmpdf2 %>% dplyr::filter(tmpdf2[,1] > thresh) #This dataframe is all gains above .3 in first biopsy
    tmpdf4 <- tmpdf2 %>% dplyr::filter(tmpdf2[,1] < -(thresh)) #This dataframe is all losses below  -.3 in first biopsy
    greatergain <- function(x) all((x >= x[1]))
    greaterloss <- function(x) all((x <= x[1]))
    plusminus <- function(x) all((abs(x - x[1]) < thresh))
    tmpdf2[, "Consensus"] <- apply(tmpdf2, 1, plusminus)
    tmpdf3[, "GreaterGains"] <- apply(tmpdf3, 1, greatergain)
    tmpdf4[, "GreaterLoss"] <- apply(tmpdf4, 1, greaterloss)
    ###add back in the idx
    tmpdf2$idx <- rownames(tmpdf2)
    tmpdf5 <- tmpdf2 %>% dplyr::filter(tmpdf2[idx,1] > thresh) %>% select(idx)
    tmpdf3$idx <-tmpdf5$idx
    tmpdf6 <- tmpdf2 %>% dplyr::filter(tmpdf2[idx,1] < -(thresh)) %>% select(idx)
    tmpdf4$idx <-tmpdf6$idx
    ####plot
    pldf <- full_join(tmpdf2,tmpdf3)
    pldf <- full_join(pldf,tmpdf4)
    pldf$color <- ifelse(pldf$Consensus==TRUE,"black",ifelse(pldf$GreaterGains==TRUE & pldf[,1] > thresh ,"black",ifelse(pldf$GreaterLoss==TRUE & pldf[,1] > -(thresh),"black","red")))
    pldf$color <- ifelse(is.na(pldf$color),"red",pldf$color)
    if (plot) {
      pdf(paste(outdir,i,"plusmins.pdf"))
    }
    if (plot) {
      plot(pldf[,1], pldf[,2],xlim = c(-2, 4),ylim = c(-2, 4),ann=FALSE)
      abline(h=0); abline(v=0)
      abline(coef=c(thresh,1), col="green", lty=2)
      abline(coef=c(-(thresh),1), col="green", lty=2)
      points(pldf[,1], pldf[,2], col=pldf[,"color"])
      title(main=paste("Case: ", i, sep=""),xlab=paste(colnames(pldf)[1]," Copy State", sep=""),ylab=paste(colnames(pldf)[2]," Copy State", sep=""))
      
    }
    if (plot) {
      dev.off()
    }
    df2[,i] <- ifelse(tmpdf2$Consensus == TRUE, tmpdf2[,1],NA)
    for (j in tmpdf3$idx){
      df2[j,i] <- ifelse(tmpdf3$GreaterGains[tmpdf3$idx==j] == TRUE, tmpdf3[as.numeric(which(tmpdf3$idx==j)),1],df2[j,i])}
    for (k in tmpdf4$idx){
      df2[k,i] <- ifelse(tmpdf4$GreaterLoss[tmpdf4$idx==k] == TRUE,  tmpdf4[as.numeric(which(tmpdf4$idx==k)),1],df2[k,i])}
    tmpdf2 <- NULL
    tmpdf3 <- NULL
    tmpdf4 <- NULL
    tmp <- NULL
  }
  write.table(df2, file = paste(outdir,"earlystableregions.txt",sep=""), append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", col.names = TRUE, row.names = FALSE)
  return(df2)
  }
  
  if(change == "mahalanobis"){
    require(tidyr)
    require(stats)
    require(rlang)
    require(stringr)
    if (is.na(plot)) {
      plot = F
    } else {
      fname = plot
      plot = T
    }
    ###############
    #####get_selected_vars and mahalanobis_Distance is written by Alboukadel Kassambra and from the rstatix package (version 0.7.0 )
    ###########
    get_selected_vars <- function(x, ..., vars = NULL){
      if(is_grouped_df(x))
        x <- x %>% dplyr::ungroup()
      dot.vars <- rlang::quos(...)
      if(length(vars) > 0){
        return(vars)
      }
      if (length(dot.vars) == 0) selected <- colnames(x)
      else selected <- tidyselect::vars_select(names(x), !!! dot.vars)
      selected %>% as.character()
    }
    mahalanobis_distance <- function(data, ...){
      if(is_grouped_df(data)){
        results <- data %>%
          doo(~mahalanobis_distance(.))
      }
      data <- data %>% select_if(is.numeric)
      data <- data %>% drop_na()
      vars <- data %>% get_selected_vars(...)
      n.vars <- length(vars)
      threshold <- stats::qchisq(0.999, n.vars)
      .data <- data %>% select(!!!syms(vars)) %>% as.matrix()
      distance <- stats::mahalanobis(.data,center = colMeans(.data),cov = cov(.data))
      results <- data %>% mutate(mahal.dist = round(distance, 3),is.outlier = .data$mahal.dist > threshold)
      results
    }
    for (i in unique(names$case)){
      tmpname <- dplyr::filter(names, names$case == i)  
      tmpname <- unique(tmpname$sample)
      tmp <- df[,tmpname]
      ##Adding in idx function back to dataframe
      tmp$idx <- df2$idx
      tmp <-mahalanobis_distance(tmp, -idx) 
      tmp <- tmp %>% dplyr::filter(is.outlier == TRUE)
      write.table(tmp, paste(outdir,as.character(i),'mahalanobisremovedegments.txt',sep=''), append = FALSE, quote = FALSE, sep = "\t",
                  eol = "\n", na = "NA", dec = ".", col.names = TRUE, row.names = FALSE) 
      tmp <- NULL
      tmpname <- NULL
    }
    dfc <- df2[1:4]
    for (i in unique(names$case)){
      tmpname <- dplyr::filter(names, names$case == i)  
      tmpname <- unique(tmpname$sample)
      tmp <- df[,tmpname]
      ##Adding in idx function back to dataframe
      tmp$idx <- df2$idx
      tmp <- mahalanobis_distance(tmp, -idx) 
      tmp$color <- ifelse(tmp$is.outlier == "TRUE", "red", "black")
      if (plot) {
        pdf(paste(outdir,i,"mahal.pdf"))
      }
      if (plot) {
        plot(tmp[,1], tmp[,2],xlim = c(-2, 4),ylim = c(-2, 4),ann=FALSE)
        abline(h=0); abline(v=0)
        points(tmp[,1], tmp[,2], col=tmp[,"color"])
        title(main=paste("Case: ", i, sep=""),xlab=paste(colnames(tmp)[1]," Copy State", sep=""),ylab=paste(colnames(tmp)[2]," Copy State", sep=""))
      }
      if (plot) {
        dev.off()
      }
      tmp <- tmp %>% dplyr::filter(is.outlier == FALSE)
      tmp <- as.data.frame(tmp)
      for (j in 1:length(dfc$idx)){
        dfc[j, i] <- ifelse(any(tmp$idx == j) != 0, tmp[which(tmp$idx == j),1], NA)
      }
      tmp <- NULL
      tmpname <- NULL
    }
    
    write.table(dfc, file = paste(outdir,"earlystableregions_mahal.txt",sep=""), append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", col.names = TRUE, row.names = FALSE)
    return(dfc)
  }
  
  if(change == "linearpc1"){
    pcdist = function(x, outlier=0.9, plot=F, id=NA) {
      xc = as.matrix(x)
      for (j in 1:ncol(x)) {
        xc[,j] = x[,j] - mean(x[,j], na.rm=T)
      }
      udv = svd(xc)
      d = rep(NA, nrow(xc))
      for (i in 1:nrow(xc)) {
        d[i] = sqrt(sum((xc[i,] - udv$v[,1]*sum(xc[i,]*udv$v[,1]))^2))
      }
      prop.var = udv$d[1]^2/sum(udv$d^2)
      if (plot) {
        plot(xc[,1], xc[,2])
        abline(h=0); abline(v=0)
        dydx = udv$v[2,1]/udv$v[1,1]
        abline(coef=c(0,dydx), col="green")
        thresh = quantile(d, outlier)
        ok = d > thresh
        d0 = thresh/(cos(atan(dydx)))
        print(d0)
        print(dydx)
        abline(coef=c(d0,dydx), col="green", lty=2)
        abline(coef=c(-d0,dydx), col="green", lty=2)
        points(xc[ok,1], xc[ok,2], col="red")
        title(paste("Case: ", id, " (n=", ncol(x),"; PC1=", round(prop.var,2), ")", sep=""))
      }
      list(d=d, prop.var=prop.var, thresh=thresh)
    }
    
    # Distance function wrapper
    run_pcdist = function(df, groups, skip=c(1,2,3,4), outlier=0.9, plot=NA) {
      if (is.na(plot)) {
        plot = F
      } else {
        fname = plot
        plot = T
      }
      annot = df[,skip,drop=F]
      data = df[,-skip]
      groups = as.character(groups)
      ugroup = unique(groups)
      y = as.data.frame(matrix(NA, nrow(data), length(ugroup)))
      z = data.frame(matrix(NA,1,length(ugroup)))
      out = data.frame(matrix(NA, nrow(data), length(ugroup)))
      names(z) = ugroup
      names(y) = ugroup
      names(out) = ugroup
      if (plot) {
        pdf(fname)
      }
      for (e in ugroup) {
        ok = apply(data[,which(groups==e)], 1, function(x) all(!is.na(x)))
        res = pcdist(data[ok,which(groups==e)], outlier=outlier, plot=plot, id=e)
        y[ok,e] = res$d
        z[e] = res$prop.var
        out[ok,e] = res$d > res$thresh
      }
      if (plot) {
        dev.off()
      }
      res = cbind(annot,y)
      list(table=res, prop.var=z, outlier=out)
    }
    groups=names$case
    res = run_pcdist(df, groups, plot=paste(outdir,"pcdis_Figure1.pdf",sep=""))
    write.table(res, file = paste(outdir,"earlystableregions_predictions_pcdist.txt",sep=""), append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", col.names = TRUE, row.names = FALSE)
    #######Here I want to resturn a table of early stable CN regions based on these results 
    res_f <- as.data.frame(df2$idx)
    colnames(res_f) <- c("idx")
    print("combining results")
    for (k in unique(groups)){
      tmpname <- dplyr::filter(names, names$case == k)  
      tmpname <- unique(tmpname$sample)
      tmp <- df[,tmpname]
      length_f <- length(tmp)
      ##Adding in idx function back to dataframe
      #tmp$idx <- df2$idx
      ###filter outlier
      o <- res$outlier[,k]
      for (j in 1:length(o)){
        for (l in 1:length_f){
      tmp[j,l] <- ifelse(o[j] == TRUE,NA,tmp[j,l])
      }}
      res_f <- cbind(tmp,res_f)
      print(k)
      }
    res_f <- left_join(df2[,c(1:4)],res_f)
    write.table(res_f, file = paste(outdir,"earlystableregions_pcdist.txt",sep=""), append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", col.names = TRUE, row.names = FALSE)
    return(res_f)
  }
  
  if(return_segmentfile == TRUE){
    #####next lets turn it back into some resemblance of segment files
    ###Each row is a segment per sample
    datalist = list()
    for (i in 5:length(df2)){
      tmp <- df2[,c(1:4,i)]
      j=paste(colnames(tmp[5]))
      tmp <- tmp %>% dplyr::filter(!is.na(tmp[j]))
      #####fill in first chr and start
      ## Find the changes in chr or P1.
      tmp$diff <- c(NA, (tmp[2:nrow(tmp), 'chr'] != tmp[1:(nrow(tmp) - 1), 'chr']) |
                      (tmp[2:nrow(tmp), j] != tmp[1:(nrow(tmp) - 1), j]))
      ## Assign a run number to the rows that increments only with the changes found
      ## above.
      tmp$runNum <- NA
      runNum <- 1
      tmp[1, 'runNum'] <- runNum
      for (ii in 2:nrow(tmp)) {
        if (tmp[ii, 'diff']) runNum <- runNum + 1
        tmp[ii, 'runNum'] <- runNum
      }
      ##get rid of factors
      tmp$chr <- as.numeric(tmp$chr)
      ## Collect data from each run number to a single row.
      tmp2 <- data.frame(chr = tapply(tmp[, 'chr'], tmp[, 'runNum'], min),
                         pos = tapply(tmp[, 'pos'], tmp[, 'runNum'], min),
                         end = tapply(tmp[, 'end'], tmp[, 'runNum'], max),
                         idx = tapply(tmp[, 'idx'], tmp[, 'runNum'], min),
                         case = tapply(tmp[, j], tmp[, 'runNum'], min),
                         runNum = tapply(tmp[, 'runNum'], tmp[, 'runNum'], min),
                         stringsAsFactors = FALSE)
      tmp2 <- tmp2[order(tmp2[, 'runNum']), , drop = FALSE]
      ##get rid of idx and run num
      tmp2$runNum <- NULL
      tmp2$idx <- NULL
      tmp2$Sample <- j
      colnames(tmp2) <- c("Chromosome", "Start", "End","SegCNa", "Sample")
      datalist[[i]] <- tmp2
      tmp <- NULL
      tmp2 <- NULL
    }
    concensus_seg = do.call(rbind, datalist)
    return(concensus_seg)
    write.table(concensus_seg, file = paste(outdir,"earlystableregions_segmentationfile.txt",sep=""), append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", col.names = TRUE, row.names = FALSE)
  }
}
  