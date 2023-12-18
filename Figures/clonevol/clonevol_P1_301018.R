library(dplyr)
library(plyr)
library(clonevol)
library(fishplot)
library(reshape)

####coverage: from reseq BAM files 
####step1: Using mpileup to extract coverage for each site. (Prepared 27/10/18), default parameter

####step2: Prepare input for pyclone, using segments from ASCAT

####step3: Run Pyclone twice. The first run is to identify founding cluster, 
####       the second run adding --tumour_content based on the cellular prevalence for each biopsy (29/10/18)

###step4: using clonevol to build evolution tree and model (in this script) (30/10/18)

###       using Fishplot to visulaize

##P1: 2 biopsies, clinic transformed

#T1:LN_right_neck
#T2:LN_left_inguinal

setwd("F:/FL_exome_final/fl_latest/11_pyclone/BED_bam/counts_new/pyclone_output_301018/")

pyclone_out<- read.table(file="ouput_pyclone/output_200X_i2/GCF0150-0001-N01_200X_i2/tables/loci.tsv",sep="\t",header=T)

mutations<-cast(pyclone_out[,1:4], mutation_id~sample_id, value="cellular_prevalence")
id<-unique(pyclone_out[, c(1,3)] )

P_case<- merge(id, mutations,by="mutation_id")

vafs = data.frame(cluster=P_case$cluster_id,
                  T1_vaf=(P_case$`GCF0150-0001-T01`/2)*100,
                  T2_vaf=(P_case$`GCF0150-0001-T02`/2)*100,
                  stringsAsFactors=F)

vafs$cluster<- vafs$cluster+1

###needs manully check density plot to identify which cluster is founding cluster
##clonevol requires founding cluster=1
##in this case founding cluster is #3

max_id<- max(vafs$cluster)+1

vafs$cluster[vafs$cluster==1]<-max_id
vafs$cluster[vafs$cluster==3]<-1
vafs$cluster[vafs$cluster==max_id]<-3



samples<-c("P1_1","P1_2")
samples1<- c("P1_1\nPrimary\nLN_right_neck","P1_2\nRelapsed\nLN_left_inguinal")

names(vafs)[2:3] = samples 
##step 2: run infer.clonal.models, run twice: 1. include all cluster. 
###########2. manual review density plot and exlcude cluster with small number of mutations


dir.create("./clonevol/P1")
##


##first: use all clusters, no consensus models

res = infer.clonal.models(variants=vafs, cluster.col.name="cluster", vaf.col.names=samples,
                          
                          subclonal.test="bootstrap", subclonal.test.model="non-parametric",
                          founding.cluster=1,
                          cluster.center="mean", num.boots=1000,
                          min.cluster.vaf=0.01, sum.p=0.01, alpha=0.01)


#Finding consensus models across samples...
#Found  0 consensus model(s)
#Found 0 consensus model(s)


##second: ignore cluster 5:7 consensus model based on manual review

##to be able to run clonevol, have to change the cellular_prevalance
vafs$P1_2[vafs$cluster=="2"]<- vafs$P1_2[vafs$cluster=="2"]-1.2

vafs_used<- subset(vafs, !cluster %in% c(5:7))

vafs_used$cluster[vafs_used$cluster=="8"]<- 5
vafs_used$cluster[vafs_used$cluster=="9"]<- 6

res = infer.clonal.models(variants=vafs_used, cluster.col.name="cluster", vaf.col.names=samples,
                          subclonal.test="bootstrap", subclonal.test.model="non-parametric",
                          founding.cluster=1,
                          cluster.center="mean", num.boots=1000,
                          min.cluster.vaf=0.01, sum.p=0.01, alpha=0.01)



res<-convert.consensus.tree.clone.to.branch(res, branch.scale = 'sqrt')

pdf("./clonevol/P1/P1_trees.pdf", useDingbats = FALSE)
plot.all.trees.clone.as.branch(res, branch.width = 0.5, node.size = 1, node.label.size = 0.5)
dev.off()



plot.clonal.models(res,
                   # box plot parameters
                   box.plot = TRUE,
                   fancy.boxplot = TRUE,
                   fancy.variant.boxplot.highlight = 'is.driver',
                   fancy.variant.boxplot.highlight.shape = 21,
                   fancy.variant.boxplot.highlight.fill.color = 'red',
                   fancy.variant.boxplot.highlight.color = 'black',
                   fancy.variant.boxplot.highlight.note.col.name = 'gene',
                   fancy.variant.boxplot.highlight.note.color = 'blue',
                   fancy.variant.boxplot.highlight.note.size = 2,
                   fancy.variant.boxplot.jitter.alpha = 1,
                   fancy.variant.boxplot.jitter.center.color = 'grey50',
                   fancy.variant.boxplot.base_size = 12,
                   fancy.variant.boxplot.plot.margin = 1,
                   fancy.variant.boxplot.vaf.suffix = '.VAF',
                   # bell plot parameters
                   clone.shape = 'bell',
                   bell.event = TRUE,
                   bell.event.label.color = 'blue',
                   bell.event.label.angle = 60,
                   clone.time.step.scale = 1,
                   bell.curve.step = 2,
                   # node-based consensus tree parameters
                   merged.tree.plot = TRUE,
                   tree.node.label.split.character = NULL,
                   tree.node.shape = 'circle',
                   tree.node.size = 30,
                   tree.node.text.size = 0.5,
                   merged.tree.node.size.scale = 1.25,
                   merged.tree.node.text.size.scale = 2,
                   merged.tree.cell.frac.ci = FALSE,
                   # branch-based consensus tree parameters
                   merged.tree.clone.as.branch = TRUE,
                   mtcab.event.sep.char = ',',
                   mtcab.branch.text.size = 1,
                   mtcab.branch.width = 0.75,
                   mtcab.node.size = 3,
                   mtcab.node.label.size = 1,
                   mtcab.node.text.size = 1.5,
                   # cellular population parameters
                   cell.plot = TRUE,
                   num.cells = 100,
                   cell.border.size = 0.25,
                   cell.border.color = 'black',
                   clone.grouping = 'horizontal',
                   #meta-parameters
                   scale.monoclonal.cell.frac = TRUE,
                   show.score = FALSE,
                   cell.frac.ci = TRUE,
                   disable.cell.frac = TRUE,
                   # output figure parameters
                   out.dir = './clonevol/P1/',
                   out.format = 'pdf',
                   overwrite.output = TRUE,
                   width = 10,
                   height = 4,
                   # vector of width scales for each panel from left to right
                   panel.widths = c(1.5,2.5,1.5,2.5,2))



plot.clonal.models(res,
                   # box plot parameters
                   box.plot = TRUE,
                   fancy.boxplot = TRUE,
                   fancy.variant.boxplot.highlight = 'is.driver',
                   fancy.variant.boxplot.highlight.shape = 21,
                   fancy.variant.boxplot.highlight.fill.color = 'red',
                   fancy.variant.boxplot.highlight.color = 'black',
                   fancy.variant.boxplot.highlight.note.col.name = 'gene',
                   fancy.variant.boxplot.highlight.note.color = 'blue',
                   fancy.variant.boxplot.highlight.note.size = 2,
                   fancy.variant.boxplot.jitter.alpha = 1,
                   fancy.variant.boxplot.jitter.center.color = 'grey50',
                   fancy.variant.boxplot.base_size = 12,
                   fancy.variant.boxplot.plot.margin = 1,
                   fancy.variant.boxplot.vaf.suffix = '.VAF',
                   # bell plot parameters
                   clone.shape = 'bell',
                   bell.event = TRUE,
                   bell.event.label.color = 'blue',
                   bell.event.label.angle = 60,
                   clone.time.step.scale = 1,
                   bell.curve.step = 2,
                   # node-based consensus tree parameters
                   merged.tree.plot = TRUE,
                   tree.node.label.split.character = NULL,
                   tree.node.shape = 'circle',
                   tree.node.size = 30,
                   tree.node.text.size = 0.5,
                   merged.tree.node.size.scale = 1.25,
                   merged.tree.node.text.size.scale = 2,
                   merged.tree.cell.frac.ci = FALSE,
                   # branch-based consensus tree parameters
                   merged.tree.clone.as.branch = TRUE,
                   mtcab.event.sep.char = ',',
                   mtcab.branch.text.size = 1,
                   mtcab.branch.width = 0.75,
                   mtcab.node.size = 3,
                   mtcab.node.label.size = 1,
                   mtcab.node.text.size = 1.5,
                   # cellular population parameters
                   cell.plot = TRUE,
                   num.cells = 100,
                   cell.border.size = 0.25,
                   cell.border.color = 'black',
                   clone.grouping = 'horizontal',
                   #meta-parameters
                   scale.monoclonal.cell.frac = TRUE,
                   show.score = FALSE,
                   cell.frac.ci = TRUE,
                   disable.cell.frac = FALSE,
                   # output figure parameters
                   out.dir = './clonevol/P1/',
                   out.format = 'pdf',
                   overwrite.output = TRUE,
                   width = 10,
                   height = 4,
                   # vector of width scales for each panel from left to right
                   panel.widths = c(1.5,2.5,1.5,2.5,2))

###removing cell.frac annotation



plot.clonal.models(res,
                   # box plot parameters
                   box.plot = TRUE,
                   fancy.boxplot = TRUE,
                   fancy.variant.boxplot.highlight = 'is.driver',
                   fancy.variant.boxplot.highlight.shape = 21,
                   fancy.variant.boxplot.highlight.fill.color = 'red',
                   fancy.variant.boxplot.highlight.color = 'black',
                   fancy.variant.boxplot.highlight.note.col.name = 'gene',
                   fancy.variant.boxplot.highlight.note.color = 'blue',
                   fancy.variant.boxplot.highlight.note.size = 2,
                   fancy.variant.boxplot.jitter.alpha = 1,
                   fancy.variant.boxplot.jitter.center.color = 'grey50',
                   fancy.variant.boxplot.base_size = 12,
                   fancy.variant.boxplot.plot.margin = 1,
                   fancy.variant.boxplot.vaf.suffix = '.VAF',
                   # bell plot parameters
                   clone.shape = 'bell',
                   bell.event = TRUE,
                   bell.event.label.color = 'blue',
                   bell.event.label.angle = 60,
                   clone.time.step.scale = 1,
                   bell.curve.step = 2,
                   # node-based consensus tree parameters
                   merged.tree.plot = TRUE,
                   tree.node.label.split.character = NULL,
                   tree.node.shape = 'circle',
                   tree.node.size = 30,
                   tree.node.text.size = 0.5,
                   merged.tree.node.size.scale = 1.25,
                   merged.tree.node.text.size.scale = 2,
                   merged.tree.cell.frac.ci = FALSE,
                   # branch-based consensus tree parameters
                   merged.tree.clone.as.branch = TRUE,
                   mtcab.event.sep.char = ',',
                   mtcab.branch.text.size = 1,
                   mtcab.branch.width = 0.75,
                   mtcab.node.size = 3,
                   mtcab.node.label.size = 1,
                   mtcab.node.text.size = 1.5,
                   # cellular population parameters
                   cell.plot = TRUE,
                   num.cells = 100,
                   cell.border.size = 0.25,
                   cell.border.color = 'black',
                   clone.grouping = 'horizontal',
                   #meta-parameters
                   scale.monoclonal.cell.frac = TRUE,
                   show.score = FALSE,
                   cell.frac.ci = TRUE,
                   disable.cell.frac = TRUE,
                   # output figure parameters
                   out.dir = './clonevol/P1/',
                   out.format = 'pdf',
                   overwrite.output = TRUE,
                   width = 10,
                   height = 4,
                   # vector of width scales for each panel from left to right
                   panel.widths = c(1.5,2.5,1.5,2.5,2))


##generating fish plot
f<- generateFishplotInputs(results = res)
fishes=createFishPlotObjects(f)



samples1<- c("P1_1\nPrimary", "P1_2\nRelapsed")




samples2<- c("P1_1\nPrimary\nLN_right_neck","P1_2\nRelapsed\nLN_left_inguinal")


pdf('./clonevol/P1/P1_fish_200x_pyclone_anno_new.pdf', width=14, height=7)
for (i in 1:length(fishes)){
  
  fish = layoutClones(fishes[[i]])
  fish = setCol(fish,f$clonevol.clone.colors)
  fishPlot(fish,shape="spline", title.btm="P1", cex.title=0.7,cex.vlab = 1.4,
           vlines=seq(1, length(samples2)), vlab=samples2, pad.left=0.5)
}

dev<-dev.off()

##other plots

pdf("./clonevol/P1/P1_box.pdf", width=3, height=3,useDingbats = FALSE, title='')
pp<-plot.variant.clusters(vafs_used,
                          cluster.col.name = 'cluster',
                          show.cluster.size = FALSE,
                          cluster.size.text.color = 'blue',
                          vaf.col.names = samples,
                          vaf.limits = 70,
                          sample.title.size = 20,
                          violin = FALSE,
                          box = FALSE,
                          jitter = TRUE,
                          jitter.shape = 1,
                          
                          jitter.size = 3,
                          jitter.alpha = 1,
                          jitter.center.method = 'median',
                          jitter.center.size = 1,
                          jitter.center.color = 'darkgray',
                          jitter.center.display.value = 'none',
                          highlight = 'is.driver',
                          highlight.shape = 21,
                          highlight.color = 'blue',
                          highlight.fill.color = 'green',
                          highlight.note.col.name = 'gene',
                          highlight.note.size = 2,
                          order.by.total.vaf = FALSE)

dev.off()



pdf('./clonevol/P1/P1_flow.pdf')
plot.cluster.flow(vafs_used, vaf.col.names = samples,
                  sample.names = c('Primary', 'Relapse'))
dev.off()

####checking coverage for cluster 2 in P1_1

c1<- pyclone_out%>%filter(cluster_id=="1")

P1_1<- read.table(file="input_pyclone_271018_newpara/200X/GCF0150-0001-N01_200X/GCF0150-0001-T01.txt",sep="\t",header=T)
P1_2<- read.table(file="input_pyclone_271018_newpara/200X/GCF0150-0001-N01_200X/GCF0150-0001-T02.txt",sep="\t",header=T)
P1<- rbind(P1_1, P1_2)

c1inP1_1<- subset(P1_1, mutation_id %in% unique(c1$mutation_id))

########
#min (dp) =min(c1inP1_1$var_counts+c1inP1_1$ref_counts)=206

#max (dp) =max(c1inP1_1$var_counts+c1inP1_1$ref_counts)=1500
#median (dp) =median(c1inP1_1$var_counts+c1inP1_1$ref_counts)=463
#mean (dp) =mean(c1inP1_1$var_counts+c1inP1_1$ref_counts)=544.9

library(ggplot2)

ggplot(P1, aes(x=var_counts+ref_counts), group_by=sample)+geom_histogram()

pdf("clonevol/P1/coverage.pdf")
ggplot(P1, aes(x=var_counts+ref_counts, fill=sample))+geom_histogram()

dev.off()

