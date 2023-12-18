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

###P59, 2 biopsies, nFL
#T1: LN_right_inguinal
#T2: LN_left_inguinal






setwd("G:/FL_resequncing/FL_exome_final/fl_latest/11_pyclone/BED_bam/counts_new/pyclone_output_301018")

pyclone_out<- read.table(file="ouput_pyclone/output_200X_i2/P59_200X_i2/tables/loci.tsv",sep="\t",header=T)

mutations<-cast(pyclone_out[,1:4], mutation_id~sample_id, value="cellular_prevalence")
id<-unique(pyclone_out[, c(1,3)] )

P_case<- merge(id, mutations,by="mutation_id")

vafs = data.frame(cluster=P_case$cluster_id,
                  T1_vaf=(P_case$`GCF0150-0059-T01`/2)*100,
                  T2_vaf=(P_case$`GCF0150-0059-T02`/2)*100,
                  stringsAsFactors=F)

vafs$cluster<- vafs$cluster+1

###needs manully check density plot to identify which cluster is founding cluster
##clonevol requires founding cluster=1

samples<-c("P59_1","P59_2")
samples1<-c("P59_1\nPretreat\nLN_right_inguinal","P59_2\nRelapsed\nLN_left_inguinal")



names(vafs)[2:3] = samples 
##step 2: run infer.clonal.models, run twice: 1. include all cluster. 
###########2. manual review density plot and exlcude cluster with small number of mutations

max_id<- max(vafs$cluster)+1

vafs$cluster[vafs$cluster==1]<-max_id
vafs$cluster[vafs$cluster==5]<-1
vafs$cluster[vafs$cluster==max_id]<-5



dir.create("./clonevol/P59_new")
##


##first: use all clusters, no consensus models

#vafs$P14_2[vafs$cluster=="4"]<-vafs$P14_2[vafs$cluster=="4"]-1
res = infer.clonal.models(variants=vafs, cluster.col.name="cluster", vaf.col.names=samples,
                          
                          subclonal.test="bootstrap", subclonal.test.model="non-parametric",
                          founding.cluster=1,
                          cluster.center="mean", num.boots=1000,
                          min.cluster.vaf=0.01, sum.p=0.01, alpha=0.01)


res = infer.clonal.models(variants=vafs, cluster.col.name="cluster", vaf.col.names=samples,
                          
                          subclonal.test="bootstrap", subclonal.test.model="non-parametric",
                          founding.cluster=1,ignore.clusters = c(8:9,11:12),
                          ,cluster.center="mean", num.boots=1000,
                          min.cluster.vaf=0.01, sum.p=0.01, alpha=0.01)

vafs$P59_1[vafs$cluster=="5"]<-vafs$P59_1[vafs$cluster=="5"]+1




vafs_used<- subset(vafs, !cluster %in% c(8:9,11:12))

vafs_used$cluster[vafs_used$cluster==10]<-8

res = infer.clonal.models(variants=vafs_used, cluster.col.name="cluster", vaf.col.names=samples,
                          subclonal.test="bootstrap", subclonal.test.model="non-parametric",
                          founding.cluster=1,
                          cluster.center="mean", num.boots=1000,
                          min.cluster.vaf=0.01, sum.p=0.01, alpha=0.01)


res<-convert.consensus.tree.clone.to.branch(res, branch.scale = 'sqrt')

pdf("./clonevol/P59_new/P59_trees.pdf", useDingbats = FALSE)
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
                   disable.cell.frac = FALSE,
                   # output figure parameters
                   out.dir = './clonevol/P59_new/',
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
                   out.dir = './clonevol/P59_new/',
                   out.format = 'pdf',
                   overwrite.output = TRUE,
                   width = 10,
                   height = 4,
                   # vector of width scales for each panel from left to right
                   panel.widths = c(1.5,2.5,1.5,2.5,2))


##generating fish plot
f<- generateFishplotInputs(results = res)
fishes=createFishPlotObjects(f)



pdf('./clonevol/P59_new/P59_fish_200x_pyclone_anno_loc.pdf', width=14, height=7)
for (i in 1:length(fishes)){
  
  fish = layoutClones(fishes[[i]])
  fish = setCol(fish,f$clonevol.clone.colors)
  fishPlot(fish,shape="bezier", title.btm="P1", cex.title=0.7,cex.vlab = 1.4,
           vlines=seq(1, length(samples1)), vlab=samples1, pad.left=0.5)
}
dev.off()



pdf("./clonevol/P59_new/P59_box.pdf", width=3, height=3,useDingbats = FALSE, title='')
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


plot.pairwise(vafs_used, col.names = samples,
              out.prefix = './clonevol/P59_new/P59_variants.pairwise.plot')


pdf('./clonevol/P59_new/P59_flow.pdf')
plot.cluster.flow(vafs_used, vaf.col.names = samples,
                  sample.names = c('Pretreat', 'Relapsed'))
dev.off()


####checking coverage f
Pcase<-do.call("rbind", lapply( list.files("input_pyclone_271018_newpara/200X/GCF0150-0059-N01_200X/",full=TRUE),
                                read.table, header=TRUE, sep="\t"))




########
#min (dp) =min(c1inP4_1$var_counts+c1inP4_1$ref_counts)=273

#max (dp) =max(c1inP4_1$var_counts+c1inP4_1$ref_counts)=648
#median (dp) =median(c1inP4_1$var_counts+c1inP4_1$ref_counts)=380
#mean (dp) =mean(c1inP4_1$var_counts+c1inP4_1$ref_counts)=417

library(ggplot2)



pdf("clonevol/P59_new/coverage.pdf")
ggplot(Pcase, aes(x=(var_counts+ref_counts)))+
  geom_histogram(position="dodge")+
  facet_grid(~sample)

dev.off()

save.image("G:/FL_resequncing/FL_exome_final/fl_latest/11_pyclone/BED_bam/counts_new/pyclone_output_301018/P59_new.RData")


