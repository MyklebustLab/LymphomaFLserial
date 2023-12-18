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

###       using Fishplot to visulaize

###P31, 4 biopsies, nFL
#T1: LN_left_neck
#T2: LN_left_neck
#T3: LN_left_inguinal
#T3: LN_right_axilla



setwd("G:/FL_resequncing/FL_exome_final/fl_latest/11_pyclone/BED_bam/counts_new/pyclone_output_301018")

pyclone_out<- read.table(file="ouput_pyclone/output_200X_i2/P31_200X_i2/tables/loci.tsv",sep="\t",header=T)

mutations<-cast(pyclone_out[,1:4], mutation_id~sample_id, value="cellular_prevalence")
id<-unique(pyclone_out[, c(1,3)] )

P_case<- merge(id, mutations,by="mutation_id")

vafs = data.frame(cluster=P_case$cluster_id,
                  T1_vaf=(P_case$`GCF0150-0031-T01`/2)*100,
                  T2_vaf=(P_case$`GCF0150-0031-T02`/2)*100,
                  T3_vaf=(P_case$`GCF0150-0031-T03`/2)*100,
                  T4_vaf=(P_case$`GCF0150-0031-T04`/2)*100,
                  stringsAsFactors=F)

vafs$cluster<- vafs$cluster+1

###needs manully check density plot to identify which cluster is founding cluster
##clonevol requires founding cluster=1


samples<-c("P31_1","P31_2","P31_3","P31_4")
samples1<-c("P31_1\nPretreat_1","P31_2\nPretreat_2", "P31_3\nRelapsed_1", "P31_4\nRelapsed_2")
samples2<-c("P31_1\nPrimary\nLN_left_neck","P31_2\nPretreat\nLN_left_neck", 
            "P31_3\nRelapsed_1\nLN_left_inguinal", "P31_4\nRelapsed_2\nLN_right_axilla")




names(vafs)[2:5] = samples 
##step 2: run infer.clonal.models, run twice: 1. include all cluster. 
###########2. manual review density plot and exlcude cluster with small number of mutations
##step 2: run infer.clonal.models, run twice: 1. include all cluster. 
###########2. manual review density plot and exlcude cluster with small number of mutations


dir.create("./clonevol/P31")
##


##first: use all clusters, no consensus models

vafs$P31_3[vafs$cluster=="1"]<-vafs$P31_3[vafs$cluster=="1"]+ 0.5
#vafs$P27_2[vafs$cluster=="1"]<-vafs$P27_2[vafs$cluster=="1"]+ 0.5
#vafs$P27_3[vafs$cluster=="3"]<-vafs$P27_3[vafs$cluster=="3"]+1

res = infer.clonal.models(variants=vafs, cluster.col.name="cluster", vaf.col.names=samples,
                          
                          subclonal.test="bootstrap", subclonal.test.model="non-parametric",
                          founding.cluster=1,ignore.clusters = c(6,12:14),
                          cluster.center="mean", num.boots=1000,
                          min.cluster.vaf=0.01, sum.p=0.01, alpha=0.01)



vafs_used<- subset(vafs, !cluster %in% c(6,12:14))

vafs_used$cluster[vafs_used$cluster==11]<-6

res = infer.clonal.models(variants=vafs_used, cluster.col.name="cluster", vaf.col.names=samples,
                          
                          subclonal.test="bootstrap", subclonal.test.model="non-parametric",
                          founding.cluster=1,
                          cluster.center="mean", num.boots=1000,
                          min.cluster.vaf=0.01, sum.p=0.01, alpha=0.01)




res<-convert.consensus.tree.clone.to.branch(res, branch.scale = 'sqrt')

pdf("./clonevol/P31/P31_trees.pdf", useDingbats = FALSE)
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
                   out.dir = './clonevol/P31/',
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
                   out.dir = './clonevol/P31/',
                   out.format = 'pdf',
                   overwrite.output = TRUE,
                   width = 10,
                   height = 4,
                   # vector of width scales for each panel from left to right
                   panel.widths = c(1.5,2.5,1.5,2.5,2))


##generating fish plot
f<- generateFishplotInputs(results = res)
fishes=createFishPlotObjects(f)



pdf('./clonevol/P31/P31_fish_200x_pyclone_anno_loc.pdf', width=14, height=7)
for (i in 1:length(fishes)){
  
  fish = layoutClones(fishes[[i]])
  fish = setCol(fish,f$clonevol.clone.colors)
  fishPlot(fish,shape="bezier", title.btm="P1", cex.title=0.7,cex.vlab = 1.4,
           vlines=seq(1, length(samples2)), vlab=samples2, pad.left=0.5)
}
dev.off()



pdf("./clonevol/P31/P31_box.pdf", width=3, height=3,useDingbats = FALSE, title='')
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
              out.prefix = './clonevol/P31/P31_variants.pairwise.plot')


pdf('./clonevol/P31/P31_flow.pdf')
plot.cluster.flow(vafs, vaf.col.names = samples,
                  sample.names = c('Primary', 'Pretreat', 'Relapsed2', "Relapsed_2"))
dev.off()


####checking coverage 
Pcase<-do.call("rbind", lapply( list.files("input_pyclone_271018_newpara/200X/GCF0150-0031-N01_200X/",full=TRUE),
                                read.table, header=TRUE, sep="\t"))




########
#min (dp) =min(c1inP4_1$var_counts+c1inP4_1$ref_counts)=273

#max (dp) =max(c1inP4_1$var_counts+c1inP4_1$ref_counts)=648
#median (dp) =median(c1inP4_1$var_counts+c1inP4_1$ref_counts)=380
#mean (dp) =mean(c1inP4_1$var_counts+c1inP4_1$ref_counts)=417

library(ggplot2)



pdf("clonevol/P31/coverage.pdf")
ggplot(Pcase, aes(x=(var_counts+ref_counts)))+
  geom_histogram(position="dodge")+
  facet_grid(~sample)

dev.off()

save.image("G:/FL_resequncing/FL_exome_final/fl_latest/11_pyclone/BED_bam/counts_new/pyclone_output_301018/P31.RData")


