library("DESeq2")
library("ggplot2")
library("pheatmap")
library("vsn")
library("RColorBrewer")
library("genefilter")

cat(sprintf(c('Working directory',getwd())))

cat(sprintf('Setting parameters'))

rds <- snakemake@input[['rds_object']]
cat(sprintf(c('RDS: ',rds)))

pca_plot <- snakemake@output[['pca']]
cat(sprintf(c('PCA plot: ',pca_plot)))

labels <- snakemake@params[['pca_labels']]
cat(sprintf(c('PCA Labels: ',labels)))

sd_mean_plot <- snakemake@output[['sd_mean_plot']]
cat(sprintf(c('SD Mean plot: ',sd_mean_plot,'\n')))

distance_plot <- snakemake@output[['distance_plot']]
cat(sprintf(c('Distance plot: ',distance_plot,'\n')))

ggplot_pca_factor <- snakemake@output[['ggplot_pca_factor']]
cat(sprintf(c('ggplot2 PCA factor plot: ', ggplot_pca_factor,'\n')))

cat(sprintf('Load rlog DESeqTransform object'))
rld <-readRDS(rds)

#####
plot_cols <- snakemake@config[['meta_columns_to_plot']]
subset_cols = names(plot_cols)
#####


pdf(pca_plot)
plotPCA(rld, intgroup=labels)
dev.off()

# SD mean plot
pdf(sd_mean_plot)
meanSdPlot(assay(rld))
dev.off()

# Heatmap of distances

pdf(distance_plot)
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)

rownames(sampleDistMatrix) <- paste(rld[[labels]], sep="-")
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix, fontsize=5, scale="row",
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()

# ggplot PCA

pcaData <- plotPCA(rld, intgroup=c(labels), returnData=TRUE)

pdf(ggplot_pca_factor)
percentVar <- round(100 * attr(pcaData, "percentVar"))
g<-ggplot(pcaData, aes_string("PC1", "PC2", color=labels)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed()
g + theme_gray()
dev.off()

