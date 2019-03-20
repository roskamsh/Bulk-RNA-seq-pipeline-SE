library("DESeq2")
library("pheatmap")

print('Setting parameters')

rds = snakemake@input[['rds']]
cat(sprintf(c('RDS object: ',rds,'\n')))

rld = snakemake@input[['rds']]
cat(sprintf(c('RDS object: ',rds,'\n')))

condition = snakemake@params[['condition']]
cat(sprintf(c('Condition: ',condition,'\n')))

ma_plot = snakemake@output[['ma_plot']]
cat(sprintf(c('MA plot', ma_plot,'\n')))

out_table = snakemake@output[['table']]
cat(sprintf(c('Summary results table', out_table,'\n')))

panel_ma = snakemake@output[['panel_ma']]
cat(sprintf(c('MA panel', panel_ma,'\n')))

heatmap_plot = snakemake@output[['heatmap_plot']]
cat(sprintf(c('Heatmap plot', heatmap_plot,'\n')))

var_heat = snakemake@output[['var_heat']]
cat(sprintf(c('Variance Heatmap plot', var_heat,'\n')))

cat(sprintf('Load dds DESeqTransform object'))
dds <- readRDS(rds)

cat(sprintf('Load rlog DESeqTransform object'))
rld <- readRDS(rld)

plot_cols <- snakemake@config[['meta_columns_to_plot']]
subset_cols = names(plot_cols)

contrast <- c(condition, snakemake@params[["contrast"]])

parallel <- FALSE
if (snakemake@threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}

res <- results(dds, contrast=contrast, parallel=parallel, independentFiltering = FALSE, cooksCutoff = Inf)
# shrink fold changes for lowly expressed genes
res <- lfcShrink(dds, contrast=contrast, res=res)

# MA plot
pdf(ma_plot)
plotMA(res, ylim=c(-2,2))
dev.off()

# P-histogram
p_hist = snakemake@output[['p_hist']]
pdf(p_hist)
hist(res$pvalue[res$baseMean > 1], breaks = 0:20/20, col = "grey50", border = "white", main='P values for genes with mean normalized count larger than 1',xlab='pvalue')
dev.off()

#panel ma plot
pdf(panel_ma)
par(mfrow=c(2,2),mar=c(2,2,1,1) +0.1)
ylim <- c(-2.5,2.5)
resGA <- results(dds, contrast=contrast, lfcThreshold=.5, altHypothesis="greaterAbs")
resLA <- results(dds, contrast=contrast, lfcThreshold=.5, altHypothesis="lessAbs")
resG <- results(dds, contrast=contrast, lfcThreshold=.5, altHypothesis="greater")
resL <- results(dds, contrast=contrast, lfcThreshold=.5, altHypothesis="less")
drawLines <- function() abline(h=c(-.5,.5),col="dodgerblue",lwd=2)
plotMA(resGA, ylim=ylim); drawLines()
plotMA(resLA, ylim=ylim); drawLines()
plotMA(resG, ylim=ylim); drawLines()
plotMA(resL, ylim=ylim); drawLines()
mtext(resG@elementMetadata$description[[2]], outer=T, cex=.6,line=-1)
dev.off()

# Heatmap of top 20 genes
topGenes <- head(order(res$padj),20)

df <- as.data.frame(colData(rld))
if (length(subset_cols)==1) {
  annot <- as.data.frame(cbind(rownames(df), paste(df[[subset_cols[1]]])))
  names(annot) <- c("SampleID", subset_cols[1])
  rownames(annot) <- annot$SampleID
  annot$SampleID <- NULL
} else {
  annot <- df[,subset_cols]
}

pdf(heatmap_plot)
pheatmap(assay(rld)[topGenes,], cluster_rows=T, scale="row", fontsize=6,fontsize_row=6,fontsize_col=6,show_rownames=T, cluster_cols=T, annotation_col=annot, labels_col=as.character(rownames(df)), main = paste("Heatmap of top 20 DE genes:", contrast[2], "vs", contrast[3]))
dev.off()

# Variance Heatmap
pdf(var_heat)
topVarGenes <- head(order(rowVars(assay(rld)), decreasing = TRUE), 20)
mat  <- assay(rld)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
pheatmap(mat, scale="row", annotation_col = annot,fontsize=6, main = paste("Heatmap of top 20 most variable genes:", contrast[2], "vs", contrast[3]))
dev.off()

# sort by p-value
res <- res[order(res$padj),]
write.table(as.data.frame(res), file=out_table, quote=FALSE, sep='\t')
