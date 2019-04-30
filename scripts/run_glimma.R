library(Glimma)
library(limma)
library(DESeq2)

condition = snakemake@params[['condition']]
cat(sprintf(c('Condition: ',condition,'\n')))

title = snakemake@params[["contrast"]]
print(title)
contrast = c(condition, snakemake@params[["contrast"]])
rds = snakemake@input[['rds']]
cat(sprintf(c('RDS object: ',rds,'\n')))

md_path = snakemake@output[['ma_plot']]
va_path = snakemake@output[['volcano_plot']]
mdout = tail(strsplit(md_path,'/')[[1]],n=1)
mdout = strsplit(mdout,'.html')[1]
vaout = tail(strsplit(va_path, '/')[[1]],n=1)
vaout = strsplit(vaout,'.html')[1]

out_path = file.path(getwd(),'results','diffexp')
dir.create(out_path)
print(out_path)
rds = readRDS(rds)
groups.df = as.data.frame(colData(rds))


#### by contrasts
#contrasts_to_plot = resultsNames(rds)
res <- results(rds, contrast=contrast)
res$padj[is.na(res$padj)] = 1

rnaseq = as.data.frame(counts(rds, normalized=T))
genes = as.data.frame(row.names(res))
colnames(genes) = 'GeneID'

status_frame = res[,c('log2FoldChange','padj')]
status_frame['status'] = 0
status_frame$padj[is.na(status_frame$padj)] = 1
status_frame[status_frame$padj<0.05 & status_frame$log2FoldChange < 0 ,'status'] = -1
status_frame[status_frame$padj<0.05 & status_frame$log2FoldChange > 0 ,'status'] = 1

title = paste(title[1],'vs',title[2],sep='-')

glMDPlot(res, anno=genes, status=status_frame$status, samples=colnames(rnaseq), 
         counts=log2(rnaseq + 0.0001), groups=groups.df$Condition, main=strsplit(res@elementMetadata$description[2],': ')[[1]][2], 
         transform=F, side.ylab='Log2-expression',launch=FALSE,side.main='GeneID', html = mdout, path=out_path)

## Volcano plot
glXYPlot(x=res$log2FoldChange, y=-log10(res$pvalue), xlab="logFC", ylab="logodds",path=out_path,
         status=status_frame$status, launch=FALSE,counts=log2(rnaseq + 0.0001), groups=groups.df$Condition, anno=genes,main=strsplit(res@elementMetadata$description[2],': ')[[1]][2], html = vaout)

