library(DESeq2)
library(ggplot2)
library(reshape2)
library(data.table)
library(plyr)
library(RColorBrewer)


rld <- snakemake@input[['rld']]
cat(sprintf(c('rld: ', rld, '\n')))


condition <- snakemake@params[['linear_model']]
cat(sprintf(c('condition: ', condition, '\n')))

project_id <- snakemake@params[['project_id']]

density_plot <- snakemake@output[['density']]
cat(sprintf(c('Density plot : ', density_plot, '\n')))


colors <- snakemake@params['colors']
discrete <- snakemake@params['discrete']


rld = readRDS(rld)
normed_values = assay(rld)
normed_t = t(normed_values)
meta = colData(rld)


if(colors!=NA & discrete ==NA){
    if(brewer.pal.info[colors]$maxcolors >= length(levels(meta[condition]))){
        pal <- colorRampPalette(brewer.pal(length(levels(meta[condition])),name=colors))
    }else if(discrete !=NA){
        pal <- as.vector(discrete)
    }
}

joined_counts = cbind(meta[condition],normed_t)

x = as.data.table(joined_counts)
mm <- melt(x,id=condition)

mu <- ddply(mm, condition, summarise, grp.mean=mean(value))
pdf(density_plot)
p<-ggplot(mm, aes_string(x='value', color=condition)) +
  geom_density()+
  geom_vline(data=mu, aes_string(xintercept='grp.mean', color=condition),
             linetype="dashed") + xlab('regularized log expression') + ggtitle(eval(project_id)) + theme(plot.title = element_text(hjust = 0.5))
p
dev.off()