
######## Snakemake header ########
library(methods)
Snakemake <- setClass(
    "Snakemake",
    slots = c(
        input = "list",
        output = "list",
        params = "list",
        wildcards = "list",
        threads = "numeric",
        log = "list",
        resources = "list",
        config = "list",
        rule = "character"
    )
)
snakemake <- Snakemake(
    input = list('data/cfRNA_9_14_counts.txt', "counts" = c('data/cfRNA_9_14_counts.txt')),
    output = list('results/diffexp/cfRNA_9_14_all.rds', 'results/tables/cfRNA_9_14_normed_counts.txt', 'results/diffexp/cfRNA_9_14_rlog_dds.rds', "rds" = 'results/diffexp/cfRNA_9_14_all.rds', "normed_counts" = 'results/tables/cfRNA_9_14_normed_counts.txt', "rld_out" = 'results/diffexp/cfRNA_9_14_rlog_dds.rds'),
    params = list('/home/exacloud/lustre1/CEDAR/roskamsh/projects/Omics-QC-pipeline/data/PP_cohort_sample_info_log_JE.2018.09.28.txt', c('Type'), 'StudyID', "samples" = '/home/exacloud/lustre1/CEDAR/roskamsh/projects/Omics-QC-pipeline/data/PP_cohort_sample_info_log_JE.2018.09.28.txt', "design" = c('Type'), "row_names" = 'StudyID'),
    wildcards = list(),
    threads = 6,
    log = list('logs/deseq2/init.log'),
    resources = list(),
    config = list("omic_meta_data" = '/home/exacloud/lustre1/CEDAR/roskamsh/projects/Omics-QC-pipeline/data/PP_cohort_sample_info_log_JE.2018.09.28.txt', "base_dir" = '/home/exacloud/lustre1/CEDAR/roskamsh/projects/Omics-QC-pipeline', "gtf_file" = '/home/exacloud/lustre1/CEDAR/anurpa/genomes/gencode.v27.annotation.gtf', "bed_file" = '/home/exacloud/lustre1/CEDAR/roskamsh/projects/Omics-QC-pipeline/data/gencode_v27.bed', "fastq_dir" = '/home/exacloud/lustre1/CEDAR/cfrna/data/LIB180515JW', "star_tool" = '/home/exacloud/lustre1/CEDAR/tools/STAR/tags/2.5.3a/bin/Linux_x86_64/STAR', "picard_tool" = '/home/exacloud/lustre1/CEDAR/tools/picard-tools-1.119/MarkDuplicates.jar', "genome" = '/home/exacloud/lustre1/CEDAR/anurpa/genomes/GRCh38.primary_assembly.genome.fa', "refflat2" = '/home/exacloud/lustre1/CEDAR/Genomes/gencode27_refflat2.txt', "star_index" = '/home/exacloud/lustre1/CEDAR/anurpa/genomes/', "project_id" = 'cfRNA_9_14', "baseline" = 'HD', "linear_model" = 'Type', "seq_layout" = 'PE', "sample_id" = 'StudyID', "dds_design" = 'SampleID', "rseqc_env" = 'source activate rseqc-omic-wf', "meta_columns_to_plot" = list("Type" = 'Type', "Notes" = 'Notes', "RNA_extraction_date" = 'RNA_extraction_date'), "diffexp" = list("contrasts" = list("HD-vs-MM" = c('HD', 'MM'))), "pca" = list("labels" = c('Type')), "samples" = '/home/exacloud/lustre1/CEDAR/roskamsh/projects/Omics-QC-pipeline/data/PP_cohort_sample_info_log_JE.2018.09.28.txt'),
    rule = 'deseq2_init'
)
######## Original script #########

library("DESeq2")

counts = snakemake@input[['counts']]

params = snakemake@params[['samples']]

output = snakemake@output[['rds']]

dds_design = snakemake@params[['design']]

row_names = snakemake@params[['row_names']]

out_table = snakemake@output[['normed_counts']]

rld_out = snakemake@output[['rld_out']]

parallel <- FALSE
if (snakemake@threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}

# colData and countData must have the same sample order, but this is ensured
# by the way we create the count matrix
cts <- read.table(counts, header=TRUE, row.names=1, sep="\t")
coldata <- read.table(params, header=TRUE, row.names=row_names,sep="\t")

coldata <- coldata[colnames(cts),]
cts <- cts[, rownames(coldata)]

dds <- DESeqDataSetFromMatrix(countData=cts,
                              colData=coldata,
                              design= as.formula(paste('~',dds_design)))

# remove uninformative columns
dds <- dds[ rowSums(counts(dds)) > 1, ]
# normalization and preprocessing
dds <- DESeq(dds, parallel=parallel)

saveRDS(dds, file=output)

normed_counts <-counts(dds,normalized=TRUE)
write.table(normed_counts,quote=F,sep='\t',file=out_table)

# obtain normalized counts
rld <- rlog(dds, blind=FALSE)
saveRDS(rld, file=rld_out)
