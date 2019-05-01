annoFile = snakemake@params[['anno']]

biotypes <- snakemake@params[['biotypes']]

countsFile <- snakemake@input[['countsFile']]

mito <- snakemake@params[['mito']]

counts <- read.delim(file=countsFile)

##----------load counts------------#
print("Loading counts table")
print(countsFile)

## must be a tsv or txt tab sep file

counts <- read.delim(file=countsFile)

##----------load anno------------#
print("Loading annotation table")
print(annoFile)

## load
anno <- get(load(file=annoFile))

if(strsplit(biotypes, split='\\,')[[1]]!=""){
    anno.sub <- anno[paste(anno$gene_biotype) %in% strsplit(biotypes, split='\\,')[[1]] ,]
    counts.sub <- counts[paste(counts$Genes) %in% unique(paste(anno.sub$external_gene_name)) , ]
}else{
    print("no biotypes provided")
    counts.sub <- counts
}

if(mito==1){
    print("tossing MT- genes")
    counts.sub <- counts.sub[grep("^MT-", paste(counts.sub$Genes), invert=TRUE), ]
}

write.table(counts.sub, file=sub(".txt", ".filt.txt", countsFile), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
