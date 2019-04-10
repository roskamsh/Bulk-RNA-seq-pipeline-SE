args <- commandArgs()

annoFile = snakemake@params[['anno']]

biotypes <- snakemake@params[['biotypes']]

countsFile <- snakemake@input[['countsFile']]

mito <- snakemake@params[['mito']]


help <- function(){
    cat("RNAseq_filterCounts.R :
- Script to filter counts table based on biotype and/or mito genes.
- Assumes that the genes are in their own column called Genes and are common names not Ensembl ids.
- Anno file must have columns: external_gene_name, gene_biotype.
- You can find biotype discriptions here:
     https://uswest.ensembl.org/info/genome/genebuild/biotypes.html
     https://www.gencodegenes.org/pages/biotypes.html
\n")
    cat("Usage: \n")
    cat("--countsFile : counts table to be filtered (.txt tab sep )     [ required ]\n")
    cat("--annoFile   : FDR adjusted p-value cutoff      [ required ]\n")
    cat("--mito       : remove mito genes (1/0)          [ default = 1 ]\n")
    cat("--biotypes   : comma seperated biotypes to keep [ required ]
           example: --biotypes=protein_coding,lincRNA,sRNA,rRNA,snRNA,snoRNA\n")
    
    cat("\n")
    q()
}

## Save values of each argument
if(!is.na(charmatch("--help",args)) || !is.na(charmatch("-h",args))){
    help()
}

## set defaults if options are not provided
if (identical( mito, character(0) )){
    mito <- 1
}else{
    mito <- as.numeric(mito)
}

if (identical( biotypes, character(0) )){
    biotypes <- ""
}

##----------load counts------------#
print("Loading counts table")
print(countsFile)

## check if an rda file or tab sep
if(grepl('rda|RData|Rdata',countsFile)){
    counts <- get(load(file=countsFile))
}
if(grepl('txt|tsv',countsFile)){
    counts <- read.delim(file=countsFile)
}

##----------load counts------------#
print("Loading counts table")
print(countsFile)

## must be a tsv or txt tab sep file

counts <- read.delim(file=countsFile)

##----------load anno------------#
print("Loading annotation table")
print(annoFile)

## check if an rda file or tab sep
if(grepl('rda|RData|Rdata',annoFile)){
    anno <- get(load(file=annoFile))
}
if(grepl('txt|tsv',annoFile)){
    anno <- read.delim(file=annoFile)
}

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
