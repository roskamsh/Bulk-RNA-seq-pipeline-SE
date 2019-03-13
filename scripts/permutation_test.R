library(DESeq2)
library(dplyr)

# Generate subdata 
counts <- snakemake@input[['counts']]

metadata <- snakemake@params[['samples']]

sampleID <- snakemake@params[['sample_id']]

hist <- snakemake@output[['histogram']]

numGenes <- snakemake@output[['numGenes']]

permList <- snakemake@output[['permList']]

Type <- snakemake@params[['linear_model']]

contrast <- snakemake@params[['contrast']]

baseline <- contrast[[1]]

target <- contrast[[2]]

md <- read.delim(file=metadata, sep = "\t", stringsAsFactors = FALSE)

# Extract only the Types that we want in further analysis & only the PP_ID and Status informative columns
md <- select(md, sampleID, Type)
md <- filter(md, !!as.name(Type) == baseline | !!as.name(Type) == target)

# Read in counts table
subdata <- read.table(counts, header=TRUE, row.names=1, sep="\t")

# Keep only the PP_IDs of the types we have chosen in the metadata table above
rownames(md) <- md[[sampleID]]
md[[sampleID]] <- NULL
keep <- colnames(subdata)[colnames(subdata) %in% rownames(md)]
subdata <- subdata[, keep]
dim(subdata)

# Get the number of Cancer samples and number of HD samples from md table
num1 = sum(md[[Type]] == baseline)
num2 = sum(md[[Type]] == target)

# Create a vector for both HD and Can, with a 1 for every HD and a 2 for every Cancer sample
One_vector = rep(c(1), times = num1)
Two_vector = rep(c(2), times = num2)

# Permutation
# Concatenate the HD and Can vector to create your "start group" vector
start_group = c(One_vector, Two_vector)
cutoff=0.01
number_of_diff_genes=c()
group_list = list()
number_of_try = 500

for (i in 1:number_of_try)
{
  print(i)
  group = data.frame(type=factor(sample(start_group)))
  
  dds = DESeqDataSetFromMatrix(countData = subdata,
                               colData = group,
                               design = ~ type)
  
  # Extract normalized counts
  dds = estimateSizeFactors(dds)
  
  # Remove genes with zero counts over all samples
  dds <- dds[ rowSums(counts(dds)) >= 1, ]
  
  # Make sure of reference, set it by rlevel
  dds$type = relevel(dds$type, ref = 1)
  
  # The standard differential expression analysis steps are wrapped into a single function, DESeq
  dds = DESeq(dds)
  
  # Extract results
  res = results(dds, contrast = c("type", "1", "2"), independentFiltering = FALSE,cooksCutoff = Inf)
  
  tmp=sum(res$padj < cutoff, na.rm=TRUE)
  number_of_diff_genes = c(number_of_diff_genes,tmp)
  group_list[[i]] <- group
  
}

# Obtain the number of genes that meet padj<0.01 for reference line in histogram
dds <- DESeqDataSetFromMatrix(countData=subdata,
                              colData=md,
                              design= as.formula(paste('~',Type)))

# Remove uninformative columns
dds <- dds[ rowSums(counts(dds)) >= 1, ]

# Normalization and pre-processing
dds <- DESeq(dds)

# Extract results and the number of significant genes with padj<0.01
results = results(dds, contrast = c(Type, target, baseline), independentFiltering = FALSE,cooksCutoff = Inf)
numSig <- sum(results$padj < cutoff, na.rm=TRUE)

pdf(hist)
hist(number_of_diff_genes, col = "red", breaks = 100, main=paste(number_of_try, "Random Permutations:", baseline, "vs.", target, sep=" "), xlab="Number of genes, padj<0.01",ylim=c(0,500),xlim=c(0,max(number_of_diff_genes, numSig)))
abline(v = numSig, col = "black", lwd = 2)
dev.off()

df <- data.frame(stringsAsFactors = FALSE)

for (i in 1:number_of_try) {
  if (i==1) {
    df = group_list[[i]]
  }
  else {
    df = cbind(df, group_list[[i]])
  }
  colnames(df)[i] = paste("perm",i, sep = "_")
}

write.csv(number_of_diff_genes, numGenes)
write.csv(df, permList)
