library("XenofilteR")

graftBam <- snakemake@input[["graft"]]
hostBam <- snakemake@input[["host"]]
sampname <- gsub("_bam", "", unlist(strsplit(graftBam,"/"))[3])
outDir <- paste("samples", "xeno", sampname, sep="/")

sample.list <- as.data.frame(cbind(graftBam, hostBam))
bp.param <- SnowParam(workers = 12, type = "SOCK")

dir.create(outDir, recursive = TRUE)
XenofilteR(sample.list, destination.folder = outDir, bp.param = bp.param)

