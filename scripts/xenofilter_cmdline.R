library("XenofilteR")
args <- commandArgs()

help <- function(){
  cat("filterXenograft.R :
      - Filter out host genome reads from graft reads\n")
  cat("Usage: \n")
  cat("--graftBam  : bam file of graft organism               [ required ]\n")
  cat("--hostBam   : bam file of host organism                [ required ]\n")
  cat("\n")
  q()
}

## Save values of each argument
if(!is.na(charmatch("--help",args)) || !is.na(charmatch("-h",args)) ){
  help()
} else {
  graftBam     <- sub('--RDS=', '', args[grep('--RDS=', args)])
  hostBam  <- sub( '--outDir=', '',args[grep('--outDir=',args)])
}
sampname <- gsub("_bam", "", unlist(strsplit(graftBam,"/"))[3])
outDir <- paste("samples", "xeno", sampname, sep="/")
sample.list <- as.data.frame(cbind(graftBam, hostBam))
bp.param <- SnowParam(workers = 1, type = "SOCK")

dir.create(outDir, recursive=TRUE)
XenofilteR(sample.list, destination.folder = outDir, bp.param = bp.param)

