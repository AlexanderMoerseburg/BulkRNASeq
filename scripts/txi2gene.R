library(tximport)
library(rjson)
library(readr)
library("edgeR")
library("gplots")

library("GenomicFeatures")

args <- commandArgs(trailingOnly = TRUE)
annotation = args[1]
outfile = paste(args[2])

txdb <- makeTxDbFromGFF(file=annotation)
saveDb(x=txdb, file = "annotation.TxDb")
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME") 
head(tx2gene)
write.csv(tx2gene, outfile)
