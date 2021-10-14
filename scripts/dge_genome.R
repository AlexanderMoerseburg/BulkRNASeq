#library(tximport)
#library(rjson)
#library(readr)
library("edgeR")
#library("gplots")
#library(biomaRt)
#library(dplyr) 
#library("AnnotationDbi")
#library("org.Hs.eg.db")
#library(magrittr)
#library(pathview)
#library(gage)
#library(gageData)
#library(data.table)

args <- commandArgs(trailingOnly = TRUE)
control = args[1]
treat = args[2]
N = args [3]

control_list = list.files(path = ".",sprintf("*.%s.txt",control) )
control_list

treat_list = list.files(path = ".",sprintf("*.%s.txt",treat) )
treat_list

group = factor(c(rep(args[1], N), rep(args[2],N)) )
files <- c(control_list, treat_list) 

dge  <- readDGE(files,header=FALSE, group =group)

countsPerMillion <- cpm(dge)
summary(countsPerMillion)
countCheck <- countsPerMillion > 1 
summary(countCheck)
keep <- which(rowSums(countCheck) >= 2)
dge <- dge[keep,]
summary(cpm(dge))
dge <- calcNormFactors(dge, method="TMM")
dge <- estimateCommonDisp(dge)
dge <- estimateTagwiseDisp(dge)
dge <- estimateTrendedDisp(dge)
et <- exactTest(dge, pair=c(args[1], args[2])) 
etp <- topTags(et, n= 100000, adjust.method="BH", sort.by="PValue", p.value = 1)

outname = paste(control,treat, sep="_")
out = paste(outname,"cpm.csv", sep="_")
write.csv(countsPerMillion, out) 
out = paste(outname,"dge.csv", sep="_")
write.csv(etp$table, out)

