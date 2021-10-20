library(tximport)
library(rjson)
library(readr)
library("edgeR")
library("gplots")


args <- commandArgs(trailingOnly = TRUE)
control = args[1]
treat = args[2]
N = args [3]

files = list.files(path = ".",sprintf("*.quant.sf") )
files



#samples <- read.table(file.path("/home/sa941/RNASeq/SLX-16108", "samplestwo2.txt"), header = FALSE)
#files <- file.path("/home/sa941/RNASeq", "SLX-16108", samples$V1, "quant.sf")

names(files) <- paste0("sample", 1:4)
print(names)
file.exists(files)
tr2gene <- read.csv("tx2gene.csv")
#txi <- tximport(files, type = "salmon",txIn = TRUE, tx2gene = tr2gene)
txi <- tximport(files, type = "salmon", txIn = TRUE, txOut = TRUE, countsFromAbundance = "scaledTPM") #, tx2gene = NULL, reader = read.delim, geneIdCol, txIdCol, abundanceCol, countsCol, lengthCol, importer, collatedFiles, ignoreTxVersion = FALSE)

names(txi)
head(txi$counts)

outname = paste(control,treat, sep="_")
out = paste(outname, "salmon.counts.csv", sep="_") 
write.csv(txi$counts, out)

cts <- txi$counts

group = factor(c(rep(control,N), rep(treat,N)))
dge = DGEList(counts=cts, genes= rownames(data), group=group)
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
#dge <- estimateTrendedDisp(dge)
et <- exactTest(dge, pair=c(control, treat)) 
etp <- topTags(et, n= 100000, adjust.method="BH", sort.by="PValue", p.value = 1)

out = paste(outname,"salmon.cpm.csv", sep="_")
write.csv(etp$table, out)

