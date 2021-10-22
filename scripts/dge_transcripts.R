library(tximport)
library(rjson)
library(readr)
library("edgeR")
library("gplots")
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
control = args[1]
treat = args[2]
N = args [3]
tximport_file = args[4] 
files = list.files(path = ".",sprintf("*.quant.sf") )
files
no_samples = as.integer(N) * 2
names(files) <- paste0("sample", 1:no_samples)
print(names)
file.exists(files)
tr2gene <- read.csv(tximport_file)
head(tr2gene)
tx2gene <- select(tr2gene,TXNAME,GENEID)
head(tx2gene)
txi <- tximport(files, type = "salmon", txIn = TRUE, txOut = FALSE, tx2gene = tx2gene, countsFromAbundance = "scaledTPM", ignoreTxVersion=TRUE, geneIdCol="GENEID", txIdCol="TXNAME") #Rest of tximport paramters include: tx2gene = NULL, reader = read.delim, geneIdCol, txIdCol, abundanceCol, countsCol, lengthCol, importer, collatedFiles, ignoreTxVersion = FALSE)

names(txi)
head(txi$counts)

outname = paste(control,treat, sep="-")
out = paste(outname, "salmon.counts.csv", sep="-") 
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
#removed for small size sample, uncomment otherwise
#dge <- estimateTrendedDisp(dge)
et <- exactTest(dge, pair=c(control, treat)) 
etp <- topTags(et, n= 100000, adjust.method="BH", sort.by="PValue", p.value = 1)

out = paste(outname,"salmon.cpm.csv", sep="-")
write.csv(etp$table, out)

#Plotting starts here 
for (i in 1:N) {
labels <- append(labels, paste(control,i, sep="-"))
}

for (i in 1:N) {
labels <- append(labels, paste(treat,i, sep="-"))
}

labels <- c(labels[2:5])
print( length(labels))  
#-------------------------------
#MA Plot
#-------------------------------
fig = paste(outname,"MA.pdf", sep="-")
title = paste(outname, "MA Plot", sep=" ")
pdf(fig)
plot(
  etp$table$logCPM,
  etp$table$logFC,
  xlim=c(-3, 20), ylim=c(-12, 12), pch=20, cex=.3,
  col = ifelse( etp$table$FDR <1.1, "black", "red" ), xlab="Log CPM", ylab ="Log FC", main =title )
dev.off()
#-------------------------------
#MDS Plot 
#-------------------------------
fig = paste(outname,"MDS.pdf", sep="-")
title = paste(outname, "MDS Plot", sep =" ")
pdf(fig)
plotMDS(dge, labels=labels, main=title)
dev.off()
#-------------------------------
#Volcano Plot
#-------------------------------
fig =paste(outname,"volcano.pdf", sep="-")
title = paste(outname, "Volcano Plot", sep= " ")
pdf(fig)
res <- etp$table
with(res, plot(logFC, -log10(PValue), pch=20, main=title, xlim=c(-12,12)))
# Add colored points: red if FDR<0.05, orange of log2FC>1, green if both)
with(subset(res, FDR<.05 ), points(logFC, -log10(PValue), pch=20, col="red"))
with(subset(res, abs(logFC)>1), points(logFC, -log10(PValue), pch=20, col="orange"))
with(subset(res, FDR<.05 & abs(logFC)>1), points(logFC, -log10(PValue), pch=20, col="green"))
# Label points with the textxy function from the calibrate plot
#library(calibrate)
with(subset(res,FDR<.05 & abs(logFC)>1), textxy(logFC, -log10(PValue), labs=Gene, cex=.8))
dev.off()
#-------------------------------
#Heatmap Plot for log CPM  
#-------------------------------
fig = paste(outname,"heatmap.pdf", sep="-")
title = paste(outname, "Heatmap", sep=" ") 
pdf(fig)
#logCPM = etp$table 
logCPM =countsPerMillion
o = rownames(etp$table[abs(etp$table$logFC)>1 & etp$table$PValue<0.05, ])
logCPM <- logCPM[o[1:25],]
colnames(logCPM) = labels
logCPM <- t(scale(t(logCPM)))
require("RColorBrewer")
require("gplots")
myCol <- colorRampPalette(c("dodgerblue", "black", "yellow"))(25)
myBreaks <- seq(-3, 3, length.out=26)
heatmap.2(logCPM, col=myCol, breaks=myBreaks,symkey=F, Rowv=TRUE,Colv=TRUE, main=title, key=T, keysize=0.7,scale="none",trace="none", dendrogram="both", cexRow=0.7, cexCol=0.9, density.info="none",margin=c(10,9), lhei=c(2,10), lwid=c(2,6),reorderfun=function(d,w) reorder(d, w, agglo.FUN=mean),  distfun=function(x) as.dist(1-cor(t(x))), hclustfun=function(x) hclust(x, method="ward.D2"))
dev.off()

