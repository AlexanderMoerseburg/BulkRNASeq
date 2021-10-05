library(readr)
library("edgeR")
library("gplots")

group = factor(c(rep(G1, 6), rep(G2,6)) )

count.control <- read.table(paste(control, '.counts.txt', sep = ''), header = TRUE, row.names = 1)
count.treat <- read.table(paste(treat, '.counts.txt', sep = ''), header = TRUE, row.names = 1)
count.table <- cbind(count.control, count.treat)  # merge the control and treat tables together

num.sample <- ncol(count.table)
num.sample.control <- ncol(count.control)
num.sample.treat <- ncol(count.treat)

sample.control <- colnames(count.control)
sample.treat <- colnames(count.treat)
  
gene.list <- rownames(count.table)
  
samples <- colnames(count.table)

subject <- factor(subject.all[c(which(group.all == control), which(group.all == treat))])
group <- relevel(factor(group.all[c(which(group.all == control), which(group.all == treat))]), ref = control)

y.control <- DGEList(counts = count.control, genes = gene.list)
y.treat <- DGEList(counts = count.treat, genes = gene.list)

y.control <- calcNormFactors(y.control, method="TMM")
count.table.control.norm <- cpm(y.control)
write.table(count.table.control.norm, paste(output.path, '/countGroup/', control, '_gene_norm.tsv', sep = ''), quote = FALSE, sep = "\t")

y.treat <- calcNormFactors(y.treat, method="TMM")
count.table.treat.norm <- cpm(y.treat)
write.table(count.table.treat.norm, paste(output.path, '/countGroup/', treat, '_gene_norm.tsv', sep = ''), quote = FALSE, sep = "\t")

dge <- DGEList(counts = count.table, genes = gene.list)

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
et <- exactTest(dge, pair=c(G1, G2)) 
etp <- topTags(et,n=100000000,adjust.method="BH", sort.by="PValue", p.value = 1)

write.csv(countsPerMillion, paste(myoutput,sep="","_cpm.csv"))  
write.csv(etp$table, paste(myoutput,sep="","_dge.csv")) 

#-----------------Raw counts heatmap--------------


pdf(paste(myoutput,sep="","_counts_heatmap.pdf"))
library(gplots)
hclustfunc <- function(x) hclust(x, method="centroid")
distfunc <- function(x) dist(x,method="euclidean")


dat = countsPerMillion
data <- na.omit(dat)
myCol <- colorRampPalette(c("dodgerblue", "black", "yellow"))(10)
colnames(data) = labels
myBreaks <- seq(0, 10, length.out=11)
data <- t(scale(t(data))) 
heatmap.2(data, col=myCol, breaks=myBreaks,Rowv=TRUE,Colv=TRUE, main=mytitle, key=T, keysize=0.7,scale="none",trace="none", reorderfun=function(d,w) reorder(d, w, agglo.FUN=mean), distfun = distfunc, hclustfun= hclustfunc, margin=c(10,9),cexCol =1.4, cexRow=1, labRow = FALSE,lhei=c(2,10), lwid=c(2,6))
dev.off()
 

#----Volcano plot
pdf(paste(myoutput,sep="","_volcano.pdf")) 
res <- etp$table
# Make a basic volcano plot
log2FoldChange = etp$logFC
Pvalue = etp$PValue
head(Pvalue)
FDR = etp$FDR 
with(res, plot(logFC, -log10(PValue), pch=20, main=mytitle, xlim=c(-12,12)))

# Add colored points: red if FDR<0.05, orange of log2FC>1, green if both)
with(subset(res, FDR<.05 ), points(logFC, -log10(PValue), pch=20, col="red"))
with(subset(res, abs(logFC)>1), points(logFC, -log10(PValue), pch=20, col="orange"))
with(subset(res, FDR<.05 & abs(logFC)>1), points(logFC, -log10(PValue), pch=20, col="green"))
dev.off()
#------

pdf(paste(myoutput,sep="","_heatmap.pdf"))
logCPM =countsPerMillion
o = rownames(etp$table[abs(etp$table$logFC)>1 & etp$table$FDR<0.05, ])
o
logCPM <- logCPM[o[1:4],]
colnames(logCPM) = labels
logCPM <- t(scale(t(logCPM)))
require("RColorBrewer")
require("gplots")
myCol <- colorRampPalette(c("dodgerblue", "black", "yellow"))(100)
myBreaks <- seq(-3, 3, length.out=101)
heatmap.2(logCPM, col=myCol, breaks=myBreaks,symkey=F, Rowv=TRUE,Colv=TRUE, main=mytitle, key=T, keysize=0.7,scale="none",trace="none", dendrogram="both", cexRow=0.7, cexCol=0.9, density.info="none",margin=c(10,9), lhei=c(2,10), lwid=c(2,6),reorderfun=function(d,w) reorder(d, w, agglo.FUN=mean),  distfun=function(x) as.dist(1-cor(t(x))), hclustfun=function(x) hclust(x, method="ward.D2"))
dev.off()

