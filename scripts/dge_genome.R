library("edgeR")
library("gplots")

args <- commandArgs(trailingOnly = TRUE)
control = args[1]
treat= args[2]
N = args [3]
ORGANISM  =args[4] 

control_list = list.files(path = ".",sprintf("*.%s.txt",control) )
control_list

treat_list = list.files(path = ".",sprintf("*.%s.txt",treat) )
treat_list

group = factor(c(rep(args[1], N), rep(args[2],N)) )

files <- c(control_list, treat_list) 
files 

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
write.csv(countsPerMillion, out,quote=FALSE) 
out = paste(outname,"dge.csv", sep="_")
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
with(subset(res, FDR<0.05 ), points(logFC, -log10(PValue), pch=20, col="red"))
with(subset(res, abs(logFC)> 2), points(logFC, -log10(PValue), pch=20, col="orange"))
with(subset(res, FDR<0.05 & abs(logFC)> 2), points(logFC, -log10(PValue), pch=20, col="green"))
dev.off()
#-------------------------------
#Heatmap Plot for log CPM  
#-------------------------------
fig = paste(outname,"heatmap.pdf", sep="-")
title = paste(outname, "Heatmap", sep=" ")
pdf(fig)
logCPM =countsPerMillion
o = rownames(etp$table[abs(etp$table$logFC)>1 & etp$table$FDR<0.05, ])
logCPM <- logCPM[o[1:50],]
colnames(logCPM) = labels
logCPM <- t(scale(t(logCPM)))
require("RColorBrewer")
require("gplots")
myCol <- colorRampPalette(c("dodgerblue", "black", "yellow"))(25)
myBreaks <- seq(-3, 3, length.out=26)
heatmap.2(logCPM, col=myCol, breaks=myBreaks,symkey=F, Rowv=TRUE,Colv=TRUE, main=title, key=T, keysize=0.7,scale="none",trace="none", dendrogram="both", cexRow=0.7, cexCol=0.9, density.info="none",margin=c(10,9), lhei=c(2,10), lwid=c(2,6),reorderfun=function(d,w) reorder(d, w, agglo.FUN=mean),  distfun=function(x) as.dist(1-cor(t(x))), hclustfun=function(x) hclust(x, method="ward.D2"))
dev.off()


