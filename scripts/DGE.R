library(readr)
library("edgeR")
library("gplots")

myoutput ="IgG_AlphaB9"
mytitle ="IgG vs AlphaB9"
gr1 ="igg"
gr2 ="b9" 
G1 = "IgG" 
G2 = "AlphaB9" 

group = factor(c(rep(G1, 6), rep(G2,6)) )

files <- c(paste(gr1,sep = "","1.counts.txt"),paste(gr1,sep = "","2.counts.txt"),paste(gr1,sep = "","3.counts.txt"),paste(gr1,sep = "","4.counts.txt"),paste(gr1,sep = "","5.counts.txt"), paste(gr1,sep = "","6.counts.txt"),paste(gr2,sep = "","1.counts.txt"), paste(gr2,sep = "","2.counts.txt"),paste(gr2,sep = "","3.counts.txt"),paste(gr2,sep = "","4.counts.txt"), paste(gr2,sep = "","5.counts.txt"), paste(gr2,sep = "","6.counts.txt") ) 

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
et <- exactTest(dge, pair=c(G1, G2)) 
etp <- topTags(et,n=100000000,adjust.method="BH", sort.by="PValue", p.value = 1)

write.csv(countsPerMillion, paste(myoutput,sep="","_cpm.csv"))  

labels=c(paste(G1,sep="","1"), paste(G1,sep="","2"), paste(G1,sep="","3"), paste(G1,sep="","4"),paste(G1,sep="","5"), paste(G1,sep="","6"), paste(G2,sep="","1"), paste(G2,sep="","2"), paste(G2,sep="","3"), paste(G2,sep="","4"), paste(G2,sep="","5"), paste(G2,sep="","6") ) 
pdf(paste(myoutput,sep="","_MA.pdf")) 
plot(
  etp$table$logCPM,
  etp$table$logFC,
  xlim=c(-3, 20), ylim=c(-12, 12), pch=20, cex=.3,
  col = ifelse( etp$table$FDR < .2, "black", "red" ) )
dev.off()
pdf(paste(myoutput,sep="","_MDS.pdf"))
plotMDS(dge, labels=labels)
dev.off()

write.csv(etp$table, paste(myoutput,sep="","_dge.csv")) 

#-----------------Raw counts heatmap--------------


#pdf(paste(myoutput,sep="","_counts_heatmap.pdf"))
#library(gplots)
#hclustfunc <- function(x) hclust(x, method="centroid")
#distfunc <- function(x) dist(x,method="euclidean")


#dat = countsPerMillion
#data <- na.omit(dat)
#myCol <- colorRampPalette(c("dodgerblue", "black", "yellow"))(10)
#colnames(data) = labels
#myBreaks <- seq(0, 10, length.out=11)
##data <- t(scale(t(data))) 
#heatmap.2(data, col=myCol, breaks=myBreaks,Rowv=TRUE,Colv=TRUE, main=mytitle, key=T, keysize=0.7,scale="none",trace="none", reorderfun=function(d,w) reorder(d, w, agglo.FUN=mean), distfun = distfunc, hclustfun= hclustfunc, margin=c(10,9),cexCol =1.4, cexRow=1, labRow = FALSE,lhei=c(2,10), lwid=c(2,6))
#dev.off()
 

#----pdf("star_volcano_BL_Age.pdf")
pdf(paste(myoutput,sep="","_volcano.pdf")) 
res <- etp$table
# Make a basic volcano plot
#log2FoldChange = etp$logFC
#Pvalue = etp$PValue
#head(Pvalue)
#FDR = etp$FDR 
with(res, plot(logFC, -log10(PValue), pch=20, main=mytitle, xlim=c(-12,12)))

# Add colored points: red if FDR<0.05, orange of log2FC>1, green if both)
with(subset(res, FDR<.05 ), points(logFC, -log10(PValue), pch=20, col="red"))
with(subset(res, abs(logFC)>1), points(logFC, -log10(PValue), pch=20, col="orange"))
with(subset(res, FDR<.05 & abs(logFC)>1), points(logFC, -log10(PValue), pch=20, col="green"))
# Label points with the textxy function from the calibrate plot
#library(calibrate)
#with(subset(res, FDR<.05 & abs(logFC)>1), textxy(logFC, -log10(PValue), labs=Gene, cex=.8))
dev.off()
#------

pdf(paste(myoutput,sep="","_heatmap.pdf"))
##logCPM = etp$table 
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

