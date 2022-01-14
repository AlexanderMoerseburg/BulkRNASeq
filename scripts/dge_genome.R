library("edgeR")
library(biomaRt)
library(dplyr)
library("AnnotationDbi")
library("org.Hs.eg.db")
library(data.table)
library("gplots")
library(magrittr)
library(pathview)
library("org.Mm.eg.db")
library(gage)
library(gageData)
library(data.table)


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
# Label points with the textxy function from the calibrate plot
#library(calibrate)
#with(subset(res,FDR<.05 & abs(logFC)>1), textxy(logFC, -log10(PValue), labs=Gene, cex=.8))
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
logCPM <- logCPM[o[1:50],]
colnames(logCPM) = labels
logCPM <- t(scale(t(logCPM)))
require("RColorBrewer")
require("gplots")
myCol <- colorRampPalette(c("dodgerblue", "black", "yellow"))(25)
myBreaks <- seq(-3, 3, length.out=26)
heatmap.2(logCPM, col=myCol, breaks=myBreaks,symkey=F, Rowv=TRUE,Colv=TRUE, main=title, key=T, keysize=0.7,scale="none",trace="none", dendrogram="both", cexRow=0.7, cexCol=0.9, density.info="none",margin=c(10,9), lhei=c(2,10), lwid=c(2,6),reorderfun=function(d,w) reorder(d, w, agglo.FUN=mean),  distfun=function(x) as.dist(1-cor(t(x))), hclustfun=function(x) hclust(x, method="ward.D2"))
dev.off()

#---------------------------
#KEGG Analysis 
#---------------------------
if (ORGANISM == "HUMAN")
{
data(kegg.sets.hs)
data(go.sets.hs)
data(carta.hs)
data(sigmet.idx.hs)
data(go.subs.hs)
kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs]
head(kegg.sets.hs,3)
foldchanges = res$logFC
names(foldchanges) = res$entrez
#---------------------------------------------------Use Kegg and gage to get upregulated and downregulated pathways
out = paste(outname, "KEGG.csv", sep="-")
outUP = paste(outname, "KEGG_UP.csv", sep ="-")
outDOWN = paste(outname, "KEGG_DOWN.csv", sep ="-")

data(kegg.gs)
keggres = gage(foldchanges, gsets =kegg.sets.hs, same.dir = TRUE, compare="paired",make.plot = TRUE)
lapply(keggres, head)
write.csv(keggres,out)
keggrespathwaysup = data.frame(id=rownames(keggres$greater), keggres$greater) %>%
  tbl_df() %>%
  filter(row_number()<=30) %>%
  .$id %>%
  as.character()
keggrespathwaysdn = data.frame(id=rownames(keggres$less), keggres$less) %>%
  tbl_df() %>%
  filter(row_number()<=30) %>%
  .$id %>%
  as.character()
write.csv(keggrespathwaysup, outUP)
write.csv(keggrespathwaysdn, outDOWN)
#-------------------------------------------------------------------------------------------------------------------------------
keggresidsup = substr(keggrespathwaysup, start=1, stop=8)
keggresidsup
keggresidsdn = substr(keggrespathwaysdn, start=1, stop=8)
#gobpres = gage(foldchanges, gsets=kegg.sets.hs, same.dir =FALSE, compare ="paired",cutoff=0.05)
#lapply(gobpres, head)
#---------------------------------------------------------Define plotting function for applying later
plot_pathway = function(pid) pathview(gene.data=foldchanges,gene.idtype="ENTREZID", pathway.id=pid, species="human", new.signature=FALSE)
#---------------------------------------------------------plot multiple pathways ups and downs 
tmpup = sapply(keggresidsup, function(pid) pathview(gene.data=foldchanges,gene.idtype="ENTREZID", pathway.id=pid, species="human"))
tmpdn = sapply(keggresidsdn, function(pid) pathview(gene.data=foldchanges,gene.idtype="ENTREZID", pathway.id=pid, species="human"))
} else if (ORGANISM == "MOUSE"){
res$symbol <- mapIds(org.Mm.eg.db,
                     keys=row.names(res),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
res$entrez <- mapIds(org.Mm.eg.db,
                     keys=row.names(res),
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")

data(kegg.sets.mm)
data(go.sets.mm)
data(carta.mm)
data(sigmet.idx.mm)
data(go.subs.mm)


kegg.sets.mm = kegg.sets.mm[sigmet.idx.mm]
head(kegg.sets.mm,3)
foldchanges = res$logFC
names(foldchanges) = res$entrez
#---------------------------------------------------Use Kegg and gage to get upregulated and downregulated pathways
data(kegg.gs)
keggres = gage(foldchanges, gsets =kegg.sets.mm, same.dir = TRUE, compare="paired",make.plot = TRUE)
lapply(keggres, head)
write.csv(keggres,out)
keggrespathwaysup = data.frame(id=rownames(keggres$greater), keggres$greater) %>%
  tbl_df() %>%
  filter(row_number()<=30) %>%
  .$id %>%
  as.character()
keggrespathwaysdn = data.frame(id=rownames(keggres$less), keggres$less) %>%
  tbl_df() %>%
  filter(row_number()<=30) %>%
  .$id %>%
  as.character()
write.csv(keggrespathwaysup, outUP)
write.csv(keggrespathwaysdn, outDOWN)
#-------------------------------------------------------------------------------------------------------------------------------
keggresidsup = substr(keggrespathwaysup, start=1, stop=8)
keggresidsup
keggresidsdn = substr(keggrespathwaysdn, start=1, stop=8)
#gobpres = gage(foldchanges, gsets=kegg.sets.mm, same.dir =FALSE, compare ="paired",cutoff=0.05)
#lapply(gobpres, head)
#---------------------------------------------------------Define plotting function for applying later
plot_pathway = function(pid) pathview(gene.data=foldchanges,gene.idtype="ENTREZID", pathway.id=pid, species="mouse", new.signature=FALSE)
#---------------------------------------------------------plot multiple pathways ups and downs 
tmpup = sapply(keggresidsup, function(pid) pathview(gene.data=foldchanges,gene.idtype="ENTREZID", pathway.id=pid, species="mouse"))
tmpdn = sapply(keggresidsdn, function(pid) pathview(gene.data=foldchanges,gene.idtype="ENTREZID", pathway.id=pid, species="mouse"))
} 

