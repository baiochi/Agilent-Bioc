# Agilent - pre-miR101 Microarray analysis

#load functions
library(limma)
source("/Users/Baiochi/Desktop/Agilent-Bioc/Scripts/Functions.R")

#read targets
targets <- read.table(file="/Users/Baiochi/Dropbox/USP/Lab/RawData/HCT/Targets.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)

#Read using ProcessedSignal raw <- readAFE(targets = targets, path = "/Users/Baiochi/Dropbox/USP/Lab/RawData/HCT", skip.lines =  9)

#Read using MedianSignal
raw = read.maimages(targets$FileName, source="agilent", green.only=TRUE, path="/Users/Baiochi/Dropbox/USP/Lab/RawData/HCT")
colnames(raw$E) <- targets$Nomeclature
rownames(raw$targets) <- targets$Nomeclature

# Start Analysis ----------------------------------------------------------

#correct artifacts
raw$Eb[which(raw$Eb[,3]>60),3] <- round(mean(raw$Eb[,3]))
raw$Eb[which(raw$Eb[,6]>60),6] <- round(mean(raw$Eb[,6]))
#background correct
bgc = backgroundCorrect(raw,method='normexp')
#normalize
norm = normalizeBetweenArrays(bgc,method='quantile')
#filter control probes
eset = norm[norm$genes$ControlType==0,]
#avg probes
eset = avereps(eset,ID=eset$genes[,"SystematicName"])

#design
f <- factor(targets$Condition, levels = unique(targets$Condition))
design = model.matrix(~0 + f)
colnames(design) <- levels(f)
#contrasts
contrast.matrix = makeContrasts(contrasts='PM101-Ctr', levels=design)
#fit model
fit <- lmFit(eset$E, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
#bayes statistics
fit2 = eBayes(fit2)

#rank genes
exprs = topTable(fit2, adjust="fdr", coef='PM101-Ctr', genelist=eset$genes, number=Inf)
rank = decideTests(fit2, method="separate", adjust.method="fdr", p.value=0.05)
summary(rank)

#diff exprs transcripts, p<0.05 fc<0.7
diff = topTable(fit2, adjust="fdr", coef='PM101-Ctr', genelist=eset$genes, p.value=0.05, lfc=0, number=Inf)




# Plotting ----------------------------------------------------------------

#Setup printer
raw$genes$Block <- 1
names(raw$genes)[2] <- "Column"
raw$printer <- getLayout(raw$genes)
r <- raw$genes$Row
c <- raw$genes$Col
nr <- max(r)
nc <- max(c)
y <- rep(NA,nr*nc)
j <- (r-1)*nc+c
#raw data visualization
boxplot(log2(raw$E), main='Raw data Boxplot', col='forestgreen')
hierclust(log2(raw$E), methdis="euclidean", methclu="complete", sel=FALSE, size=100)
cor.matrix(log2(raw$E), title="Pearson Correlation Matrix")
plotDensities(raw, legend = 'topright', main='Raw Foreground densities')
par(mfrow=c(2,3))
for(i in 1:6)
  plotMA(raw, array=i, main=paste('Array', i, '- Raw'))
par(mfrow=c(2,3),oma = c(0, 0, 2, 0))
for(i in 1:6)
  plotFB(raw, array=i)
mtext('Foreground vs Background', outer = TRUE)

#post-norm visualization
boxplot(eset$E, main='Normalized Boxplot', col='deepskyblue2')
plotDensities(eset, legend = 'topright', main='Normalized Foreground Densities')
MAPlot(eset, 2, 3, title = "Post-Normalization MA Plot")

#differential expression visualization
plotSA(fit2)
plotMA(fit2)
volcanoplot(fit2, 'PM101-Ctr',title='pre-miR101 vs Control VolcanoPlot')
diffexprs.hist(exprs, coef='PM101-Ctr')



#g <- diff[grep('example', diff$GeneName),]
#text(x=g$logFC, y=-log2(g$adj.P.Val), label=g$GeneName, pos=3, offset=0.2)
#points(x=g$logFC, y=-log2(g$adj.P.Val), pch=20, col='deepskyblue2')
#pval hist

# Database Cross ----------------------------------------------------------

#read data from miRanda
load(file='/Users/Baiochi/Dropbox/USP/Lab/RawData/miRanda/predictGood.rda')
#read data from miRTarBase
tarBaseMirs <- read.delim(file='/Users/Baiochi/Dropbox/USP/Lab/RawData/miRTarBase/hsa_MTI.txt',header=TRUE,sep='\t', stringsAsFactors = FALSE)


#get only hsa-miR-101-*
mir101_tar <- tarBaseMirs[grep(pattern = 'hsa-miR-101-', tarBaseMirs$miRNA),]
#unique genes Array vs miRTarBase = 82
table(unique(diff$GeneName) %in% mir101_tar$Target.Gene)
#all genes Array vs miRTarBase = 87
table(diff$GeneName %in% mir101_tar$Target.Gene)
#Result Table: GeneName, ProbeName, Systematic Name, p.val, adj.P.Val, FC, Start, Sequence, Description
validatedGene <- diff[which((diff$GeneName %in% mir101_tar$Target.Gene)),c(9,8,7,14,15,11,3,4,10)]
validatedGene <- validatedGene[with(validatedGene, order(adj.P.Val, logFC)), ]
rownames(validatedGene) <- NULL

#get only hsa-miR-101*
mir101_miranda <- predictGood[grep('hsa-miR-101', predictGood$mirna_name),]
#all transcripts Array vs miRanda = 808
table(diff$SystematicName %in% mir101_miranda$ext_transcript_id)
#all genes Array vs miRanda = 1252
table(diff$GeneName %in% mir101_miranda$gene_symbol)
#Result Table for transcripts
predictTranscript <- diff[which((diff$SystematicName %in% mir101_miranda$ext_transcript_id)),c(9,8,7,14,15,11,3,4,10)]
predictTranscript <- predictTranscript[with(predictTranscript, order(adj.P.Val, logFC)), ]
rownames(predictTranscript) <- NULL
#Result Table for genes
predictGene <- diff[which((diff$GeneName %in% mir101_miranda$gene_symbol)),c(9,8,7,14,15,11,3,4,10)]
predictGene <- predictGene[with(predictGene, order(adj.P.Val, logFC)), ]
rownames(predictGene) <- NULL

#VennDiagram
venn.diagram(x=list(Microarray=diff$SystematicName, 
                    miRTarBase=validatedGene$SystematicName, 
                    miRanda=predictTranscript$SystematicName),file='Venn-AllTranscripts.tiff')
#Down/UP reg [which(exprs$logFC<0)]
#venn.diagram(x=list(Microarray=diff$GeneName[which(diff$logFC<0)], 
#                    miRTarBase=validatedGene$GeneName[which(validatedGene$logFC<0)], 
#                   miRanda=predictTranscript$GeneName[which(predictTranscript$logFC<0)]),file='Venn-DownGenestiff')
#complete diagram
venn.diagram(x=list(Microarray=diff$GeneName[which(diff$logFC<0)],
                    miRanda=mir101_miranda$gene_symbol,
                    miRTarBase=mir101_tar$Target.Gene),file='down_regulated_genes_veendiagram.tiff')

mirandagenes = mir101_miranda$gene_symbol
targenes = mir101_tar$Target.Gene
arraygenes = diff$GeneName[which(diff$logFC<0)]
intergenes <- Reduce(intersect, list(mirandagenes,targenes,arraygenes))
array_miranda <- intersect(mir101_miranda$ext_transcript_id, diff$SystematicName[which(diff$logFC<0)])
write.table(array_miranda, file='Array_miranda_transcript_list.txt', quote=FALSE, sep = '\t', col.names=FALSE, row.names=FALSE)
write.table(intergenes, file='intersect_genes_list.txt', quote=FALSE, sep = '\t', col.names=FALSE, row.names=FALSE)




# Save Plots --------------------------------------------------------------

RawPlots(raw)
NormPlots(eset, sn=2, rep=3)
DiffExprsPlots(fit2, 'PM101-Ctr', rank, eset)

# Results -----------------------------------------------------------------

#annotation
WriteResults(fit2, 'PM101-Ctr', eset, 'miR-101', lists = TRUE)
#Save R Data
saveData(path='/Users/Baiochi/Dropbox/USP/Lab/Results/pre-miR101/RData/',raw,bgc,norm,eset,design,fit2,exprs,rank)
#Load Data
loadData(path='/Users/Baiochi/Dropbox/USP/Lab/Results/pre-miR101/RData')


#-------------------NOT RUN-------------------#
###############################################

#mirtarbase target genes in volcanoplot
points(x=exprs$logFC[which(exprs$GeneName %in% TARgenes)], 
       y=-log(exprs$adj.P.Val[which(exprs$GeneName %in% TARgenes)]), 
       pch=20, cex=.5, col='deepskyblue2')


#foldchange check
diff = topTable(fit2, adjust="fdr", coef='PM101-Ctr', genelist=eset$genes, p.value=0.05, lfc=2, number=Inf)
diff <- diff[,c(3,4,7,8,9,11,14,15)]
rownames(diff) <- NULL
crc <- read.delim(file='/Users/Baiochi/Dropbox/USP/Lab/Results/pre-miR101/miRTarBase/genelist_UP_mir101_TarBase.txt', header=FALSE, sep='\n', stringsAsFactors = FALSE)
d <- diff[diff$GeneName %in% crc$V1,]
d.int <- eset$E[eset$genes$SystematicName %in% d$SystematicName,]
d <- d[with(d, order(SystematicName)),]
d.int <- d.int[order(rownames(d.int)),]
d.int <- cbind(d.int, log=d$logFC)
aux <- numeric()
for(i in 1:dim(d.int)[1]){
  aux[i] <- mean(d.int[i,1:3]) - mean(d.int[i,4:6])
}
d.int <- cbind(d.int, hand.log2=aux)

highlight(d$SystematicName, diff)

require(ggplot2)
gene_list <- exprs
##Highlight genes that have an absolute fold change > 2 and a p-value < Bonferroni cut-off
gene_list$threshold = as.factor(gene_list$adj.P.Val < 0.05)

##Construct the plot object
g = ggplot(data=gene_list, aes(x=logFC, y=-log10(P.Value), colour=threshold)) +
  geom_point(alpha=0.8, size=1.75) +
  xlim(c(-7, 5)) + ylim(c(0, 8)) +
  xlab("log2 fold change") + ylab("-log10 p-value")
g







#4 last columns
outlierProbes <- as.data.frame(which(rawdata$genes[,9:12]==1, arr.ind = TRUE))
outlierProbes <- unique(outlierProbes[,1])
outlierGenes <- new('character')
for(i in 9:12){
  aux <- rawdata$genes$GeneName[which(rawdata$genes[,i]==1)]
  outlierGenes <- c(outlierGenes, aux)
}

