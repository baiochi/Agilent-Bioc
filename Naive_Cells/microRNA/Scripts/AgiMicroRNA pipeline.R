#AgiMicroRna Pipeline



source("/Users/Baiochi/Desktop/Agilent-Bioc/Scripts/Functions.R")
library(AgiMicroRna)
library(limma)
library(marray)
library(gplots)

#------------AgiMicroRna analysis--------#
#     iTreg-Naive iTreg-Teff Naive-Teff
# -1          46         14         14
# 0          736        801        820
# 1           69         36         17

# 23a, 30a, 636, 1299, - 422a, 610


#Create Target File
targets <- read.table(file="/Users/Baiochi/Dropbox/USP/Lab/RawData/Naive microRNA/RawData/Targets.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
targets.micro <- data.frame(FileName=targets$FileName, 
                            Treatment=rep(c("iTreg", "Naive", "Teff"), each=3) , 
                            GErep=rep(c(1,2,3), each=3), 
                            Subejct=rep(c(1,2,3),3), 
                            stringsAsFactors=FALSE)
tempdir <- getwd()
setwd('/Users/Baiochi/Dropbox/USP/Lab/RawData/Naive microRNA/RawData')
dd.micro=readMicroRnaAFE(targets.micro,verbose=TRUE)
setwd(tempdir)
colnames(dd.micro$TGS) <- targets$Nomeclature
rownames(dd.micro$targets) <- targets$Nomeclature

#Extract TotalGeneSignal
tgs <- tgsMicroRna(dd.micro, half=FALSE, makePLOT=FALSE, verbose=TRUE)
tgs <- tgs[grep('hsa',tgs$genes$GeneName),]

#Normalization
norm <- tgsNormalization(tgs, "quantile", makePLOTpre=FALSE, 
                        makePLOTpost=FALSE, targets.micro, verbose=FALSE)

#Create Expression Set
esetPROC <- esetMicroRna(norm, targets.micro, makePLOT=TRUE, verbose=TRUE)

#Compute Differential Expression
f <- factor(targets$Condition, levels = unique(targets$Condition))
design <- model.matrix(~0 + f)
colnames(design) <- levels(f)
fit <- lmFit(esetPROC, design)
contrast.matrix <- makeContrasts('CD4tgf-Naive','CD4tgf-CD4med','CD4med-Naive', levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
fit2$genes <- norm$genes

#Rank Genes
rank <- decideTests(fit2, method="separate", adjust.method="fdr", p.value=0.05, lfc=0)
DE <- getDecideTests(fit2, DEmethod='separate', MTestmethod='BH', PVcut=0.05,verbose=FALSE)
summary(DE)

#Statistical Results
exprs = topTable(fit2, adjust="fdr", coef='CD4tgf-Naive', number=Inf
                #, p.value = 0.05
                 )
exprs2 = topTable(fit2, adjust="fdr", coef='CD4tgf-CD4med', number=Inf)
exprs3 = topTable(fit2, adjust="fdr", coef='CD4med-Naive', number=Inf)
length(Reduce(union, list(unique(exprs$GeneName),
                          unique(exprs2$GeneName),
                          unique(exprs3$GeneName))))

#Write diffexprs list
write(names(DE[DE[,1]==1,1]), file='CD4tgf-Naive_up.txt', sep='\n')
write(names(DE[DE[,1]==-1,1]), file='CD4tgf-Naive_down.txt', sep='\n')
write(names(DE[DE[,2]==1,2]), file='CD4tgf-CD4med_up.txt', sep='\n')
write(names(DE[DE[,2]==-1,2]), file='CD4tgf-CD4med_down.txt', sep='\n')
write(names(DE[DE[,3]==1,3]), file='CD4med-Naive_up.txt', sep='\n')
write(names(DE[DE[,3]==-1,3]), file='CD4med-Naive_down.txt', sep='\n')


#Generate Complete Results Table
tempdir <- getwd()
setwd('/Users/Baiochi/Dropbox/USP/Lab/RawData/Naive microRNA/RawData/GeneView/')
data <- data.frame(SystematicName=norm$genes$GeneName)
files <- dir()
for (i in 1:9) {
  cat(paste('Reading', files[i],'\n'))
  aux <- read.delim(file=files[i], header=TRUE, skip=1, stringsAsFactors=FALSE)
  aux <- aux[grep('hsa',aux$SystematicName),] #keep only hsa mirs
  aux <- aux[,-c(2,4)] #remove control and tgerror columns
  colnames(aux)[2:3] <- c(paste(targets$Nomeclature[i],'TGS'),paste('gDect',i,sep=''))
  data <- merge(data, aux, by='SystematicName')
}
setwd(tempdir)
norma <- as.data.frame(cbind(norm$TGS, norm$genes$GeneName), stringsAsFactors=FALSE)
colnames(norma)[10] <- 'SystematicName'
tabela <- merge(data, norma, by='SystematicName')
k=20 #concatenate columns raw to norm data
for(i in seq(2,18,by=2)){
  tabela[,i] <- round(as.numeric(tabela[,k]), digits=6)
  k = k + 1
}
tabela <- tabela[,-c(20:28)]
colnames(tabela)[1] <- 'GeneName'
completa <- merge(tabela, exprs[,c(3,7,8,4)], by='GeneName')
colnames(completa)[20:22] <- c('P.Value_iTreg_Naive','Adjusted.P.Value_iTreg_Naive','logFC_iTreg_Naive2')
completa <- merge(completa, exprs2[,c(3,7,8,4)], by='GeneName')
colnames(completa)[23:25] <- c('P.Value_iTreg_Teff','Adjusted.P.Value_iTreg_Teff','logFC_iTreg_Teff')
completa <- merge(completa, exprs3[,c(3,7,8,4)], by='GeneName')
colnames(completa)[26:28] <- c('P.Value_Teff_Naive','Adjusted.P.Value_Teff_Naive','logFC_Teff_Naive')
for(i in seq(3,19,by=2)) colnames(completa)[i] <- 'gDetected'
colnames(completa)[1] <- 'miR_Name'
#head(completa)
write.table(completa, file='Complete_Table_Results.txt', sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

#Only diffexpressed
diffexprs <- subset(completa, Adjusted.P.Value_iTreg_Naive < 0.05 | Adjusted.P.Value_iTreg_Teff < 0.05 | Adjusted.P.Value_Teff_Naive < 0.05)
#diffexprs <- diffexprs[with(diffexprs, order(miR_Name)), ]
#head(diffexprs)
write.table(completa, file='Differential_Expressed_miRs.txt', sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)








# Plots -------------------------------------------------------------------
lograw <- log2(tgs$TGS)
hierclust(lograw, methdis="spearman", methclu="complete", sel=FALSE) 
cor.matrix(lograw, title="Correlation Matrix")
Boxplot(lograw, title="Raw log Boxplot", color="forestgreen")
densityplot(lograw, dim(lograw)[2], title="Raw log Density")
Boxplot(norm$TGS, title="Normalized Boxplot", color="deepskyblue3")
densityplot(norm$TGS, dim(norm$TGS)[2], title="Normalized Density")

volcanoplot(fit2, 'CD4tgf-Naive', title='CD4tgf/atRA-Naive VolcanoPlot')
diffexprs.hist(exprs, coef='CD4tgf-Naive')
volcanoplot(fit2, 'CD4tgf-CD4med',title='CD4tgf/atRA-CD4med VolcanoPlot')
diffexprs.hist(exprs2, coef='CD4tgf-CD4med')
volcanoplot(fit2, 'CD4med-Naive',title='CD4med-Naive VolcanoPlot')
diffexprs.hist(exprs3, coef='CD4med-Naive')

#vennDiagram
colnames(rank)[1:2] <- c('CD4tgf/atRA-Naive','CD4tgf/atRA-CD4med')
rank2 <- rank[,c(1,3,2)]
vennDiagram(rank2, include=c('up', 'down'))

#hierarchical clustering
colnames(norm$TGS)[1:3] <- c('CD4tgf/atRA_1','CD4tgf/atRA_2','CD4tgf/atRA_4')
hierclust(norm$TGS, methdis="spearman", methclu="complete", sel=FALSE) 

#HEATMAP
data <- diffexprs[,seq(2,18,by=2)]
colnames(data) <- targets$Nomeclature
rownames(data) <- diffexprs$miR_Name
data <- as.matrix(data)
rbg = maPalette(low = "springgreen", high = "firebrick2", mid = "grey10", k = 100)
par(oma = c(0, 0, 2, 0))
pdf(file = 'microRNA_Clustering.pdf', width = 8, height = 12)
heatmap.2(data, 
          distfun = function(x) dist(x,method = 'euclidean'),
          hclustfun = function(x) hclust(x,method = 'complete'),
          labCol = colnames(data), labRow = rownames(data), 
          scale = "none", symkey=FALSE, symbreaks=FALSE, 
          keysize = 1.3, col = rbg, trace = 'none', margins = c(6,10),
          key.title = 'Color Key', main='', lhei=c(1.5, 10))
title(main='Differentially Expressed miRs\nHeatmap')
dev.off()
mtext(text='Differentially Expressed miRs Heatmap', outer=TRUE)


# Save objet
save(tgs, file = 'tgs.rda')
save(norm, file = 'norm.rda')
save(data, file = 'diffexprs.rda')


# References
citation("AgiMicroRna")


