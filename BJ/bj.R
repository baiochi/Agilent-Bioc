#Agilent - pre-miR101 Microarray analysis
#Using Raw MeanSignal
#load functions
library(limma)
source("/Users/Baiochi/Desktop/Agilent-Bioc/Scripts/Functions.R")

#read targets
targets <- read.table(file="/Users/Baiochi/Dropbox/USP/Lab/RawData/BJ/Targets.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)

#read raw
raw = read.maimages(targets$FileName, source="agilent", green.only=TRUE, path="/Users/Baiochi/Dropbox/USP/Lab/RawData/BJ")
colnames(raw$E) <- targets$Nomeclature
rownames(raw$targets) <- targets$Nomeclature
#correct artifacts
raw$Eb[which(raw$Eb[,6]>68),6] <- round(mean(raw$Eb[,6]))
#setup printer
raw$genes$Block <- 1
names(raw$genes)[2] <- "Column"
raw$printer <- getLayout(raw$genes)
r <- raw$genes$Row
c <- raw$genes$Col
nr <- max(r)
nc <- max(c)
y <- rep(NA,nr*nc)
j <- (r-1)*nc+c

# Array Visualization
par(mfrow=c(2,3),oma = c(0, 0, 2, 0))
for(i in 1:6){
  y[j] <- log2(raw$Eb[,i])
  imageplot(y, raw$printer)
}
mtext('Array Images', outer = TRUE)
par(mfrow=c(2,3),oma = c(0, 0, 2, 0))
for(i in 1:6)
  plotFB(raw, array=i)
mtext('Foreground vs Background', outer = TRUE)

#bgc correct
bgc = backgroundCorrect(raw,method='normexp')


# ProcessedSignal
ps <- readAFE(targets = targets, path = "/Users/Baiochi/Dropbox/USP/Lab/RawData/BJ/", skip.lines =  9)
colnames(ps$E) <- targets$Nomeclature
rownames(ps$targets) <- targets$Nomeclature

#normalize
norm = normalizeBetweenArrays(ps,method='quantile')

#filter control probes
eset = norm[norm$genes$ControlType==0,]

#avg probes
eset = avereps(eset,ID=eset$genes[,"SystematicName"])

#design
f <- factor(targets$Condition, levels = unique(targets$Condition))
design = model.matrix(~0 + f)
colnames(design) <- levels(f)
#contrasts
contrast.matrix = makeContrasts(contrasts='PM29a-Ctr', levels=design)
#fit model
fit <- lmFit(eset$E, design)
fit2 <- contrasts.fit(fit, contrast.matrix)

#bayes statistics
fit2 = eBayes(fit2)

#rank genes
exprs = topTable(fit2, adjust="fdr", coef='PM29a-Ctr', genelist=eset$genes, number=Inf)
rank = decideTests(fit2, method="separate", adjust.method="fdr", p.value=0.05, lfc=0)
summary(rank)
as.numeric(table(exprs$P.Value<0.05 & exprs$logFC>0)[2])
as.numeric(table(exprs$P.Value<0.05 & exprs$logFC<0)[2])


#-------------------Plots-------------------#
raw <- ps

boxplot(log2(raw$E), main='Raw data Boxplot', col='forestgreen')
hierclust(log2(raw$E), methdis="spearman", methclu="complete", sel=FALSE, size=100)
cor.matrix(log2(raw$E), title="Pearson Correlation Matrix")
plotDensities(raw, legend = 'topright', main='Raw Foreground densities')
volcanoplot(fit2, title='PM29a vs Control VolcanoPlot')
densityplot(exprs[,14:15],2,title='P-values density')

# VolcanoPlot
Pval = -log(exprs$P.Value)
FC = exprs$logFC
pval <- -log(0.05)
fc <- 0
plot(x = FC, y = Pval,
     ylab = "-log Adjusted P.value", xlab = "log2 Fold Change", 
     main = 'PM29a vs Control VolcanoPlot', pch = 20, col = "black", 
     cex=.3, #xlim=(c(-6, 4)), ylim=(c(0, 15)), 
     cex.axis=0.9, cex.lab=0.9)
points(FC[(Pval>pval & FC>fc)],
       Pval[(Pval>pval & FC>fc)],
       col = "springgreen1",pch = 20, cex=.3)
points(FC[(Pval>pval & FC< -fc)],
       Pval[(Pval>pval & FC< -fc)],
       col = "firebrick1",pch = 20, cex=.3)
legend("bottomright",
       legend = c(paste('Up-regulated:',as.numeric(table(exprs$P.Value<0.05 & FC>0)[2]), sep=' '),
                  paste('Down-regulated:',as.numeric(table(exprs$P.Value<0.05 & FC<0)[2]), sep=' ') ),
       lty= 0,# lwd=c(2.5,2.5),
       pch = 20, cex=0.75, col=c('springgreen1','firebrick1'))


#-------------------Results-------------------#

# Systematic, NormSignal, P.value, FC, adj.Pval, Gene, Probe, Sequence
results <- exprs[,c(9, 14,15,11,12,13,16)]
norm <- norm[norm$genes$ControlType==0,]
normalized <- as.data.frame(norm$E)
normalized$SystematicName <- norm$genes$SystematicName
normalized$GeneName <- norm$genes$GeneName
normalized$ProbeName <- norm$genes$ProbeName
normalized$Sequence <- norm$genes$Sequence
results <- merge(normalized, results)

result.table <- results[,c(7:16,1:6)]

only.pval <- result.table[result.table$P.Value<0.05,]

write.table(only.pval, file='BJ-PM29a_Statistical_Results.txt', 
            sep="\t", quote=FALSE, row.names=FALSE)
write.table(result.table, file='BJ-PM29a_Complete_Statistical_Results.txt', 
            sep="\t", quote=FALSE, row.names=FALSE)


#annotation
WriteResults(fit = fit2, coefs = 'PM29a-Ctr', pval = 1,eset = eset, filename = 'miR-29a')
#Save R Data
saveData(path='/Users/Baiochi/Dropbox/USP/Lab/Results/pre-miR101/RData/',raw,bgc,norm,eset,design,fit2,exprs,rank)
#Load Data
loadData(path='/Users/Baiochi/Dropbox/USP/Lab/Results/pre-miR101/RData')

