# Minimal - Naive Cells Microarray analysis
#load functions
library(limma)
source("/Users/Baiochi/Desktop/Agilent-Bioc/Scripts/Functions.R")

#read targets
targets <- read.table(file="/Users/Baiochi/Dropbox/USP/Lab/RawData/Naive mRNA//Targets.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)

#Read using ProcessedSignal
raw <- readAFE(targets = targets, path = "/Users/Baiochi/Dropbox/USP/Lab/RawData/Naive mRNA/", skip.lines =  0)

#Read using MedianSignal
#raw = read.maimages(targets$FileName, source="agilent", green.only=TRUE, path="/Users/Baiochi/Dropbox/USP/Lab/RawData/Naive mRNA/")
#colnames(raw$E) <- targets$Nomeclature
#colnames(raw$Eb) <- targets$Nomeclature
#rownames(raw$targets) <- targets$Nomeclature
#bgc correct bgc = backgroundCorrect(raw,method='normexp')

#normalize
norm = normalizeBetweenArrays(raw,method='quantile')
#filter control probes
eset = norm[norm$genes$ControlType==0,]
#avg probes
eset = avereps(eset,ID=eset$genes[,"SystematicName"])

#design
f <- factor(targets$Condition, levels = unique(targets$Condition))
design = model.matrix(~0 + f)
colnames(design) <- levels(f)
#contrasts
contrast.matrix = makeContrasts('iTreg-Naive','iTreg-Teff','Naive-Teff',levels=design)
#fit model
fit <- lmFit(eset$E, design)
fit2 <- contrasts.fit(fit, contrast.matrix)

#bayes statistics
fit2 = eBayes(fit2)

# iTreg vs Naive
exprs = topTable(fit2, adjust="fdr", coef='iTreg-Naive', genelist=eset$genes, number=Inf)
# iTreg vs Teff
exprs2 = topTable(fit2, adjust="fdr", coef='iTreg-Teff', genelist=eset$genes, number=Inf)
# Naive vs Teff
exprs3 = topTable(fit2, adjust="fdr", coef='Naive-Teff', genelist=eset$genes, number=Inf)

#rank genes
rank = decideTests(fit2, method="separate", adjust.method="fdr", p.value=0.05, lfc=0)
summary(rank)


volcanoplot(fit2, 'iTreg-Naive',title='iTreg-Naive vs Control VolcanoPlot')
volcanoplot(fit2, 'iTreg-Teff',title='iTreg-Teff vs Control VolcanoPlot')
volcanoplot(fit2, 'Naive-Teff',title='Naive-Teff vs Control VolcanoPlot')

#ProcessedSignal analysis
#iTreg-Naive    iTreg-Teff Naive-Teff
#-1        3229        630       2314
#0        26006      30516      28273
#1         3242       1331       1890

#MedianSignal analysis
#iTreg-Naive      #iTreg-Teff       #Teff-Naive
#-1        2813   #-1        496    #-1       1816
#0        26773   #0       30725    #0       28610
#1         2891   #1        1256    #1        2051


# Plots -------------------------------------------------------------------

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
par(mfrow=c(3,3))
for(i in 1:9)
  plotMA(raw, array=i, main=paste('Array', i, '- Raw'))

#post-norm visualization
boxplot(eset$E, main='Normalized Boxplot', col='deepskyblue2')
plotDensities(eset, legend = 'topright', main='Normalized Foreground Densities')
MAPlot(eset, 3, 3, title = "Post-Normalization MA Plot")

#differential expression visualization
plotSA(fit2)
plotMA(fit2)
volcanoplot(fit2, 'iTreg-Naive',title='iTreg-Naive vs Control VolcanoPlot')
diffexprs.hist(exprs, coef='iTreg-Naive')
volcanoplot(fit2, 'iTreg-Teff',title='iTreg-Teff vs Control VolcanoPlot')
diffexprs.hist(exprs2, coef='iTreg-Teff')
volcanoplot(fit2, 'Naive-Teff',title='Naive-Teff vs Control VolcanoPlot')
diffexprs.hist(exprs3, coef='Naive-Teff')









#-------------------Plots-------------------#
#MvsA Plots
par(mfrow=c(2,3))
for(i in 1:6)
  plotMA(raw, array=i, main=paste('Array', i, '- Raw'))
RawPlots(raw=raw)
NormPlots(eset=eset, sn=3, rep=3)
DiffExprsPlots(fit=fit2, coefs=c('iTreg-Naive', 'iTreg-Teff', 'Naive-Teff'), rank=rank, eset=eset)

#-------------------Results-------------------#
#annotation
WriteResults(fit2, 'R', eset, filename=)
#Save R Data
saveData(path='',raw,bgc,norm,eset,design,fit2,exprs)
#Load Data
loadData(path='')