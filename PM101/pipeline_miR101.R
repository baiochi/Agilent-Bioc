#Agilent - pre-miR101 Microarray analysis
#Using Raw MeanSignal
#load functions
library(limma)
source("/Users/Baiochi/Desktop/Agilent-Bioc/Scripts/Functions.R")

#read targets
targets <- read.table(file="/Users/Baiochi/Dropbox/USP/Lab/RawData/HCT/Targets.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)

#read raw
raw = read.maimages(targets$FileName, source="agilent", green.only=TRUE, path="/Users/Baiochi/Dropbox/USP/Lab/RawData/HCT")
colnames(raw$E) <- targets$Nomeclature
colnames(raw$Eb) <- targets$Nomeclature
rownames(raw$targets) <- targets$Nomeclature
#correct artifacts
raw$Eb[which(raw$Eb[,3]>60),3] <- round(mean(raw$Eb[,3]))
raw$Eb[which(raw$Eb[,6]>60),6] <- round(mean(raw$Eb[,6]))
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

#bgc correct
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
rank = decideTests(fit2, method="separate", adjust.method="fdr", p.value=0.05, lfc=0)
summary(rank)






#-------------------Plots-------------------#
#MvsA Plots
par(mfrow=c(2,3))
for(i in 1:6)
  plotMA(raw, array=i, main=paste('Array', i, '- Raw'))
RawPlots(raw)
NormPlots(eset, sn=2, rep=3)
DiffExprsPlots(fit2, 'PM101-Ctr', rank, eset)

#-------------------Results-------------------#
#annotation
WriteResults(fit2, 'PM101-Ctr', eset, 'miR-101')
#Save R Data
saveData(path='/Users/Baiochi/Dropbox/USP/Lab/Results/pre-miR101/RData/',raw,bgc,norm,eset,design,fit2,exprs,rank)
#Load Data
loadData(path='/Users/Baiochi/Dropbox/USP/Lab/Results/pre-miR101/RData')

