# Agilent - pre-miR101 Microarray analysis
#Using ProcessedSignal from AFE
#load functions
library(limma)
source("/Users/Baiochi/Desktop/Agilent-Bioc/Scripts/Functions.R")

tempdir <- getwd()
setwd("/Users/Baiochi/Dropbox/USP/Lab/RawData/HCT")
targets <- read.table(file="/Users/Baiochi/Dropbox/USP/Lab/RawData/HCT/Targets.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)

#Read data
data <- new('list')
for (i in 1:6) {
  cat(paste('Reading', targets$FileName[i],'\n'))
  data[[i]] <- read.delim(file=targets$FileName[i], header=TRUE, skip=9, stringsAsFactors=FALSE)
}

#Create Expression set
signal <- data.frame(numeric(45015))
mean <- data.frame(numeric(45015))
bg <- data.frame(numeric(45015))
for (i in 1:6){
  signal <- cbind(signal, data[[i]]$gProcessedSignal)
  mean <- cbind(mean, data[[i]]$gMeanSignal)
  bg <- cbind(bg, data[[i]]$gBGMeanSignal)
}
aux <- data[[1]]
genes <- aux[,c(9,10,11,12,13,14,15,16,37,38,39,40)]
rm(aux)
signal <- signal[,-1]
mean <- mean[,-1]
bg <- bg[,-1]
colnames(signal) <- targets$Nomeclature
colnames(mean) <- targets$Nomeclature
colnames(bg) <- targets$Nomeclature

rawdata <- new('EListRaw')
rawdata$E <- signal
rawdata$Eb <- bg
rawdata$Raw <- mean
rawdata$targets <- as.data.frame(targets$FileName, row.names=targets$Nomeclature)
rawdata$genes <- genes
rawdata$source <- 'agilent'

setwd(tempdir)

#4 last columns
outlierProbes <- as.data.frame(which(rawdata$genes[,9:12]==1, arr.ind = TRUE))
outlierProbes <- unique(outlierProbes[,1])
outlierGenes <- new('character')
for(i in 9:12){
  aux <- rawdata$genes$GeneName[which(rawdata$genes[,i]==1)]
  outlierGenes <- c(outlierGenes, aux)
}



#--------Start Analysis----------#

#normalize
norm = normalizeBetweenArrays(rawdata,method='quantile')
#filter control probes
eset = rawdata[rawdata$genes$ControlType==0,]
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


