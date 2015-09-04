# Minimal - Naive Cells Microarray analysis
#load functions
library(limma)
source("/Users/Baiochi/Desktop/Agilent-Bioc/Scripts/Functions.R")

#read targets
targets <- read.table(file="/Users/Baiochi/Dropbox/USP/Lab/RawData/Naive mRNA//Targets.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)

#read raw
raw = read.maimages(targets$FileName, source="agilent", green.only=TRUE, path="/Users/Baiochi/Dropbox/USP/Lab/RawData/Naive mRNA/")
colnames(raw$E) <- targets$Nomeclature
colnames(raw$Eb) <- targets$Nomeclature
rownames(raw$targets) <- targets$Nomeclature
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

#read ProcessedSignal from AFE
tempdir <- getwd()
setwd("/Users/Baiochi/Dropbox/USP/Lab/RawData/Naive mRNA/")
data <- new('list')
for (i in 1:9) {
  cat(paste('Reading', targets$FileName[i],'\n'))
  data[[i]] <- read.delim(file=targets$FileName[i], header=TRUE, stringsAsFactors=FALSE)
}
signal <- data.frame(numeric(45015))
for (i in 1:9) signal <- cbind(signal, data[[i]]$gProcessedSignal)
aux <- data[[1]]
flag <- aux[,c(37,38,39,40)]
rm(aux)
signal <- signal[,-1]
signal <- as.matrix(signal)
colnames(signal) <- targets$Nomeclature
raw$E <- signal
raw$flag <- flag
setwd(tempdir)


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





#normal analytsis
#iTreg-Naive
#-1        2813
#0        26773
#1         2891
#iTreg-Teff
#-1        496
#0       30725
#1        1256
#Teff-Naive
#-1       1816
#0       28610
#1        2051




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