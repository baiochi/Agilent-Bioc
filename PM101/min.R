# Minimal - PM101 Microarray analysis
#load functions
library(limma)
source("/Users/Baiochi/Dropbox/USP/Lab/Agilent Bioc/Functions.R")
#read targets
targets <- read.table(file="/Users/Baiochi/Dropbox/USP/Lab/RawData/HCT/Targets.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
#read raw
raw = read.maimages(targets$FileName, source="agilent", green.only=TRUE, path="/Users/Baiochi/Dropbox/USP/Lab/RawData/HCT")
colnames(raw$E) <- targets$Nomeclature
rownames(raw$targets) <- targets$Nomeclature
#bgc correct
bgc = backgroundCorrect(raw,method='normexp')
#normalize
norm = normalizeBetweenArrays(bgc,method='quantile')
#filter
neg95 <- apply(norm$E[norm$genes$ControlType==-1,],2,function(x) quantile(x,p=0.95))
cutoff <- matrix(1.1*neg95,nrow(norm),ncol(norm),byrow=TRUE)
isexpr <- rowSums(norm$E > cutoff) >= length(unique(targets$Sample))
eset = norm[norm$genes$ControlType==0 & isexpr,]
#avg probes
eset = avereps(eset,ID=eset$genes[,"SystematicName"])
#design
f <- factor(targets$Condition, levels = unique(targets$Condition))
design = model.matrix(~0 + f)
colnames(design) <- levels(f)
#fit model
fit <- lmFit(eset$E, design)
#contrasts
contrast.matrix = makeContrasts(contrasts='PM101-Ctr', levels=design)
#bayes statistics
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 = eBayes(fit2)
#rank genes
rank = decideTests(fit2, method="separate", adjust.method="fdr", p.value=0.05, lfc=0)
summary(rank)
exprs = topTable(fit2, adjust="fdr", coef='PM101-Ctr', genelist=eset$genes, number=Inf)
#plots
RawPlots(raw)
NormPlots(eset)
DiffExprsPlots(fit2, 'PM101-Ctr', eset)
#annotation
WriteResults(fit2, 'PM101-Ctr', eset)
#Save R Data
saveData(path='',raw,bgc,norm,eset,design,fit2,exprs)
#Load Data
loadData(path='')