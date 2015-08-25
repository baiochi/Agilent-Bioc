#Naive-iTreg mRNA Analysis

eset <- read.delim('/Users/Baiochi/Dropbox/USP/Lab/Analises/Naive.iTreg.Teff.mRNA/RawData/Naive_1.txt',
					header=TRUE, sep='\t', stringsAsFactors=FALSE)



# Load libraries/functions
library(limma)
source("/Users/Baiochi/Dropbox/USP/Lab/Analises/Functions.R")

# Read targets
setwd('/Users/Baiochi/Dropbox/USP/Lab/Analises/Naive.iTreg.Teff.mRNA')
targets <- read.table(file="Targets.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
targets <- cbind(targets,Nomeclature=paste(targets$Condition,targets$Sample,sep="_"))

# Read Images
raw <- read.maimages(files = targets$FileName,
					path = "RawData",
					source = "agilent", 
					green.only = TRUE, 
					verbose = TRUE)
colnames(raw$E) <- targets$Nomeclature
colnames(raw$Eb) <- targets$Nomeclature
rownames(raw$targets) <- targets$Nomeclature

# Log transform for plot visualization
rawlog <- raw
rawlog$E <- log2(raw$E)

# Correct Background
bgc <- backgroundCorrect(raw, method="subtract", offset=16)

# Normalization
#norm <- normalizeBetweenArrays(bgc, method="cyclicloess")
norm <- normalizeBetweenArrays(raw, method='quantile')

# Remove control probes
filter <- norm[norm$genes$ControlType==0,]
# Average Probes
exprs <- avereps(filter, ID=filter$genes$ProbeName)

# DiffExprs with ANOVA
f <- factor(targets$Condition, levels = unique(targets$Condition))
design <- model.matrix(~0 + f)
colnames(design) <- levels(f)
fit <- lmFit(exprs$E, design)
contrast.matrix <- makeContrasts(Teff-Naive, iTreg-Teff, iTreg-Naive, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# TopTable
coef1 <- topTable(fit2, adjust="fdr", coef=1, genelist=exprs$genes, number=Inf)
coef2 <- topTable(fit2, adjust="fdr", coef=2, genelist=exprs$genes, number=Inf)
coef3 <- topTable(fit2, adjust="fdr", coef=3, genelist=exprs$genes, number=Inf)

# DecideTests
results <- decideTests(fit2, method="separate", adjust.method="fdr", p.value=0.05, lfc=0)
summary(results)
vennDiagram(results, include='up')
vennDiagram(results, include='down')
vennDiagram(results, include=c('up', 'down'))

par(mfrow=c(1,3))
volcanoplot(coef1$logFC, coef1$adj.P.Val, results[,1], 
			title="Volcano Plot\nTeff-Naive", y="-log(Adj.P.Val)")
volcanoplot(coef2$logFC, coef2$adj.P.Val, results[,2], 
			title="Volcano Plot\niTreg-Teff", y="-log(Adj.P.Val)")
volcanoplot(coef3$logFC, coef3$adj.P.Val, results[,3], 
			title="Volcano Plot\niTreg-Naive", y="-log(Adj.P.Val)")


# Save R Objects
tempdir <- getwd()
setwd('/Users/Baiochi/Dropbox/USP/Lab/Analises/Naive.iTreg.Teff.mRNA/Scripts/RData/')
save(raw, file='raw.rda')
save(exprs, file='exprs.rda')
save(fit2, file='fit2.rda')
save(diffexprs, file='diffexprs.rda')
setwd(tempdir)

############# AFTER PIPELINE ##############

# Load R Objects
tempdir <- getwd()
setwd('/Users/Baiochi/Dropbox/USP/Lab/Analises/Naive.iTreg.Teff.mRNA/Scripts/RData/')
load(file='raw.rda')
load(file='exprs.rda')
load(file='fit2.rda')
load(file='diffexprs.rda')
setwd(tempdir)

# Array Q&A
library(arrayQualityMetrics)
arrayQualityMetrics(raw, outdir = "./qa-raw/")
arrayQualityMetrics(exprs, outdir = "./qa-norm/")

# Save plots
raw.plot <- raw$E
exprs.plot <- exprs
diffexprs.plot <- diffexprs
saveAllPlots(raw = raw.plot, eset = exprs.plot, diffexprs = diffexprs.plot, testresult)



