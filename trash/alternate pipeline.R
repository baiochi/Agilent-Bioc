# limma Studie Case Analysis - HCT
# Modify:
# Line 30: number of replicates
# Line 53: coef for statistical test

library(limma)

targets <- read.table(file="Targets.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)

# Reading raw data
x <- read.maimages(files = targets$FileName,
					path = "RawData",
					source = "agilent", 
					green.only = TRUE, 
					verbose = TRUE)

# Backround Correction
y <- backgroundCorrect(x,method="normexp")

# Normalization
y <- normalizeBetweenArrays(y,method="quantile")

# Filtering probes
# 95% percentile of the neg-crtl probes on each array
# Keep probes that are at least 10% brigter than negative control at 
# least N arrays(N = number of replicates)
neg95 <- apply(y$E[y$genes$ControlType==-1,],2,function(x) quantile(x,p=0.95))
cutoff <- matrix(1.1*neg95,nrow(y),ncol(y),byrow=TRUE)
isexpr <- rowSums(y$E > cutoff) >= 3
table(isexpr)

# Remove control probes
y0 <- y[y$genes$ControlType==0 & isexpr,]

# Differential Expression
Treatment <- targets$Condition
levels <- unique(targets$Condition)
Treatment <- factor(Treatment,levels=levels)
design <- model.matrix(~Treatment)
fit <- lmFit(y0$E,design)
fit <- eBayes(fit,trend=TRUE)
#plotSA(fit, main="Probe-level")
summary(decideTests(fit[,-1]))

# DiffExprs with averaging probes for each gene
yave <- avereps(y0,ID=y0$genes[,"SystematicName"])
fit2 <- lmFit(yave,design)
fit2 <- eBayes(fit2,trend=TRUE)
#plotSA(fit, main="Gene-level")
summary(decideTests(fit2[,-1]))

# Extracing test results
st.tests <- topTable(fit, adjust="fdr", coef="TreatmentPM101", genelist=y0$genes, number=Inf)
st.tests2 <- topTable(fit2, adjust="fdr", coef="TreatmentPM101", genelist=y0$genes, number=Inf)

# Evaluating Diff Exprs
diffexprs <- decideTests(fit, method="separate", adjust.method="fdr", p.value=0.05, lfc=0)
diffexprs2 <- decideTests(fit2, method="separate", adjust.method="fdr", p.value=0.05, lfc=0)
table(diffexprs)
table(diffexprs2)

# Diff Exprs Plot
volcano(st.tests$logFC, st.tests$adj.P.Val, title="Volcano Plot", y="-log(Adj.P.Val)")

# Down fit:  5509
table(st.tests$logFC<0 & st.tests$adj.P.Val<0.05)
# Up fit: 6106
table(st.tests$logFC>0 & st.tests$adj.P.Val<0.05)
# Down fit2: 4318
table(st.tests2$logFC<0 & st.tests2$adj.P.Val<0.05)
# Up fit2: 4248
table(st.tests2$logFC>0 & st.tests2$adj.P.Val<0.05)

table(st.tests2$logFC<0 & st.tests2$adj.P.Val<0.05)
down <- which(st.tests$logFC<0 & st.tests$adj.P.Val<0.05)
down <- st.tests2[down,]
down <- down[,c(5,7,8,11,14,15)]
rownames(down) <- NULL
alt.down <- down$GeneName

table(st.tests2$logFC>0 & st.tests2$adj.P.Val<0.05)
up <- which(st.tests2$logFC<0 & st.tests2$adj.P.Val<0.05)
up <- st.tests2[up,]
up <- up[,c(5,7,8,11,14,15)]
rownames(up) <- NULL
alt.up <- up$GeneName

# diffexprs - fit
#    -1     0     1 
#  5509  9301 27022 

#  diffexprs2 - fit2
#    -1     0     1 
#  4318  6602 19416 

