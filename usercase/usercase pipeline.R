#limma UserGuide Agilent Example

# Read Targets File
SDRF <- read.delim("E-GEOD-33005.sdrf.txt",check.names=FALSE,stringsAsFactors=FALSE)

# Read raw data
x <- read.maimages(SDRF[,"Array Data File"],source="agilent", path='RawData', green.only=TRUE)

# Background subtraction
y <- backgroundCorrect(x,method="normexp")

# Normalization
y <- normalizeBetweenArrays(y,method="quantile")

# Filtering
neg95 <- apply(y$E[y$genes$ControlType==-1,],2,function(x) quantile(x,p=0.95))
cutoff <- matrix(1.1*neg95,nrow(y),ncol(y),byrow=TRUE)
isexpr <- rowSums(y$E > cutoff) >= 4
table(isexpr)
y0 <- y[y$genes$ControlType==0 & isexpr,]

# Creating design
Treatment <- SDRF[,"Characteristics[treatment]"]
levels <- c("10 ml/kg saline","2 ml/kg corn oil","5 ml/kg corn oil","10 ml/kg corn oil")
Treatment <- factor(Treatment,levels=levels)
design <- model.matrix(~Treatment)

# Probe Level
fit <- lmFit(y0,design)
fit <- eBayes(fit,trend=TRUE)
plotSA(fit, main="Probe-level")
summary(decideTests(fit[,-1]))

# Gene Level
yave <- avereps(y0,ID=y0$genes[,"SystematicName"])
fit2 <- lmFit(yave,design)
fit2 <- eBayes(fit2,trend=TRUE)
plotSA(fit2, main="Gene-level")
summary(decideTests(fit2[,-1]))

f <- factor(SDRF[,"Characteristics[treatment]"], levels = unique(SDRF[,"Characteristics[treatment]"]))
design <- model.matrix(~0 + f)
colnames(design) <- levels(f)
fit <- lmFit(y0$E, design)
coef <- colnames(design)
contrast.matrix <- makeContrasts(contrasts=coef, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# TopTable
diffexprs <- topTable(fit[,-1], adjust="fdr", coef="Treatment10 ml/kg corn oil", genelist=y0$genes, number=Inf)
rownames(diffexprs) <- NULL
diffexprs2 <- topTable(fit2[,-1], adjust="fdr", coef="Treatment10 ml/kg corn oil", genelist=y0$genes, number=Inf)
rownames(diffexprs2) <- NULL
# DecideTests
testresult <- decideTests(fit[,-1], method="separate", adjust.method="fdr", p.value=0.05, lfc=0)
summary(testresult)
testresult2 <- decideTests(fit2[,-1], method="separate", adjust.method="fdr", p.value=0.05, lfc=0)
summary(testresult2)


#Plots
cor.matrix(log2(x$E), title="Correlation Matrix", save=TRUE)
Boxplot(log2(x$E), title="Raw log Boxplot", color="forestgreen", save=TRUE)
densityplot(log2(x$E), 19, title="Raw log Density", save=TRUE)

Boxplot(y$E, title="Normalized Boxplot", color="deepskyblue1", save=TRUE)
Boxplot(y0$E, title="Filtered Boxplot", color="deepskyblue3", save=TRUE)
densityplot(y0$E, 19, title="Normalized Density", save=TRUE)

#Probe Level
DiffExprsPlots(fit=fit,coef="Treatment10 ml/kg corn oil",eset=y0, title='Probe Level', save=TRUE, path='Probe Level')

#Gene Level
DiffExprsPlots(fit=fit2,coef="Treatment10 ml/kg corn oil",eset=y0, title='Gene Level', save=TRUE, path='Gene Level')




