# Array reproducibility
cvArray(dd.micro,"MeanSignal", targets.micro, verbose=TRUE)
cvArray(dd.micro,"ProcessedSignal", targets.micro, verbose=TRUE)
# Total Gene Signal
ddTGS <- ddTGS[grep('hsa',ddTGS$genes$GeneName),]
ddTGS=tgsMicroRna(dd.micro, half=FALSE, makePLOT=FALSE, verbose=FALSE)
# Normalization between arrays
ddNORM=tgsNormalization(ddTGS, "quantile", makePLOTpre=FALSE, makePLOTpost=FALSE, targets.micro, verbose=TRUE)
# Filtering probes 
ddPROC=filterMicroRna(ddNORM,
                      dd.micro,
                      control=TRUE,
                      IsGeneDetected=FALSE,
                      wellaboveNEG=FALSE,
                      limIsGeneDetected=0,
                      limNEG=0,
                      makePLOT=FALSE,
                      targets.micro,
                      verbose=TRUE,
                      writeout=FALSE)
# Creating Expression Set
esetPROC=esetMicroRna(ddNORM, targets.micro, makePLOT=FALSE, verbose=FALSE)
#writeEset(esetPROC, ddPROC, targets.micro, verbose=TRUE)
# Differential Expression
f <- factor(targets$Condition[-c(3,7,11)], levels = unique(targets$Condition))
design <- model.matrix(~0 + f)
colnames(design) <- levels(f)
fit <- lmFit(esetPROC, design)
contrast.matrix <- makeContrasts('iTreg-Naive', 'iTreg-Teff', 'Teff-Naive', levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
#Selecting significant miRNAs
rank = decideTests(fit2, method="separate", adjust.method="fdr", p.value=0.05, lfc=0)
summary(rank)
DE <- getDecideTests(fit2, DEmethod='separate', MTestmethod='BH', PVcut=0.05,verbose=FALSE)
# iTreg vs Naive
exprs = topTable(fit2, adjust="fdr", coef='iTreg-Naive', genelist=esetPROC$genes, number=Inf)
# iTreg vs Teff
exprs2 = topTable(fit2, adjust="fdr", coef='iTreg-Teff', genelist=eset$genes, number=Inf)
# Naive vs Teff
exprs3 = topTable(fit2, adjust="fdr", coef='Teff-Naive', genelist=eset$genes, number=Inf)
rank = decideTests(fit2, method="separate", adjust.method="fdr", p.value=0.05, lfc=0)
summary(rank)
volcanoplot(fit2, 'iTreg-Naive',title='iTreg-Naive vs Control VolcanoPlot')
diffexprs.hist(exprs, coef='iTreg-Naive')
volcanoplot(fit2, 'iTreg-Teff',title='iTreg-Teff vs Control VolcanoPlot')
diffexprs.hist(exprs2, coef='iTreg-Teff')
volcanoplot(fit2, 'Naive-Teff',title='Naive-Teff vs Control VolcanoPlot')
diffexprs.hist(exprs3, coef='Naive-Teff')
# iTreg-Naive (up/down)
write(names(DE[DE[,1]==1,1]), file='iTreg-Naive_up.txt', sep='\n')
write(names(DE[DE[,1]==-1,1]), file='iTreg-Naive_down.txt', sep='\n')
# iTreg-Teff 
write(names(DE[DE[,2]==1,2]), file='iTreg-Teff_up.txt', sep='\n')
write(names(DE[DE[,2]==-1,2]), file='iTreg-Teff_down.txt', sep='\n')
# Naive-Teff
write(names(DE[DE[,3]==1,3]), file='Naive-Teff_up.txt', sep='\n')
write(names(DE[DE[,3]==-1,3]), file='Naive-Teff_down.txt', sep='\n')


#------------limma analysis--------#
#Using Median Signal
#     iTreg-Naive iTreg-Teff Naive-Teff
# -1          80         60         62
# 0          752        814        805
# 1          107         65         72
# Using Processed Signal
#     iTreg-Naive iTreg-Teff Naive-Teff
# -1          59         12         27
# 0          781        894        887
# 1           80         14          6
# Using TotalGene Signal (warning: More than half of residual variances are exactly zero: eBayes unreliable)
#     iTreg-Naive iTreg-Teff Naive-Teff
# -1           8          2         12
# 0          873        918        926
# 1           58         19          1


targets <- read.table(file="/Users/Baiochi/Dropbox/USP/Lab/RawData/Naive microRNA/RawData/Targets.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
targets <- cbind(targets,Nomeclature=paste(targets$Condition,targets$Sample,sep="_"))
raw = read.maimages(targets$FileName, source="agilent", green.only=TRUE, path="/Users/Baiochi/Dropbox/USP/Lab/RawData/Naive microRNA/RawData/")
colnames(raw$E) <- targets$Nomeclature
rownames(raw$targets) <- targets$Nomeclature

#raw <- readAFE(targets = targets, value='TGS',path = "/Users/Baiochi/Dropbox/USP/Lab/RawData/Naive microRNA/RawData/", skip.lines =  9)


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
contrast.matrix = makeContrasts('iTreg-Naive', 'iTreg-Teff', 'Naive-Teff', levels=design)
#fit model
fit <- lmFit(eset$E, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
#bayes statistics
fit2 = eBayes(fit2)
#rank genes
# iTreg vs Naive
exprs = topTable(fit2, adjust="fdr", coef='iTreg-Naive', genelist=eset$genes, number=Inf)
# iTreg vs Teff
exprs2 = topTable(fit2, adjust="fdr", coef='iTreg-Teff', genelist=eset$genes, number=Inf)
# Naive vs Teff
exprs3 = topTable(fit2, adjust="fdr", coef='Naive-Teff', genelist=eset$genes, number=Inf)
rank = decideTests(fit2, method="separate", adjust.method="fdr", p.value=0.05, lfc=0)
summary(rank)
volcanoplot(fit2, 'iTreg-Naive',title='iTreg-Naive vs Control VolcanoPlot')
diffexprs.hist(exprs, coef='iTreg-Naive')
volcanoplot(fit2, 'iTreg-Teff',title='iTreg-Teff vs Control VolcanoPlot')
diffexprs.hist(exprs2, coef='iTreg-Teff')
volcanoplot(fit2, 'Naive-Teff',title='Naive-Teff vs Control VolcanoPlot')
diffexprs.hist(exprs3, coef='Naive-Teff')

RawPlots(raw)
NormPlots(eset, sn=3, rep=4)
DiffExprsPlots(fit2, 'iTreg-Naive', rank, eset)
DiffExprsPlots(fit2, 'iTreg-Teff', rank, eset)
DiffExprsPlots(fit2, 'Naive-Teff', rank, eset)



cor.matrix <- function(dat, title="Correlation Matrix", save=FALSE)
  
  
  
  
  
  
  
  #targets<- data.frame(FileName=list.files(),
  #                     Condition=c(rep('iTreg',4), rep('Naive',4), rep('Teff',4)),
  #                     Sample=rep(c(1,2,3,4),3)
  #                     )
#targets <- cbind(targets,Nomeclature=paste(targets$Condition,targets$Sample,sep="_"))
#write.table(targets, file='Targets.txt', row.names = FALSE, sep='\t', quote=FALSE)



























########################################################################### PLOT

#Modified function to save all plots (JPEG files)
dd <- dd.micro
MMMfunction <- function(MMM, offset=5){
  min = min(MMM)
  for (i in 1:dim(MMM)[2]) {
    MMM[, i] = MMM[, i] + (abs(min) + offset)
  }
  MMM = log2(MMM)
  return(MMM)
}

# MeanSignal
MMM <- MMMfunction(dd$meanS)
maintitle = "MeanSignal"
colorfill = "orange"
boxplotMicroRna(MMM, maintitle, colorfill)
plotDensityMicroRna(MMM, maintitle)

ddaux = dd.micro
ddaux$meanS = MMM
mvaMicroRna(ddaux, maintitle, verbose = FALSE)
rm(ddaux)

# ProcessedSignal
MMM <- MMMfunction(dd$procS)
maintitle = "ProcessedSignal"
colorfill = "blue"
boxplotMicroRna(MMM, maintitle, colorfill)
plotDensityMicroRna(MMM, maintitle)

ddaux = dd
ddaux$TPS = MMM
mvaMicroRna(ddaux, maintitle, verbose = FALSE)
rm(ddaux)
maintitle = "ProcessedSignal - RLE "
RleMicroRna(MMM, maintitle, colorfill)   

# TotalProbeSignal
up = which(duplicated(dd$genes$ProbeName) == FALSE)
ddaux = dd[up, ]
MMM <- MMMfunction(ddaux$TPS)
maintitle = "TotalProbeSignal"
colorfill = "red"

boxplotMicroRna(MMM, maintitle, colorfill)
plotDensityMicroRna(MMM, maintitle)
RleMicroRna(MMM, maintitle="TotalProbeSignal - RLE ", colorfill)    

# TotalGeneSignal
ddTGS = tgsMicroRna(dd, offset, half = FALSE, makePLOT = FALSE, 
                    verbose = FALSE)
MMM = log2(ddTGS$TGS)
maintitle = "TotalGeneSignal"
colorfill = "green"

boxplotMicroRna(MMM, maintitle, colorfill)
plotDensityMicroRna(MMM, maintitle)
maintitle = " TotalGeneSignal - RLE "
RleMicroRna(MMM, maintitle, colorfill)    

# BGMedianSignal
MMM = log2(dd$other$gBGMedianSignal)
maintitle = "BGMedianSignal"
colorfill = "yellow"
boxplotMicroRna(MMM, maintitle, colorfill)

# BGUsed
MMM = log2(dd$other$gBGUsed)
maintitle = "BGused"
colorfill = "cyan"
boxplotMicroRna(MMM, maintitle, colorfill)

mvaMicroRna2 <- function (uRNAList, maintitle, verbose = FALSE) 
{
  if (verbose) {
    cat("\n")
    cat("------------------------------------------------------", 
        "\n")
    cat("mvaMicroRna info:", "\n")
  }
  nGEN = c(1:dim(uRNAList$meanS)[1])
  indexrep = which(uRNAList$genes$ControlType == 0)
  LL = length(indexrep)
  if (verbose) {
    cat(" FEATURES :\t", LL, "\n")
  }
  index1 = which(uRNAList$genes$ControlType == 1)
  index2 = which(uRNAList$genes$GeneName[index1] != "DarkCorner")
  index3 = which(uRNAList$genes$GeneName[index1[index2]] != 
                   "miRNABrightCorner30")
  index4 = which(uRNAList$genes$GeneName[index1[index2[index3]]] != 
                   "SCorner3")
  indexpos = nGEN[index1[index2[index3[index4]]]]
  LL = length(indexpos)
  if (verbose) {
    cat(" POSITIVE CTRL:\t\t", LL, "\n")
  }
  indexneg = which(uRNAList$genes$ControlType == -1)
  LL = length(indexneg)
  if (verbose) {
    cat(" NEGATIVE CTRL:\t\t", LL, "\n")
  }
  indexDC = which(uRNAList$genes$GeneName == "DarkCorner")
  indexGE = which(uRNAList$genes$GeneName == "miRNABrightCorner30")
  indexSC = which(uRNAList$genes$GeneName == "SCorner3")
  LL = length(indexDC) + length(indexGE) + length(indexSC)
  if (verbose) {
    cat(" STRUCTURAL:\t\t", LL, "\n")
  }
  nARR = dim(uRNAList$meanS)[2]
  y = apply(uRNAList$meanS, 1, median)
  for (i in 1:nARR) {
    dev.new()
    if (!missing(maintitle)) {
      what = paste(maintitle, colnames(uRNAList$meanS)[i], 
                   " - ", "genewise Median")
    }
    else {
      what = paste(colnames(uRNAList$meanS)[i], " - ", 
                   "genewise Median")
    }
    x = uRNAList$meanS[, i]
    A = (x + y)/2
    M = (x - y)
    mva = plot(A, M, type = "p", cex = 0.7, col = "blue", 
               xlab = "A", ylab = "M")
    points(A[indexrep], M[indexrep], cex = 0.3, col = "cyan3", 
           pch = 19)
    points(A[indexDC], M[indexDC], cex = 0.5, col = "orange", 
           pch = 19)
    points(A[indexpos], M[indexpos], cex = 0.5, col = "red", 
           pch = 19)
    points(A[indexneg], M[indexneg], cex = 0.5, col = "green", 
           pch = 19)
    points(A[indexGE], M[indexGE], cex = 1, col = "yellow", 
           pch = 19)
    points(A[indexSC], M[indexSC], cex = 1, col = "pink", 
           pch = 19)
    smooth.fit = fitted(loess(M ~ A))
    points(A, smooth.fit, col = "black", cex = 0.5, pch = 19)
    title(main = what)
    abline(0, 0, col = "black", lty = 2, lwd = 1)
    abline(2, 0, col = "navyblue")
    abline(-2, 0, col = "navyblue")
    colors = c("cyan", "red", "green", "cyan3", "orange", 
               "yellow", "pink", "black")
    samples = c("features", "pos Ctrl", "neg Ctrl", "Rep NonCtrl", 
                "DarkCorner", "miRNABrightCorner", "SCorner3", "M~A smooth fit")
    legend(x = "topright", legend = samples, cex = 0.8, fill = colors, 
           inset = 0.05)
  }
}


qcPlotsJPEG <- function(dd, offset = 5, MeanSignal = TRUE, ProcessedSignal = FALSE, 
                        TotalProbeSignal = FALSE, TotalGeneSignal = FALSE, BGMedianSignal = FALSE, BGUsed = FALSE, targets) 
{
  if (!is(dd, "uRNAList")) {
    stop("'input' must be a uRNAList")
    if (is.null(dim(dd)[1])) {
      stop("'input' is empty")
    }
  }
  if (MeanSignal) {
    MMM = dd$meanS
    min = min(MMM)
    for (i in 1:dim(MMM)[2]) {
      MMM[, i] = MMM[, i] + (abs(min) + offset)
    }
    MMM = log2(MMM)
    maintitle = "MeanSignal"
    colorfill = "orange"
    dev.new()
    jpeg(filename = "MeanSignalBoxPlot.jpg", quality = 100, width = 653, height = 653)
    boxplotMicroRna(MMM, maintitle, colorfill)
    dev.off()
    dev.new()
    jpeg(filename = "MeanSignalDensityPlot.jpg", quality = 100, width = 653, height = 653)
    plotDensityMicroRna(MMM, maintitle)
    dev.off()       
    dev.new()
    jpeg(filename = "MeanSignalMVAPlot.jpg", quality = 100, width = 653, height = 597)
    ddaux = dd
    ddaux$meanS = MMM
    mvaMicroRna(ddaux, maintitle, verbose = FALSE)
    rm(ddaux)
    dev.off()
    maintitle = "MeanSignal- RLE "
    dev.new()
    jpeg(filename = "MeanSignalRLE.jpg", quality = 100, width = 653, height = 653)
    RleMicroRna(MMM, maintitle, colorfill)
    dev.off()
    dev.new()
    jpeg(filename = "ClusterEuclidean.jpg", quality = 100, width = 653, height = 653)
    hierclusMicroRna(MMM, targets$GErep, methdis = "euclidean", 
                     methclu = "complete", sel = FALSE, 100)
    dev.off()
  }
  if (ProcessedSignal) {
    jpeg(filename = "ProcessedSignal.jpg", quality = 100, width = 653, height = 653)
    MMM = dd$procS
    min = min(MMM)
    for (i in 1:dim(MMM)[2]) {
      MMM[, i] = MMM[, i] + (abs(min) + offset)
    }
    MMM = log2(MMM)
    maintitle = "ProcessedSignal"
    colorfill = "blue"
    dev.new()
    boxplotMicroRna(MMM, maintitle, colorfill)
    dev.new()
    plotDensityMicroRna(MMM, maintitle)
    dev.new()
    ddaux = dd
    ddaux$TPS = MMM
    mvaMicroRna(ddaux, maintitle, verbose = FALSE)
    rm(ddaux)
    maintitle = "ProcessedSignal - RLE "
    dev.new()
    RleMicroRna(MMM, maintitle, colorfill)
    dev.off()
  }
  if (TotalProbeSignal) {
    jpeg(filename = "TotalProbeSignal.jpg", quality = 100, width = 653, height = 653)
    up = which(duplicated(dd$genes$ProbeName) == FALSE)
    ddaux = dd[up, ]
    MMM = ddaux$TPS
    min = min(MMM)
    for (i in 1:dim(MMM)[2]) {
      MMM[, i] = MMM[, i] + (abs(min) + offset)
    }
    MMM = log2(MMM)
    maintitle = "TotalProbeSignal"
    colorfill = "red"
    dev.new()
    boxplotMicroRna(MMM, maintitle, colorfill)
    dev.new()
    plotDensityMicroRna(MMM, maintitle)
    maintitle = " TotalProbeSignal - RLE "
    dev.new()
    RleMicroRna(MMM, maintitle, colorfill)
    dev.off()
  }
  if (TotalGeneSignal) {
    jpeg(filename = "TotalGeneSignal.jpg", quality = 100, width = 653, height = 653)
    ddTGS = tgsMicroRna(dd, offset, half = FALSE, makePLOT = FALSE, 
                        verbose = FALSE)
    MMM = log2(ddTGS$TGS)
    maintitle = "TotalGeneSignal"
    colorfill = "green"
    dev.new()
    boxplotMicroRna(MMM, maintitle, colorfill)
    dev.new()
    plotDensityMicroRna(MMM, maintitle)
    maintitle = " TotalGeneSignal - RLE "
    dev.new()
    RleMicroRna(MMM, maintitle, colorfill)
    dev.off()
  }
  if (BGMedianSignal) {
    jpeg(filename = "BGMedianSignal.jpg", quality = 100, width = 653, height = 653)
    MMM = log2(dd$other$gBGMedianSignal)
    maintitle = "BGMedianSignal"
    colorfill = "yellow"
    dev.new()
    boxplotMicroRna(MMM, maintitle, colorfill)
    dev.off()
  }
  if (BGUsed) {
    jpeg(filename = "BGUsed.jpg", quality = 100, width = 653, height = 653)
    MMM = log2(dd$other$gBGUsed)
    maintitle = "BGused"
    colorfill = "cyan"
    dev.new()
    boxplotMicroRna(MMM, maintitle, colorfill)
    dev.off()
  }
}

#Modified function to save all plots (PDF file)
qcPlotsPDF <- function (dd, offset = 5, MeanSignal = TRUE, ProcessedSignal = FALSE, 
                        TotalProbeSignal = FALSE, TotalGeneSignal = FALSE, BGMedianSignal = FALSE, BGUsed = FALSE, targets) 
{
  if (!is(dd, "uRNAList")) {
    stop("'input' must be a uRNAList")
    if (is.null(dim(dd)[1])) {
      stop("'input' is empty")
    }
  }
  pdf(file="qcPlotsReport.pdf")
  if (MeanSignal) {
    MMM = dd$meanS
    min = min(MMM)
    for (i in 1:dim(MMM)[2]) {
      MMM[, i] = MMM[, i] + (abs(min) + offset)
    }
    MMM = log2(MMM)
    maintitle = "MeanSignal"
    colorfill = "orange"
    boxplotMicroRna(MMM, maintitle, colorfill)
    plotDensityMicroRna(MMM, maintitle)
    ddaux = dd
    ddaux$meanS = MMM
    mvaMicroRna(ddaux, maintitle, verbose = FALSE)
    rm(ddaux)
    maintitle = "MeanSignal- RLE "
    RleMicroRna(MMM, maintitle, colorfill)
    hierclusMicroRna(MMM, targets$GErep, methdis = "euclidean", 
                     methclu = "complete", sel = FALSE, 100)
  }
  if (ProcessedSignal) {
    MMM = dd$procS
    min = min(MMM)
    for (i in 1:dim(MMM)[2]) {
      MMM[, i] = MMM[, i] + (abs(min) + offset)
    }
    MMM = log2(MMM)
    maintitle = "ProcessedSignal"
    colorfill = "blue"        
    boxplotMicroRna(MMM, maintitle, colorfill)       
    plotDensityMicroRna(MMM, maintitle)       
    ddaux = dd
    ddaux$TPS = MMM
    mvaMicroRna(ddaux, maintitle, verbose = FALSE)
    rm(ddaux)
    maintitle = "ProcessedSignal - RLE " 
    RleMicroRna(MMM, maintitle, colorfill)
  }
  if (TotalProbeSignal) {
    up = which(duplicated(dd$genes$ProbeName) == FALSE)
    ddaux = dd[up, ]
    MMM = ddaux$TPS
    min = min(MMM)
    for (i in 1:dim(MMM)[2]) {
      MMM[, i] = MMM[, i] + (abs(min) + offset)
    }
    MMM = log2(MMM)
    maintitle = "TotalProbeSignal"
    colorfill = "red"
    boxplotMicroRna(MMM, maintitle, colorfill)
    plotDensityMicroRna(MMM, maintitle)
    maintitle = " TotalProbeSignal - RLE "
    RleMicroRna(MMM, maintitle, colorfill)
  }
  if (TotalGeneSignal) {
    ddTGS = tgsMicroRna(dd, offset, half = FALSE, makePLOT = FALSE, 
                        verbose = FALSE)
    MMM = log2(ddTGS$TGS)
    maintitle = "TotalGeneSignal"
    colorfill = "green"
    boxplotMicroRna(MMM, maintitle, colorfill)
    plotDensityMicroRna(MMM, maintitle)
    maintitle = " TotalGeneSignal - RLE "
    RleMicroRna(MMM, maintitle, colorfill)
  }
  if (BGMedianSignal) {
    MMM = log2(dd$other$gBGMedianSignal)
    maintitle = "BGMedianSignal"
    colorfill = "yellow"
    boxplotMicroRna(MMM, maintitle, colorfill)
  }
  if (BGUsed) {
    MMM = log2(dd$other$gBGUsed)
    maintitle = "BGused"
    colorfill = "cyan"
    boxplotMicroRna(MMM, maintitle, colorfill)
  }
  dev.off()
}

#All Plots in PDF
qcPlotsPDF(dd.micro,offset=5,
           MeanSignal=TRUE,
           ProcessedSignal=FALSE,
           TotalProbeSignal=FALSE,
           TotalGeneSignal=FALSE,
           BGMedianSignal=FALSE,
           BGUsed=FALSE,
           targets.micro)

#BoxPlot
boxplotMicroRna(log2(dd.micro$meanS), maintitle='log2 Mean Signal', colorfill= 'orange')
#Density
plotDensityMicroRna(log2(dd.micro$meanS), maintitle='log2 Mean Signal')
#log2 Mean Signal
ddaux=dd.micro
ddaux$G=log2(dd.micro$meanS)
mvaMicroRna(ddaux, maintitle='log2 Mean Signal', verbose=FALSE)

rm(ddaux)
#Mean Signal RLE
RleMicroRna(log2(dd.micro$meanS), maintitle='log2 Mean Signal - RLE')
#Clustering - Euclidean
hierclusMicroRna(log2(dd.micro$meanS),targets.micro$GErep,
                 methdis="euclidean",
                 methclu="complete",
                 sel=TRUE,100)
#Save to file template
jpeg(filename = "name.jpg", quality = 100, width = 653, height = 653)
pdf(file = "PlotsReport.pdf")
dev.off()

#Linear Model
# levels.treatment=levels(factor(targets.micro$Treatment))
# treatment=factor(as.character(targets.micro$Treatment),levels=levels.treatment)
# levels.subject=levels(factor(targets.micro$Subject))
# subject=factor(as.character(targets.micro$Subject),levels=levels.subject)
# #Fitting the Model
# design=model.matrix(~ -1 + treatment)
# print(design)
# fit=lmFit(esetPROC,design)
# names(fit)
# print(head(fit$coeff))
# CM=cbind(MSC_AvsMSC_B=c(1,-1,0),MSC_AvsMSC_C=c(1,0,-1))
# print(CM)
# fit2=contrasts.fit(fit,CM)
# names(fit2)
# print(head(fit2$coeff))
# fit2=eBayes(fit2)
# names(fit2)
# fit2=basicLimma(esetPROC,design,CM,verbose=TRUE)

# MMM = dd.micro$meanS
# min = min(MMM)
# for (i in 1:dim(MMM)[2]) 
# {
#     MMM[, i] = MMM[, i] + (abs(min) + 5)
# }
# MMM = log2(MMM)
# colnames(MMM) <- c("Naive1", "Naive2",  "Naive3",  "Naive4", "Teff1", "Teff2", "Teff3", "Teff4", "iTreg1", "iTreg2", "iTreg3", "iTreg4")
# hierclusMicroRna(MMM, targets.micro$GErep, methdis = "euclidean", methclu = "complete", sel = FALSE, 100)
