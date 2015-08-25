# Somple CRC Analysis




# Example for Treatment(n=3) and Control(n=3)
targets2 <- data.frame(FileName = list.files()[c(-7,-8)],
						Condition = c(rep("PM101",3), rep("Ctr", 3)),
						Sample = rep(c(1,2,3),2),
						stringsAsFactors=FALSE)
targets2 <- cbind(targets2,Nomeclature=paste(targets2$Condition,targets2$Sample,sep="_"))
# write.table(x=targets, file='Targets.txt', quote=FALSE, row.names=FALSE, sep='\t')

# Read Targets File
targets <- read.table(file="Targets.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)

# targets: Targets file
# data.path: Path to read RawData
# coefs: Coefficients of interest
runAnalysis <- function(targets, data.path='RawData', coefs='Treat1-Treat2'){

	cat('Dinamic Agilent Single-Channel Analysis.\n')
	cat('Pipeline steps: > Read raw data\n> Background correction\n> Normalize arrays\n> Filter probes\n> Average replicated probes\n> Create experiment design\n> Fit linear models\n> Make contrats of interest\n> Compute Bayesian Statistics\n> Save Plots\n> Save results\n> Save R data\n\n')
	if(!file.exists(data.path))
		stop('Invalid path.')

	cat('Loading libraries and functions...\n')
	library(limma)
	source("/Users/Baiochi/Dropbox/USP/Lab/Agilent Bioc/Functions.R")

	# Read raw data
	cat('\n\nReading raw data\n')
	raw <- read.maimages(targets$FileName, source="agilent", green.only=TRUE, path=data.path)

	#Bgc correct
	cat('Background correction:\n')
	input <- readBgc()
	bgc <- backgroundCorrect(raw,method=input)
	rm(input)

	#Normalization
	cat('Normalization between arrays:\n')
	input <- readNorm()
	norm <- normalizeBetweenArrays(bgc,method=input)
	rm(input)

	#Filtering
	cat('Filter probes? (y/n)\n')
	input <- ask()
	if(input=='y'){
		neg95 = apply(norm$E[norm$genes$ControlType==-1,],2,function(x) quantile(x,p=0.95))
		cutoff = matrix(1.1*neg95,nrow(norm),ncol(norm),byrow=TRUE)
		isexpr = rowSums(norm$E > cutoff) >= length(unique(targets$Sample))
		#table(isexpr)
		eset = norm[norm$genes$ControlType==0 & isexpr,]
		cat('Probes filtered.\n')
	}
	else{
		eset = norm
	}
	rm(input)

	#Averaging replicated spots
	eset <- avereps(eset,ID=eset$genes[,"SystematicName"])

	#Design
	cat('Creating design matrix...\n')
	f <- factor(targets$Condition, levels = unique(targets$Condition))
	design <- model.matrix(~0 + f)
	colnames(design) <- levels(f)

	#Fit linear models
	cat('Fitting linear model...\n')
	fit <- lmFit(eset$E, design)

	#Contrast Matrix
	#coefs = colnames(design)
	cat('Creating contrast matrix...\n')
	contrast.matrix <- makeContrasts(contrasts=coefs, levels=design)

	#Compute Bayesian Statistics
	cat('Computing Bayesian Statistics...\n')
	fit2 <- contrasts.fit(fit, contrast.matrix)
	fit2 <- eBayes(fit2)

	#Rank Genes
	rank <- decideTests(fit2, method="separate", adjust.method="fdr", p.value=0.05, lfc=0)
	cat('Differential Expression summary:\n')
	summary(rank)

	#Saving Plots
	input = readline('Enter a directory to save plots: ')
	if(!file.exists(input)){
		cat('Directory doesn\'t exist, creating a new one...\n')
		dir.create(input)
	}
	tempdir = getwd()
	setwd(input)

	RawPlots(raw)
	NormPlots(eset)
	DiffExprsPlots(fit2, coefs, eset)

	setwd(tempdir)
	rm(input)

	#Write down Results
	input = readline('Enter a directory to save results: ')
	if(!file.exists(input)){
		cat('Directory doesn\'t exist, creating a new one...\n')
		dir.create(input)
	}
	tempdir = getwd()
	setwd(input)

	WriteResults(fit2, coefs, eset)

	setwd(tempdir)
	rm(input)
}


runAnalysis(targets= targets, data.path="/Users/Baiochi/Dropbox/USP/Lab/RawData/HCT", coefs='PM101-Ctr')




#profile plot with the bcl
n <- which(eset$genes$GeneName=='BCL2')
plot(eset$E[n,], ylab='Expression', xlab='Arrays')

#BCL2 probe: ATCAGAGTTGTTGCTTCCCGGCGTCCCTACCTCCTCCTCTGGACAAAGCGTTCACTCCCA

BCL2_F 1815 GAAGTCTGGGAATCGATCTGG
BCL2_R 1816 TCCCATCAATCTTCAGCACTC














