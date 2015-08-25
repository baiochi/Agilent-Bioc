

# MODIFY CONTRASTS OF INTEREST

# Experiment
single.run("/Targets.txt",
			file.path="/Experiment",
			rawdata.path="/RawData",
			save.plot=TRUE,
			save.data=TRUE)

single.run <- function(target.file,
	file.path,
	rawdata.path,
	bgc=TRUE,
	norm=TRUE,
	avereps=TRUE,
	save.plot=FALSE,
	save.data=FALSE,
	plot.path=NULL,
	data.path=NULL,
	markdown){

	#Targets format example:
	#	FileName 	Condition 	Sample	
	# file1.txt 	Ctr 		1 		
	# file2.txt 	Ctr 		2 		
	# file3.txt 	Ctr 	 	3 		
	# file4.txt 	Treatment 	1 		
	# file5.txt 	Treatment 	2 		
	# file6.txt 	Treatment 	3 		

	# Data Integrity Check
	cat("Checking data...\n")
	if(missing(target.file))
		stop("Please input target.file name.")
	if(missing(file.path)){
		file.path="." 
		cat("Warning: file.path not defined, setting default to current directory.\n")
	}
	default.path <- getwd()
	if(missing(rawdata.path)){
		rawdata.path="."
		cat("Warning: rawdata.path not defined, setting default to current directory.\n")
	}

	targets <- read.table(file=target.file, header=TRUE, sep="\t", stringsAsFactors=FALSE)

	if(!is(targets, "data.frame"))
		stop("Targets file must be a Data Frame")
	if(!colnames(targets)[1]=="FileName")
		stop("Invalid file format")
	if(!colnames(targets)[2]=="Condition")
		stop("Invalid file format")
	if(!colnames(targets)[3]=="Sample")
		stop("Invalid file format")

	targets <- cbind(targets,Nomeclature=paste(targets$Condition,targets$Sample,sep="_"))

	cat("Ok.\n")

	cat("Loading built-in functions...\n")
	# Pre-load functions
	source("/Users/Baiochi/Dropbox/USP/Lab/Functions.R")
	cat("Ok.\n")

	cat("Loading packages...\n")
	library(limma)
	cat("Ok.\n")

	cat("Reading raw data...\n")
	raw <- read.maimages(files = targets$FileName,
					path = rawdata.path,
					source = "agilent", 
					green.only = TRUE, 
					verbose = TRUE)
	colnames(raw$E) <- targets$Nomeclature
	colnames(raw$Eb) <- targets$Nomeclature
	rownames(raw$targets) <- targets$Nomeclature

	cat("Normalizing data...\n")

	dataset <- raw

	# Background Correction
	if(bgc){
		dataset <- backgroundCorrect(dataset, method="normexp", offset=1)
		bgc.plot <- dataset
	}

	# Normalization
	if(norm){
		dataset <- normalizeBetweenArrays(dataset, method="quantile")
		norm.plot <- dataset
	}

	# Create ExpressionSet
	exprs <- new("MAList", list(targets=dataset$targets, 
								genes=dataset$genes, 
								source=dataset$source, 
								M=dataset$E, 
								A=dataset$E))

	# Filter Probes
	if(avereps){
		exprs <- avereps(exprs, ID=exprs$genes$ProbeName)
	}

	cat("Computing statistics for differential expression...\n")

	# DiffExprs
	f <- factor(targets$Condition, levels = unique(targets$Condition))

	cat("design matrix...\n")
	design <- model.matrix(~0 + f)
	colnames(design) <- levels(f)

	cat("fitting model...\n")
	fit <- lmFit(exprs$M, design)

	cat("contrast matrix...\n")
	contrast <- c("*********************CONTRATS OF INTEREST HERE********************************")
	contrast.matrix <- makeContrasts("*CONTRATS OF INTEREST HERE**", levels=design)

	cat("finishing analysis...\n")
	fit2 <- contrasts.fit(fit, contrast.matrix)
	fit2 <- eBayes(fit2)

	cat("preparing data...\n")
	tests<- new("list")
	for (i in 1:length(colnames(contrast.matrix))) {
		tests[[i]] <- topTable(fit2, adjust="fdr", coef=colnames(contrast.matrix)[i], genelist=exprs$genes, number=Inf)
		rownames(tests[[i]]) <- NULL
	}

	cat("calculating differential expression...\n")
	diff.exprs <- decideTests(fit2, method="separate", adjust.method="fdr", p.value=0.05, lfc=0)
	print(summary(diff.exprs))
	
	# Save Plots
	if(save.plot){
		cat("Generating plots...\n")

		if(is.null(plot.path)){
			plot.path=paste(default.path, "/Plots", sep="")
			dir.create(path=plot.path)
		}
		else{
			dir.create(path=plot.path)
		}
		cat(paste("Saving plot files at '", plot.path, "/'\n", sep=""))
		setwd(plot.path)

		# Raw Data Visualization
		cluster.tree(raw$E, save=TRUE)
		cor.matrix(raw$E, save=TRUE)
		my.boxplot(log2(raw$E), 
			title="RawData Boxplot",
			color="forestgreen",
			save=TRUE)
		my.densityplot(log2(raw$E), 
			samples=dim(targets)[1],
			title="RawData Density",
			save=TRUE)
		my.distribuitionplot(log2(raw$E), 
			x=length(unique(targets$Condition)),
			y=length(unique(targets$Sample)),
			title="RawData Distribuition",
			save=TRUE)
		my.MAPlot(raw, 
			x=length(unique(targets$Condition)),
			y=length(unique(targets$Sample)),
			title="RawData M vs A Plot",
			save=TRUE)

		# if(bgc){
		# 	my.boxplot(log2(bgc.plot$E), 
		# 		title="BackgroundCorrection Boxplot",
		# 		color="palevioletred1",
		# 		save=TRUE)
		# 	my.densityplot(log2(bgc.plot$E), 
		# 		samples=dim(targets)[1],
		# 		title="BackgroundCorrection Density",
		# 		save=TRUE)
		# 	my.distribuitionplot(log2(bgc.plot$E), 
		# 		x=length(unique(targets$Condition)),
		# 		y=length(unique(targets$Sample)),
		# 		title="BackgroundCorrection Distribuition",
		# 		save=TRUE)
		# 	my.MAPlot(bgc.plot, 
		# 		x=length(unique(targets$Condition)),
		# 		y=length(unique(targets$Sample)),
		# 		title="BackgroundCorrection M vs A Plot",
		# 		save=TRUE)
		# }
		if(norm){
			my.boxplot(norm.plot$E, 
				title="Normalization Between Arrays Boxplot",
				color="steelblue3",
				save=TRUE)
			my.densityplot(norm.plot$E, 
				samples=dim(targets)[1],
				title="Normalization Between Arrays Density",
				save=TRUE)
			my.distribuitionplot(norm.plot$E, 
				x=length(unique(targets$Condition)),
				y=length(unique(targets$Sample)),
				title="Normalization Between Arrays Distribuition",
				save=TRUE)
			my.MAPlot(norm.plot, 
				x=length(unique(targets$Condition)),
				y=length(unique(targets$Sample)),
				title="Normalization Between Arrays M vs A Plot",
				save=TRUE)

		}
		if(avereps){
			my.boxplot(exprs$M, 
				title="Normalized with Average Probes Boxplot",
				color="gold",
				save=TRUE)
			my.densityplot(exprs$M, 
				samples=dim(targets)[1],
				title="Normalized with Average Probes Density",
				save=TRUE)
			my.distribuitionplot(exprs$M, 
				x=length(unique(targets$Condition)),
				y=length(unique(targets$Sample)),
				title="Normalized with Average Probes Distribuition",
				save=TRUE)
		}
		for (i in 1:length(contrast)) {
			my.volcanoplot(tests[[i]]$logFC, tests[[i]]$adj.P.Val, 
						title=paste(contrast[i], " Volcano Plot"),
						save=TRUE)
		}

		setwd(default.path)
		cat("Done.\n")
	}
	
	# Save Data
	if(save.data){
		cat("Saving Data...\n")
		if(is.null(data.path)){
			data.path=paste(default.path, "/DiffExprs", sep="")
			dir.create(path=data.path)
		}
		else{
			dir.create(path=data.path)
		}
		cat(paste("Saving data files at '", data.path, "/'\n", sep=""))
		setwd(data.path)

		for (i in 1:length(colnames(contrast.matrix))) {
			cat(paste("Saving ", contrast[i], " DiffExprs.txt\n"))
			write.results(tests[[i]], filename=paste(contrast[i], "DiffExprs.txt"))
		}
		write.table(as.data.frame(diff.exprs), file="Differential Expression Summary.txt", sep="\t", quote=FALSE, row.names=FALSE)

		setwd(default.path)
		cat("Done.\n")
	}

	# Generate Analysis Report
	if(!missing(markdown)){
		library(knitr)
		cat("Generating Analysis Report...\n")
		knit2html(markdown, output = NULL,
			envir = parent.frame(), 
			text = NULL, 
			quiet = FALSE, 
			encoding = getOption("encoding"))
		print("Done.\n")
	}
	
	cat("Analysis complete.\n")

}


