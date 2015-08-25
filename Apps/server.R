#Server

library(shiny)
library(limma)

source("/Users/Baiochi/Dropbox/USP/Lab/Functions.R")
targets <- read.table(file='/Users/Baiochi/Dropbox/USP/Lab/Naive.iTreg.miRNA/Targets.txt', header=TRUE, sep="\t", stringsAsFactors=FALSE)
targets <- cbind(targets,Nomeclature=paste(targets$Condition,targets$Sample,sep="_"))

raw <- read.maimages(files = targets$FileName,
				path = '/Users/Baiochi/Dropbox/USP/Lab/Naive.iTreg.miRNA/RawData',
				source = "agilent", 
				green.only = TRUE, 
				verbose = TRUE)
colnames(raw$E) <- targets$Nomeclature
colnames(raw$Eb) <- targets$Nomeclature
rownames(raw$targets) <- targets$Nomeclature

dataset <- raw
# Background Correction
dataset <- backgroundCorrect(dataset, method="normexp", offset=1)
bgc.plot <- dataset

# Normalization
dataset <- normalizeBetweenArrays(dataset, method="quantile")
norm.plot <- dataset

# Create ExpressionSet
exprs <- new("MAList", list(targets=dataset$targets, 
							genes=dataset$genes, 
							source=dataset$source, 
							M=dataset$E, 
							A=dataset$E))
exprs <- avereps(exprs, ID=exprs$genes$ProbeName)
print('Finished')

shinyServer(function(input, output) {

	stepInput <- reactive({
		if(input$step=='raw')
			step <- raw$E
		if(input$step=='norm')
			step <- norm.plot$E
		if(input$step=='exprs')
			step <- exprs$M
		return(step)
  	})

	output$plot <- renderPlot({
		cluster.tree(stepInput(), title="Cluster Dendogram", 
			dist=input$dist, clust=input$clust)
	})

})
