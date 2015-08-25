# ChIP ChIP Bioconductor Tutorial Pipeline

# Load packages
source("/Users/Baiochi/Dropbox/USP/Lab/Analises/Functions.R")
packages <- c("ggplot2", "limma", "Ringo", "biomaRt", "topGO", "ccTutorial")
ipak(packages)

pairDir <- system.file("PairData",package="ccTutorial")
list.files(pairDir, pattern="pair$")

read.delim(file.path(pairDir,"files_array1.txt"), header=TRUE)

RGs <- lapply(sprintf("files_array%d.txt",1:4), readNimblegen, "spottypes.txt", path=pairDir)
# save(RGs, file='ccRGs.rda')
load('ccRGs.rda')

# Visualize data
head(RGs[[1]]$R)
head(RGs[[1]]$G)
tail(RGs[[1]]$genes)
table(RGs[[1]]$genes$Status)

# Read files into RGList
RG1breaks <- c(0,quantile(RGs[[1]]$G, probs=seq(0,1,by=0.1)),2^16)

jpeg("ccTutorialArrayImages.jpg", quality=100, height=1074*1.5, width=768*1.5)
par(mar=c(0.01,0.01,2.2,0.01))
layout(matrix(c(1,2,5,6,3,4,7,8,9,10,13,14,11,12,15,16), ncol=4,byrow=TRUE))
for (this.set in 1:4){
	thisRG <- RGs[[this.set]]
	for (this.channel in c("green","red")){
		my.colors <- colorRampPalette(c("black",paste(this.channel,c(4,1),sep="")))(length(RG1breaks)-1)
		for (arrayno in 1:2){
			image(thisRG, arrayno, channel=this.channel,mybreaks=RG1breaks, mycols=my.colors)
			mtext(side=3, line=0.2, font=2, text=colnames(thisRG[[toupper(substr(this.channel,1,1))]])[arrayno])
		}
	}
}
dev.off()

corPlot(log2(RGs[[2]]$G))
corPlot(log2(RGs[[2]]$R))


probeAnno <- posToProbeAnno(file.path(system.file("exonerateData", package="ccTutorial"), "allChromExonerateOut.txt"))
allChrs <- chromosomeNames(probeAnno)
genome(probeAnno) <- "M. musculus (mm9)"
arrayName(probeAnno) <- "2005-06-17_Ren_MM5Tiling"
show(probeAnno)

