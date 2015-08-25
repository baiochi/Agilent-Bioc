


foo <- function(bar, bah){
	if(!is(bar, "numeric"))
		stop("Abandon Ship!")
	else
		print("Yaaar!")
	if(missing(bah))
		stop("bah is missing")
	print("42")
}

foo <- function(pt=NULL){
	if(is.null(pt)){
		pt="DefaultName"
		dir.create(path=pt)
	}
	else{
		dir.create(path=pt)
	}
	print(pt)
}

targets <- data.frame(Filename = list.files(),
						Condition = c(rep("Naive",3), rep("Teff", 3), rep("iTreg", 3)),
						Sample = rep(c(1,2,3),3))




##########################################################################################
##########################################################################################

#Filter control probes
neg <- which(exprs$genes$ControlType==-1)
ctr.probes <- which(exprs$genes$ControlType==1)
ctr.probes[93] <- neg
exprs$M <- exprs$M[-ctr.probes,]

#T-test parameters
ctr <- colnames(exprs$M)[1:3]
pm29 <- colnames(exprs$M)[4:6]

#Students T-Test
test <- apply(exprs$M,1,t.test.all.genes,s1=ctr,s2=pm29)

#Table of significant genes(TRUE)
table(test<0.05) 

#Ploting
hist(test, col="yellow")

#Fold Change
ctr.mean <- apply(exprs[,ctr],1,mean,na.rm=TRUE)
pm29.mean <- apply(exprs[,pm29],1,mean,na.rm=TRUE)

#Get fold changes
fc <- log2(ctr.mean/pm29.mean)

#Histograms
hist(fc, col="yellow")

#Multiple Testing Correction
adj.Pval <- p.adjust(test, method = "fdr")
hist(adj.Pval, col="green")
table(adj.Pval<0.05)

##########################################################################################
##########################################################################################
# 										KNIT 											 #
##########################################################################################
##########################################################################################

knit2html(input, output = NULL,
		envir = parent.frame(), 
		text = NULL, 
		quiet = FALSE, 
		encoding = getOption("encoding"))

@chunck options: results = "asis"
kable(x, format, digits = getOption("digits"), row.names = NA, col.names = colnames(x), 
    align, caption = NULL, escape = TRUE, ...)

##########################################################################################
##########################################################################################
# 									 ArrayQuality 										 #
##########################################################################################
##########################################################################################

#ArrayQualityReport
RawData <- new("MAList", list(targets=raw$targets, 
						genes=raw$genes, 
						source=raw$source, 
						M=raw$E, 
						A=raw$E))
library(arrayQualityMetrics)
arrayQualityMetrics(expressionset = RawData,
					outdir = "Array Quality Report",
					force = TRUE)





####### Running SAM

	install.packages(c("samr", "matrixStats", "GSA", "shiny", "openxlsx"))
	source("http://bioconductor.org/biocLite.R")
	biocLite("impute")

	library(shiny)
	runGitHub("SAM", "MikeJSeo")

##################
# Houtan
##################

	# scatter plot matrix
	dat <- read.table("C:\\gecolon.dat",header=T)
	dimnames(dat)[[1]] <- as.character(dat[,1])
	dat <- dat[,-1]; dat <- as.data.frame(dat)

	# other data sets in R to use
	library(Biobase)
	library(annotate)
	data(geneData)
	dat <- geneData

	# Pearson's correlation matrix
	dat.cor <- cor(dat)
	image(dat.cor,axes=F)
	axis(2,at=seq(0,1,length=ncol(dat.cor)),label=dimnames(dat.cor)[[2]])
	axis(3,at=seq(0,1,length=ncol(dat.cor)),label=dimnames(dat.cor)[[2]])


	# Profile plot
	rand.genes <- sample(dimnames(dat)[[1]],5,replace=F)
	plot(c(1,ncol(dat)),range(dat[rand.genes,]),type='n',main="Profile plot of 5 random
			genes",xlab="Samples",ylab="Expression")
	for(i in 1:length(rand.genes)) {
		dat.y <- as.numeric(dat[rand.genes[i],])
		lines(c(1:ncol(dat)),dat.y,col=i)
	}

	# load the yeast cell cycle data set
	dat <- read.table(“C:\\spellman.txt”,header=T)
	dimnames(dat)[[1]] <- as.character(dat[,1])
	dat <- dat[,-1]
	dat <- dat[,23:46]
	dat[is.na(dat)] <- 0

	# pca biplot
	biplot(prcomp(t(dat[450:500,])),cex=0.6)

	# k-means cluster profiles
	dd <- dat[names(f.p)[f.p<0.001],]
	d.k <- kmeans(dd,9)
	par(mfrow=c(3,3))
	for(i in 1:9) {
		tmp <- scale(dd[d.k$cluster==i,])
		matplot(c(1:ncol(dat)),t(tmp),type='l',col=i,xlab='Time',ylab='Expression')
	}

	# cv vs. mean plot
	dat.mean <- apply(dat,1,mean) # calculate mean for each gene
	dat.sd <- sqrt(apply(dat,1,var)) # calculate st.deviation for each gene
	dat.cv <- dat.sd/dat.mean #calculate cv
	plot(dat.mean,dat.cv,main="Sample CV vs. Mean",xlab="Mean",ylab="CV",col='blue',cex=1.5)

	# 2D sample pca plot
	dat.pca <- prcomp(t(dat))
	dat.loads <- dat.pca$x[,1:2]
	plot(dat.loads[,1],dat.loads[,2],main="Sample PCA plot",xlab="p1",ylab="p2",col='red',cex=1.5,pch=16)


	# k-means clustering for missing value imputation
	dat <- dat[2:30,] # only use 29 genes for example
	cl <- kmeans(dat[,-1],centers=5, iter.max=20) # cluster into 5 groups
	# we pretend to be missing a value at sample#1 gene #2
	groups <- cl$cluster # get cluster membership for each gene
	groups # look at groups to see where gene 2 is
	group.2 <- groups==2 # since gene 2 is in group 2, get all other members
	genes.cluster <- dimnames(dat)[[1]][group.2]
	genes.cluster # look at all other genes in cluster #2
	gene.dist <- dist(dat[genes.cluster,-1],method="euclidean") # get distances from genes in cluster 2 to
	# gene #2
	gene.dist <- as.matrix(gene.dist)
	gene.dist <- gene.dist[2:5,1]
	gene.weight <- as.numeric(gene.dist/sum(gene.dist)) # get weights for each gene
	weight.mean <- weighted.mean(dat[genes.cluster[-1],1], gene.weight) # calculate weighted mean for
	# gene #2

	# perspective plot
	data(volcano) # load volcano data set
	persp(volcano, theta=45, phi=30, col="red")

	# MvA plot
	library(sma)
	data(MouseArray)
	mouse.lratio <- stat.ma(mouse.data, mouse.setup)
	plot.mva(mouse.data, mouse.setup, norm="l", 2, extra.type="pci", plot.type="n",main="MvA plot")


	# calculate mean for some genes, with respect to class
	library(multtest)
	data(golub)
	dat <- as.data.frame(golub)
	ann <- golub.cl
	dat.aml <- apply(dat[,ann==1],1,mean)
	dat.all <- apply(dat[,ann==0],1,mean)
	tab <- data.frame(rbind(dat.aml[1:20],dat.all[1:20]))
	dimnames(tab)[[1]] <- c("AML","ALL")
	names(tab) <- dimnames(dat)[[1]][1:20]
	mp <- barplot(tab)
	tot <- colMeans(tab)
	text(mp, tot + 3, format(tot), xpd = TRUE, col = "blue")
	barplot(as.matrix(tab),beside=T,col=c("red","yellow"),legend=rownames(as.matrix(tab)),ylim=c(-
	5,5),ylab="Expression")
	title(main = "Mean Expression Levels of first 20 genes")
	# cluster tree
	dat <- t(dat) #transpose dat
	dat.dist <- dist(dat,method="euclidean") # calculate distance
	dat.clust <- hclust(dat.dist,method="single") # calculate clusters
	plot(dat.clust,labels=names(dat),cex=0.75) # plot cluster tree

	# calculate mean for some genes, with respect to class
	library(multtest)
	data(golub)
	dat <- as.data.frame(golub)
	ann <- golub.cl
	dat.aml <- apply(dat[,ann==1],1,mean)
	dat.all <- apply(dat[,ann==0],1,mean)
	tab <- data.frame(rbind(dat.aml[1:20],dat.all[1:20]))
	dimnames(tab)[[1]] <- c("AML","ALL")
	names(tab) <- dimnames(dat)[[1]][1:20]
	mp <- barplot(tab)
	tot <- colMeans(tab)
	text(mp, tot + 3, format(tot), xpd = TRUE, col = "blue")
	barplot(as.matrix(tab),beside=T,col=c("red","yellow"),
		legend=rownames(as.matrix(tab)),ylim=c(-5,5),ylab="Expression")
	title(main = "Mean Expression Levels of first 20 genes")

	# cluster tree
	dat <- t(dat) #transpose dat
	dat.dist <- dist(dat,method="euclidean") # calculate distance
	dat.clust <- hclust(dat.dist,method="single") # calculate clusters
	plot(dat.clust,labels=names(dat),cex=0.75) # plot cluster tree
##########################################################################################
##########################################################################################
#										Plot Functions  								 #
##########################################################################################
##########################################################################################

my.boxplot <- function(exprs, title="Boxplot", color="white", save=FALSE){
	if(save){
		png(filename= paste(title, ".png", sep = ""))
		boxplot(exprs,
			names = colnames(exprs),
			las=2,
			ylab="Spot Intensity",
			col = color,
			main = title)
	dev.off()
	}
	else{
		boxplot(exprs,
		names = colnames(exprs),
		las=2,
		ylab="Spot Intensity",
		col = color,
		main = title)
	}
}

my.densityplot <- function(exprs, samples, title="Density", save=FALSE){
	colors <- c("green", "red", "yellow", "blue", "purple", "orange",
			"green3", "red3", "yellow3", "blue3", "purple3", "orange3")
	if(save){
		png(filename= paste(title, ".png", sep = ""))
		plot(density(exprs[,1]), col=colors[1], main=title)
		for (i in 1:samples) {
			lines(density(exprs[,i]), col=colors[i])
		}
		legend("topright",
			legend = colnames(exprs),
			lty=c(1,1),
			lwd=c(2.5,2.5),col=colors)
		dev.off()
	}
	else{
		plot(density(exprs[,1]), col=colors[1], main=title)
		for (i in samples:j) {
			lines(density(exprs[,i]), col=colors[i])
		}
		legend("topright",
			legend = colnames(exprs),
			lty=c(1,1),
			lwd=c(2.5,2.5),col=colors)
	}
}

my.distribuitionplot <- function(exprs, x, y, title="Distribuition", save=FALSE){
	if(missing(x) | missing(y))
		stop("Please input the plot dimensions(x = rows, y = cols.")
	if(save){
		png(filename=paste(title, ".png", sep = ""))
		par(mfrow = c(x,y), oma = c(0, 0, 3, 0))
		dim <- x*y
		for (i in 1:dim) {
			plot(exprs[,i], pch=16, cex=0.3, main=colnames(exprs)[i])
		}
		mtext(title, outer = TRUE, cex = 1.5)
		dev.off()
	}
	else{
		par(mfrow = c(x,y), oma = c(0, 0, 3, 0))
		dim <- x*y
		for (i in 1:dim) {
			plot(exprs[,i], pch=16, cex=0.3, main=colnames(exprs)[i])
		}
		mtext(title, outer = TRUE, cex = 1.5)
	}
}

my.MAPlot <- function(obj, x, y, title="MA Plot", save=FALSE){
	if(missing(x) | missing(y))
		stop("Please input the plot dimensions(x = rows, y = cols.")
	
	if(save){
		png(filename=paste(title, ".png", sep = ""))
		par(mfrow = c(x,y), oma = c(0, 0, 3, 0))
		dim <- x*y
		for (i in 1:dim) {
			plotMA(object = obj, array = i)
		}
		mtext(title, outer = TRUE, cex = 1.5)
		dev.off()
	}
	else{
		par(mfrow = c(x,y), oma = c(0, 0, 3, 0))
		dim <- x*y
		for (i in 1:dim) {
			plotMA(object = obj, array = i)
		}
		mtext(title, outer = TRUE, cex = 1.5)
	}
}

my.volcanoplot <- function(fc, pvalue, title="Volcano Plot", x="log Fold Change", y="-log(Adjusted P.value)", upcol="springgreen1", downcol="firebrick1", save=FALSE){
	pvalue <- -log(pvalue)
	if(save){
		png(filename=paste(title, ".png", sep = ""))
		plot(x = fc, 
			y = pvalue,
		    ylab = y,
		    xlab = x,
		    main = title,
		    pch = 20,
		    col = "black")
		#Add lines to visualize more significant genes
		abline(v = 0)
		abline(h = 3)
		#Add green to up regulated genes and red to down regulated genes
		points(fc[(pvalue>3 & fc>0)],
	        pvalue[(pvalue>3 & fc>0)],
	        col = upcol,
	        pch = 20)
		points(fc[(pvalue>3 & fc<0)],
		    pvalue[(pvalue>3 & fc<0)],
		    col = downcol,
		    pch = 20)
		dev.off()
	}
	else{
		plot(x = fc, 
			y = pvalue,
		    ylab = y,
		    xlab = x,
		    main = title,
		    pch = 20,
		    col = "black")
		#Add lines to visualize more significant genes
		abline(v = 0)
		abline(h = 3)
		#Add green to up regulated genes and red to down regulated genes
		points(fc[(pvalue>3 & fc>0)],
	        pvalue[(pvalue>3 & fc>0)],
	        col = upcol,
	        pch = 20)
		points(fc[(pvalue>3 & fc<0)],
		    pvalue[(pvalue>3 & fc<0)],
		    col = downcol,
		    pch = 20)
	}
}

##########################################################################################
##########################################################################################
# 									   Other Functions 	 								 #
##########################################################################################
##########################################################################################

t.test.all.genes <- function(x,s1,s2){
  x1 <- x[s1]
  x2 <- x[s2]
  x1 <- as.numeric(x1)
  x2 <- as.numeric(x2)
  t.out <- t.test(x1,x2, 
                  alternative="two.sided",
                  var.equal=TRUE, 
                  paired = FALSE)
  out <- as.numeric(t.out$p.value)
  return(out)
}

# Classify Gene Regulation
evaluate.exprs <- function(data){
  up.down <- " "
  for (i in 1:length(data)) {
    if(data[i] < 0)
      up.down <- c(up.down, "Down")
    if(data[i] > 0)
      up.down <- c(up.down, "Up")
    if(data[i]==0)
      up.down <- c(up.down, "Equal")
  }
  return(up.down[-1])
}

# Write Results
write.results <- function(diff, filename){
	diff <- subset(diff, diff$P.Value<0.05)
	diff$Regulation <- evaluate.exprs(diff$logFC)
	df <- diff[, c(4,6,5,8,9,11,7,10,1,2,3)]
	write.table(diff, file=filename, sep="\t", quote=FALSE, row.names=FALSE)
}

#params: fn=function
extract_help <- function(pkg, fn = NULL, to = c("txt", "html", "latex", "ex"))
{
  to <- match.arg(to)
  rdbfile <- file.path(find.package(pkg), "help", pkg)
  rdb <- tools:::fetchRdDB(rdbfile, key = fn)
  convertor <- switch(to, 
      txt   = tools::Rd2txt, 
      html  = tools::Rd2HTML, 
      latex = tools::Rd2latex, 
      ex    = tools::Rd2ex
  )
  f <- function(x) capture.output(convertor(x))
  if(is.null(fn)) lapply(rdb, f) else f(rdb)
}


