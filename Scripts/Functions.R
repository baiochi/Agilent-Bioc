# Bulti-in Functions for Microarray Analysis
# contact: baiochi@usp.br

# =================================================
#          Save Essential Plots
#
#   Hierclust, Boxplot, Density, PCA and MvsA Plots
# =================================================

cat('Loading functions for Agilent Microarray Analysis...')

RawPlots <- function(raw, x=2, y=3){
	cat('Plotting raw data.\n')
	# if(!is(raw ,"EList"))
	# 	stop('fit must be a \'MArrayLM\' ')

	lograw = log2(raw$E)

	hierclust(lograw, methdis="euclidean", methclu="complete", sel=FALSE, size=100, save=TRUE) 
	cor.matrix(lograw, title="Correlation Matrix", save=TRUE)
	Boxplot(lograw, title="Raw log Boxplot", color="forestgreen", save=TRUE)
	#pca.plot(raw, title="Sample PCA plot", save=TRUE)
	densityplot(lograw, dim(lograw)[2], title="Raw log Density", save=TRUE)
	MAPlot(raw, x, y, title = "Raw MA Plot", save=TRUE)
}

NormPlots <- function(eset, sn=1, rep=1){
	cat('Plotting processed data.\n')
	# if(!is(eset ,"EList"))
	# 	stop('fit must be a \'MArrayLM\' ')
	Boxplot(eset$E, title="Normalized Boxplot", color="deepskyblue3", save=TRUE)
	if(anyNA(eset$E)){
		cat('\nExpression-set containg NA-values, unable to plot densities and MAPlot.\n')
	}
	else{
		densityplot(eset$E, dim(eset$E)[2], title="Normalized Density", save=TRUE)
		MAPlot(eset, sn, rep, title = "Post-Normalization MvsA Plot", save=TRUE)
	}
}

#Diff Exprs Plots
DiffExprsPlots <- function(fit, coefs, rank, eset){
	cat('Plotting statistical results.\n')
	# if(!is(eset ,"EList"))
	# 	stop('fit must be a \'MArrayLM\' ')
	# if(!is(fit ,"MArrayLM"))
	# 	stop('eset must be a \'EList\' ')
	# result = decideTests(fit, method="separate", adjust.method="fdr", p.value=0.05, lfc=0)
	for (i in 1:length(coefs)) {
		exprs = topTable(fit, adjust="fdr", coef=coefs[i], genelist=eset$genes, number=Inf)
		pval = exprs$P.Value
		fc = exprs$logFC
		adj.pval = exprs$adj.P.Val
		aveexprs = exprs$AveExpr

		#Plotting
		histogram(pval, title = paste(coefs[i],' P-Value Histogram'), save=TRUE, col="yellow", xlab='p-value')
		histogram(fc, title = paste(coefs[i], ' Fold Change Histogram'), save=TRUE, col="firebrick2", xlab='log Fold Change')
		histogram(adj.pval, title = paste(coefs[i], ' Adjusted P-value Histogram'), save=TRUE, col="royalblue", xlab='adjusted p-value')
		histogram(aveexprs, title = paste(coefs[i], ' Average Expression Histogram'), save=TRUE, col="mediumorchid4", xlab='average expression')
		densityplot(exprs[,15:16],2, title=paste(coefs[i], ' P-values density'), save=TRUE)
		volcanoplot(fc, adj.pval, rank, title=paste(coefs[i], ' Volcanoplot'), save=TRUE)
	}
}


saveAllPlots <- function(raw, eset, diffexprs, result){

	if(!sum(dir()=='Plots')){
		dir.create('Plots')
	}
	tempdir <- getwd()
	setwd('Plots')
	#set vars
	lograw <- log2(raw)
	exprs <- eset$E
	pval <- diffexprs$P.Value
	fc <- diffexprs$logFC
	adj.pval <- diffexprs$adj.P.Val
	aveexpr <- diffexprs$AveExpr
	size <- dim(exprs)[2]

	# Raw outliers
	hierclust(lograw, methdis="euclidean", methclu="complete", sel=FALSE, size=100, save=TRUE) 
	cor.matrix(lograw, title="Correlation Matrix", save=TRUE)
	Boxplot(lograw, title="Raw log Boxplot", color="forestgreen", save=TRUE)
	pca.plot(raw, title="Sample PCA plot", save=TRUE)
	densityplot(lograw, size, title="Raw log Density", save=TRUE)
	MAPlot(raw, 2, 3, title = "Raw MvsA Plot", save=TRUE)

	# Normalized
	Boxplot(exprs, title="Normalized Boxplot", color="deepskyblue3", save=TRUE)
	densityplot(exprs, size, title="Normalized Density", save=TRUE)
	MAPlot(eset, 2, 3, title = "Post-Normalization MvsA Plot", save=TRUE)

	# DiffExprs

	histogram(pval, title = 'P-Value Histogram', save=TRUE, col="yellow", xlab='p-value')
	histogram(fc, title = 'Fold Change Histogram', save=TRUE, col="firebrick2", xlab='log Fold Change')
	histogram(adj.pval, title = 'Adjusted P-value Histogram', save=TRUE, col="royalblue", xlab='adjusted p-value')
	histogram(aveexpr, title = 'Average Expression Histogram', save=TRUE, col="mediumorchid4", xlab='average expression')
	densityplot(diffexprs[,14:15],2, title='P-values density', save=TRUE)
	volcanoplot(fc, adj.pval, result,
				title="Volcano Plot\nControl vc PM101", y="-log(Adj.P.Val)",
				save=TRUE)

	setwd(tempdir)
}


Boxplot <- function(exprs, title="Boxplot", color="white", save=FALSE){
  if(save){
    png(filename= paste(title, ".png", sep = ""))
    Boxplot <- function(exprs, title, color, save=FALSE)
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

densityplot <- function(exprs, samples, title="Density", save=FALSE){
	colors <- c("green", "red", "yellow", "blue", "purple", "orange",
			"green3", "red3", "yellow3", "blue3", "purple3", "orange3")
	if(save){
		png(filename= paste(title, ".png", sep = ""))
		densityplot(exprs, samples, title, save=FALSE)
		dev.off()
	}
	else{
		plot(density(exprs[,1]), col=colors[1], main=title)
		for (i in 1:samples) {
			lines(density(exprs[,i]), col=colors[i])
		}
		legend("topright",
			legend = colnames(exprs),
			lty=c(1,1),
			lwd=c(2.5,2.5),col=colors)
	}
}

distribuitionplot <- function(exprs, x, y, title="Distribuition", save=FALSE){
	if(missing(x) | missing(y))
		stop("Please input the plot dimensions(x = rows, y = cols.")
	if(save){
		png(filename=paste(title, ".png", sep = ""))
		distribuitionplot(exprs, x, y, title, save=FALSE)
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

MAPlot <- function(obj, x, y, title="MA Plot", save=FALSE){
  if(missing(x) | missing(y))
    stop("Please input the plot dimensions(x = rows, y = cols.")
	if(save){
		png(filename=paste(title, ".png", sep = ""))
		MAplot(obj, x, y, title, save=FALSE)
		dev.off()
	}
	else{
		par(mfrow = c(x,y), oma = c(0, 0, 3, 0))
		dim <- x*y
		for (i in 1:dim) {
			limma::plotMA(object = obj, array = i)
		}
		mtext(title, outer = TRUE, cex = 1.5)
	}
}

histogram <- function (dat, title='Histogram' ,save=FALSE, ...){
	if(save){
		png(filename=paste(title, ".png", sep = ""))
		hist(dat, main=title, ...)
		dev.off()
	}
	else{
		hist(dat, main=title, ...)
	}
}

diffexprs.hist <- function(exprs, coef='', save=FALSE){
  if(save){
    png(filename=paste(title, ".png", sep = ""))
    diffexprs.hist(exprs, coef, save=FALSE)
    dev.off()
  }
  par(mfrow=c(2,2), oma = c(0, 0, 2, 0))
  histogram(exprs$P.Value, title = 'P-Value Histogram', col="yellow", xlab='p-value')
  histogram(exprs$adj.P.Val, title = 'Adjusted P-value Histogram', col="royalblue", xlab='Adjusted p-value')
  histogram(exprs$logFC, title = ' Fold Change Histogram', col="firebrick2", xlab='log Fold Change', axes=FALSE)
  axis(side=1, at=seq(round(min(exprs$logFC)),round(max(exprs$logFC)), by=2))
  densityplot(cbind(P.value=exprs$P.Value, Adjusted.P.value=exprs$adj.P.Val),2,title='P-values density')
  mtext(coef, outer = TRUE)
}


volcanoplot <- function(fit, coeff=NULL, hl =NULL, pval=0.05, fc=0, title="Volcano Plot", 
                        x="log2 Fold Change", y="-log Adjusted P.value", 
                        upcol="springgreen1", downcol="firebrick1", ab=FALSE, save=FALSE){
  #pvalue <- -log(pvalue)
  #plot(x = fc, y = pvalue,ylab = y,xlab = x,main = title,pch = 20,col = "black")
  #points(fc[(pvalue>3 & fc>0.7)],pvalue[(pvalue>3 & fc>0.7)],col = upcol,pch = 20)
  #points(fc[(pvalue>3 & fc< -0.7)],pvalue[(pvalue>3 & fc< -0.7)],col = downcol,pch = 20)
  if(save){
    png(file=paste(title, ".png", sep = ""))
    volcanoplot(fit, coeff, hl, pval, fc, title, x, y, upcol, downcol, ab, save=FALSE)
    dev.off()
  }
  else{
    exprs = topTable(fit, adjust="fdr", coef=coeff, number=Inf)
    diff = topTable(fit, adjust="fdr", coef=coeff, p.value=pval, lfc=fc,number=Inf)
    rank = decideTests(fit, method="separate", adjust.method="fdr", p.value=pval, lfc=fc)
    rank = rank[,colnames(rank)==coeff]
    rank = table(rank)
    Pval = -log(exprs$adj.P.Val)
    FC = exprs$logFC
    pval = -log(pval)
    plot(x = FC, y = Pval,
         ylab = y, xlab = x, main = title, pch = 20, col = "black", 
         cex=.3, #xlim=(c(-6, 4)), ylim=(c(0, 15)), 
         cex.axis=0.9, cex.lab=0.9)
    if(ab){
      #Add lines to visualize more significant genes
      abline(v = log(fc))
      abline(v = -log(fc))
      abline(h = pval)
    }
    #Add green to up regulated genes and red to down regulated genes
    points(FC[(Pval>pval & FC>fc)],
           Pval[(Pval>pval & FC>fc)],
           col = upcol,pch = 20, cex=.3)
    points(FC[(Pval>pval & FC< -fc)],
           Pval[(Pval>pval & FC< -fc)],
           col = downcol,pch = 20, cex=.3)
    legend("bottomright",
           legend = c(paste('Up-regulated:',rank[3], sep=' '),
                      paste('Down-regulated:',rank[1], sep=' ') ),
           lty= 0,# lwd=c(2.5,2.5),
           pch = 20, cex=0.75, col=c('springgreen1','firebrick1'))
    if(!is.null(hl))
      highlight(hl, diff)
  }
}

#highlith transcripts from Volcanoplot
highlight <- function(genes, exprs, by.gene=TRUE){
  if(!by.gene)
    listname <- exprs$SystematicName
  else 
    listname <- exprs$GeneName
  for(i in 1:length(genes)){
    g <- exprs[grep(paste('^',genes[i],'$',sep=''), listname),]
    if(nrow(g)==0 | nrow(g)>1){
      
      cat(paste(genes[i],'transcript has duplicated instances, or none.\n'))
    }
    else{
      if(g$logFC>0) p=4
      else  p=2
      text(x=g$logFC, y=-log(g$adj.P.Val), label=g$GeneName, pos=p, offset=0.2, cex=.75)
      points(x=g$logFC, y=-log(g$adj.P.Val), pch=20, cex=.5, col='deepskyblue2')
    }
  }
}



Heatmap <- function (object, size, maintitle='HeatMap') {
  require(marray)
  require(gplots)
  require(gtools)
  require(gdata)
  if(missing(size)) 
    size = dim(object)[1]
  samples = colnames(object)
  names = rownames(object)
  genes.var = apply(object, 1, var, na.rm = TRUE)
  genes.var.select = order(genes.var, decreasing = TRUE)[1:size]
  DD.s = object[genes.var.select, ]
  samples = colnames(DD.s)
  names = rownames(DD.s)
  c <- rainbow(ncol(DD.s), start = 0, end = 0.3)
  rc <- rainbow(nrow(DD.s), start = 0, end = 0.3)
  rbg = maPalette(low = "green", high = "red", mid = "black", 
                  k = 50)
  heatmap.2(DD.s, labCol = samples, labRow = names, scale = "none", 
            col = rbg, margin = c(10, 10), tracecol = "cyan")
  subtitle = paste(as.character(size), " high variance genes")
  title(main = maintitle, sub = subtitle)
}

# function (fit, coef = 1, highlight = 0, names = fit$genes$ID, 
#     xlab = "Log Fold Change", ylab = "Log Odds", pch = 16, cex = 0.35, 
#     ...) 
# {
#     if (!is(fit, "MArrayLM")) 
#         stop("fit must be an MArrayLM")
#     if (is.null(fit$lods)) 
#         stop("No B-statistics found, perhaps eBayes() not yet run")
#     x <- as.matrix(fit$coef)[, coef]
#     y <- as.matrix(fit$lods)[, coef]
#     plot(x, y, xlab = xlab, ylab = ylab, pch = pch, cex = cex, 
#         ...)
#     if (highlight > 0) {
#         if (is.null(names)) 
#             names <- 1:length(x)
#         names <- as.character(names)
#         o <- order(y, decreasing = TRUE)
#         i <- o[1:highlight]
#         text(x[i], y[i], labels = substring(names[i], 1, 8), 
#             cex = 0.8, col = "blue")
#     }
#     invisible()
# }


# Cluster Dendogram
cluster.tree <- function(dat, title="Cluster Dendogram", save=FALSE, dist='euclidean', clust='single'){
	dat <- t(dat) #transpose dat
	dat.dist <- dist(dat, method=dist) # calculate distance
	dat.clust <- hclust(dat.dist, method=clust) # calculate clusters
	if(save){
		png(filename=paste(title, ".png", sep = ""))
		plot(dat.clust, main=title, labels=names(dat), cex=0.75) # plot cluster tree
		dev.off()
	}
	else{
		plot(dat.clust, main=title, labels=names(dat), cex=0.75) # plot cluster tree
	}
}

# Pearson's Correlation Matrix
cor.matrix <- function(dat, title="Correlation Matrix", save=FALSE, ...){
  dat.cor <- cor(dat)
  if(save){
    png(filename=paste(title, ".png", sep = ""))
    cor.matrix(dat, title, save=FALSE, ...)
    dev.off()
  }
  else{
    par(oma = c(2, 2, 2, 0))
    image(dat.cor,axes=F, ...)
    axis(2,at=seq(0,1,length=ncol(dat.cor)),label=dimnames(dat.cor)[[2]],las=2)
    axis(1,at=seq(0,1,length=ncol(dat.cor)),label=dimnames(dat.cor)[[2]], las=2)
    mtext(title, outer = TRUE, cex = 1.5)
  }
}

# PCA Plot
pca.plot <- function(dat, title="Sample PCA plot", save=FALSE){
	dat.pca <- prcomp(t(dat))
	dat.loads <- dat.pca$x[,1:2]
	if(save){
		png(filename=paste(title, ".png", sep = ""))
		plot(dat.loads[,1], dat.loads[,2], main=title, xlab="p1", ylab="p2", col='red', cex=1.5, pch=16, xlim=c(min(dat.loads),max(dat.loads)+150000))
		text(dat.loads, labels = colnames(dat), pos = 4)
		dev.off()
	}
	else{
		plot(dat.loads[,1], dat.loads[,2], main=title, xlab="p1", ylab="p2", col='red', cex=1.5, pch=16, xlim=c(min(dat.loads),max(dat.loads)+150000))
		text(dat.loads, labels = colnames(dat), pos = 4)
	}


	# png(filename="pca.png",1280,800)
	# k <- 1
	# par(mfrow = c(3,3), oma = c(0, 0, 3, 0))
	# for (i in 1:3) {
	# 	for (j in 1:3) {
	# 		dat.loads <- dat.pca$x[,c(i,j)]
	# 		plot(dat.loads[,1], dat.loads[,2], col=k, cex=1.5, pch=16, xlim=c(min(dat.loads),max(dat.loads)+150000))
	# 		text(dat.loads, labels = colnames(dat), pos = 4)
	# 	}
	# 	k <- k+1
	# }
	# dev.off()
}

# Coefficient of Variation Plot
cv.plot <- function(dat){
	dat.mean <- apply(dat,1,mean) # calculate mean for each gene
	dat.sd <- sqrt(apply(dat,1,var)) # calculate st.deviation for each gene
	dat.cv <- dat.sd/dat.mean #calculate cv
	plot(dat.mean,dat.cv,main="Sample CV vs. Mean",xlab="Mean",ylab="CV",col='blue',cex=1.5)
}

# # K-Means Clustering
# kmean.clust <- function(dat){
#   dd <- dat[names(f.p)[f.p<0.001],]
#   d.k <- kmeans(mydata,6)
#   par(mfrow=c(3,3))
#   for(i in 1:9) {
#   	tmp <- scale(mydata[d.k$cluster==i,])
#   	matplot(c(1:ncol(mydata)),t(tmp),type='l',col=i,xlab='Time',ylab='Expression')
#   }
# }

# methclu = c("complete","median", "centroid", "average")

# hierclust(raw$E, methdis="euclidean", methclu="complete", sel=FALSE, size=100)

hierclust <- function (object, methdis, methclu, sel, size, save=FALSE) 
{
    samples = colnames(object)
    if(missing(size))
      size = dim(object)[1]
    if (sel == "TRUE") {
        genes.var = apply(object, 1, var)
        genes.var.select = order(genes.var, decreasing = TRUE)[1:size]
        object = object[genes.var.select, ]
    }
    if (methdis == "euclidean") {
        d <- dist(t(object), method = methdis)
        title <- "Hierarchical Clustering - Euclidean"
    }
    if (methdis == "pearson") {
        d <- as.dist(1 - cor(object, use = "complete.obs", method = methdis))
        title <- "Hierarchical Clustering - Pearson"
    }
    if (methdis == "spearman") {
      d <- as.dist(1 - cor(object, use = "complete.obs", method = methdis))
      title <- "Hierarchical Clustering - Spearman"
    }
    if (sel == "TRUE") {
        title <- paste(title, "\nhigh variance genes")
    }
    else {
        title <- paste(title, "\nall genes")
    }
    dim <- dim(as.matrix(d))
    hc <- hclust(d, method = methclu)
    if(save){
    	png(filename=paste(title, ".png", sep = ""))
    	plot(hc, labels = samples, main = title)
    	dev.off()
	  }
	  else{
		  plot(hc, labels = samples, main = title)
	  }
    
}

# ==============================================
#          Read line functions
#
#   user inputs for dinamic analysis
# ==============================================


readBgc <- function(){ 
  n <- readline("Background correction method: (normexp, subtract, none)")
  n <- as.character(n)
  if (!n %in% c('normexp','subtract','none')){
  	cat('Invalid method, please input again: ')
    n <- readBgc()
  }
  return(n)
}

readNorm <- function(){ 
  n <- readline("Normalization method: (quantile, cyclicloess)")
  n <- as.character(n)
  if (!n %in% c('quantile', 'cyclicloess')){
  	cat('Invalid method, please input again: ')
    n <- readNorm()
  }
  return(n)
}

readCoef <- function(design){
	coefs <- paste(colnames(design), collapse=', ')
	n <- readline(paste("Coefficient of interest: ", coefs))
	n <- as.character(n)
	if (!n %in% colnames(design)){
		cat('Invalid coefficient, please input again: ')
	n <- readCoef()
	}
	return(n)
}

ask <- function(){
	n <- readline()
	n <- as.character(n)
	if (!n %in% c('y','n')){
		n <- ask()
	}
	return(n)
}

ask.dir <- function(){
	n <- readline('Enter directory: ')
	n <- as.character(n)
	if(!sum(dir()==n)){
		dir.create(n)
	}
	return(n)
}



##########################################################################################
# 									   Other Functions 	 								 #
##########################################################################################

# ==============================================
#          Ipak Function
#
#   install and load multiple R packages
# ==============================================
 
ipak <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
        install.packages(new.pkg, dependencies = TRUE)
    sapply(pkg, require, character.only = TRUE)
}

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

WriteResults <- function(fit, coefs, pval = 0.05, log=0, eset, filename, lists=FALSE){
  
	comp <- topTable(fit, adjust="fdr", coef=coefs, genelist=eset$genes, number=Inf, sort.by='logFC')
	exprs <- topTable(fit, adjust="fdr", coef=coefs, p.value = pval, lfc=log, genelist=eset$genes, number=Inf, sort.by='logFC')
	rownames(exprs) <- NULL
	cat('Ja foi')
	exprs$Regulation <- evaluate.exprs(exprs$logFC)
	cat('batata')
	up = subset(exprs, exprs$Regulation=='Up')
	up = up[with(up, order(adj.P.Val)), ]
	down = subset(exprs, exprs$Regulation=='Down')
	down = down[with(down, order(adj.P.Val)), ]
	up <- up[,c(7,8,9,15,16,12,4,10)]
	down <- down[,c(7,8,9,15,12,16,4,10)]
  
	if(lists){
	  #Up regulated probes
	  write.table(up$ProbeName, file=paste(filename,'Up_ProbeName.txt',sep='_'), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
	  write.table(unique(up$GeneName), file=paste(filename, 'Up_GeneName.txt',sep='_'), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
	  write.table(unique(up$SystematicName), file=paste(filename, 'Up_SystematicName.txt',sep='_'), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
	  
	  #Down regulated probes
	  write.table(down$ProbeName, file=paste(filename,'Down_ProbeName.txt',sep='_'), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
	  write.table(unique(down$GeneName), file=paste(filename, 'Down_GeneName.txt',sep='_'), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
	  write.table(unique(down$SystematicName), file=paste(filename, 'Down_SystematicName.txt',sep='_'), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
	}

	#Up, Down and All
	write.table(up, file=paste(filename, 'Up_Regulated.txt',sep='_') , sep="\t", quote=FALSE, row.names=FALSE)
	write.table(down, file=paste(filename, 'Down_Regulated.txt',sep='_'), sep="\t", quote=FALSE, row.names=FALSE)
	write.table(exprs[,c(7,8,9,15,16,12,4,10)], file=paste(filename, 'All_genes.txt',sep='_'), sep="\t", quote=FALSE, row.names=FALSE)

	#Complete
	write.table(comp, file=paste(filename, 'Complete_Statistical_Results.txt',sep='_'), sep="\t", quote=FALSE, row.names=FALSE)
	
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


saveData <- function(path, raw, bgc, norm, eset, design, fit2, exprs,rank){
	tempdir <- getwd()
	if(!file.exists(path))
		stop('Invalid path\n.')
	setwd(path)
	save(raw, file='raw.rda')
	save(bgc, file='bgc.rda')
	save(norm, file='norm.rda')
	save(eset, file='eset.rda')
	save(design, file='design.rda')
	save(fit2, file='fit2.rda')
	save(exprs, file='exprs.rda')
	save(rank, file='rank.rda')
	setwd(tempdir)
}

loadData <- function(path){
	if(!file.exists(path))
		stop('Invalid path\n.')
	load(file=paste(path,'raw.rda',sep='/'))
	load(file=paste(path,'bgc.rda',sep='/'))
	load(file=paste(path,'norm.rda',sep='/'))
	load(file=paste(path,'eset.rda',sep='/'))
	load(file=paste(path,'design.rda',sep='/'))
	load(file=paste(path,'fit2.rda',sep='/'))
	load(file=paste(path,'exprs.rda',sep='/'))
	load(file=paste(path,'rank.rda',sep='/'))

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


#Create Fasta Format
createFasta <- function(sequence, identifier){

	#odd which(length(sequence) %% 2 == 1)

	file <- new('character', '')
	k <- 1 	#constant index
	i <- 1 	#id index
	j <- 2 	#seq index
	for (l in 1:length(sequence) ) {
		file[i] <- paste('>',identifier[k], sep='')		#id
		file[j] <- sequence[k]						#sequence
		k <- k + 1
		i <- i + 2
		j <- j + 2
	}
	file
}

createFasta2 <- function(sequence, sn, probe, gene, dsc){

	#odd which(length(sequence) %% 2 == 1)

	file <- new('character', '')
	k <- 1 	#constant index
	i <- 1 	#id index
	j <- 2 	#seq index
	for (l in 1:length(sequence) ) {
		file[i] <- paste('>', sn[k], '|', probe[k], '|', gene[k], '| ', dsc[k], sep='')		#id
		file[j] <- sequence[k]						#sequence
		k <- k + 1
		i <- i + 2
		j <- j + 2
	}
	file
}

getFasta <- function(exprs){
	#arrange data
	seq <- exprs$Sequence
	sn <- exprs$SystematicName
	probe <- exprs$ProbeName
	gene <- exprs$GeneName

	file <- new('character', '')
	k <- 1 	#constant index
	i <- 1 	#id index
	j <- 2 	#seq index
	for (l in 1:length(seq) ) {
		file[i] <- paste('>', sn[k], '|', probe[k], '|', gene[k], sep='')		#id
		file[j] <- seq[k]						#sequence
		k <- k + 1
		i <- i + 2
		j <- j + 2
	}
	file
}

createTargetFile <- function(filename,conditions, replicates){
	targets<- data.frame(FileName=filename,
						Condition=conditions,
						Sample=)
	targets <- cbind(targets,Nomeclature=paste(targets$Condition,targets$Sample,sep="_"))
	return(targets)
}

readAFE <- function(targets, path, value='PS',skip.lines=0){
  
  #Read Files
  tempdir <- getwd()
  setwd(path)
  data <- new('list')
  for (i in 1:length(targets$FileName)) {
    cat(paste('Reading', targets$FileName[i],'\n'))
    data[[i]] <- read.delim(file=targets$FileName[i], header=TRUE, skip=skip.lines, stringsAsFactors=FALSE)
  }
  setwd(tempdir)
  
  #Extract ProcessedSignal and Flags
  n <- dim(data[[1]])[1]
  signal <- data.frame(numeric(n))
  if(value=='PS'){
    for (i in 1:length(targets$FileName)) signal <- cbind(signal, data[[i]]$gProcessedSignal)
  }
  if(value=='TGS'){
    for (i in 1:length(targets$FileName)) signal <- cbind(signal, data[[i]]$gTotalGeneSignal)
  }
  signal <- signal[,-1]
  signal <- as.matrix(signal)
  flags <- data[[1]][,37:40]
  
  #Read using limma (easier to get the genes dataset)
  raw = read.maimages(targets$FileName, source="agilent.mean", green.only=TRUE, path=path)
  
  #Adjust values
  raw$mean <- raw$E
  raw$flags <- flags
  raw$E <- signal
  colnames(raw$E) <- targets$Nomeclature
  colnames(raw$Eb) <- targets$Nomeclature
  colnames(raw$mean) <- targets$Nomeclature
  
  return(raw)
}










