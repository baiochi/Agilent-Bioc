#Functions for DAVID analysis

chart = read.delim(file='/Users/Baiochi/Desktop/PM101/Annotation/DAVID/down chart.txt', header=TRUE, sep='\t', stringsAsFactors=FALSE)
chart = chart[,c(1,2,6)]
down.genes <- getGeneName(chart, exprs)
chart = read.delim(file='/Users/Baiochi/Desktop/PM101/Annotation/DAVID/up chart.txt', header=TRUE, sep='\t', stringsAsFactors=FALSE)
chart = chart[,c(1,2,6)]
up.genes <- getGeneName(chart, exprs)

all.up <- allGenes(up.genes)
all.down <- allGenes(down.genes)

write(all.up, file='Up - All pathway genes.txt', sep='\t')
write(all.down, file='Down - All pathway genes.txt', sep='\t')

writePathways(down.genes, path='Pathways - down regulated')
writePathways(up.genes, path = 'Pathways - up regulated')


allGenes <- function(list){
	genes <- new('character')
	for (i in 1:length(list)){
		genes <- c(genes, as.character(list[[i]]$GeneName))
	}
		
	return(unique(genes))
}


#Given a DAVID pathway by probe ID, retrieve the gene names
getGeneName <- function(chart, exprs){

	ids <- exprs[with(exprs, order(ProbeName)), ]
	rownames(ids) <- NULL
	l <- new('list')
	lost <- new('character')
	for (i in 1:dim(chart)[1]) {
		q <- sort(unlist(strsplit(chart$Genes[i],', '))) #sort genes in chart file, one pathway
		p <- ids[which(ids$ProbeName %in% q),]			 #eset with the genes containing in q
		cat(paste('q: ', length(q)))
		cat(paste('p: ', length(p$ProbeName),'\n'))
		if(length(q)>length(p$ProbeName)){
			lost <- paste(lost, ids$GeneName[which(!q %in% p$ProbeName)])
		} 
		l[[i]] <- data.frame(ProbeName=p$ProbeName,
			   				GeneName=p$GeneName, 
			    			adj.PVal=p$adj.P.Val,
			    			FC=p$logFC)
	}
	names(l) <- chart$Term
	if(length(lost)>0)
		cat(paste('Genes not found:'), unique(unlist(strsplit(lost,' ')))[-1])

	return(l)

}

writePathways <- function(genes, path='Pathway'){
	tempdir <- getwd()
	if(!file.exists(path))
		dir.create(path)
	setwd(path)
	for (i in 1:length(genes)) 
		write.xlsx(genes[[i]], file=paste(names(genes)[i],'.xlsx',sep=""), row.names=FALSE)
	setwd(tempdir)
}
