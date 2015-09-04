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




#-------------miRTarBase---------------#

tarbase <- read.delim(file='Dropbox/USP/Lab/RawData/miRTarBase/hsa_MTI.txt', sep='\t',header=TRUE,stringsAsFactors=FALSE)
tar_mir101 <- tarbase[grep('hsa-miR-101-',tarbase$miRNA),]

#unique genes found on array matching tar database
table(unique(exprs$GeneName) %in% tar_mir101$Target.Gene)
FALSE  TRUE 
28916   174 

#all genes found on array matching tar database
table(exprs$GeneName %in% tar_mir101$Target.Gene)
FALSE  TRUE 
32172   207 

#Result Table: GeneName, ProbeName, Systematic Name, adj.P.Val, FC, Start, Sequence
exprs$Regulation <- evaluate.exprs(exprs$logFC)
genesFound <- exprs[which((exprs$GeneName %in% tar_mir101$Target.Gene)),c(8,7,9,15,16,18,3,4, 10)]
genesFound <- genesFound[with(genesFound, order(adj.P.Val, Regulation)), ]
rownames(genesFound) <- NULL

table(genesFound$adj.P.Val<0.05)
FALSE  TRUE 
  120    87 

tarUpGenes <- genesFound[which(genesFound$adj.P.Val<0.05 & genesFound$Regulation=='Up'),1]
tarDownGenes <- genesFound[which(genesFound$adj.P.Val<0.05 & genesFound$Regulation=='Down'),1]

write.table(genesFound , file='targets_mir101_TarBase.txt', sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(tarDownGenes , file='genelist_DOWN_mir101_TarBase.txt', sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(tarUpGenes , file='genelist_UP_mir101_TarBase.txt', sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)



