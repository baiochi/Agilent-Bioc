# Miranda Codes
setwd('/Users/Baiochi/Dropbox/USP/Lab/Analises/miRanda')

# ===============================================
#          DAVID vs Predicted
#
#   GIven a Pathway list, check predicted targets
# ===============================================


# ===========================================
#          Validated miRNA - miRTarBase
#
#   Read, save and load all predicted regions
# ===========================================

load(file='/Users/Baiochi/Dropbox/USP/Lab/Results/pre-miR101/RData/exprs.rda')
tarBaseMirs <- read.delim(file='/Users/Baiochi/Dropbox/USP/Lab/RawData/miRTarBase/hsa_MTI.txt',header=TRUE,sep='\t', stringsAsFactors = FALSE)
mir101 <- tarBaseMirs[grep(pattern = 'hsa-miR-101-', tarBaseMirs$miRNA),]

#unique genes found on array matching tar database
table(unique(exprs$GeneName) %in% mir101$Target.Gene)
#FALSE  TRUE 
#28916   174 

#all genes found on array matching tar database
table(exprs$GeneName %in% mir101$Target.Gene)
#FALSE  TRUE 
#32172   207 

#Result Table: GeneName, ProbeName, Systematic Name, p.val, adj.P.Val, FC, Start, Sequence, Description
genesFound <- exprs[which((exprs$GeneName %in% mir101$Target.Gene)),c(9,8,7,15,16,12,3,4,10)]
genesFound <- genesFound[with(genesFound, order(adj.P.Val, logFC)), ]
rownames(genesFound) <- NULL

table(genesFound$adj.P.Val<0.05)
#FALSE  TRUE 
#20    87 

tarUpGenes <- genesFound[which(genesFound$adj.P.Val<0.05 & genesFound$Regulation=='Up'),1]
tarDownGenes <- genesFound[which(genesFound$adj.P.Val<0.05 & genesFound$Regulation=='Down'),1]

write.table(genesFound , file='targets_mir101_TarBase.txt', sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(tarDownGenes , file='genelist_DOWN_mir101_TarBase.txt', sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(tarUpGenes , file='genelist_UP_mir101_TarBase.txt', sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

# ===========================================
#          Predicted miRNA - miRanda
#
#   Read, save and load all predicted regions
# ===========================================

# Read
predictGood <- read.table(file='human_predictions_good_conserved.txt', header=TRUE, sep="\t", stringsAsFactors=FALSE)
predictBad <- read.table(file='human_predictions_bad_conserved.txt', header=TRUE, sep="\t", stringsAsFactors=FALSE)
# Save
save(predictGood, file='predictGood.rda')
save(predictBad, file='predictBad.rda')
# Load
load(file='/Users/Baiochi/Dropbox/USP/Lab/RawData/miRanda/predictGood.rda')
load(file='/Users/Baiochi/Dropbox/USP/Lab/RawData/miRanda/predictBad.rda')



# ===========================================
#          Cross miRNA data
#
#   Predicted miRs vs miRS found in analysis
# ===========================================

# Selete mir101 predicted targets
aux <- predictGood
grep(pattern='hsa-miR-101',unique(aux$mirna_name))
mir101G <- aux[which(aux$mirna_name=='hsa-miR-101-'),]
aux <- predictBad
grep(pattern='hsa-miR-101',unique(aux$mirna_name))
mir101B <- aux[which(aux$mirna_name=='hsa-miR-101'),]
rm(aux)


#test
mirGene <- read.delim(file='/Users/Baiochi/Dropbox/USP/Lab/Results/pre-miR101/List/miR-101_Down_GeneName.txt', header=FALSE, quote='', stringsAsFactors=FALSE)
colnames(mirGene) <- 'GeneName'
mirSN <- read.delim(file='/Users/Baiochi/Dropbox/USP/Lab/Results/pre-miR101/List/miR-101_Down_SystematicName.txt', header=FALSE, quote='', stringsAsFactors=FALSE)
colnames(mirSN) <- 'SystematicName'
targetsGene <- subset(mir101G, mir101G$gene_symbol %in% mirGene$GeneName)
targetsSN <- subset(mir101G, mir101G$ext_transcript_id %in% mirSN$SystematicName)
dim(targetsGene)[1]
dim(targetsSN)[1]

#write
write.table(targetsSN, file='targets_mir101_miRanda_bySN.txt', sep="\t", quote=FALSE, row.names=FALSE)
write.table(targetsGene, file='targets_mir101_miRanda_byGene.txt', sep="\t", quote=FALSE, row.names=FALSE)
write.table(targetsSN$ext_transcript_id, file='mir101_miRanda_SN_list.txt', sep="\t", quote=FALSE, row.names=FALSE)
write.table(targetsGene$gene_symbol, file='mir101_miRanda_GeneName_list.txt', sep="\t", quote=FALSE, row.names=FALSE)

# Load down regulated miRs from analysis
microRNA <- read.delim(file='HCT.Down - SystematicName.txt', header=FALSE, quote='', stringsAsFactors=FALSE)
colnames(microRNA) <- 'SystematicName'
# sys 5164
# probe 6352
# gene 4924

# Find similar targets
targetsG <- merge(microRNA, mir101G, by.x='SystematicName', by.y='ext_transcript_id')
targetsB <- merge(microRNA, mir101B, by.x='SystematicName', by.y='ext_transcript_id')
dim(targetsG)[1]    #625
dim(targetsB)[1]    #1505

# Subset unique values
utargetsG <- unique(targetsG$SystematicName)    #526
utargetsB <- unique(targetsB$SystematicName)    #959

# Write data
write.table(targetsG$SystematicName, file='Down_vs_Predicted_Good.txt', sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(targetsB$SystematicName, file='Down_vs_Predicted_Bad.txt', sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(utargetsG, file='unique_Down_vs_Predicted_Good.txt', sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(utargetsB, file='unique_Down_vs_Predicted_Bad.txt', sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)






# ========================================
#          Read RefSeq
#
#   Read data with all RefSeqs
# ========================================

# Read file
i <- 1
x <- new('character')
conn <- file("sequence.txt", "r")
while(length(line <- readLines(conn, 1)) != 0) {
    x[i] <- line
    i = i + 1
    cat(paste(i,'\n'))
}
rm(i)

# Begin here
setwd("/Users/Baiochi/Desktop")
# Load object
load('allRefSqe.rda')

# Remove blank spaces
y <- subset(x, x!='')

# Subset data
id <- new('character')
seq <- new('character')
aux <- new('character')
k <- 1
l <- 0
for (i in 1:length(y)) {
    if(substring(y[i], 1, 1)=='>'){
        id[k] <- y[i]
        aux <- new('character')
        k = k+1
        l = l+1
    }
    else{
        seq[l] <- paste(aux, y[i], sep='')
        aux <- seq[l]
    }
    #cat(paste(i,'\n'))
}
rm(list=c('i','l','k', 'aux'))



















# ========================================
#          Title
#
#   Description
# ========================================





