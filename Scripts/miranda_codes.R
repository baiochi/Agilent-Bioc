# Miranda Codes
setwd('/Users/Baiochi/Dropbox/USP/Lab/Analises/miRanda')

# ===============================================
#          DAVID vs Predicted
#
#   GIven a Pathway list, check predicted targets
# ===============================================






# ===========================================
#          Predicted miRNA
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
tempdir <- getwd()
setwd('/Users/Baiochi/Dropbox/USP/Lab/Analises/miRanda')
load(file='predictGood.rda')
load(file='predictBad.rda')
setwd(tempdir)
rm(tempdir)


# ===========================================
#          Cross miRNA data
#
#   Predicted miRs vs miRS found in analysis
# ===========================================

# Selete mir101 predicted targets
aux <- predictGood
grep(pattern='hsa-miR-101',unique(aux$mirna_name))
mir101G <- aux[which(aux$mirna_name=='hsa-miR-101'),]
aux <- predictBad
grep(pattern='hsa-miR-101',unique(aux$mirna_name))
mir101B <- aux[which(aux$mirna_name=='hsa-miR-101'),]
rm(aux)

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





