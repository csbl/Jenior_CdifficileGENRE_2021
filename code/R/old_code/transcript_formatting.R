
source('~/Desktop/repos/Cdiff_modeling/src/init.R')

cef <- read.delim('data/transcript/cefoperazone_630.mapped.norm.tsv', sep='\t', header=TRUE)
cef$contigLen <- NULL
clinda <- read.delim('data/transcript/clindamycin_630.mapped.norm.tsv', sep='\t', header=TRUE)
clinda$contigLen <- NULL
strep <- read.delim('data/transcript/streptomycin_630.mapped.norm.tsv', sep='\t', header=TRUE)
strep$contigLen <- NULL
gnoto <- read.delim('data/transcript/gnotobiotic_630.mapped.norm.tsv', sep='\t', header=TRUE)
gnoto$contigLen <- NULL

# Merge tables
transcript <- merge(x=cef, y=clinda, by='contigName')
transcript <- merge(x=transcript, y=strep, by='contigName')
transcript <- merge(x=transcript, y=gnoto, by='contigName')
colnames(transcript) <- c('gene','cefoperazone','clindamycin','streptomycin','gnotobiotic')
rm(cef, clinda, strep, gnoto)

# Remove rRNA genes (dude to rRNA depletion, their abundances are not reliable)
rRNA_test <- as.vector(strsplit(as.character(transcript$gene), 'ribosomal'))
x <- 1
omit <- c()
for (index in rRNA_test) {
  if (length(index) > 1) {
    omit <- c(omit, x)
  }
  x <- x + 1
}
transcript <- transcript[-omit,]
rm(omit, x, index, rRNA_test)

# Reformat values to integers
transcript$cefoperazone <- as.integer(round(transcript$cefoperazone))
transcript$clindamycin <- as.integer(round(transcript$clindamycin))
transcript$streptomycin <- as.integer(round(transcript$streptomycin))
transcript$gnotobiotic <- as.integer(round(transcript$gnotobiotic))

# Calculate optimal subsample size
subSize <- round(as.numeric(min(c(sum(transcript$cefoperazone), sum(transcript$clindamycin), 
                                  sum(transcript$streptomycin), sum(transcript$gnotobiotic)))) * 0.90)

# Subsample table
transcript$cefoperazone <- as.vector(rrarefy(transcript$cefoperazone, sample=subSize))
transcript$clindamycin <- as.vector(rrarefy(transcript$clindamycin, sample=subSize))
transcript$streptomycin <- as.vector(rrarefy(transcript$streptomycin, sample=subSize))
transcript$gnotobiotic <- as.vector(rrarefy(transcript$gnotobiotic, sample=subSize))

# Write to a new file
write.table(transcript, file='data/transcript/cdf_transcription.sub.tsv', sep='\t', row.names=FALSE, quote=FALSE)
rm(transcript, subSize)


