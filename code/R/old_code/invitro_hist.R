
library(vegan)

# Read in data
m9 <- read.delim('/home/mjenior/Desktop/repos/Cdiff_modeling/data/transcript/Monk_et_al_2016/normalized.tsv', header=TRUE, row.names=1)
m9 <- m9[,1:3]
m9 <- as.vector(apply(m9, 1, median))
lb <- read.delim('/home/mjenior/Desktop/repos/Cdiff_modeling/data/transcript/SRR941894.mapped.norm.tsv', header=TRUE, row.names=1)
lb <- as.vector(lb$normalized_abundance)
lb <- round(lb * 10.0)

# Format equally
subsize <- round(min(c(sum(m9), sum(lb))) * 0.75)
m9 <- as.vector(rrarefy(m9, sample=subsize))
lb <- as.vector(rrarefy(lb, sample=subsize))
rm(subsize)

# Generate figures
png(filename='~/Desktop/repos/Cdiff_modeling/results/m9_hist.png', 
    units='in', width=4, height=3, res=300)
par(mar=c(3,3,0.5,1), las=1, mgp=c(2,0.7,0), xaxs='i', yaxs='i')
hist(m9, main='', xlim=c(0,200), ylim=c(0,2000), breaks=500, col='#4145ba',
     xlab='Transcript Density', ylab='Gene Frequency', cex.lab=1.1, cex.axis=0.7, lwd=2)
abline(v=as.vector(quantile(m9, probs=c(0.5,0.625,0.75,0.875))), lty=3, lwd=2, col='firebrick')
box(lwd=2)
text(x=c(10,26,50,110,180), y=1900, labels=c('1704','425','426','425','426'), cex=0.7)
dev.off()

png(filename='~/Desktop/repos/Cdiff_modeling/results/lb_hist.png', 
    units='in', width=4, height=3, res=300)
par(mar=c(3,3,0.5,1), las=1, mgp=c(2,0.7,0), xaxs='i', yaxs='i')
hist(lb, main='', xlim=c(0,300), ylim=c(0,1400), breaks=20, col='#ffa05d',
     xlab='Transcript Density', ylab='Gene Frequency', cex.lab=1.1, cex.axis=0.7, lwd=2)
abline(v=as.vector(quantile(lb, probs=c(0.5,0.625,0.75,0.875))), lty=3, lwd=2, col='firebrick')
box(lwd=2)
text(x=c(26,65,105,165,250), y=1350, labels=c('2070','517','518','517','518'), cex=0.7)
dev.off()

# Clean up
rm(list=ls())
gc()

