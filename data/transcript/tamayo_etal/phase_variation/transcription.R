
# Start with clean environment
rm(list=ls())
gc()

# Flux sampling files
rough1 <- '~/Desktop/tamayo_analysis/rough_1.tsv'
rough2 <- '~/Desktop/tamayo_analysis/rough_2.tsv'
smooth2 <- '~/Desktop/tamayo_analysis/smooth_2.tsv'
smooth3 <- '~/Desktop/tamayo_analysis/smooth_3.tsv'

# Read in data
rough1 <- read.delim(rough1, sep='\t', header=FALSE, row.names=1)
rough2 <- read.delim(rough2, sep='\t', header=FALSE, row.names=1)
smooth2 <- read.delim(smooth2, sep='\t', header=FALSE, row.names=1)
smooth3 <- read.delim(smooth3, sep='\t', header=FALSE, row.names=1)

# Merge data
rough <- merge(rough1, rough2, by='row.names')
rownames(rough) <- rough$Row.names
rough$Row.names <- NULL
colnames(rough) <- c('rough1','rough2')
smooth <- merge(smooth2, smooth3, by='row.names')
rownames(smooth) <- smooth$Row.names
smooth$Row.names <- NULL
colnames(smooth) <- c('smooth2','smooth3')
rm(rough1, rough2, smooth2, smooth3)

# Subsample reads
library(vegan)
sub_lvl <- round(min(c(sum(rough$rough1), sum(rough$rough2), sum(smooth$smooth2), sum(smooth$smooth3)))*0.9)
rough$rough1 <- as.vector(rrarefy(rough$rough1, sample=sub_lvl))
rough$rough2 <- as.vector(rrarefy(rough$rough2, sample=sub_lvl))
smooth$smooth2 <- as.vector(rrarefy(smooth$smooth2, sample=sub_lvl))
smooth$smooth3 <- as.vector(rrarefy(smooth$smooth3, sample=sub_lvl))
rm(sub_lvl)

# Find median within each group
library(matrixStats)
rough$rough_median <- rowMedians(as.matrix(rough))
smooth$smooth_median <- rowMedians(as.matrix(smooth))

# Merge data
transcription <- merge(rough, smooth, by='row.names')
rownames(transcription) <- transcription$Row.names
transcription$Row.names <- NULL
rm(rough, smooth)

# Test for correlations
cor.test(x=transcription$rough1, y=transcription$rough2, method='spearman', exact=FALSE)
png(filename='~/Desktop/tamayo_analysis/rough_transcript_corr.png', units='in', width=6, height=6, res=300)
par(mar=c(4,4,1,1), las=1, mgp=c(2.6,1,0), lwd=2, xpd=FALSE)
plot(x=transcription$rough1, y=transcription$rough2, pch=21, bg='firebrick',
     xlab='Rough 1', ylab='Rough 2', cex=2,
     xlim=c(0,7000), ylim=c(0,7000), cex.axis=0.8)
abline(lm(transcription$rough2 ~ transcription$rough1), lwd=3)
legend('bottomright', legend=c('R = 0.994','p << 0.001'), pt.cex=0, bty='n', cex=1.5)
dev.off()
# Highest = 645463.3.peg.1595, Cysteine synthase (EC 2.5.1.47) | cysM

cor.test(x=transcription$smooth2, y=transcription$smooth3, method='spearman', exact=FALSE)
png(filename='~/Desktop/tamayo_analysis/smooth_transcript_corr.png', units='in', width=6, height=6, res=300)
par(mar=c(4,4,1,1), las=1, mgp=c(2.6,1,0), lwd=2, xpd=FALSE)
plot(x=transcription$smooth2, y=transcription$smooth3, pch=21, bg='blue2',
     xlab='Smooth 2', ylab='Smooth 3', cex=2,
     xlim=c(0,7000), ylim=c(0,7000), cex.axis=0.8)
abline(lm(transcription$smooth3 ~ transcription$smooth2), lwd=3)
legend('bottomright', legend=c('R = 0.995','p << 0.001'), pt.cex=0, bty='n', cex=1.5)
dev.off()
# Highest = NAD-specific glutamate dehydrogenase (EC 1.4.1.2) | gluD

cor.test(x=transcription$rough_median, y=transcription$smooth_median, method='spearman', exact=FALSE)
png(filename='~/Desktop/tamayo_analysis/all_transcript_corr.png', units='in', width=6, height=6, res=300)
par(mar=c(4,4,1,1), las=1, mgp=c(2.6,1,0), lwd=2, xpd=FALSE)
plot(x=transcription$rough_median, y=transcription$smooth_median, pch=21, bg='darkorchid3',
     xlab='Rough Median', ylab='Smooth Median', cex=2,
     xlim=c(0,7000), ylim=c(0,7000), cex.axis=0.8)
abline(lm(transcription$smooth_median ~ transcription$rough_median), lwd=3)
legend('bottomright', legend=c('R = 0.98','p << 0.001'), pt.cex=0, bty='n', cex=1.5)
dev.off()

# Find largest differences
transcription$diff <- transcription$rough_median - transcription$smooth_median
transcription$abs_diff <- abs(transcription$rough_median - transcription$smooth_median)
transcription <- transcription[order(-transcription$abs_diff), ]
transcription <- subset(transcription, abs_diff > sd(transcription$abs_diff)*3) # >3 standard dev

# Add gene names
genes_names <- read.delim('~/Desktop/tamayo_analysis/r20291_patric_genes.tsv', sep='\t', header=TRUE)
transcription <- merge(transcription, genes_names, by.x='row.names', by.y='patric')
transcription$product <- gsub('_', ' ', transcription$product)
transcription <- transcription[order(transcription$abs_diff), ]
rownames(transcription) <- transcription$Row.names
transcription$Row.names <- NULL
rm(genes_names)

# Generate figure
png(filename='~/Desktop/tamayo_analysis/diff_transcript.png', units='in', width=8, height=5, res=300)
par(mar=c(4,24,1,1), las=1, mgp=c(2.6,1,0), lwd=2, xpd=FALSE)
barplot(transcription$diff, horiz=TRUE, xlim=c(-3000,3000), col='darkorchid3', cex.axis=0.75, 
        xlab='Differential Transcription', names.arg=transcription$product, cex.names=0.8)
abline(v=0, lty=3)
text(x=-2200, y=0.25, labels='Smooth')
text(x=2200, y=0.25, labels='Rough')
box()
dev.off()
