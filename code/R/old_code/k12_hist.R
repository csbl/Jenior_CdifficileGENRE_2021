
# Load in data
transcription <- read.delim('~/Desktop/repos/Cdiff_modeling/data/transcript/clinda_k12.mapped.norm.tsv', header=TRUE, row.names=1)

# Format data
transcription <- transcription$normDepth
k12_quantiles <- as.vector(quantile(transcription, c(0.5, 0.625, 0.75, 0.875)))

# Transcript abundance histogram
pdf(file='~/Desktop/repos/Cdiff_modeling/results/k12_invivo_transcript.pdf', width=3.5, height=5)

par(mar=c(4,4,1,1), las=1, mgp=c(2.5,1,0), xpd=FALSE, xaxs='i', yaxs='i')
hist(transcription, main='', xlim=c(0,400), ylim=c(0,800), breaks=200, col='dimgray',
     xlab='Transcript Density', ylab='Number of Genes', cex.lab=1.2, cex.axis=0.9, lwd=2)
abline(v=k12_quantiles, lty=5, lwd=2, col=c('firebrick4','firebrick3','firebrick2','firebrick1'))
legend('topright', legend=c('Q50', 'Q62.5','Q75','Q87.5'), pt.cex=0, lwd=2, lty=5,
       col=c('firebrick4','firebrick3','firebrick2','firebrick1'), bty='n')
box(lwd=2)

dev.off()

# Clean up
rm(list=ls())
gc()
