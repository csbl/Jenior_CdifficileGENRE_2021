#------------Read in data--------------

library(vioplot)

# Metabolomes
metabolome <- read.delim('/home/mjenior/Desktop/repos/Cdiff_modeling/data/metabolome/scaled_intensities.log10.tsv', sep='\t', header=TRUE)

# Metadata
metadata <- read.delim('/home/mjenior/Desktop/repos/Cdiff_modeling/data/metadata.tsv', sep='\t', header=T, row.names=1)

# Growth rates
invivo_rates <- as.data.frame(read.delim(file='~/Desktop/repos/Cdiff_modeling/data/invivo_growth_rates.tsv', 
                                         sep='\t', header=FALSE))[,1][1:8304]

#------------Format data----------------

# Metadata
metadata$type <- NULL
metadata$cage <- NULL
metadata$mouse <- NULL
metadata$gender <- NULL
metadata$susceptibility <- NULL
metadata$clearance <- NULL

# Metabolomes
metabolome$BIOCHEMICAL <- gsub('_', ' ', metabolome$BIOCHEMICAL)
rownames(metabolome) <- metabolome$BIOCHEMICAL
metabolome$BIOCHEMICAL <- NULL
metabolome$PUBCHEM <- NULL
metabolome$KEGG <- NULL
metabolome$SUB_PATHWAY <- NULL
metabolome$SUPER_PATHWAY <- NULL
metabolome <- as.data.frame(t(metabolome))
colnames(metabolome) <- make.names(colnames(metabolome))

# Combine tables and subset data
metabolome <- merge(metadata, metabolome, by='row.names')
rownames(metabolome) <- metabolome$Row.names
metabolome$Row.names <- NULL
metabolome <- subset(metabolome, infection == 'mock')
metabolome <- subset(metabolome, abx %in% c('clindamycin','none'))
metabolome$infection <- NULL
clindamycin <- subset(metabolome, abx == 'clindamycin')
clindamycin$abx <- NULL
clindamycin <- clindamycin[, 'glucose.6.phosphate']
no_abx <- subset(metabolome, abx == 'none')
no_abx$abx <- NULL
no_abx <- no_abx[, 'glucose.6.phosphate']
rm(metabolome, metadata)

# Convert to numeric matrices
clindamycin <- as.numeric(clindamycin)
no_abx <- as.numeric(no_abx)

#------------Analyze data-----------------

# Test for significant differences
pval <- wilcox.test(clindamycin, no_abx, exact=FALSE)$p.value


invivo_pvals <- p.adjust(c(round(wilcox.test(invivo_rates, m9_gluc_aerobic, exact=FALSE)$p.value,5),
                           round(wilcox.test(invivo_rates, m9_gluc_anaerobic, exact=FALSE)$p.value,5),
                           round(wilcox.test(invivo_rates, lb_aerobic, exact=FALSE)$p.value,5)), method='BH')

#------------Generate figure--------------

pdf(file='/home/mjenior/Desktop/repos/Cdiff_modeling/results/figure_S4.pdf', width=6, height=4)
layout(matrix(c(1,2), nrow=1, ncol=2, byrow=TRUE))

par(mar=c(3,3.3,1,1), xpd=FALSE, las=1, mgp=c(2,0.75,0), lwd=2)
plot(0,0,type="n",xlim=c(0,1.2), ylim=c(0,1.5),  xaxt = 'n', xlab='', xaxt='n', 
     yaxt='n', ylab=expression(paste('Growth Rate (hr'^'-1',')')))
axis(2, at=c(0,0.5,1,1.5), labels=c('0.0','0.5','0.1','1.5'), lwd=2) 
vioplot(invivo_rates, at=0.5, col='firebrick', lwd=2, drawRect=FALSE, add=TRUE)
text(x=0.5, y=1.25, '***', font=2, cex=3)
mtext('in vivo\nBiomass Flux Samples', side=1, at=0.5, padj=1, cex=0.9)
# https://www.ncbi.nlm.nih.gov/pubmed/7592332?dopt=Abstract
segments(x0=1, x1=1.1, y0=0.67, lwd=3, col='firebrick2')
segments(x0=1.1, y0=0.67, y1=1.33, lwd=3, col='firebrick2')
segments(x0=1, x1=1.1, y0=1.33, lwd=3, col='firebrick2')
#text(x=0.9, y=1.33, labels='17', cex=0.8, col='firebrick2')
mtext('A',side=3, padj=0.5, cex=1.2, font=2, at=-0.4)

par(mar=c(3,4,1,1), las=1, mgp=c(2.6,1,0), lwd=2)
plot(1, type='n', xlab='', ylab='Relative Intesity (Log10)', axes=FALSE, xlim=c(0.75,3.25), ylim=c(0,6))
axis(side=2, at=c(0:6), labels=c('0.0','1.0','2.0','3.0','4.0','5.0','6.0'), cex.axis=0.9, lwd=2)
mtext(text=c('Clindamycin:','-','+'), side=1, at=c(0.2,1.5,2.5), cex=0.7)
mtext(text='Glucose-6-phosphate', side=1, at=2, padj=2)
abline(v=c(3.5,6.5,9.5,12.5))
segments(x0=1.5, x1=2.5, y0=5, lwd=3)
if (pval <= 0.05) {text(x=2, y=5.3, labels='*', cex=2)} else {text(x=2, y=5.1, labels='n.s.')}
box()
stripchart(no_abx, vertical=T, pch=21, xaxt='n', yaxt='n', bg='navy', at=1.5,
           cex=1.5, method='jitter', jitter=0.2, cex.lab=1.2, add=TRUE)
segments(1.1, x1=1.9, y0=median(no_abx), lwd=3)
stripchart(clindamycin, vertical=T, pch=21, xaxt='n', yaxt='n', bg='gold2', at=2.5,
           cex=1.5, method='jitter', jitter=0.2, cex.lab=1.2, add=TRUE)
segments(2.9, x1=2.1, y0=median(clindamycin), lwd=3)
mtext('B',side=3, padj=0.5, cex=1.2, font=2, at=-0.2)

dev.off()

#------------Clean up----------------

# Clean up
rm(list=ls())
gc()
