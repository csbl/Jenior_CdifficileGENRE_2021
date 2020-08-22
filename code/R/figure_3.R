
# Read in data
rough_fluxes <- read.delim('~/Desktop/repos/Jenior_Cdifficile_2019/data/transcript/tamayo_etal/phase_variation/riptide_rough/flux_samples.tsv', sep='\t', header=TRUE)
smooth_fluxes <- read.delim('~/Desktop/repos/Jenior_Cdifficile_2019/data/transcript/tamayo_etal/phase_variation/riptide_smooth/flux_samples.tsv', sep='\t', header=TRUE)

# Get reactions of interest
min_range <- c()
# input
smooth_glucose <- as.vector(smooth_fluxes[,'EX_cpd00027_e'])
min_range <- c(min_range, length(smooth_glucose))
smooth_glytyr <- as.vector(smooth_fluxes[,'EX_cpd15606_e'])
min_range <- c(min_range, length(smooth_glytyr))
rough_glytyr <- as.vector(rough_fluxes[,'EX_cpd15606_e'])
min_range <- c(min_range, length(rough_glytyr))
# output
rough_hydroxymandelate <- as.vector(rough_fluxes[,'EX_cpd03170_e'])
min_range <- c(min_range, length(rough_hydroxymandelate))
rough_ammonia <- as.vector(rough_fluxes[,'EX_cpd00013_e'])
min_range <- c(min_range, length(rough_ammonia))
smooth_aspartate <- as.vector(smooth_fluxes[,'EX_cpd00320_e'])
min_range <- c(min_range, length(smooth_aspartate))
rough_aspartate <- as.vector(rough_fluxes[,'EX_cpd00320_e'])
min_range <- c(min_range, length(rough_aspartate))
smooth_methylbutyrate <- as.vector(smooth_fluxes[,'EX_cpd19585_e'])
min_range <- c(min_range, length(smooth_methylbutyrate))
rough_methylbutyrate <- as.vector(rough_fluxes[,'EX_cpd19585_e'])
min_range <- c(min_range, length(rough_methylbutyrate))
rm(rough_fluxes, smooth_fluxes)

# Test differences in predicted flux distributions
glytyr_pval <- round(wilcox.test(rough_glytyr, smooth_glytyr, exact=FALSE)$p.value, 3)
aspartate_pval <- round(wilcox.test(rough_aspartate, smooth_aspartate, exact=FALSE)$p.value, 3)
methylbutyrate_pval <- round(wilcox.test(rough_methylbutyrate, smooth_methylbutyrate, exact=FALSE)$p.value, 3)

# Generate figures
smooth_col <- 'white'
rough_col <- 'cornflowerblue'

library(vioplot)
png(filename='~/Desktop/repos/Jenior_Cdifficile_2019/results/figures/Figure_3B.png', 
    units='in', width=4, height=5, res=300)
layout(matrix(c(1,2,
                3,4), nrow=2, ncol=2, byrow=TRUE))
par(mar=c(2,3,1.5,0.5), xpd=FALSE, las=1, mgp=c(1.8,0.7,0), lwd=1.7)

vioplot(smooth_glucose, 0, col=c(smooth_col,rough_col), main='D-Glucose', cex.main=1,
        ylim=c(-900,900), ylab='Sampled Flux', lwd=1.7, drawRect=FALSE, yaxt='n')
abline(h=0, lty=5, col='gray70')
vioplot(smooth_glucose, 0, col=c(smooth_col,rough_col), add=TRUE,
        ylim=c(-900,900), ylab='Sampled Flux', lwd=1.7, drawRect=FALSE, yaxt='n')
axis(side=2, at=seq(-900,900,300), cex.axis=0.7)
text(x=1.5, y=900, labels='Net Production', cex=0.8)
text(x=1.5, y=-900, labels='Net Consumption', cex=0.8)
mtext(c('Smooth','Rough'), side=1, padj=0.5, adj=c(0.1,0.9))
text(x=2, y=100, labels='inactive', cex=0.9)

vioplot(smooth_glytyr, rough_glytyr, col=c(smooth_col,rough_col), main='Gly-Tyr', cex.main=1,
        ylim=c(-300,300), ylab='Sampled Flux', lwd=1.7, drawRect=FALSE, yaxt='n')
abline(h=0, lty=5, col='gray70')
vioplot(smooth_glytyr, rough_glytyr, col=c(smooth_col,rough_col), add=TRUE,
        ylim=c(-300,300), ylab='Sampled Flux', lwd=1.7, drawRect=FALSE, yaxt='n')
axis(side=2, at=seq(-300,300,100), cex.axis=0.7)
text(x=1.5, y=300, labels='Net Production', cex=0.8)
text(x=1.5, y=-300, labels='Net Consumption', cex=0.8)
mtext(c('Smooth','Rough'), side=1, padj=0.5, adj=c(0.1,0.9))
segments(x0=1, x1=2, y0=50)
text(x=1.5, y=80, labels='***', cex=1.5, font=2)

vioplot(0, rough_hydroxymandelate, col=c(smooth_col,rough_col), main='4-Hydroxymandelate', cex.main=1,
        ylim=c(-300,300), ylab='Sampled Flux', lwd=1.7, drawRect=FALSE, yaxt='n')
abline(h=0, lty=5, col='gray70')
vioplot(0, rough_hydroxymandelate, col=c(smooth_col,rough_col), add=TRUE,
        ylim=c(-300,300), ylab='Sampled Flux', lwd=1.7, drawRect=FALSE, yaxt='n')
axis(side=2, at=seq(-300,300,100), cex.axis=0.7)
text(x=1.5, y=300, labels='Net Production', cex=0.8)
text(x=1.5, y=-300, labels='Net Consumption', cex=0.8)
mtext(c('Smooth','Rough'), side=1, padj=0.5, adj=c(0.1,0.9))
text(x=1, y=33, labels='inactive', cex=0.9)

vioplot(0, rough_ammonia, col=c(smooth_col,rough_col), main='Ammonia', cex.main=1,
        ylim=c(-1000,1000), ylab='Sampled Flux', lwd=1.7, drawRect=FALSE, yaxt='n')
abline(h=0, lty=5, col='gray70')
vioplot(0, rough_ammonia, col=c(smooth_col,rough_col), add=TRUE,
        ylim=c(-1000,1000), ylab='Sampled Flux', lwd=1.7, drawRect=FALSE, yaxt='n')
axis(side=2, at=seq(-1000,1000,500), cex.axis=0.7)
text(x=1.5, y=1000, labels='Net Production', cex=0.8)
text(x=1.5, y=-1000, labels='Net Consumption', cex=0.8)
mtext(c('Smooth','Rough'), side=1, padj=0.5, adj=c(0.1,0.9))
text(x=1, y=100, labels='inactive', cex=0.9)

dev.off()

#--------------------------------------------------#

# in vitro data
glucose <- read.delim('~/Desktop/repos/Jenior_Cdifficile_2019/data/invitro/glucose.tsv', sep='\t', header=TRUE)

# Subset data
wt <- as.numeric(glucose$wt)
mutant <- as.numeric(glucose$mutant)
rm(glucose)

# Calculate significance
p_val <- round(wilcox.test(wt, mutant, exact=F)$p.value, 3)

# Generate figure
par(mar=c(0.6,5,1,1))
stripchart(cfu_vegetative~treatment, data=vegetative_cfu, vertical=T, pch=19, lwd=2.2,
           ylim=c(1,9), xaxt='n', cex=0, col=select_palette,
           ylab='Glucose concentration', method='jitter', jitter=0.15, cex.lab=1.2)
stripchart(cfu_vegetative~treatment, data=vegetative_cfu, vertical=T, pch=19, lwd=2.5,
           ylim=c(1,9), xaxt='n', yaxt='n', cex=2, col=select_palette,
           ylab='Vegetative cfu/g content', method='jitter', jitter=0.15, cex.lab=1.2, add=TRUE)

# Draw axis break
axis.break(2, 1.5, style='slash')

# Draw median
segments(0.7, median(wt), 1.3, median(wt), lwd=3)
segments(1.7, median(mutant), 2.3, median(mutant), lwd=3)


# Adding significance to plot
segments(1, 5, 2, 5, lwd=3)
text(1.5, 5.4, labels='*', cex=3, font=2)



