
# Read in data
rough_fluxes <- read.delim('~/Desktop/repos/Jenior_Cdifficile_2019/data/transcript/tamayo_etal/phase_variation/riptide_rough/flux_samples.tsv', sep='\t', header=TRUE)
smooth_fluxes <- read.delim('~/Desktop/repos/Jenior_Cdifficile_2019/data/transcript/tamayo_etal/phase_variation/riptide_smooth/flux_samples.tsv', sep='\t', header=TRUE)

# Get reactions of interest
smooth_biomass <- as.vector(smooth_fluxes[,'biomass'])
rough_biomass <- as.vector(rough_fluxes[,'biomass'])
smooth_glucose <- as.vector(smooth_fluxes[,'EX_cpd00027_e'])
smooth_proline <- as.vector(smooth_fluxes[,'EX_cpd00567_e'])
rough_proline <- as.vector(rough_fluxes[,'EX_cpd00567_e'])
smooth_glcnac <- as.vector(smooth_fluxes[,'EX_cpd00122_e'])
rough_glcnac <- as.vector(rough_fluxes[,'EX_cpd00122_e'])
smooth_cys <- as.vector(smooth_fluxes[,'EX_cpd00424_e']) 
rough_glutamate <- as.vector(rough_fluxes[,'EX_cpd00023_e'])
rm(rough_fluxes, smooth_fluxes)

# Do math
pvals <- c()
for (x in c(1:1000)) {
    test_1 <- sample(smooth_proline, size=10)
    test_2 <- sample(rough_proline, size=10)
    pvals[x] <- round(wilcox.test(test_1, test_2, exact=FALSE)$p.value, 3)
}
proline_pval <- median(pvals)
pvals <- c()
for (x in c(1:1000)) {
    test_1 <- sample(smooth_glcnac, size=10)
    test_2 <- sample(rough_glcnac, size=10)
    pvals[x] <- round(wilcox.test(test_1, test_2, exact=FALSE)$p.value, 3)
}
glcnac_pval <- median(pvals)
pvals <- c()
for (x in c(1:1000)) {
    test_1 <- sample(smooth_biomass, size=10)
    test_2 <- sample(rough_biomass, size=10)
    pvals[x] <- round(wilcox.test(test_1, test_2, exact=FALSE)$p.value, 3)
}
biomass_pval <- median(pvals)
rm(pvals, test_1, test_2)
smooth_biomass <- (1/smooth_biomass) * 3600.0
rough_biomass <- (1/rough_biomass) * 3600.0
smooth_biomass <- smooth_biomass - 40
rough_biomass <- rough_biomass - 40
smooth_glucose <- smooth_glucose * -1.0
smooth_proline <- smooth_proline * -1.0
rough_proline <- rough_proline * -1.0
smooth_cys <- smooth_cys * -1.0
rough_glutamate <- subset(rough_glutamate, rough_glutamate <= 0)
rough_glutamate <- rough_glutamate * -1.0

# Generate figures
smooth_col <- 'white'
rough_col <- 'cornflowerblue'
library(vioplot)
library(plotrix)

png(filename='~/Desktop/repos/Jenior_Cdifficile_2019/results/figures/Figure_3A.png', 
    units='in', width=1.5, height=3, res=300)
par(mar=c(3,2.5,0.5,0.5), xpd=FALSE, las=1, mgp=c(1.5,0.7,0), lwd=1.7)
boxplot(smooth_biomass, at=0.5, xlim=c(0,2), ylab='Predicted Doubling Time (min)', 
        ylim=c(0,40), yaxt='n', col=smooth_col, outline=FALSE, cex.lab=0.9,
        boxlwd=2, medlwd=2, staplelwd=2, whisklwd=2, whisklty=1, width=0.9, cex.main=0.9)
boxplot(rough_biomass, at=1.5, xlim=c(0,2), add=TRUE, yaxt='n', col=rough_col,
        boxlwd=2, medlwd=2, staplelwd=2, whisklwd=2, whisklty=1, width=0.9, outline=FALSE)
axis(side=2, at=seq(0,40,5), labels=c(0,seq(35,70,5)), cex.axis=0.7, lwd=1.7)
segments(x0=0.5, y0=35, x1=1.5, lwd=2)
text(x=1, y=37, 'n.s.', cex=0.7)
axis.break(2, 2.5, style='slash', brw=0.04)
par(xpd=TRUE)
text(x=c(0.5,1.5), y=-8, labels=c('Smooth','Rough'), srt=55, cex=1)
par(xpd=FALSE)
dev.off()

png(filename='~/Desktop/repos/Jenior_Cdifficile_2019/results/figures/Figure_3B.png', 
    units='in', width=2.5, height=3, res=300)
par(mar=c(1.5,3,1.5,0.5), xpd=FALSE, las=1, mgp=c(1.8,0.7,0), lwd=2, yaxs='i')
vioplot(smooth_glcnac, rough_glcnac, col=c(smooth_col,rough_col), main='N-Acetylglucosamine', cex.main=1,
        ylim=c(0,1100), ylab='Predicted Secretion Flux', lwd=2, drawRect=FALSE, yaxt='n')
axis(side=2, at=seq(0,1000,200), cex.axis=0.7, lwd=2)
box(lwd=2)
mtext(c('Smooth','Rough'), side=1, padj=0.5, adj=c(0.15,0.85), cex=1.1)
segments(x0=1, x1=2, y0=1000, lwd=2)
text(x=1.5, y=1040, labels='*', cex=1.6, font=2)
dev.off()

png(filename='~/Desktop/repos/Jenior_Cdifficile_2019/results/figures/Figure_3FG.png', 
    units='in', width=2, height=4, res=300)
layout(matrix(c(1,2), nrow=2, ncol=1, byrow=TRUE))
par(mar=c(1.5,2.5,1,0.5), xpd=FALSE, las=1, mgp=c(1.8,0.7,0), lwd=1.7, yaxs='i')

vioplot(smooth_glucose, 0, col=c(smooth_col,rough_col), 
        ylim=c(0,1000), lwd=2, drawRect=FALSE, yaxt='n')
axis(side=2, at=seq(0,1000,200), cex.axis=0.5, lwd=2)
box(lwd=2)
text(x=2, y=75, labels='inactive', cex=0.6)
par(xpd=TRUE)
text(x=c(1,2), y=-100, labels=c('Smooth','Rough'), cex=0.8)
text(x=1.5, y=1070, labels='D-Glucose', cex=0.8, font=2)
text(x=-0.2, y=500, labels='Predicted Uptake Flux', cex=0.8, srt=90)
par(xpd=FALSE)
vioplot(smooth_proline, rough_proline, col=c(smooth_col,rough_col),
        ylim=c(0,1000), lwd=2, drawRect=FALSE, yaxt='n')
axis(side=2, at=seq(0,1000,200), cex.axis=0.5, lwd=2)
box(lwd=2)
segments(x0=1, x1=2, y0=900, lwd=2)
text(x=1.5, y=940, labels='**', cex=1.2, font=2)
par(xpd=TRUE)
text(x=c(1,2), y=-100, labels=c('Smooth','Rough'), cex=0.8)
text(x=1.5, y=1070, labels='Proline', cex=0.8, font=2)
text(x=-0.2, y=500, labels='Predicted Uptake Flux', cex=0.8, srt=90)
par(xpd=FALSE)
dev.off()

