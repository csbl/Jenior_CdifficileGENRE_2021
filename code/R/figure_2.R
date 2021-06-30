
# Read in data
rough_fluxes <- read.delim('~/Desktop/repos/Jenior_Cdifficile_2019/data/transcript/tamayo_etal/phase_variation/riptide_rough_maxfit/flux_samples.tsv', sep='\t', header=TRUE, row.names=1)
smooth_fluxes <- read.delim('~/Desktop/repos/Jenior_Cdifficile_2019/data/transcript/tamayo_etal/phase_variation/riptide_smooth_maxfit/flux_samples.tsv', sep='\t', header=TRUE, row.names=1)

# Get reactions of interest
smooth_biomass <- as.vector(smooth_fluxes[,'biomass'])
rough_biomass <- as.vector(rough_fluxes[,'biomass'])
rough_proline_uptake <- as.vector(rough_fluxes[,'EX_cpd00567_e']) * -1.0
smooth_glcnac_efflux <- as.vector(smooth_fluxes[,'EX_cpd00122_e'])
smooth_glucosamine_efflux <- as.vector(smooth_fluxes[,'EX_cpd00276_e'])
rough_glucosamine_efflux <- as.vector(rough_fluxes[,'EX_cpd00276_e'])
smooth_orn_uptake <- as.vector(smooth_fluxes[,'EX_cpd00064_e']) * -1.0
smooth_glc6p_uptake <- as.vector(smooth_fluxes[,'EX_cpd00079_e']) * -1.0
rough_glc6p_uptake <- as.vector(rough_fluxes[,'EX_cpd00079_e']) * -1.0
smooth_co2_efflux <- as.vector(smooth_fluxes[,'EX_cpd00011_e'])
rough_co2_efflux <- as.vector(rough_fluxes[,'EX_cpd00011_e'])
smooth_isoval_efflux <- as.vector(smooth_fluxes[,'EX_cpd05178_e'])
rough_isoval_efflux <- as.vector(rough_fluxes[,'EX_cpd05178_e'])
smooth_acetate_efflux <- as.vector(smooth_fluxes[,'EX_cpd00029_e'])
rough_acetate_efflux <- as.vector(rough_fluxes[,'EX_cpd00029_e'])

rough_alanine_efflux <- as.vector(rough_fluxes[,'EX_cpd00035_e'])
smooth_alanine_uptake <- as.vector(smooth_fluxes[,'EX_cpd00035_e'])
alanine_pval <- round(wilcox.test(rough_alanine_efflux, smooth_alanine_uptake, exact=FALSE)$p.value, 5) * 100.0
rough_alanine_efflux <- subset(rough_alanine_efflux, rough_alanine_efflux > 0.0)
smooth_alanine_uptake <- subset(smooth_alanine_uptake, smooth_alanine_uptake < 0.0) * -1.0

# Test differences
shapiro.test(smooth_glucosamine_efflux)$p.value
shapiro.test(rough_glucosamine_efflux)$p.value
shapiro.test(smooth_co2_efflux)$p.value
shapiro.test(rough_co2_efflux)$p.value
co2_pval <- round(wilcox.test(smooth_co2_efflux, rough_co2_efflux, exact=FALSE)$p.value, 5) * 100.0
shapiro.test(smooth_isoval_efflux)$p.value
shapiro.test(rough_isoval_efflux)$p.value
isoval_pval <- round(wilcox.test(smooth_isoval_efflux, rough_isoval_efflux, exact=FALSE)$p.value, 5) * 100.0
shapiro.test(smooth_acetate_efflux)$p.value
shapiro.test(rough_acetate_efflux)$p.value
acetate_pval <- round(wilcox.test(smooth_acetate_efflux, rough_acetate_efflux, exact=FALSE)$p.value, 5) * 100.0
shapiro.test(smooth_glc6p_uptake)$p.value
shapiro.test(rough_glc6p_uptake)$p.value
glc6p_pval <- round(wilcox.test(smooth_glc6p_uptake, rough_glc6p_uptake, exact=FALSE)$p.value, 5) * 100.0
shapiro.test(smooth_biomass)$p.value
shapiro.test(rough_biomass)$p.value
biomass_pval <- round(wilcox.test(smooth_biomass, rough_biomass, exact=FALSE)$p.value, 5) * 100.0
smooth_doubling <- (1/smooth_biomass) * 3600.0
rough_doubling <- (1/rough_biomass) * 3600.0

plot_efflux <- function(smooth, rough, substrate='test', ymax=1000) {
    smooth_flux <- subset(smooth, smooth >= 0)
    rough_flux <- subset(rough, rough >= 0)
    vioplot(smooth_flux, rough_flux, col=c(smooth_col,rough_col), main=substrate, cex.main=1.1,
            ylim=c(0, ymax), ylab='Predicted Efflux', lwd=1.7, drawRect=FALSE, yaxt='n', yaxs='i')
    axis(side=2, at=seq(0, ymax, ymax/4), cex.axis=0.7, lwd=2)
    vioplot(smooth_flux, rough_flux, col=c(smooth_col,rough_col), main=substrate, cex.main=1.1, add=TRUE,
            ylim=c(0, ymax), ylab='Predicted Efflux', lwd=1.7, drawRect=FALSE, yaxt='n', yaxs='i')
    abline(h=0, lwd=2, lty=2, col='gray25')
    box(lwd=2)
    par(xpd=TRUE)
    text(x=1, y=-(ymax*0.11), labels='Smooth', cex=1.2)
    text(x=2, y=-(ymax*0.12), labels='Rough', cex=1.2)
    par(xpd=FALSE)
    if (length(smooth) == 1) {text(x=1, y=ymax*0.1, 'inactive', cex=1.1)}
    if (length(rough) == 1) {text(x=2, y=ymax*0.1, 'inactive', cex=1.1)}}

plot_uptake <- function(smooth, rough, substrate='test', ymax=1000) {
    smooth_flux <- subset(smooth, smooth >= 0)
    rough_flux <- subset(rough, rough >= 0)
    vioplot(smooth_flux, rough_flux, col=c(smooth_col,rough_col), main=substrate, cex.main=1.1,
            ylim=c(0, ymax), ylab='Predicted Uptake', lwd=1.7, drawRect=FALSE, yaxt='n', yaxs='i')
    axis(side=2, at=seq(0, ymax, ymax/4), cex.axis=0.7, lwd=2)
    vioplot(smooth_flux, rough_flux, col=c(smooth_col,rough_col), main=substrate, cex.main=1.1, add=TRUE,
            ylim=c(0, ymax), ylab='Predicted Uptake', lwd=1.7, drawRect=FALSE, yaxt='n', yaxs='i')
    abline(h=0, lwd=2, lty=2, col='gray25')
    box(lwd=2)
    par(xpd=TRUE)
    text(x=1, y=-(ymax*0.11), labels='Smooth', cex=1.2)
    text(x=2, y=-(ymax*0.12), labels='Rough', cex=1.2)
    par(xpd=FALSE)
    if (length(smooth) == 1) {text(x=1, y=ymax*0.1, 'inactive', cex=1.1)}
    if (length(rough) == 1) {text(x=2, y=ymax*0.1, 'inactive', cex=1.1)}}

# Machine learning results
rf_mda <- read.delim('~/Desktop/repos/Jenior_Cdifficile_2019/data/phase_flux_mda.tsv', sep='\t', header=TRUE)
rf_mda$name <- gsub('_', ' ', rf_mda$name)
rf_mda <- rf_mda[order(-rf_mda$mda),]
rf_mda <- rf_mda[c(1:15),]
rf_mda <- rf_mda[order(rf_mda$mda),]

#--------------------------------------------------------------------------------------------------#

# Generate figures
smooth_col <- 'lightsteelblue2'
rough_col <- 'red3'
library(vioplot)
library(plotfunctions)

pdf(file='~/Desktop/repos/Jenior_Cdifficile_2019/results/figures/Figure_2.pdf', width=6, height=4.5)
layout(matrix(c(1,2,2,
                1,2,2,
                3,4,5), nrow=3, ncol=3, byrow=TRUE))

plot(0,type='n',axes=FALSE,ann=FALSE)

par(mar=c(2.5, 1.5, 1, 1), mgp=c(1.3, 0.4, 0), xpd=FALSE, lwd=2)
dotchart(rf_mda$mda,  xlab='Mean Decrease Accuracy (%)', xlim=c(0,22),  
         pch=16, lwd=1.7, xaxs='i', pt.cex=0.1, cex=0.9, xaxt='n')
axis(1, at=seq(0,20,5), lwd=2, cex=0.7)
text(x=-0.5, y=seq(1.4,15.4,1), labels=rf_mda$name, cex=0.9, pos=4)
points(x=rf_mda$mda, y=c(1:15), pch=21, cex=1.5, bg='chartreuse3')
par(xpd=TRUE)
text(x=-1.8, y=15.5, 'B', cex=1.2, font=2)
text(x=-1.8, y=13, 'A', cex=1.2, font=2)
par(xpd=FALSE)

par(mar=c(2,3.7,1.5,1.5), xpd=FALSE, las=1, mgp=c(2.1,0.7,0), lwd=2)
plot_uptake(smooth_alanine_uptake, c(0,0,0), substrate='Alanine', ymax=400)
segments(x0=1, y0=300, x1=2, lwd=2)
text(x=1.5, y=330, '***', cex=1.4, font=2)
par(xpd=TRUE)
text(x=-0.2, y=425, 'C', cex=1.4, font=2)
par(xpd=FALSE)
plot_efflux(smooth_glcnac_efflux, 0, substrate='N-Acetylglucosamine', ymax=200)
par(xpd=TRUE)
text(x=-0.2, y=212.5, 'D', cex=1.4, font=2)
par(xpd=FALSE)
plot_efflux(smooth_glucosamine_efflux, rough_glucosamine_efflux, substrate='D-Glucosamine')
segments(x0=1, y0=840, x1=2, lwd=2)
text(x=1.5, y=900, '***', cex=1.4, font=2)
par(xpd=TRUE)
text(x=-0.2, y=1062.5, 'E', cex=1.4, font=2)
par(xpd=FALSE)

dev.off()

