
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
    par(mar=c(2,3,1.5,1), xpd=FALSE, las=1, mgp=c(2.1,0.7,0), lwd=2)
    vioplot(smooth_flux, rough_flux, col=c(smooth_col,rough_col), main=substrate, cex.main=1.1,
            ylim=c(0, ymax), ylab='Predicted Efflux', lwd=1.7, drawRect=FALSE, yaxt='n', yaxs='i')
    axis(side=2, at=seq(0, ymax, ymax/4), cex.axis=0.7, lwd=2)
    vioplot(smooth_flux, rough_flux, col=c(smooth_col,rough_col), main=substrate, cex.main=1.1, add=TRUE,
            ylim=c(0, ymax), ylab='Predicted Efflux', lwd=1.7, drawRect=FALSE, yaxt='n', yaxs='i')
    abline(h=0, lwd=2, lty=2, col='gray25')
    box(lwd=2)
    par(xpd=TRUE)
    text(x=1, y=-(ymax*0.1), labels='Smooth', cex=1.2)
    text(x=2, y=-(ymax*0.1), labels='Rough', cex=1.2)
    par(xpd=FALSE)
    if (length(smooth) == 1) {text(x=1, y=ymax*0.1, 'inactive', cex=1.1)}
    if (length(rough) == 1) {text(x=2, y=ymax*0.1, 'inactive', cex=1.1)}}

plot_uptake <- function(smooth, rough, substrate='test', ymax=1000) {
    smooth_flux <- subset(smooth, smooth >= 0)
    rough_flux <- subset(rough, rough >= 0)
    par(mar=c(2,3,1.5,1), xpd=FALSE, las=1, mgp=c(2.1,0.7,0), lwd=2)
    vioplot(smooth_flux, rough_flux, col=c(smooth_col,rough_col), main=substrate, cex.main=1.1,
            ylim=c(0, ymax), ylab='Predicted Uptake', lwd=1.7, drawRect=FALSE, yaxt='n', yaxs='i')
    axis(side=2, at=seq(0, ymax, ymax/4), cex.axis=0.7, lwd=2)
    vioplot(smooth_flux, rough_flux, col=c(smooth_col,rough_col), main=substrate, cex.main=1.1, add=TRUE,
            ylim=c(0, ymax), ylab='Predicted Uptake', lwd=1.7, drawRect=FALSE, yaxt='n', yaxs='i')
    abline(h=0, lwd=2, lty=2, col='gray25')
    box(lwd=2)
    par(xpd=TRUE)
    text(x=1, y=-(ymax*0.1), labels='Smooth', cex=1.2)
    text(x=2, y=-(ymax*0.1), labels='Rough', cex=1.2)
    par(xpd=FALSE)
    if (length(smooth) == 1) {text(x=1, y=ymax*0.1, 'inactive', cex=1.1)}
    if (length(rough) == 1) {text(x=2, y=ymax*0.1, 'inactive', cex=1.1)}}


# Format data
shared_rxns <- intersect(colnames(rough_fluxes), colnames(smooth_fluxes))
rough_fluxes <- rough_fluxes[,shared_rxns]
smooth_fluxes <- smooth_fluxes[,shared_rxns]
rm(shared_rxns)

# Subsample data
sample_size <- min(c(nrow(rough_fluxes), nrow(smooth_fluxes), 250))
sub_sample <- sample(1:min(c(nrow(rough_fluxes), nrow(smooth_fluxes))), sample_size, replace=FALSE)
rough_fluxes <- rough_fluxes[sub_sample,]
smooth_fluxes <- smooth_fluxes[sub_sample,]
rm(sample_size, sub_sample)

# Add column and merge
rough_fluxes$phase <- rep('rough', nrow(rough_fluxes))
smooth_fluxes$phase <- rep('smooth', nrow(smooth_fluxes))
phase_fluxes <- as.data.frame(rbind(rough_fluxes, smooth_fluxes))
rownames(phase_fluxes) <- c(paste0('rough_', c(1:nrow(rough_fluxes))), paste0('smooth_', c(1:nrow(smooth_fluxes))))
rm(rough_fluxes, smooth_fluxes)

# Create phase variant groups
phase_groups <- as.data.frame(cbind(rownames(phase_fluxes), phase_fluxes$phase))
colnames(phase_groups) <- c('sample', 'phase')
phase_groups <- list(rough=which(phase_groups$phase == 'rough'),
                     smooth=which(phase_groups$phase == 'smooth'))

# Find most informative metabolites with RF
library(randomForest)
condition <- as.factor(phase_fluxes$phase)
phase_fluxes$phase <- NULL
phase_fluxes <- droplevels(phase_fluxes)
rf_obj <- randomForest(condition ~ ., data=phase_fluxes, importance=TRUE, err.rate=TRUE, na.action=na.roughfix, ntree=1500, mtry=15)
rf_obj <- importance(rf_obj, type=1, scale=TRUE)
rf_mda <- as.data.frame(subset(rf_obj, rf_obj > (abs(min(rf_obj))*5))) # significance
phase_fluxes <- phase_fluxes[,rownames(rf_mda)]
rm(condition, rf_obj, rf_mda)

# Calculate correlation matrix
library(Hmisc)
corr_data <- rcorr(as.matrix(t(phase_fluxes)*10.0), type='spearman')
r_vals <- corr_data$r # correlation coefficients
p_vals <- corr_data$P # p-values
p_vals[is.na(p_vals)] <- 1
rm(corr_data)

# Prune correlations based on p-values
r_vals[which(p_vals > 0.05)] <- 0.0

#--------------------------------------------------------------------------------------------------#

# Generate figures
smooth_col <- 'lightblue1'
rough_col <- 'sienna2'
library(vioplot)
library(qgraph)
library(plotfunctions)

pdf(file='~/Desktop/repos/Jenior_Cdifficile_2019/results/figures/Figure_2.pdf', width=6.5, height=5.5)
layout(matrix(c(1,2,3,4,
                1,5,6,7,
                8,9,9,10), nrow=3, ncol=4, byrow=TRUE))

par(mar=c(5,3,0.5,0.5), xpd=FALSE, las=1, mgp=c(1.6,0.7,0), lwd=2)
boxplot(smooth_biomass, at=0.5, xlim=c(0,2), ylab='Simulated Biomass Flux', 
        ylim=c(0,60), yaxt='n', col=smooth_col, outline=FALSE, cex.lab=1.2,
        boxlwd=2, medlwd=2, staplelwd=2, whisklwd=2, whisklty=1, width=0.5)
boxplot(rough_biomass, at=1.5, add=TRUE, yaxt='n', col=rough_col,
        boxlwd=2, medlwd=2, staplelwd=2, whisklwd=2, whisklty=1, width=0.5, outline=FALSE)
axis(side=2, at=seq(0,60,10), cex.axis=0.9, lwd=1.7)
segments(x0=0.5, y0=56, x1=1.5, lwd=2)
text(x=1, y=57.5, '***', cex=1.5, font=2)
par(xpd=TRUE)
text(x=c(0.5,1.5), y=-8, labels=c('Smooth','Rough'), srt=55, cex=1.3)
text(x=-0.7, y=61, 'A', cex=1.5, font=2)
text(x=-0.6, y=-8, 'H', cex=1.5, font=2)
par(xpd=FALSE)

plot_uptake(0, rough_proline_uptake, substrate='Proline', ymax=300)
par(xpd=TRUE)
text(x=-0.24, y=300, 'B', cex=1.5, font=2)
par(xpd=FALSE)
plot_uptake(smooth_alanine_uptake, c(0,0,0), substrate='Alanine', ymax=400)
segments(x0=1, y0=300, x1=2, lwd=2)
text(x=1.5, y=320, '***', cex=1.4, font=2)
par(xpd=TRUE)
text(x=-0.24, y=400, 'C', cex=1.5, font=2)
par(xpd=FALSE)
plot_uptake(smooth_glc6p_uptake, rough_glc6p_uptake, substrate='Glucose-6P', ymax=500)
segments(x0=1, y0=400, x1=2, lwd=2)
text(x=1.5, y=425, '*', cex=1.4, font=2)
par(xpd=TRUE)
text(x=-0.24, y=500, 'D', cex=1.5, font=2)
par(xpd=FALSE)

plot_efflux(smooth_glcnac_efflux, 0, substrate='GlcNAc', ymax=200)
par(xpd=TRUE)
text(x=-0.24, y=200, 'E', cex=1.5, font=2)
par(xpd=FALSE)
plot_efflux(smooth_glucosamine_efflux, rough_glucosamine_efflux, substrate='D-Glucosamine')
segments(x0=1, y0=825, x1=2, lwd=2)
text(x=1.5, y=870, '***', cex=1.4, font=2)
par(xpd=TRUE)
text(x=-0.24, y=1000, 'F', cex=1.5, font=2)
par(xpd=FALSE)
plot_efflux(c(0,0,0), rough_alanine_efflux, substrate='Alanine', ymax=500)
segments(x0=1, y0=375, x1=2, lwd=2)
text(x=1.5, y=400, '***', cex=1.4, font=2)
par(xpd=TRUE)
text(x=-0.24, y=500, 'G', cex=1.5, font=2)
par(xpd=FALSE)

plot(0, type='n', ylim=c(0,10), xlim=c(0,10), ylab='', xlab='', xaxt='n', yaxt='n', axes=FALSE)

par(mar=c(0,0,0,0), xpd=NA, lwd=1.5)
qgraph(r_vals, minimum=0.3, vsize=3, layout='spring', weighted=TRUE,
       groups=phase_groups, borders=TRUE, labels=FALSE, legend=FALSE,
       posCol='royalblue3', negCol='red3', color=c(rough_col, smooth_col))
legend('topleft', legend=c('Smooth', 'Rough'), bty='n',
       pt.bg=c(smooth_col, rough_col), pch=21, pt.cex=2, cex=0.9)
plot(0, type='n', ylim=c(0,10), xlim=c(0,10), ylab='', xlab='', xaxt='n', yaxt='n', axes=FALSE)
corr_gradient <- colorRampPalette(c('red3','white','royalblue3'))(50)
points(x=seq(3,6.92,0.08), y=rep(7,50), pch=15, cex=4, col=corr_gradient)
points(x=c(7,7,7), y=c(6.9,7,7.1), pch=15, cex=4, col='white')
text(x=4, y=8.5, 'Correlation coefficient', cex=0.9)
text(x=c(2,4,6), y=5.7, c('-1','0','+1'), cex=0.8)

dev.off()
