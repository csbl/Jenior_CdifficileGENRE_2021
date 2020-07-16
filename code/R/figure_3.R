
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
# biomass
rough_biomass <- unique(as.vector(rough_fluxes[,'biomass']))
min_range <- c(min_range, length(rough_biomass))
smooth_biomass <- unique(as.vector(smooth_fluxes[,'biomass']))
min_range <- c(min_range, length(smooth_biomass))
rm(rough_fluxes, smooth_fluxes)

# Subsample data
min_range <- min(min_range)
sample_size <- round(min_range * 0.8)
sub_sample <- sample(1:min_range, sample_size, replace=FALSE)
smooth_glucose <- smooth_glucose[sub_sample]
smooth_glytyr <- smooth_glytyr[sub_sample]
rough_glytyr <- rough_glytyr[sub_sample]
rough_hydroxymandelate <- rough_hydroxymandelate[sub_sample]
rough_ammonia <- rough_ammonia[sub_sample]
smooth_aspartate <- smooth_aspartate[sub_sample]
rough_aspartate <- rough_aspartate[sub_sample]
smooth_methylbutyrate <- smooth_methylbutyrate[sub_sample]
rough_methylbutyrate <- rough_methylbutyrate[sub_sample]
rough_biomass <- rough_biomass[sub_sample]
smooth_biomass <- smooth_biomass[sub_sample]
rm(sub_sample, sample_size, min_range)

# Test differences in predicted flux distributions
biomass_pval <- round(wilcox.test(rough_biomass, smooth_biomass, exact=FALSE)$p.value, 3)
glytyr_pval <- round(wilcox.test(rough_glytyr, smooth_glytyr, exact=FALSE)$p.value, 3)
aspartate_pval <- round(wilcox.test(rough_aspartate, smooth_aspartate, exact=FALSE)$p.value, 3)
methylbutyrate_pval <- round(wilcox.test(rough_methylbutyrate, smooth_methylbutyrate, exact=FALSE)$p.value, 3)

# Generate figures
smooth_col <- 'chocolate2'
rough_col <- 'cyan3'

library(plotrix)
png(filename='~/Desktop/repos/Jenior_Cdifficile_2019/results/figures/figure_3A.png', 
    units='in', width=1.5, height=4, res=300)
par(mar=c(4,3,1.5,0.5), xpd=FALSE, las=1, mgp=c(2,0.8,0), lwd=1.7)
barplot(c(15.7, 16.0), col=c(smooth_col,rough_col), 
        ylim=c(0,20), ylab='Optimal Doubling Time (min)', lwd=1.7, main='iCdR700', yaxt='n')
axis(side=2, at=seq(0,20,10), labels=c(0,30,40), cex.axis=0.9, lwd=1.7)
box()
axis.break(2, 4, style='slash')
par(xpd=TRUE)
text(x=c(0.035,0.095), y=-2.5, labels=c('Smooth','Rough'), srt=55, cex=1.1)
par(xpd=FALSE)
dev.off()

library(vioplot)
flux_plot <- function(smooth, rough, cpd_name, pval) {
  #plot_title <- paste0(cpd_name, '\nExchange')
  par(mar=c(2,3,1.5,0.5), xpd=FALSE, las=1, mgp=c(1.8,0.7,0), lwd=1.7)
  vioplot(smooth, rough, col=c(smooth_col,rough_col), main=cpd_name, cex.main=1,
          ylim=c(-1100,1100), ylab='Sampled Flux', yaxt='n', lwd=1.7, drawRect=FALSE)
  abline(h=0, lty=5, col='gray70')
  vioplot(smooth, rough, col=c(smooth_col,rough_col), add=TRUE,
          ylim=c(-1100,1100), ylab='Sampled Flux', yaxt='n', lwd=1.7, drawRect=FALSE)
  axis(side=2, at=seq(-1000,1000,500), cex.axis=0.7)
  text(x=1.5, y=1100, labels='Net Production', cex=0.9)
  text(x=1.5, y=-1100, labels='Net Consumption', cex=0.9)
  mtext(c('Smooth','Rough'), side=1, padj=0.5, adj=c(0.1,0.9))
  if (max(rough) > 800) {
    sigy_1 <- -500
    sigy_2 <- -430
  } else { 
    sigy_1 <- 800
    sigy_2 <- 870
  }
  if (max(smooth) > 800) {
    sigy_1 <- -500
    sigy_2 <- -430
  } else { 
    sigy_1 <- 800
    sigy_2 <- 870
  }
  if (length(smooth) != 1) {
    if (length(rough) != 1) {
      if (pval <= 0.001) {
        segments(x0=1, x1=2, y0=sigy_1)
        text(x=1.5, y=sigy_2, labels='***', cex=1.5, font=2)
      }
    }
  }
  if (length(smooth) == 1) {text(x=1, y=100, labels='inactive', cex=0.9)}
  if (length(rough) == 1) {text(x=2, y=100, labels='inactive', cex=0.9)}
}


png(filename='~/Desktop/repos/Jenior_Cdifficile_2019/results/figures/figure_3C.png', 
    units='in', width=4, height=6, res=300)
layout(matrix(c(1,2,
                3,4), nrow=2, ncol=2, byrow=TRUE))

flux_plot(smooth_glucose, 0, 'D-Glucose', 1)
flux_plot(smooth_glytyr, rough_glytyr, 'Gly-Tyr', glytyr_pval)
flux_plot(0, rough_hydroxymandelate, '4-Hydroxymandelate', 1)
flux_plot(0, rough_ammonia, 'Ammonia', 1)

dev.off()


