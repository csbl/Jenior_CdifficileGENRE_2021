
#install.packages('vegan')
library(vegan)

stepRarefy <- function(abundVect, rareStep, coverage){
  options(warn=-1)
  rareVect <- c()
  subVect <- seq(0, 25000, rareStep)
  for (x in 1:length(subVect)) {
    rareVect[x] <- as.numeric(sum(round(as.vector(rrarefy(abundVect, sample=subVect[x])) / coverage) != 0))
  }
  options(warn=0)
  return(rareVect)
}

cdf_transcription <- read.delim('~/Desktop/repos/Cdiff_modeling/data/transcript/cdf_transcription.sub.format.csv',
                                sep=',', header=TRUE, row.names=1)

cef <- cdf_transcription$cefoperazone
clinda <- cdf_transcription$clindamycin
strep <- cdf_transcription$streptomycin
gnoto <- cdf_transcription$gnotobiotic
rm(cdf_transcription)

#------------------------#

# Generate curves
cef_collectors <- stepRarefy(cef, 100, 10)
clinda_collectors <- stepRarefy(clinda, 100, 10)
strep_collectors <- stepRarefy(strep, 100, 10)
rm(stepRarefy)

# Find terminal slopes
cef_slope <- (cef_collectors[251] - cef_collectors[250]) / 100
clinda_slope <- (clinda_collectors[251] - clinda_collectors[250]) / 100
strep_slope <- (strep_collectors[251] - strep_collectors[250]) / 100

# Find score densities
cef_density <- density(cef)
clinda_density <- density(clinda)
strep_density <- density(strep)
gnoto_density <- density(gnoto)

# Find quantiles
cef_quantile <- quantile(cef, c(0.5, 0.625, 0.75, 0.875)) 
clinda_quantile <- quantile(clinda, c(0.5, 0.625, 0.75, 0.875)) 
strep_quantile <- quantile(strep, c(0.5, 0.625, 0.75, 0.875)) 
gnoto_quantile <- quantile(gnoto, c(0.5, 0.625, 0.75, 0.875)) 
quantile_1 <- median(c(cef_quantile[1],clinda_quantile[1],strep_quantile[1],gnoto_quantile[1]))
quantile_2 <- median(c(cef_quantile[2],clinda_quantile[2],strep_quantile[2],gnoto_quantile[2]))
quantile_3 <- median(c(cef_quantile[3],clinda_quantile[3],strep_quantile[3],gnoto_quantile[3]))
quantile_4 <- median(c(cef_quantile[3],clinda_quantile[3],strep_quantile[3],gnoto_quantile[3]))

#------------------------#

# Collector's curve
pdf(file='~/Desktop/repos/Cdiff_modeling/results/figures/cdf_transcript_collectors.pdf', width=5, height=4)
par(mar=c(4,4,1,1), mgp=c(2.5, 0.75, 0), las=1)
plot(1, type='n', xlim=c(0,260), ylim=c(0,700), xaxt='n', xlab='Sampling Depth', ylab='Detected Genes')
axis(1, at=seq(0,250,50), label=seq(0,25000,5000))
legend('bottomright', legend=c('Streptomycin','Cefoperazone','Clindamycin'), 
       col=c('red', 'blue', 'green'), pch=22, pt.cex=0, lwd=3)
lines(cef_collectors, col='blue', lwd=3)
lines(strep_collectors, col='red', lwd=3)
lines(clinda_collectors, col='green', lwd=3)
legend('topleft', legend='x10 coverage', pt.cex=0, bty='n', cex=1.5)
box()
dev.off()

#------------------------#

# Transcript abundance histogram

pdf(file='~/Desktop/repos/Cdiff_modeling/results/figures/cdf_transcript_density.pdf', width=5, height=4)
par(mar=c(4,4,1,1), las=1, mgp=c(2.5,1,0), xaxs='i', yaxs='i', xpd=FALSE)
plot(gnoto_density, xlim=c(-5,100), ylim=c(0,0.2), main='', yaxt='n', cex.lab=1.5, 
     xlab='Transcript Density', ylab='Number of Genes', lwd=2.5, col='darkorange2')
par(new=TRUE)
plot(clinda_density, xlim=c(-5,100), ylim=c(0,0.2), main='', yaxt='n', cex.lab=1.5, 
     xlab='', ylab='', lwd=2.5, col='firebrick3')
par(new=TRUE)
plot(strep_density, xlim=c(-5,100), ylim=c(0,0.2), main='', yaxt='n', cex.lab=1.5, 
     xlab='', ylab='', lwd=2.5, col='blue3')
par(new=TRUE)
plot(cef_density, xlim=c(-5,100), ylim=c(0,0.2), main='', yaxt='n', cex.lab=1.5, 
     xlab='', ylab='', lwd=2.5, col='chartreuse4')
axis(2, at=seq(0,0.2,0.04), label=seq(0,200,40))
legend('topright', legend=c('Cefoperazone', 'Streptomycin','Clindamycin','Gnotobiotic'), pt.cex=2, pch=15,
       col=c('chartreuse4','blue3','firebrick3','darkorange2'))
abline(v=c(0,quantile_1,quantile_2,quantile_3,quantile_4), 
       lty=5, lwd=2, col=c('darkred','gray75','gray50','gray25','black'))

dev.off()

#------------------------#


