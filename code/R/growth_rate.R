
# Campos M et al. (2014) A constant size extension drives bacterial cell size homeostasis. 
# Cell 159:1433–46. doi: 10.1016/j.cell.2014.11.022. p.1439 right column top paragraph
# bionumbers.hms.harvard.edu/bionumber.aspx?id=111767

library(vioplot)
library(scales)
library(vegan)
library(plotrix)

#------------------------------------------------------------------------------

# Load in data and format
growth_rates <- as.data.frame(t(read.delim(file='~/Desktop/repos/Cdiff_modeling/data/growth_rates.tsv', sep='\t', header=FALSE, row.names=1)))
base_pfba <- as.numeric(growth_rates$base_pfba)
base_pfba  <- base_pfba[!is.na(base_pfba)]
lb_aerobic <- as.numeric(growth_rates$lb_aerobic)
lb_aerobic  <- lb_aerobic[!is.na(lb_aerobic)]
m9_gluc_aerobic <- as.numeric(growth_rates$m9_gluc_aerobic)
m9_gluc_aerobic  <- m9_gluc_aerobic[!is.na(m9_gluc_aerobic)]
m9_gluc_anaerobic <- as.numeric(growth_rates$m9_gluc_anaerobic)
m9_gluc_anaerobic  <- m9_gluc_anaerobic[!is.na(m9_gluc_anaerobic)]
full_range  <- as.numeric(growth_rates$base_range)
full_range  <- full_range[!is.na(full_range)]
rm(growth_rates)

# Subsample data
full_range <- sample(full_range, 5000, replace=FALSE)
base_pfba <- sample(base_pfba, 5000, replace=FALSE)
lb_aerobic <- sample(lb_aerobic, 5000, replace=FALSE)
m9_gluc_aerobic <- sample(m9_gluc_aerobic, 5000, replace=FALSE)
m9_gluc_anaerobic <- sample(m9_gluc_anaerobic, 5000, replace=FALSE)

# Find ceiling for data
max(c(max(full_range),max(base_pfba),max(lb_aerobic),max(m9_gluc_aerobic),max(m9_gluc_anaerobic)))

# Calculate significant differences
pvals <- p.adjust(c(round(wilcox.test(base_pfba, m9_gluc_aerobic, exact=FALSE)$p.value,5),
           round(wilcox.test(base_pfba, m9_gluc_anaerobic, exact=FALSE)$p.value,5),
           round(wilcox.test(base_pfba, lb_aerobic, exact=FALSE)$p.value,5),
           round(wilcox.test(m9_gluc_aerobic, m9_gluc_anaerobic, exact=FALSE)$p.value,5),
           round(wilcox.test(m9_gluc_aerobic, lb_aerobic, exact=FALSE)$p.value,5),
           round(wilcox.test(m9_gluc_anaerobic, lb_aerobic, exact=FALSE)$p.value,5)), method='BH')

# Generate figure
#pdf(file='~/Desktop/repos/Cdiff_modeling/results/growth_rates.pdf', width=6, height=4)
png(filename='~/Desktop/repos/Cdiff_modeling/results/figure_3c.png', units='in', width=6, height=4, res=300)

par(mar=c(3,3.3,1,1), xpd=FALSE, las=1, mgp=c(2,0.75,0), lwd=2)
plot(0,0,type="n",xlim=c(0.5,5), ylim=c(0,1.5),  xaxt = 'n', xlab='', yaxt='n', 
     ylab=expression(paste('Calculated Growth Rate (hr'^'-1',')')))
axis(2, at=c(0,0.5,1,1.5), labels=c('0.0','0.5','0.1','3.0'), lwd=2) 
abline(v=1.75, lwd=1.5, lty=2)
text(x=3.5, y=0.01, labels='RIPTiDe-contextualized', cex=1.2)

# LB
segments(x0=2.4, y0=1, y1=1.5, lty=3, lwd=3, col='firebrick2')
segments(x0=2.3, x1=2.5, y0=1, lwd=3, col='firebrick2')
segments(x0=2.3, x1=2.5, y0=1.5, lwd=3, col='firebrick2')
#text(x=2.65, y=1.5, labels='S3', cex=0.8, col='firebrick2', font=2)

# M9 - aerobic
segments(x0=3.5, y0=0.63, y1=0.73, lty=3, lwd=3, col='firebrick2')
segments(x0=3.4, x1=3.6, y0=0.63, lwd=3, col='firebrick2')
segments(x0=3.4, x1=3.6, y0=0.73, lwd=3, col='firebrick2')
#text(x=3.75, y=0.73, labels='30', cex=0.8, col='firebrick2', font=2)

# M9 - anaerobic
segments(x0=4.5, y0=0.17, y1=0.6, lty=3, lwd=3, col='firebrick2')
segments(x0=4.4, x1=4.6, y0=0.17, lwd=3, col='firebrick2')
segments(x0=4.4, x1=4.6, y0=0.6, lwd=3, col='firebrick2')
#text(x=4.75, y=0.6, labels='30', cex=0.8, col='firebrick2', font=2)

legend('topright', legend='Experimentally measured', 
       pch=1, cex=0.8, pt.cex=0, box.lwd=2, lwd=3, col='firebrick2')

#vioplot(full_range, at=1, col='gray25', lwd=2, drawRect=FALSE, add=TRUE)
vioplot(base_pfba, at=1, col='#b2b2b1', lwd=2, drawRect=FALSE, add=TRUE)
vioplot(lb_aerobic, at=2.5, col=alpha('chocolate2', 0.8), lwd=2, drawRect=FALSE, add=TRUE)
vioplot(m9_gluc_aerobic, at=3.5, col=alpha('blue3', 0.8), lwd=2, drawRect=FALSE, add=TRUE)
vioplot(m9_gluc_anaerobic, at=4.5, col=alpha('white', 0.8), lwd=2, drawRect=FALSE, add=TRUE)
mtext(c('Unweighted\npFBA','LB media\nAerobic','M9 media\nAerobic','M9 media\nAnaerobic'), 
      side=1, at=c(1,2.5,3.5,4.5), padj=1)

segments(x0=2.5, x1=3.5, y0=1.3, lwd=2)
text(x=3, y=1.325, '***', font=2, cex=1.5)
segments(x0=3.5, x1=4.5, y0=1.2, lwd=2)
text(x=4, y=1.225, '***', font=2, cex=1.5)
segments(x0=2.5, x1=4.5, y0=1.1, lwd=2)
text(x=3.5, y=1.125, '***', font=2, cex=1.5)

axis.break(2, 1.25, style='slash') # axis break

dev.off()

# Clean up
rm(list=ls())
gc()

#------------------------------------------------------------------------------

invivo_rates <- as.data.frame(read.delim(file='~/Desktop/repos/Cdiff_modeling/data/invivo_growth_rates.tsv', 
                                           sep='\t', header=FALSE))[,1][1:8304]
invivo_pvals <- p.adjust(c(round(wilcox.test(invivo_rates, m9_gluc_aerobic, exact=FALSE)$p.value,5),
                    round(wilcox.test(invivo_rates, m9_gluc_anaerobic, exact=FALSE)$p.value,5),
                    round(wilcox.test(invivo_rates, lb_aerobic, exact=FALSE)$p.value,5)), method='BH')

png(filename='~/Desktop/repos/Cdiff_modeling/results/figure_S3.png', 
    units='in', width=2, height=4, res=300)
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
text(x=0.9, y=1.33, labels='17', cex=0.8, col='firebrick2')

dev.off()

# Clean up
rm(list=ls())
gc()

#------------------------------------------------------------------------------

# Plot measure E. coli growth in LB
lb_growth <- as.data.frame(read.delim(file='~/Desktop/repos/Cdiff_modeling/data/tom_96well.ecoli_lb.15-30-46.tsv', sep=',', header=TRUE))
lb_growth$Board_temp <- NULL
lb_growth$Air_temp <- NULL

# Save times
read_times <- lb_growth$UNIX_Timestamp
lb_growth$UNIX_Timestamp <- NULL

# Calculate means of replicates
lb_growth <- data.matrix(lb_growth)
lb_growth[is.na(lb_growth)] <- 0

# Calculate summary statistics
time_means <- as.vector(apply(lb_growth, 1, mean))
time_sds <- as.vector(apply(lb_growth, 1, sd))
time_sd_upper <- time_means+time_sds
time_sd_lower <- time_means-time_sds

# Subset before diauxic growth
time_means <- time_means[1:80]
time_sd_upper <- time_sd_upper[1:80]
time_sd_lower <- time_sd_lower[1:80]

# Caculate slope during exponential phase
slope <- round(diff(time_means)[which.max(diff(time_means))], digits=3)
slope <- paste('m = ', as.character(slope), sep='')

# Generate figure
png(filename='~/Desktop/repos/Cdiff_modeling/results/figure_S4.png', units='in', width=4, height=4, res=300)
par(mar=c(3,3.3,1,1), las=1, mgp=c(2,0.75,0), lwd=2, xaxs='i', yaxs='i')
plot(x=0:(length(time_means)-1), y=time_means, type='l', xlim=c(0,80), xlab='Time (minutes)', xaxt='n',
     ylim=c(0,0.3), ylab='OD600', yaxt='n')
segments(x0=0:(length(time_means)-1), y0=time_sd_upper, y1=time_sd_lower, lwd=10, col='gray')
lines(x=0:(length(time_means)-1), y=time_means, col='firebrick3', lwd=4)
axis(1, at=c(0,20,40,60,80), labels=c(0,60,120,180,240), lwd=2) 
axis(2, at=c(0,0.1,0.2,0.3), lwd=2) 
segments(x0=which.max(diff(time_means))-1, y0=time_means[which.max(diff(time_means))], 
         x1=12, y1=0.15)
text(x=12, y=0.16, labels=slope, cex=1.1)
legend('topleft', legend='E. coli K-12 MG1655 Growth in LB 37C', cex=0.7, pt.cex=0, bty='n')
box()
dev.off()

# Clean up
rm(list=ls())
gc()


