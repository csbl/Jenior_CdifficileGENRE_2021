
# Read in data
dimensions <- read.delim('~/Desktop/repos/Jenior_Cdifficile_2019/data/microscopy/quantification.tsv', sep='\t', header=TRUE)

# Subset perimeter
smooth_bhi_perimeter <- as.numeric(as.character(subset(dimensions, variant =='smooth' & condition == 'bhi')$perimeter))
rough_bhi_perimeter <- as.numeric(as.character(subset(dimensions, variant =='rough' & condition == 'bhi')$perimeter))
smooth_gluc_perimeter <- as.numeric(as.character(subset(dimensions, variant =='smooth' & condition == 'gluc')$perimeter))
rough_gluc_perimeter <- as.numeric(as.character(subset(dimensions, variant =='rough' & condition == 'gluc')$perimeter))
smooth_nogluc_perimeter <- as.numeric(as.character(subset(dimensions, variant =='smooth' & condition =='no_gluc')$perimeter))
rough_nogluc_perimeter <- as.numeric(as.character(subset(dimensions, variant =='rough' & condition == 'no_gluc')$perimeter))
rm(dimensions)

# Test differences
bhi_pval <- round(wilcox.test(smooth_bhi_perimeter, rough_bhi_perimeter, exact=FALSE)$p.value, 4)
smooth_pval <- round(wilcox.test(smooth_gluc_perimeter, smooth_nogluc_perimeter, exact=FALSE)$p.value, 4)
rough_pval <- round(wilcox.test(rough_gluc_perimeter, rough_nogluc_perimeter, exact=FALSE)$p.value, 4)
smooth_bhi_gluc_pval <- round(wilcox.test(smooth_bhi_perimeter, smooth_gluc_perimeter, exact=FALSE)$p.value, 4)
smooth_bhi_nogluc_pval <- round(wilcox.test(smooth_bhi_perimeter, smooth_nogluc_perimeter, exact=FALSE)$p.value, 4)
rough_bhi_gluc_pval <- round(wilcox.test(rough_bhi_perimeter, rough_gluc_perimeter, exact=FALSE)$p.value, 4)
rough_bhi_nogluc_pval <- round(wilcox.test(rough_bhi_perimeter, rough_nogluc_perimeter, exact=FALSE)$p.value, 4)

#----------------------------------------------------------------------------------------------------------#

# Generate figure
smooth_col <- 'lightsteelblue2'
rough_col <- 'red3'
library(plotrix)

pdf(file='~/Desktop/repos/Jenior_Cdifficile_2019/results/figures/Figure_3DEF.pdf', width=5.4, height=1.67)
layout(matrix(c(1,2,3), nrow=1, ncol=3, byrow=TRUE))

par(mar=c(2,3,1,1), xpd=FALSE, mgp=c(1.8,0.7,0), lwd=1.5, xaxt='n', las=1)
stripchart(log10(smooth_bhi_perimeter), at=0.25, xlim=c(0,1), ylim=c(3,4.5), bg=smooth_col, vertical=T, cex=1.3,
           cex.lab=0.9, ylab='Colony Perimeter (Log10)', method='jitter', jitter=0.1, pch=21, yaxt='n')
stripchart(log10(rough_bhi_perimeter), at=0.75, bg=rough_col, vertical=T, cex=1.3, 
           method='jitter', jitter=0.1, pch=21, add=TRUE)
axis(side=2, at=c(3,3.5,4,4.5), labels=c(0,3.5,4,4.5), cex.axis=0.8, lwd=1.5)
axis.break(2, 3.1, style='slash')
segments(x0=c(0.1,0.6), x1=c(0.4,0.9), 
         y0=c(median(log10(smooth_bhi_perimeter)),median(log10(rough_bhi_perimeter))))
segments(x0=c(0.15,0.65), x1=c(0.35,0.85), 
         y0=c(quantile(log10(smooth_bhi_perimeter), 0.75),quantile(log10(rough_bhi_perimeter), 0.75)))
segments(x0=c(0.15,0.65), x1=c(0.35,0.85), 
         y0=c(quantile(log10(smooth_bhi_perimeter), 0.25),quantile(log10(rough_bhi_perimeter), 0.25)))
segments(x0=c(0.25,0.75), 
         y0=c(quantile(log10(smooth_bhi_perimeter), 0.25),quantile(log10(rough_bhi_perimeter), 0.25)),
         y1=c(quantile(log10(smooth_bhi_perimeter), 0.75),quantile(log10(rough_bhi_perimeter), 0.75)))
segments(x0=0.25, x1=0.75, y0=4.25, lwd=1.5)
text(x=0.5, y=4.33, '***', font=2, cex=1.3)
par(xpd=TRUE)
text(x=c(0.25,0.75), y=2.77, labels=c('Smooth\nBHIS','Rough\nBHIS'), cex=1)
text(x=-0.3, y=4.61, 'D', cex=1.2, font=2)
par(xpd=FALSE)
text(x=0.75, y=3.25, 'A    B    C', cex=1.2, font=2)

par(mar=c(2,3,1,1), xpd=FALSE, mgp=c(2.1,0.7,0), lwd=1.5, xaxt='n', las=1)
stripchart(smooth_bhi_perimeter, at=0.25, xlim=c(0,1.5), ylim=c(0,5050), bg=smooth_col, vertical=T, cex=1.3,
           cex.lab=1, xlab='', ylab='Colony Perimeter', method='jitter', jitter=0.1, pch=21, yaxt='n')
axis(side=2, at=seq(0,5000,1000), cex.axis=0.6, lwd=1.5)
legend('bottomright', legend='Smooth', cex=1.2, pt.cex=0, bty='n')
box()
stripchart(smooth_gluc_perimeter, at=0.75, bg=smooth_col, vertical=T, cex=1.3, 
           method='jitter', jitter=0.15, pch=21, add=TRUE)
stripchart(smooth_nogluc_perimeter, at=1.25, bg=smooth_col, vertical=T, cex=1.3, 
           method='jitter', jitter=0.15, pch=21, add=TRUE)
segments(x0=c(0.1,0.6,1.1), x1=c(0.4,0.9,1.4), 
         y0=c(median(smooth_bhi_perimeter),median(smooth_gluc_perimeter),median(smooth_nogluc_perimeter)))
segments(x0=c(0.15,0.65,1.15), x1=c(0.35,0.85,1.35), 
         y0=c(quantile(smooth_bhi_perimeter, 0.75),quantile(smooth_gluc_perimeter, 0.75),quantile(smooth_nogluc_perimeter, 0.75)))
segments(x0=c(0.15,0.65,1.15), x1=c(0.35,0.85,1.35), 
         y0=c(quantile(smooth_bhi_perimeter, 0.25),quantile(smooth_gluc_perimeter, 0.25),quantile(smooth_nogluc_perimeter, 0.25)))
segments(x0=c(0.25,0.75,1.25), 
         y0=c(quantile(smooth_bhi_perimeter, 0.25),quantile(smooth_gluc_perimeter, 0.25),quantile(smooth_nogluc_perimeter, 0.25)),
         y1=c(quantile(smooth_bhi_perimeter, 0.75),quantile(smooth_gluc_perimeter, 0.75),quantile(smooth_nogluc_perimeter, 0.75)))
par(xpd=TRUE)
text(x=c(0.25,0.75,1.25), y=-700, labels=c('BHIS','BDM\n+ glucose','BDM\n- glucose'), cex=c(0.9,0.75,0.75))
text(x=-0.45, y=5450, 'E', cex=1.2, font=2)
par(xpd=FALSE)
segments(x0=c(0.25, 0.25, 0.75), x1=c(1.25, 0.75, 1.25), 
         y0=c(4800, 4400, 4000), lwd=1.5)
text(x=c(0.75,0.5,1), y=c(5000, 4600, 4200), '*', cex=1.3, font=2)

par(mar=c(2,3,1,1), xpd=FALSE, mgp=c(2.1,0.7,0), lwd=1.5, xaxt='n', las=1)
stripchart(rough_bhi_perimeter, at=0.25, xlim=c(0,1.5), ylim=c(0,16500), bg=rough_col, vertical=T, cex=1.3,
           cex.lab=1, xlab='', ylab='Colony Perimeter', method='jitter', jitter=0.1, pch=21, yaxt='n')
axis(side=2, at=seq(0,16000,4000), cex.axis=0.6, lwd=1.5)
legend('bottomright', legend='Rough', cex=1.2, pt.cex=0, bty='n')
box()
stripchart(rough_gluc_perimeter, at=0.75, bg=rough_col, vertical=T, cex=1.3, 
           method='jitter', jitter=0.15, pch=21, add=TRUE)
stripchart(rough_nogluc_perimeter, at=1.25, bg=rough_col, vertical=T, cex=1.3, 
           method='jitter', jitter=0.15, pch=21, add=TRUE)
segments(x0=c(0.1,0.6,1.1), x1=c(0.4,0.9,1.4), 
         y0=c(median(rough_bhi_perimeter),median(rough_gluc_perimeter),median(rough_nogluc_perimeter)))
segments(x0=c(0.15,0.65,1.15), x1=c(0.35,0.85,1.35), 
         y0=c(quantile(rough_bhi_perimeter, 0.75),quantile(rough_gluc_perimeter, 0.75),quantile(rough_nogluc_perimeter, 0.75)))
segments(x0=c(0.15,0.65,1.15), x1=c(0.35,0.85,1.35), 
         y0=c(quantile(rough_bhi_perimeter, 0.25),quantile(rough_gluc_perimeter, 0.25),quantile(rough_nogluc_perimeter, 0.25)))
segments(x0=c(0.25,0.75,1.25),
         y0=c(quantile(rough_bhi_perimeter, 0.25),quantile(rough_gluc_perimeter, 0.25),quantile(rough_nogluc_perimeter, 0.25)),
         y1=c(quantile(rough_bhi_perimeter, 0.75),quantile(rough_gluc_perimeter, 0.75),quantile(rough_nogluc_perimeter, 0.75)))
par(xpd=TRUE)
text(x=c(0.25,0.75,1.25), y=-2200, labels=c('BHIS','BDM\n+ glucose','BDM\n- glucose'), cex=c(0.9,0.75,0.75))
text(x=-0.45, y=17800, 'F', cex=1.2, font=2)
par(xpd=FALSE)
segments(x0=c(0.25, 0.25, 0.75), x1=c(1.25, 0.75, 1.25), 
         y0=c(15500, 14250, 13000), lwd=1.5)
text(x=c(0.75,0.5,1), y=c(16150, 14850, 13600), c('n.s','*','*'), cex=c(0.7,1.3,1.3), font=c(1,2,2))

dev.off()
