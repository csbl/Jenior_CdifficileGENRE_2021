
# Read in data
sporulation <- read.delim('~/Desktop/repos/Jenior_Cdifficile_2019/data/630_sporulation.tsv', sep='\t', header=TRUE)

# Subset to experimental groups
bdm_veg <- as.vector(sporulation[,'bdm_veg'])
bdm_gluc_veg <- as.vector(sporulation[,'bdm_gluc_veg'])
bdm_neu5Ac_veg <- as.vector(sporulation[,'bdm_neu5Ac_veg'])
bdm_cyt_veg <- as.vector(sporulation[,'bdm_cyt_veg'])
bdm_neu5Ac_cyt_veg <- as.vector(sporulation[,'bdm_neu5Ac_cyt_veg'])
bdm_spore <- as.vector(sporulation[,'bdm_spore'])
bdm_gluc_spore <- as.vector(sporulation[,'bdm_gluc_spore'])
bdm_neu5Ac_spore <- as.vector(sporulation[,'bdm_neu5Ac_spore'])
bdm_cyt_spore <- as.vector(sporulation[,'bdm_cyt_spore'])
bdm_neu5Ac_cyt_spore <- as.vector(sporulation[,'bdm_neu5Ac_cyt_spore'])
rm(sporulation)

# Test differences
bdm_spore_pvals <- p.adjust(c(wilcox.test(bdm_spore, bdm_gluc_spore, exact=FALSE)$p.value, 
                        wilcox.test(bdm_spore, bdm_neu5Ac_spore, exact=FALSE)$p.value,
                        wilcox.test(bdm_spore, bdm_cyt_spore, exact=FALSE)$p.value,
                        wilcox.test(bdm_spore, bdm_neu5Ac_cyt_spore, exact=FALSE)$p.value, ), method='BH')
gluc_spore_pvals <- p.adjust(c(wilcox.test(bdm_gluc_spore, bdm_spore, exact=FALSE)$p.value, 
                         wilcox.test(bdm_gluc_spore, bdm_neu5Ac_spore, exact=FALSE)$p.value,
                         wilcox.test(bdm_gluc_spore, bdm_cyt_spore, exact=FALSE)$p.value,
                         wilcox.test(bdm_gluc_spore, bdm_neu5Ac_cyt_spore, exact=FALSE)$p.value, ), method='BH')
neu5Ac_spore_pvals <- p.adjust(c(wilcox.test(bdm_neu5Ac_spore, bdm_spore, exact=FALSE)$p.value, 
                           wilcox.test(bdm_neu5Ac_spore, bdm_cyt_spore, exact=FALSE)$p.value,
                           wilcox.test(bdm_neu5Ac_spore, bdm_neu5Ac_cyt_spore, exact=FALSE)$p.value, ), method='BH')
cyt_spore_pvals <- p.adjust(c(wilcox.test(bdm_cyt_spore, bdm_spore, exact=FALSE)$p.value, 
                        wilcox.test(bdm_cyt_spore, bdm_neu5Ac_spore, exact=FALSE)$p.value,
                        wilcox.test(bdm_cyt_spore, bdm_neu5Ac_cyt_spore, exact=FALSE)$p.value, ), method='BH')
neu5Ac_cyt_spore_pvals <- p.adjust(c(wilcox.test(bdm_neu5Ac_cyt_spore, bdm_spore, exact=FALSE)$p.value, 
                        wilcox.test(bdm_neu5Ac_cyt_spore, bdm_neu5Ac_spore, exact=FALSE)$p.value,
                        wilcox.test(bdm_neu5Ac_cyt_spore, bdm_cyt_spore, exact=FALSE)$p.value, ), method='BH')

# Generate figures
veg_col <- 'green3'
spore_col <- 'goldenrod3'



pdf(file='~/Desktop/repos/Jenior_Cdifficile_2019/results/figures/Figure_5.pdf', width=6.5, height=5.5)
layout(matrix(c(1,2), nrow=1, ncol=4, byrow=TRUE))

par(mar=c(5,3,0.5,0.5), xpd=FALSE, las=1, mgp=c(1.6,0.7,0), lwd=2)
boxplot(smooth_biomass, at=0.5, xlim=c(0,2), ylab='Vegetative CFU', 
        ylim=c(0,60), yaxt='n', col=smooth_col, outline=FALSE, cex.lab=1.2,
        boxlwd=2, medlwd=2, staplelwd=2, whisklwd=2, whisklty=1, width=0.5)
boxplot(rough_biomass, at=1.5, add=TRUE, yaxt='n', col=rough_col,
        boxlwd=2, medlwd=2, staplelwd=2, whisklwd=2, whisklty=1, width=0.5, outline=FALSE)
axis(side=2, at=seq(0,60,10), cex.axis=0.9, lwd=1.7)



segments(x0=0.5, y0=56, x1=1.5, lwd=2)
text(x=1, y=57.5, '***', cex=1.5, font=2)
par(xpd=TRUE)
text(x=c(0.5,1.5), y=-8, labels=c('Smooth','Rough'), cex=1.3)

par(xpd=FALSE)





par(mar=c(5,3,0.5,0.5), xpd=FALSE, las=1, mgp=c(1.6,0.7,0), lwd=2)
boxplot(smooth_biomass, at=0.5, xlim=c(0,2), ylab='Spore CFU', 
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
par(xpd=FALSE)





dev.off()
