

# Read in data
genres <- read.delim('~/Desktop/repos/Jenior_Cdifficile_2019/data/genre_comparisons.tsv', sep='\t', header=TRUE, row.names=1)

# Subset to previous GENRE quality metrics
genres <- subset(genres, !rownames(genres) %in% c('cd630_PATRIC_draft','iCdG791','iCdR758'))
quality <- genres[,c('noGPR','stoichiometric_inconsistency','free','imbalance')]
quality[is.na(quality)] <- 0
quality <- as.matrix(quality)

#--------------------------------------------------------------------------#

# Generate figure
col_palette <- c('#7C5B76','#31688EFF','#35B779FF','#DC7C40')
png(filename='~/Desktop/repos/Jenior_Cdifficile_2019/results/figures/figure_1.png', width=7, height=4, units='in', res=300)
layout(matrix(c(1,1,2,
                1,1,3), nrow=2, ncol=3, byrow=TRUE))

# GENRE quality metrics
par(mar=c(3,7,0.5,1), las=1, mgp=c(1.9,0.75,0), xpd=FALSE, lwd=1.5)
barplot(quality, xlab='Fraction of Total (%)', beside=TRUE, horiz=TRUE, xlim=c(0,100), 
        xaxt='n', yaxt='n', col=col_palette, cex.lab=1.2, ylim=c(1.2,19.8), lwd=1.5)
axis(side=1, at=seq(0,100,20), lwd=2)
abline(h=c(5.5,10.5,15.5), col='gray40', lty=5)
par(lwd=2)
legend('topright', legend=rev(rownames(quality)), pch=22, pt.bg=rev(col_palette), pt.cex=2.2, cex=1.2, pt.lwd=1.5)
box()
# mass imbalance
text(x=18, y=19.5, labels='10.74%', cex=1.1)
text(x=12, y=18.5, labels='No formula data', col='firebrick')
text(x=20, y=17.5, labels='12.33%', cex=1.1)
text(x=22, y=16.5, labels='15.3%', cex=1.1)
# free metabolites
#text(x=5, y=14.5, labels='0.0%', cex=1.1)
#text(x=81, y=13.5, labels='74.12%', cex=1.1)
#text(x=6, y=12.5, labels='0.25%', cex=1.1)
#text(x=7, y=11.5, labels='0.68%', cex=1.1)
# stoichiometry
text(x=90, y=9.5, labels='83.1%', cex=1.1)
text(x=71, y=8.5, labels='64.3%', cex=1.1)
text(x=76, y=7.5, labels='69.3%', cex=1.1)
text(x=7, y=6.5, labels='1.0%', cex=1.1)
# no gpr
text(x=31, y=4.5, labels='23.81%', cex=1.1)
text(x=94, y=3.5, labels='100%', cex=1.1, col='white')
text(x=32, y=2.5, labels='24.71%', cex=1.1)
text(x=94, y=1.5, labels='100%', cex=1.1, col='white')
par(xpd=TRUE)
text(x=-13, y=18, labels='Mass\nImbalanced\nReactions', cex=1.2, font=2)
#text(x=-13, y=13, labels='Freely\nProduced\nMetabolites', cex=1.2, font=2)
text(x=-13, y=8, labels='Overall\nStoichiometric\nInconsistency', cex=1.2, font=2)
text(x=-13, y=3, labels='Metabolic\nReactions\nWithout\nGPRs', cex=1.2, font=2)
text(x=-22, y=20.2, labels='A', font=2, cex=1.5)
par(xpd=FALSE)

# Imputed doubling time in complete media
par(mar=c(5.5,4,1,1), las=1, mgp=c(2.5,0.75,0), lwd=1.5, xpd=FALSE)
barplot(rev(genres$growth), ylab='Doubling Time (min)', ylim=c(0,600), col=rev(col_palette), yaxt='n')
abline(h=39, col='azure4', lty=3, lwd=1.2)
barplot(rev(genres$growth), ylab='Doubling Time (min)', ylim=c(0,600), col=rev(col_palette), yaxt='n', add=TRUE)
axis(side=2, at=seq(0,600,100), lwd=2)
box(lwd=2)
text(x=0.7, y=530, labels='576', font=2, cex=1.1, col='white')
text(x=1.9, y=70, labels='3.6', font=2, cex=1.1)
text(x=3.1, y=100, labels='24.57', font=2, cex=1.1)
text(x=4.3, y=70, labels='3.6', font=2, cex=1.1)
par(xpd=TRUE)
text(x=0.6, y=-110, labels='iCN900', srt=55, cex=1.1)
text(x=1.7, y=-110, labels='iHD992', srt=55, cex=1.1)
text(x=2.9, y=-110, labels='icdf834', srt=55, cex=1.1)
text(x=4.1, y=-180, labels='iMLTC806cdf', srt=55, cex=1.1)
text(x=-1.35, y=605, labels='B', font=2, cex=1.5)
par(xpd=FALSE)

# MEMOTE score
par(mar=c(5.5,4,1,1), las=1, mgp=c(2.5,0.75,0), lwd=1.5, xpd=FALSE)
barplot(rev(genres$memote_score), ylab='MEMOTE Score (%)', ylim=c(0,100), col=rev(col_palette), yaxt='n')
axis(side=2, at=seq(0,100,20), lwd=2)
box(lwd=2)
text(x=0.7, y=81, labels='74', font=2, cex=1.1)
text(x=1.9, y=23, labels='16', font=2, cex=1.1)
text(x=3.1, y=34, labels='26', font=2, cex=1.1)
text(x=4.3, y=48, labels='41', font=2, cex=1.1)
par(xpd=TRUE)
text(x=0.6, y=-19, labels='iCN900', srt=55, cex=1.1)
text(x=1.7, y=-19, labels='iHD992', srt=55, cex=1.1)
text(x=2.9, y=-19, labels='icdf834', srt=55, cex=1.1)
text(x=4.1, y=-30, labels='iMLTC806cdf', srt=55, cex=1.1)
text(x=-1.35, y=102, labels='C', font=2, cex=1.5)
par(xpd=FALSE)

dev.off()
