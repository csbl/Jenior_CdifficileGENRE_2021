

# Read in data
genres <- read.delim('~/Desktop/repos/Jenior_Cdifficile_2019/data/genre_comparisons.tsv', sep='\t', header=TRUE, row.names=1)

# Subset to new GENRE quality metrics
genres <- subset(genres, rownames(genres) %in% c('iCdR758','iCdG791'))
genres[is.na(genres)] <- 0
quality <- genres[,c('noGPR','stoichiometric_inconsistency','free','imbalance')]
quality <- as.matrix(quality)

#--------------------------------------------------------------------------#

# Generate figure
col_palette <- c('#9E2F28','#E7CCAF')
png(filename='~/Desktop/repos/Jenior_Cdifficile_2019/results/figures/figure_S1.png', width=5, height=5, units='in', res=300)
layout(matrix(c(1,2,
                1,3), nrow=2, ncol=2, byrow=TRUE))

# GENRE quality metrics
par(mar=c(3,5,0.5,1), las=1, mgp=c(1.9,0.75,0), xpd=FALSE, lwd=1.5)
barplot(quality, xlab='Fraction of Total (%)', beside=TRUE, horiz=TRUE, xlim=c(0,100), 
        xaxt='n', yaxt='n', col=col_palette, ylim=c(1,12), lwd=1.5)
axis(side=1, at=seq(0,100,50), lwd=2)
abline(h=c(3.5,6.5,9.5), col='gray40', lty=5)
par(lwd=2)
legend('topright', legend=rownames(quality), pch=22, pt.bg=rev(col_palette), pt.cex=2.2, cex=1.2, pt.lwd=1.5)
box()
# mass imbalance
text(x=14, y=11.5, labels='0.0%', font=2, cex=1.1)
text(x=14, y=10.5, labels='0.0%', font=2, cex=1.1)
# free metabolites
text(x=14, y=8.5, labels='0.0%', font=2, cex=1.1)
text(x=14, y=7.5, labels='0.0%', font=2, cex=1.1)
# stoichiometry
text(x=14, y=5.5, labels='0.0%', font=2, cex=1.1)
text(x=14, y=4.5, labels='0.0%', font=2, cex=1.1)
# no gpr
text(x=39.5, y=2.5, labels='20.48%', font=2, cex=1.1)
text(x=38, y=1.5, labels='20.7%', font=2, cex=1.1)
par(xpd=TRUE)
text(x=-28, y=11, labels='Mass\nImbalanced\nReactions', cex=0.9)
text(x=-28, y=8, labels='Freely\nProduced\nMetabolites', cex=0.9)
text(x=-28, y=5, labels='Stoichiometric\nInconsistency', cex=0.9)
text(x=-28, y=2, labels='Metabolic\nReactions\nWithout\nGPRs', cex=0.9)
text(x=-45, y=12.4, labels='A', font=2, cex=1.5)
par(xpd=FALSE)

# Imputed doubling time in complete media
par(mar=c(5.5,4,1,1), las=1, mgp=c(2.5,0.75,0), lwd=1.5, xpd=FALSE)
barplot(genres$growth, ylab='Doubling Time (min)', ylim=c(0,60), col=rev(col_palette), yaxt='n')
abline(h=39, col='azure4', lty=3, lwd=1.2)
barplot(rev(genres$growth), ylab='Doubling Time (min)', ylim=c(0,60), col=rev(col_palette), yaxt='n', add=TRUE)
axis(side=2, at=seq(0,60,10), lwd=2)
box(lwd=2)
text(x=0.7, y=44, labels='40.10', font=2, cex=1.1)
text(x=1.9, y=44, labels='40.02', font=2, cex=1.1)
par(xpd=TRUE)
text(x=0.6, y=-14, labels='iCdG791', srt=55, cex=1.1)
text(x=1.8, y=-14, labels='iCdR758', srt=55, cex=1.1)
text(x=-0.7, y=63, labels='B', font=2, cex=1.5)
par(xpd=FALSE)

# MEMOTE score
par(mar=c(5.5,4,1,1), las=1, mgp=c(2.5,0.75,0), lwd=1.5, xpd=FALSE)
barplot(genres$memote_score, ylab='MEMOTE Score (%)', ylim=c(0,100), col=rev(col_palette), yaxt='n')
axis(side=2, at=seq(0,100,20), lwd=2)
box(lwd=2)
text(x=0.7, y=90.5, labels='85', font=2, cex=1.1)
text(x=1.9, y=89.5, labels='84', font=2, cex=1.1)
par(xpd=TRUE)
text(x=0.6, y=-23, labels='iCdG791', srt=55, cex=1.1)
text(x=1.8, y=-23, labels='iCdR758', srt=55, cex=1.1)
text(x=-0.7, y=105, labels='C', font=2, cex=1.5)
par(xpd=FALSE)

dev.off()
