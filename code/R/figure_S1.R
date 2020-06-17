
# Read in data
all_genres <- read.delim('~/Desktop/repos/Jenior_Cdifficile_2019/data/old_cdifficile_genres/genre_comparisons.tsv', sep='\t', header=TRUE, row.names=1)

# Subset to new GENRE quality metrics
genres <- subset(all_genres, rownames(all_genres) %in% c('iCdG698','iCdR700',"agora_NAP07","agora_NAP08","agora_CD196","agora_R20291"))
genres[is.na(genres)] <- 0
quality <- genres[,c('noGPR','stoichiometric_inconsistency','free','imbalance')]
quality <- as.matrix(quality)
rownames(quality) <- c('iCdG698 (630)','iCdR700 (R20291)',"AGORA (NAP07)","AGORA (NAP08)","AGORA (CD196)","AGORA (R20291)")
#--------------------------------------------------------------------------#

# Generate figure
library(viridis)
col_palette <- viridis(6)
png(filename='~/Desktop/repos/Jenior_Cdifficile_2019/results/figures/figure_S1.png', width=6, height=5, units='in', res=300)
layout(matrix(c(1,1,2,
                1,1,3), nrow=2, ncol=3, byrow=TRUE))

# GENRE quality metrics
par(mar=c(3,6,0.5,1), las=1, mgp=c(1.9,0.75,0), xpd=FALSE, lwd=1.5)
barplot(quality, xlab='Fraction of Total (%)', beside=TRUE, horiz=TRUE, xlim=c(0,100), 
        xaxt='n', yaxt='n', col=col_palette, ylim=c(1.5,27.5), lwd=1.5)
axis(side=1, at=seq(0,100,50), lwd=2)
abline(h=c(7.5,14.5,21.5), col='gray40', lty=5)
par(lwd=2)
legend('topright', legend=rownames(quality), pch=22, pt.bg=rev(col_palette), pt.cex=1.9, pt.lwd=1.2)
box()
par(xpd=TRUE)
text(x=-13, y=24.5, labels='Mass\nImbalanced\nReactions')
text(x=-13, y=17.5, labels='Freely\nProduced\nMetabolites')
text(x=-13, y=10.5, labels='Stoichiometric\nInconsistency')
text(x=-13, y=3.5, labels='Metabolic\nReactions\nWithout\nGPRs')
par(xpd=FALSE)
# mass imbalance
text(x=12, y=seq(22.5,27.5,1), labels=c('0.0%','0.0%','0.24%','0.24%','0.44%','0.57%'), cex=0.8)
# free metabolites
text(x=9, y=seq(15.5,20.5,1), labels=c('0.0%','0.0%','0.0%','0.0%','0.0%','0.0%'), cex=0.8)
# stoichiometry
text(x=12, y=seq(8.5,13.5,1), labels=c('0.0%','0.0%','0.2%','0.1%','0.1%','0.2%'), cex=0.8)
# no gpr
text(x=35, y=seq(1.5,6.5,1), labels=c('15.78%','15.98%','11.22%','11.19%','16.58%','21.74%'), cex=0.8)
mtext('A', side=2, font=2, cex=1, adj=5, padj=-18.5)

# MEMOTE score
par(mar=c(5.5,4,0.5,0.5), las=1, mgp=c(2.5,0.75,0), lwd=1.5, xpd=FALSE)
barplot(genres$memote_score, ylab='MEMOTE Score (%)', ylim=c(0,100), col=rev(col_palette), yaxt='n')
axis(side=2, at=seq(0,100,20), lwd=2)
box(lwd=2)
text(x=seq(0.7,6.7,1.2), y=c(91,91,45,45,46,48), labels=c('86%','86%','40%','40%','41%','43%'), cex=0.8)
par(xpd=TRUE)
text(x=seq(0.6,6.6,1.2), y=-12, labels=c('iCdG698','iCdR700',"NAP07","NAP08","CD196","R20291"), srt=55, cex=0.9)
text(x=4.75, y=-30, labels='AGORA', font=2)
segments(x0=2.6,x1=7, y0=-25, lwd=2)
par(xpd=FALSE)
mtext('B', side=2, font=2, cex=1, adj=4, padj=-6.5)

# All GENRE doubling times
col_palette <- viridis(nrow(growth_genres))
growth_genres <- subset(all_genres, rownames(all_genres) %in% c("iCdG698","iCdR700","iMLTC806cdf","icdf834","iHD992","iCN900","agora_NAP07","agora_NAP08","agora_CD196","agora_R20291"))
growth_genres <- all_genres[c("iCdG698","iCdR700","iMLTC806cdf","icdf834","iHD992","iCN900","agora_NAP07","agora_NAP08","agora_CD196","agora_R20291"),]
rownames(growth_genres) <- c("iCdG698","iCdR700","iMLTC806cdf","icdf834","iHD992","iCN900","AGORA NAP07","AGORA NAP08","AGORA CD196","AGORA R20291")
par(mar=c(5.5,4,0.5,0.5), las=1, mgp=c(2.5,0.75,0), lwd=1.7, xpd=FALSE)
barplot(growth_genres$growth, ylab='Imputed Doubling Time (min)', ylim=c(0,675), col=rev(col_palette))
text(x=seq(0.7,11.5,1.2), y=c(100,100,50,100,50,620,100,100,100,100), 
     labels=growth_genres$growth, srt=90, cex=0.8)
box()
par(xpd=TRUE)
text(x=seq(1.4,12.2,1.2), y=-30, labels=rownames(growth_genres), srt=55, cex=0.9, pos=2)
par(xpd=FALSE)
mtext('C', side=2, font=2, cex=1, adj=4, padj=-6.5)


dev.off()
