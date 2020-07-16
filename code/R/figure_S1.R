
# Read in data
genres <- read.delim('~/Desktop/repos/Jenior_Cdifficile_2019/data/old_cdifficile_genres/genre_comparisons.tsv', sep='\t', header=TRUE, row.names=1)

# Format quality metrics
genres$annotation_genes <- NULL
genres$annotation_reactions <- NULL
genres$annotation_metabolites <- NULL
genres$yesGPR <- round(100.0 - genres$noGPR,1)
genres$consistency <- round(100.0 - genres$stoichiometric_inconsistency,1)

# Subset data
growth_genres <- genres[order(-genres$growth),]
growth_names <- rownames(growth_genres)
growth_genres <- growth_genres$growth
genres <- genres[,c('memote_score','yesGPR','consistency')]
old_genres <- subset(genres, rownames(genres) %in% c("iMLTC806cdf","icdf834","iHD992","iCN900"))
agora_genres <- subset(genres, rownames(genres) %in% c("agora_NAP07","agora_NAP08","agora_CD196","agora_R20291"))
new_genres <- subset(genres, rownames(genres) %in% c('iCdG698','iCdR700'))
rm(genres)

#--------------------------------------------------------------------------#

# Generate figure
library(wesanderson)
y_pos <- c(4.5,3.5,2.5,1.5)

col_palette <- wes_palette("Chevalier1")
png(filename='~/Desktop/repos/Jenior_Cdifficile_2019/results/figures/figure_S1A.png', width=3, height=4, units='in', res=300)
par(mar=c(2.5,3.5,1.3,0.5), las=1, mgp=c(1.5,0.5,0), xpd=FALSE, lwd=1.7)
barplot(as.matrix(old_genres), xlab='Percentage (%)', beside=TRUE, horiz=TRUE, xlim=c(0,120), 
        xaxt='n', yaxt='n', col=col_palette, cex.lab=1.2, ylim=c(1.4,15.2), lwd=1.7, cex.lab=0.9)
axis(side=1, at=seq(0,100,20), lwd=1.7, cex.axis=0.7)
abline(h=c(5.8,10.8), col='gray40', lty=5)
text(x=60, y=c(15.35,10.35,5.35), cex=0.7, font=2,
     labels=c('Stoichiometric Matrix Consistency','Metabolic Reactions With GPRs','Cumulative MEMOTE Score'))
text(x=rev(old_genres$consistency)+10, y=y_pos+10, labels=rev(paste0(as.character(old_genres$consistency),'%')), cex=0.6) # stoichiometry
text(x=rev(old_genres$yesGPR)+10, y=y_pos+5, labels=rev(paste0(as.character(old_genres$yesGPR),'%')), cex=0.6) # gpr
text(x=rev(old_genres$memote_score)+10, y=y_pos, labels=rev(paste0(as.character(old_genres$memote_score),'%')), cex=0.6) # MEMOTE score
box()
par(xpd=TRUE)
for (x in c(0,5,10)) {text(x=2, y=y_pos+x, labels=c('iCN900','iHD992','icdf834','iMLTC806cdf'), cex=0.7, pos=2)}
text(x=60, y=16.3, labels='Previous Curated GENREs', font=2, cex=0.9)
par(xpd=FALSE)
dev.off()

col_palette <- wes_palette('GrandBudapest1')
png(filename='~/Desktop/repos/Jenior_Cdifficile_2019/results/figures/figure_S1B.png', width=3, height=4, units='in', res=300)
par(mar=c(2.5,3.5,1.3,0.5), las=1, mgp=c(1.5,0.5,0), xpd=FALSE, lwd=1.7)
barplot(as.matrix(agora_genres), xlab='Percentage (%)', beside=TRUE, horiz=TRUE, xlim=c(0,120),
        xaxt='n', yaxt='n', col=col_palette, cex.lab=1.2, ylim=c(1.4,15.2), lwd=1.7, cex.lab=0.9)
axis(side=1, at=seq(0,100,20), lwd=1.7, cex.axis=0.7)
abline(h=c(5.8,10.8), col='gray40', lty=5)
text(x=60, y=c(15.35,10.35,5.35), cex=0.7, font=2,
     labels=c('Stoichiometric Matrix Consistency','Metabolic Reactions With GPRs','Cumulative MEMOTE Score'))
text(x=rev(agora_genres$consistency)+10, y=y_pos+10, labels=rev(paste0(as.character(agora_genres$consistency),'%')), cex=0.6) # stoichiometry
text(x=rev(agora_genres$yesGPR)+10, y=y_pos+5, labels=rev(paste0(as.character(agora_genres$yesGPR),'%')), cex=0.6) # gpr
text(x=rev(agora_genres$memote_score)+10, y=y_pos, labels=rev(paste0(as.character(agora_genres$memote_score),'%')), cex=0.6) # MEMOTE score
box()
par(xpd=TRUE)
for (x in c(0,5,10)) {text(x=2, y=y_pos+x, labels=c('str. R20291','str. CD196','str. NAP08','str. NAP07'), cex=0.7, pos=2)}
text(x=60, y=16.3, labels='AGORA GENREs', font=2, cex=0.9)
par(xpd=FALSE)
dev.off()

col_palette <- c('chartreuse3','blue3')
png(filename='~/Desktop/repos/Jenior_Cdifficile_2019/results/figures/figure_S1C.png', width=3, height=3, units='in', res=300)
par(mar=c(2.5,3.3,1.,0.5), las=1, mgp=c(1.5,0.5,0), xpd=FALSE, lwd=1.7)
barplot(as.matrix(new_genres), xlab='Percentage (%)', beside=TRUE, horiz=TRUE, xlim=c(0,120),
        xaxt='n', yaxt='n', col=col_palette, cex.lab=1.2, ylim=c(0.8,9.4), lwd=1.7, cex.lab=0.9)
axis(side=1, at=seq(0,100,20), lwd=1.7, cex.axis=0.7)
abline(h=c(6.6,3.6), col='gray40', lty=5)
text(x=60, y=c(9.3,6.3,3.3), cex=0.7, font=2,
     labels=c('Stoichiometric Matrix Consistency','Metabolic Reactions With GPRs','Cumulative MEMOTE Score'))
text(x=rev(new_genres$consistency)+10, y=c(7.5,8.5), labels=rev(paste0(as.character(new_genres$consistency),'%')), cex=0.6) # stoichiometry
text(x=rev(new_genres$yesGPR)+10, y=c(4.5,5.5), labels=rev(paste0(as.character(new_genres$yesGPR),'%')), cex=0.6) # gpr
text(x=rev(new_genres$memote_score)+10, y=c(1.5,2.5), labels=rev(paste0(as.character(new_genres$memote_score),'%')), cex=0.6) # MEMOTE score
box()
par(xpd=TRUE)
text(x=-2, y=c(1.5,2.5,4.5,5.5,7.5,8.5), labels=c('iCdR700','iCdG698'), cex=0.7, pos=2)
text(x=60, y=10.2, labels='New GENREs', font=2, cex=0.9)
par(xpd=FALSE)
dev.off()

library(viridis)
col_palette <- viridis(length(growth_genres))
png(filename='~/Desktop/repos/Jenior_Cdifficile_2019/results/figures/figure_S1D.png', width=3, height=3, units='in', res=300)
par(mar=c(4,4,0.5,0.5), las=1, mgp=c(2.5,0.75,0), lwd=1.7, xpd=FALSE)
barplot(growth_genres, ylab='Doubling Time (min)', ylim=c(0,675), col=rev(col_palette))
text(x=seq(0.7,13.7,1.2), y=growth_genres+60, labels=growth_genres, srt=90, cex=0.6)
box()
par(xpd=TRUE)
text(x=seq(1.4,13.4,1.2), y=-30, srt=55, cex=0.7, pos=2,
     labels=c("iCN900","PATRIC draft","AGORA NAP08","AGORA NAP07","iCdG698","icdf834","AGORA R20291","AGORA CD196","iCdR700","iMLTC806cdf","iHD992"))
par(xpd=FALSE)
dev.off()

