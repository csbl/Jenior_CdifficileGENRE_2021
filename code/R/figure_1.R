
# Read in data
genres <- read.delim('~/Desktop/repos/Jenior_Cdifficile_2019/data/old_cdifficile_genres/genre_comparisons.tsv', sep='\t', header=TRUE, row.names=1)

# Subset to previous GENRE quality metrics
genres$annotation_genes <- NULL
genres$annotation_reactions <- NULL
genres$annotation_metabolites <- NULL

genres <- subset(genres, rownames(genres) %in% c("iMLTC806cdf","icdf834","iHD992","iCN900"))
quality <- genres[,c('memote_score','noGPR','stoichiometric_inconsistency','free','imbalance')]
quality[is.na(quality)] <- 0
quality <- as.matrix(quality)

#--------------------------------------------------------------------------#

# Generate figure
#col_palette <- c('#7C5B76','#31688EFF','#35B779FF','#DC7C40')
library(wesanderson)
col_palette <- wes_palette("Chevalier1")
y_pos <- c(4.5,3.5,2.5,1.5)
png(filename='~/Desktop/repos/Jenior_Cdifficile_2019/results/figures/figure_1.png', width=3, height=6, units='in', res=300)
par(mar=c(2.5,3.8,0.5,0.5), las=1, mgp=c(1.5,0.5,0), xpd=FALSE, lwd=1.7)
barplot(quality, xlab='Ratio (%)', beside=TRUE, horiz=TRUE, xlim=c(0,105), 
        xaxt='n', yaxt='n', col=col_palette, cex.lab=1.2, ylim=c(1.4,24.9), lwd=1.7, cex.lab=0.9)
axis(side=1, at=seq(0,100,20), lwd=1.7, cex.axis=0.7)
abline(h=c(5.8,10.8,15.8,20.8), col='gray40', lty=5)
text(x=52.5, y=c(25.35,20.35,15.35,10.35,5.35), cex=0.7, font=2,
     labels=c('Mass Imbalanced Reactions','Freely Produced Metabolites','Stoichiometric Inconsistency','No GPR Metabolic Reactions','Cumulative MEMOTE Score'))
text(x=c(25,30,26,28), y=y_pos+20, labels=c('10.74%','No formula data','12.33%','15.3%'), cex=0.6) # mass imbalance
text(x=c(10,59,12,13), y=y_pos+15, labels=c('0.0%','74.12%','0.25%','0.68%'), cex=0.6) # free metabolites
text(x=c(70,77,82,11), y=y_pos+10, labels=c('83.1%','64.3%','69.3%','1.0%'), cex=0.6) # stoichiometry
text(x=c(39,88,40,88), y=y_pos+5, labels=c('23.81%','100%','24.71%','100%'), cex=0.6) # no gpr
text(x=c(83,25,35,50), y=y_pos, labels=c('74%','16%','26%','41%'), cex=0.6) # MEMOTE score
box()
par(xpd=TRUE)
for (x in c(0,5,10,15,20)) {text(x=2, y=y_pos+x, labels=c('iCN900','iHD992','icdf834','iMLTC806cdf'), cex=0.7, pos=2)}
par(xpd=FALSE)
dev.off()
