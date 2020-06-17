# Read in data
growth <- read.delim('~/Desktop/repos/Jenior_Cdifficile_2019/data/biolog_sim.values.tsv', sep='\t', header=TRUE)
growth$name <- gsub('_',' ',growth$name)
#binary <- read.delim('~/Desktop/repos/Jenior_Cdifficile_2019/data/biolog_sim.binary.tsv', sep='\t', header=TRUE)
#binary$name <- gsub('_',' ',binary$name)

# Separate strains
cd630 <- growth[,c(1:5)]
cdR20291 <- growth[,c(1:3,6:7)]
rm(growth)

# Restructure
colnames(cd630) <- c('modelseed','name','type','fba','biolog')
cd630$diff <- abs(cd630$fba - cd630$biolog)
cd630_carb <- subset(cd630, type == 'Carbohydrate')
cd630_carb <- cd630_carb[order(cd630_carb$fba),]
cd630_carb$type <- rep('a_Carbohydrate', nrow(cd630_carb))
cd630_nuc <- subset(cd630, type == 'Nucleotide')
cd630_nuc <- cd630_nuc[order(cd630_nuc$fba),]
cd630_nuc$type <- rep('b_Nucleotide', nrow(cd630_nuc))
cd630_aa <- subset(cd630, type == 'Amino_acid')
cd630_aa <- cd630_aa[order(cd630_aa$fba),]
cd630_aa$type <- rep('c_Amino_acid', nrow(cd630_aa))
cd630_carboxy <- subset(cd630, type == 'Carboxylic_acid')
cd630_carboxy <- cd630_carboxy[order(cd630_carboxy$fba),]
cd630_carboxy$type <- rep('d_Carboxylic_acid', nrow(cd630_carboxy))
cd630_other <- subset(cd630, !type %in% c('Carbohydrate','Amino_acid','Nucleotide','Carboxylic_acid'))
cd630_other <- cd630_other[order(cd630_other$fba),]
cd630_other$type <- rep('e_Other', nrow(cd630_other))
#cd630 <- as.data.frame(rbind(cd630_other,cd630_carboxy,cd630_nuc,cd630_aa,cd630_carb))
cd630 <- as.data.frame(rbind(cd630_carb,cd630_aa,cd630_nuc,cd630_carboxy,cd630_other))
cd630$type <- as.factor(cd630$type)

colnames(cdR20291) <- c('modelseed','name','type','fba','biolog')
cdR20291$diff <- abs(cdR20291$fba - cdR20291$biolog)
cdR20291_carb <- subset(cdR20291, type == 'Carbohydrate')
cdR20291_carb <- cdR20291_carb[order(cdR20291_carb$fba),]
cdR20291_carb$type <- rep('a_Carbohydrate', nrow(cdR20291_carb))
cdR20291_nuc <- subset(cdR20291, type == 'Nucleotide')
cdR20291_nuc <- cdR20291_nuc[order(cdR20291_nuc$fba),]
cdR20291_nuc$type <- rep('b_Nucleotide', nrow(cdR20291_nuc))
cdR20291_aa <- subset(cdR20291, type == 'Amino_acid')
cdR20291_aa <- cdR20291_aa[order(cdR20291_aa$fba),]
cdR20291_aa$type <- rep('c_Amino_acid', nrow(cdR20291_aa))
cdR20291_carboxy <- subset(cdR20291, type == 'Carboxylic_acid')
cdR20291_carboxy <- cdR20291_carboxy[order(cdR20291_carboxy$fba),]
cdR20291_carboxy$type <- rep('d_Carboxylic_acid', nrow(cdR20291_carboxy))
cdR20291_other <- subset(cdR20291, !type %in% c('Carbohydrate','Amino_acid','Nucleotide','Carboxylic_acid'))
cdR20291_other <- cdR20291_other[order(cdR20291_other$fba),]
cdR20291_other$type <- rep('e_Other', nrow(cdR20291_other))
cdR20291 <- as.data.frame(rbind(cdR20291_carb,cdR20291_aa,cdR20291_nuc,cdR20291_carboxy,cdR20291_other))
cdR20291$type <- as.factor(cdR20291$type)


# Generate figure
png(filename='~/Desktop/repos/Jenior_Cdifficile_2019/results/figures/figure_S2.png', units='in', width=14, height=15, res=300)
layout(matrix(c(1,2), nrow=1, ncol=2, byrow=TRUE))

par(mar=c(3,4,2,0.5), las=1, mgp=c(2,0.75,0), lwd=2)
dotchart(cd630$fba, labels=cd630$name, groups=cd630$type, gcolor='white', cex=0.8, cex.lab=1.2, cex.axis=1.1,
         xlab='Growth Enhancement Ratio', xlim=c(0.9,2.8), pch=21, pt.cex=1.5, bg='blue3', main='iCdG692')
abline(v=1)
points(x=cd630_carb$biolog, y=c((124-nrow(cd630_carb)):123), pch=1, lwd=2.5, col='chocolate2', cex=1.5)
points(x=cd630_nuc$biolog, y=c((96-nrow(cd630_nuc)):95), pch=1, lwd=2.5, col='chocolate2', cex=1.5)
points(x=cd630_aa$biolog, y=c((80-nrow(cd630_aa)):79), pch=1, lwd=2.5, col='chocolate2', cex=1.5)
points(x=cd630_carboxy$biolog, y=c((32-nrow(cd630_carboxy)):31), pch=1, lwd=2.5, col='chocolate2', cex=1.5)
points(x=cd630_other$biolog, y=c((7-nrow(cd630_other)):6), pch=1, lwd=2.5, col='chocolate2', cex=1.5)
legend('bottomright', legend=c('Predicted','Measured'), pch=21,
       pt.bg=c('blue3','white'), col=c('black','chocolate2'),
       pt.cex=2, pt.lwd=c(2,2.5), lwd=0, cex=1.5)
par(xpd=TRUE)
text(x=-0.075, y=110.5, 'Carbohydrates', cex=1.5, srt=90)
text(x=-0.075, y=88.5, 'Nucleotides', cex=1.5, srt=90)
text(x=-0.075, y=56.5, 'Amino acids', cex=1.5, srt=90)
text(x=-0.075, y=20, 'Carboxylic acids', cex=1.5, srt=90)
text(x=-0.075, y=3.5, 'Other', cex=1.5, srt=90)
par(xpd=FALSE)

dotchart(cdR20291$fba, labels=cdR20291$name, groups=cdR20291$type, gcolor='white', cex=0.8, cex.lab=1.2, cex.axis=1.1,
         xlab='Growth Enhancement Ratio', xlim=c(0.9,2.8), pch=21, pt.cex=1.5, bg='blue3', main='iCdR700')
abline(v=1)
points(x=cdR20291_carb$biolog, y=c((124-nrow(cdR20291_carb)):123), pch=1, lwd=2.5, col='chocolate2', cex=1.5)
points(x=cdR20291_nuc$biolog, y=c((96-nrow(cdR20291_nuc)):95), pch=1, lwd=2.5, col='chocolate2', cex=1.5)
points(x=cdR20291_aa$biolog, y=c((80-nrow(cdR20291_aa)):79), pch=1, lwd=2.5, col='chocolate2', cex=1.5)
points(x=cdR20291_carboxy$biolog, y=c((32-nrow(cdR20291_carboxy)):31), pch=1, lwd=2.5, col='chocolate2', cex=1.5)
points(x=cdR20291_other$biolog, y=c((7-nrow(cdR20291_other)):6), pch=1, lwd=2.5, col='chocolate2', cex=1.5)
legend('bottomright', legend=c('Predicted','Measured'), pch=21,
       pt.bg=c('blue3','white'), col=c('black','chocolate2'),
       pt.cex=2, pt.lwd=c(2,2.5), lwd=0, cex=1.5)
par(xpd=TRUE)
text(x=-0.075, y=110.5, 'Carbohydrates', cex=1.5, srt=90)
text(x=-0.075, y=88.5, 'Nucleotides', cex=1.5, srt=90)
text(x=-0.075, y=56.5, 'Amino acids', cex=1.5, srt=90)
text(x=-0.075, y=20, 'Carboxylic acids', cex=1.5, srt=90)
text(x=-0.075, y=3.5, 'Other', cex=1.5, srt=90)
par(xpd=FALSE)

dev.off()

#------------------------------------------------------------------------------------------------------------#

# Clean environment
rm(list=ls())
gc()

# Read in data
binary <- read.delim('~/Desktop/repos/Jenior_Cdifficile_2019/data/biolog_sim.binary.tsv', sep='\t', header=TRUE)
binary$name <- gsub('_',' ',binary$name)

# Separate strains
cd630 <- binary[,c(1:5)]
cdR20291 <- binary[,c(1:3,6:7)]
rm(binary)








