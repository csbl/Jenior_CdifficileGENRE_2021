
# Read in data

#growth <- read.delim('~/Desktop/repos/Jenior_Cdifficile_2019/data/biolog_sim.values.tsv', sep='\t', header=TRUE)
growth <- read.delim('~/Desktop/repos/Jenior_Cdifficile_2019/data/new_biolog_sim.values.tsv', sep='\t', header=TRUE)
growth$name <- gsub('_',' ',growth$name)

cd630 <- growth[,c(1:5)]
cdR20291 <- growth[,c(1:3,6:7)]
rm(growth)

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

# Normalize data
cd630$fba_norm <- cd630$fba / min(cd630$biolog)
cd630$biolog_norm <- cd630$biolog / min(cd630$fba)
cd630$fba <- cd630$fba_norm
cd630$fba_norm <- NULL
cd630$biolog <- cd630$biolog
cd630$biolog_norm <- NULL
cdR20291$fba_norm <- cdR20291$fba / min(cdR20291$biolog)
cdR20291$biolog_norm <- cdR20291$biolog / min(cdR20291$fba)
cdR20291$fba <- cdR20291$fba_norm
cdR20291$fba_norm <- NULL
cdR20291$biolog <- cdR20291$biolog
cdR20291$biolog_norm <- NULL

# Calculate correlations
# Overall
cor.test(x=cd630$fba, y=cd630$biolog, method='spearman', exact=FALSE) # R = 0.418 p-value = 0.000003 ***
cor.test(x=cdR20291$fba, y=cdR20291$biolog, method='spearman', exact=FALSE) # R = 0.326, p-value = 0.0004 ***
# Metabolite classes
cor.test(x=cdR20291_aa$fba, y=cdR20291_aa$biolog, method='spearman', exact=FALSE) # R = -0.198, p-value = 0.1876
cor.test(x=cdR20291_carb$fba, y=cdR20291_carb$biolog, method='spearman', exact=FALSE) # R = 0.558, p-value = 0.003 **
cor.test(x=cdR20291_carboxy$fba, y=cdR20291_carboxy$biolog, method='spearman', exact=FALSE) # R = -0.044, p-value = 0.842
cor.test(x=cdR20291_nuc$fba, y=cdR20291_nuc$biolog, method='spearman', exact=FALSE) # R = 0.237, p-value = 0.414
cor.test(x=cdR20291_other$fba, y=cdR20291_other$biolog, method='spearman', exact=FALSE) # R = 0.393, p-value = 0.441
cor.test(x=cdR20291_aa$fba, y=cdR20291_aa$biolog, method='spearman', exact=FALSE) # R = -0.198, p-value = 0.188
cor.test(x=cdR20291_carb$fba, y=cdR20291_carb$biolog, method='spearman', exact=FALSE) # R = 0.557, p-value = 0.003 **
cor.test(x=cdR20291_carboxy$fba, y=cdR20291_carboxy$biolog, method='spearman', exact=FALSE) # R = -0.044, p-value = 0.842
cor.test(x=cdR20291_nuc$fba, y=cdR20291_nuc$biolog, method='spearman', exact=FALSE) # R = 0.237, p-value = 0.414
cor.test(x=cdR20291_other$fba, y=cdR20291_other$biolog, method='spearman', exact=FALSE) # R = 0.393, p-value = 0.441

# Generate figure
library(viridis)
pall <- viridis(5)

pdf(file='~/Desktop/repos/Jenior_Cdifficile_2019/results/figures/figure_1.pdf', width=8, height=4)
layout(matrix(c(1,2), nrow=1, ncol=2, byrow=TRUE))

par(mar=c(3,3,0.5,0.5), las=1, mgp=c(1.9,0.75,0), lwd=1.5)
plot(x=cd630$fba, y=cd630$biolog, cex=0, xlab='Simulated Growth Enhancement Ratio', ylab='Experimental Growth Enhancement Ratio', 
     xlim=c(0.9,4.9), ylim=c(0.8,2.5), xaxt='n', yaxt='n')
axis(1, at=seq(1,5,1), label=c('1.0','2.0','3.0','4.0','5.0'), cex.axis=0.85)
axis(2, at=c(1,1.5,2,2.5), label=c('1.0','1.5','2.0','2.5'), cex.axis=0.85)
points(x=cd630_carb$fba, y=cd630_carb$biolog, pch=21, cex=1.5, bg=pall[1])
points(x=cd630_aa$fba, y=cd630_aa$biolog, pch=21, cex=1.5, bg=pall[2])
points(x=cd630_carboxy$fba, y=cd630_carboxy$biolog, pch=21, cex=1.5, bg=pall[3])
points(x=cd630_nuc$fba, y=cd630_nuc$biolog, pch=21, cex=1.5, bg=pall[4])
points(x=cd630_other$fba, y=cd630_other$biolog, pch=21, cex=1.5, bg=pall[5])
legend('topleft', legend=c('Carbohydrates','Amino acids','Carboxylic acids','Nucleotides','Other'), pch=21,
       pt.bg=pall, pt.cex=1.5, cex=0.9, bg='white')
legend('topright', legend='iCdG709 (str. 630)', pt.cex=0, bty='n', text.font=2, cex=0.9)
legend('bottomright', legend=c('R = 0.418', as.expression(bquote(paste(italic('p'),'-value < 0.001***')))), 
       cex=0.9, col='white', pch=15, pt.cex=0.5, bty='n')
abline(lm(cd630_carb$biolog ~ cd630_carb$fba), lwd=2)
box(lwd=2)

par(mar=c(3,3,0.5,0.5), las=1, mgp=c(1.9,0.75,0), lwd=1.5)
plot(x=cdR20291$fba, y=cdR20291$biolog, cex=0, xlab='Simulated Growth Enhancement Ratio', ylab='Experimental Growth Enhancement Ratio', 
     xlim=c(0.9,3.5), ylim=c(0.8,2.5), xaxt='n', yaxt='n')
axis(1, at=seq(1,5,1), label=c('1.0','2.0','3.0','4.0','5.0'), cex.axis=0.85)
axis(2, at=c(1,1.5,2,2.5), label=c('1.0','1.5','2.0','2.5'), cex.axis=0.85)
points(x=cdR20291_carb$fba, y=cdR20291_carb$biolog, pch=21, cex=1.5, bg=pall[1])
points(x=cdR20291_aa$fba, y=cdR20291_aa$biolog, pch=21, cex=1.5, bg=pall[2])
points(x=cdR20291_carboxy$fba, y=cdR20291_carboxy$biolog, pch=21, cex=1.5, bg=pall[3])
points(x=cdR20291_nuc$fba, y=cdR20291_nuc$biolog, pch=21, cex=1.5, bg=pall[4])
points(x=cdR20291_other$fba, y=cdR20291_other$biolog, pch=21, cex=1.5, bg=pall[5])
legend('topleft', legend=c('Carbohydrates','Amino acids','Carboxylic acids','Nucleotides','Other'), pch=21,
       pt.bg=pall, pt.cex=1.5, cex=0.9, bg='white')
legend('topright', legend='iCdR703 (str. R20291) ', pt.cex=0, bty='n', text.font=2, cex=0.9)
legend('bottomright', legend=c('R = 0.326', as.expression(bquote(paste(italic('p'),'-value < 0.001***')))), 
       cex=0.9, col='white', pch=15, pt.cex=0.5, bty='n')
abline(lm(cdR20291_carb$biolog ~ cdR20291_carb$fba), lwd=2)
box(lwd=2)
dev.off()

#-------------------------------------------------------------------------------------------------------------------#

# Read in data
#binary <- read.delim('~/Desktop/repos/Jenior_Cdifficile_2019/data/biolog_sim.binary.tsv', sep='\t', header=TRUE)
binary <- read.delim('~/Desktop/repos/Jenior_Cdifficile_2019/data/new_biolog_sim.binary.tsv', sep='\t', header=TRUE)
binary$name <- gsub('_',' ',binary$name)
binary_630 <- binary[,c(1:5)]
colnames(binary_630) <- c('metabolite','name','group','fba','biolog')
binary_R20291 <- binary[,c(1:3,6,7)]
colnames(binary_R20291) <- c('metabolite','name','group','fba','biolog')
rm(binary)

# Calculate hits
true_pos_630 <- subset(binary_630, fba == 1)
true_pos_630 <- subset(true_pos_630, biolog == 1)
true_neg_630 <- subset(binary_630, fba == 0)
true_neg_630 <- subset(true_neg_630, biolog == 0)
false_pos_630 <- subset(binary_630, fba == 1)
false_pos_630 <- subset(false_pos_630, biolog == 0)
false_neg_630 <- subset(binary_630, fba == 0)
false_neg_630 <- subset(false_neg_630, biolog == 1)
true_pos_630_perc <- round((nrow(true_pos_630)/(nrow(true_pos_630)+nrow(false_pos_630)))*100.0, 2)
true_pos_630_perc<- paste0('(',true_pos_630_perc,'%)')
false_pos_630_perc <- round((nrow(false_pos_630)/(nrow(true_pos_630)+nrow(false_pos_630)))*100.0, 2)
false_pos_630_perc <- paste0('(',false_pos_630_perc,'%)')
true_neg_630_perc <- round((nrow(true_neg_630)/(nrow(true_neg_630)+nrow(false_neg_630)))*100.0, 2)
true_neg_630_perc <- paste0('(',true_neg_630_perc,'%)')
false_neg_630_perc <- round((nrow(false_neg_630)/(nrow(true_pos_630)+nrow(false_neg_630)))*100.0, 2)
false_neg_630_perc <- paste0('(',false_neg_630_perc,'%)')
accuracy_630 <- round(((nrow(true_pos_630)+nrow(true_neg_630))/nrow(binary_630))*100.0, 2)
accuracy_630 <- paste0('(',accuracy_630,'%)')
true_neg_630 <- as.character(nrow(true_neg_630))
true_pos_630 <- as.character(nrow(true_pos_630))
false_pos_630 <- as.character(nrow(false_pos_630))
false_neg_630 <- as.character(nrow(false_neg_630))

# Print results
print('str. 630')
print(paste('True positive:', true_pos_630, true_pos_630_perc))
print(paste('True negative:', true_neg_630, true_neg_630_perc))
print(paste('False positive:', false_pos_630, false_pos_630_perc))
print(paste('False negative:', false_neg_630, false_neg_630_perc))
print(paste('Overall accuracy:', accuracy_630))


# Calculate hits
true_pos_R20291 <- subset(binary_R20291, fba == 1)
true_pos_R20291 <- subset(true_pos_R20291, biolog == 1)
true_neg_R20291 <- subset(binary_R20291, fba == 0)
true_neg_R20291 <- subset(true_neg_R20291, biolog == 0)
false_pos_R20291 <- subset(binary_R20291, fba == 1)
false_pos_R20291 <- subset(false_pos_R20291, biolog == 0)
false_neg_R20291 <- subset(binary_R20291, fba == 0)
false_neg_R20291 <- subset(false_neg_R20291, biolog == 1)
true_pos_R20291_perc <- round((nrow(true_pos_R20291)/(nrow(true_pos_R20291)+nrow(false_pos_R20291)))*100.0, 2)
true_pos_R20291_perc<- paste0('(',true_pos_R20291_perc,'%)')
false_pos_R20291_perc <- round((nrow(false_pos_R20291)/(nrow(true_pos_R20291)+nrow(false_pos_R20291)))*100.0, 2)
false_pos_R20291_perc <- paste0('(',false_pos_R20291_perc,'%)')
true_neg_R20291_perc <- round((nrow(true_neg_R20291)/(nrow(true_neg_R20291)+nrow(false_neg_R20291)))*100.0, 2)
true_neg_R20291_perc <- paste0('(',true_neg_R20291_perc,'%)')
false_neg_R20291_perc <- round((nrow(false_neg_R20291)/(nrow(true_pos_R20291)+nrow(false_neg_R20291)))*100.0, 2)
false_neg_R20291_perc <- paste0('(',false_neg_R20291_perc,'%)')
accuracy_R20291 <- round(((nrow(true_pos_R20291)+nrow(true_neg_R20291))/nrow(binary_R20291))*100.0, 2)
accuracy_R20291 <- paste0(accuracy_R20291,'%')
true_neg_R20291 <- as.character(nrow(true_neg_R20291))
true_pos_R20291 <- as.character(nrow(true_pos_R20291))
false_pos_R20291 <- as.character(nrow(false_pos_R20291))
false_neg_R20291 <- as.character(nrow(false_neg_R20291))

# Print results
print('str. R20291')
print(paste('True positive:', true_pos_R20291, true_pos_R20291_perc))
print(paste('True negative:', true_neg_R20291, true_neg_R20291_perc))
print(paste('False positive:', false_pos_R20291, false_pos_R20291_perc))
print(paste('False negative:', false_neg_R20291, false_neg_R20291_perc))
print(paste('Overall accuracy:', accuracy_R20291))
