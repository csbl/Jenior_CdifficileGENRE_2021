
# Flux sampling files
highspore_fluxes <- '~/Desktop/repos/Jenior_Cdifficile_2019/data/contextualized_iCdG709/riptide_highspore_maxfit/flux_samples.tsv'
lowspore_fluxes <- '~/Desktop/repos/Jenior_Cdifficile_2019/data/contextualized_iCdG709/riptide_lowspore_maxfit/flux_samples.tsv'

# Read in data
highspore_fluxes <- read.delim(highspore_fluxes, sep='\t', header=TRUE)
highspore_fluxes$X <- NULL
lowspore_fluxes <- read.delim(lowspore_fluxes, sep='\t', header=TRUE)
lowspore_fluxes$X <- NULL

# Subsample
sample_size <- min(c(nrow(highspore_fluxes), nrow(lowspore_fluxes), 250))
sub_sample <- sample(1:min(c(nrow(highspore_fluxes), nrow(lowspore_fluxes))), sample_size, replace=FALSE)
highspore_fluxes <- highspore_fluxes[sub_sample,]
lowspore_fluxes <- lowspore_fluxes[sub_sample,]
rm(sample_size, sub_sample)

# Subset data
highspore_biomass <- as.vector(highspore_fluxes[,'biomass'])
lowspore_biomass <- as.vector(lowspore_fluxes[,'biomass'])

# Test difference
pvals <- c()
for (x in c(1:1000)) {
  test_1 <- sample(highspore_biomass, size=10)
  test_2 <- sample(lowspore_biomass, size=10)
  pvals[x] <- wilcox.test(test_1, test_2, exact=FALSE)$p.value}
biomass_pval <- round(median(pvals), 4)

# Format data for random forest
shared_rxns <- intersect(colnames(highspore_fluxes), colnames(lowspore_fluxes))
highspore_fluxes <- highspore_fluxes[,shared_rxns]
lowspore_fluxes <- lowspore_fluxes[,shared_rxns]
highspore_fluxes[,c('biomass','dna_rxn','rna_rxn','peptidoglycan_rxn','protein_rxn','lipid_rxn','cofactor_rxn','teichoicacid_rxn','cellwall_rxn')] <- NULL
lowspore_fluxes[,c('biomass','dna_rxn','rna_rxn','peptidoglycan_rxn','protein_rxn','lipid_rxn','cofactor_rxn','teichoicacid_rxn','cellwall_rxn')] <- NULL
rm(shared_rxns)

# Add column and merge
highspore_fluxes$spore <- rep('highspore', nrow(highspore_fluxes))
lowspore_fluxes$spore <- rep('lowspore', nrow(lowspore_fluxes))
spore_fluxes <- as.data.frame(rbind(highspore_fluxes, lowspore_fluxes))
rownames(spore_fluxes) <- c(paste0('highspore_', c(1:nrow(highspore_fluxes))), paste0('lowspore_', c(1:nrow(lowspore_fluxes))))
rm(highspore_fluxes, lowspore_fluxes)

# Find most informative metabolites with RF
library(randomForest)
condition <- as.factor(spore_fluxes$spore)
spore_fluxes$spore <- NULL
spore_fluxes <- droplevels(spore_fluxes)
#rf_obj <- randomForest(condition ~ ., data=spore_fluxes, importance=TRUE, err.rate=TRUE, ntree=1500, mtry=15)
#rf_obj <- importance(rf_obj, type=1, scale=TRUE)
#rf_mda <- as.data.frame(subset(rf_obj, rf_obj > (abs(min(rf_obj)))))
#rm(condition)

# Subset to most informative features
#rf_mda <- as.data.frame(subset(rf_obj, rf_obj > (abs(min(rf_obj)))))
#rf_mda$reaction <- rownames(rf_mda)
#rf_mda <- rf_mda[order(-rf_mda$MeanDecreaseAccuracy),]
#rf_mda <- rf_mda[c(1:20),]
#rm(rf_obj)
#write.table(rf_mda, file='~/Desktop/repos/Jenior_Cdifficile_2019/data/spore_flux_mda.tsv', 
#            quote=FALSE, sep='\t', row.names=FALSE, col.names=TRUE)
rf_mda <- read.delim('~/Desktop/repos/Jenior_Cdifficile_2019/data/spore_flux_mda.tsv', sep='\t', header=TRUE)
rf_mda <- rf_mda[order(rf_mda$mda),]
rf_mda$name <- gsub('_', ' ', rf_mda$name)

# Metabolome files
metabolome <- '~/Desktop/repos/Jenior_Cdifficile_2019/data/metabolome/scaled_intensities.tsv'
metadata <- '~/Desktop/repos/Jenior_Cdifficile_2019/data/metadata.tsv'

# Read in data
metabolome <- read.delim(metabolome, sep='\t', header=T, row.names=1)
metadata <- read.delim(metadata, sep='\t', header=T, row.names=1)

# Merge metabolomics with metadata
metadata$cage <- NULL
metadata$mouse <- NULL
metadata$gender <- NULL
metabolome$SUPER_PATHWAY <- NULL
metabolome$SUB_PATHWAY <- NULL
metabolome$PUBCHEM <- NULL
metabolome$KEGG <- NULL
metabolome <- as.data.frame(t(metabolome))
metabolome <- merge(metadata, metabolome, by='row.names')
rownames(metabolome) <- metabolome$Row.names
metabolome$Row.names <- NULL
metabolome <- subset(metabolome, type == 'conventional')
metabolome$type <- NULL
rm(metadata)

select_and_plot <- function(metabolite_name, best_ylim=0, correction=FALSE, title=NA, across=FALSE) {
  metabolite <- metabolome[, c(1,2,which(colnames(metabolome) %in% c(metabolite_name)))]
  metabolite <- subset(metabolite, abx %in% c('cefoperazone','clindamycin'))
  colnames(metabolite) <- c('abx','infection','intensity')
  
  if (correction == TRUE) {metabolite$intensity <- metabolite$intensity - min(metabolite$intensity)}
  if (is.na(title)) {title <- metabolite_name}
  
  a_mock <- subset(subset(metabolite, abx=='clindamycin'), infection=='mock')[,3]
  a_630 <- subset(subset(metabolite, abx=='clindamycin'), infection=='630')[,3]
  b_mock <- subset(subset(metabolite, abx=='cefoperazone'), infection=='mock')[,3]
  b_630 <- subset(subset(metabolite, abx=='cefoperazone'), infection=='630')[,3]
  colnames(metabolite) <- c('abx','infection','intensity')
  
  if (best_ylim == 0) {best_ylim <- ceiling(as.numeric(quantile(metabolite[,3], 0.95)))}
  if (across == TRUE) {best_ylim <- best_ylim + (best_ylim * 0.5)}
  
  par(mar=c(2.7,2.8,1.5,0.5), xpd=FALSE, mgp=c(1.6,0.7,0), lwd=2, xaxt='n', las=1)
  boxplot(a_mock, at=1, xlim=c(0.6,3.4), ylim=c(0,best_ylim), col=hispore_col, cex.axis=0.7, cex.lab=1, 
          xlab="", ylab="Scaled Intensity (Log10)", outcex=0, whisklty=1, medlwd=2,  main=title, cex.main=1.1) 
  boxplot(a_630, at=1.5, xlim=c(0.6,3.4), ylim=c(0,best_ylim), col=hispore_col, 
          xlab="", ylab="", outcex=0, whisklty=1, medlwd=2, yaxt='n', add=TRUE) 
  if (across == FALSE) {abline(v=2)}
  boxplot(b_mock, at=2.5, xlim=c(0.6,3.4), ylim=c(0,best_ylim), col=lospore_col, 
          xlab="", ylab="", outcex=0, whisklty=1, medlwd=2, yaxt='n', add=TRUE) 
  boxplot(b_630, at=3, xlim=c(0.6,3.4), ylim=c(0,best_ylim), col=lospore_col, 
          xlab="", ylab="", outcex=0, whisklty=1, medlwd=2, yaxt='n', add=TRUE) 
  mtext(c('High Spores','Low Spores'), side=1, at=c(1.25,2.75), padj=1.6, cex=1)
  mtext(c('CDI:','-','+','-','+'), side=1, at=c(0.5,1,1.5,2.5,3), padj=0.2, cex=c(0.75,1,1,1,1))
  
  a_pval <- round(wilcox.test(a_mock, a_630, exact=FALSE)$p.value,3)
  b_pval <- round(wilcox.test(b_mock, b_630, exact=FALSE)$p.value,3)
  ab_pval <- round(wilcox.test(a_mock, b_mock, exact=FALSE)$p.value,3)
  ba_pval <- round(wilcox.test(a_630, b_630, exact=FALSE)$p.value,3)
  if (across == TRUE) {fctr <- 0.7} else {fctr <- 0.96}
  if (a_pval <= 0.01) {
    segments(x0=1, y0=best_ylim*fctr, x1=1.5, y1=best_ylim*fctr, lwd=1.5)
    text(x=1.25, y=best_ylim*(fctr+0.04), '**', font=2, cex=1.5)}
  if (b_pval <= 0.001) {
    segments(x0=2.5, y0=best_ylim*fctr, x1=3, y1=best_ylim*fctr, lwd=1.5)
    text(x=2.75, y=best_ylim*(fctr+0.04), '***', font=2, cex=1.5)}
  if (across == TRUE) {
    if (ab_pval <= 0.001) {
      segments(x0=1, y0=best_ylim*0.96, x1=2.5, y1=best_ylim*0.96, lwd=1.5)
      text(x=1.75, y=best_ylim, '***', font=2, cex=1.5)}
    if (ba_pval <= 0.001) {
      segments(x0=1.5, y0=best_ylim*0.86, x1=3, y1=best_ylim*0.86, lwd=1.5)
      text(x=2.25, y=best_ylim*0.9, '***', font=2, cex=1.5)}
  }
}

leucine <- metabolome[, c(1,2,which(colnames(metabolome) %in% c('leucine')))]
leucine <- subset(leucine, abx %in% c('clindamycin'))
colnames(leucine) <- c('abx','infection','intensity')
leucine_mock <- subset(subset(leucine, abx=='clindamycin'), infection=='mock')[,3]
leucine_630 <- subset(subset(leucine, abx=='clindamycin'), infection=='630')[,3]
round(wilcox.test(leucine_mock, leucine_630, exact=FALSE)$p.value,3)

# Additional CFU spores vs vegetative
# Read in data
sporulation <- read.delim('~/Desktop/repos/Jenior_Cdifficile_2019/data/supp_veg_spore_cfu.tsv', sep='\t', header=TRUE)
bdm_neu5Ac_cyt_veg <- as.vector(sporulation[,'bdm_neu5Ac_cyt_veg'])
bdm_neu5Ac_cyt_spore <- as.vector(sporulation[,'bdm_neu5Ac_cyt_spore'])
bdm_neu5Ac_veg <- as.vector(sporulation[,'bdm_neu5Ac_veg'])
bdm_neu5Ac_spore <- as.vector(sporulation[,'bdm_neu5Ac_spore'])
bdm_cyt_veg <- as.vector(sporulation[,'bdm_cyt_veg'])
bdm_cyt_spore <- as.vector(sporulation[,'bdm_cyt_spore'])
rm(sporulation)

# Assemble summary statistics
veg_med <- c(median(bdm_neu5Ac_cyt_veg), median(bdm_neu5Ac_veg), median(bdm_cyt_veg))
veg_q25 <- as.numeric(c(quantile(bdm_neu5Ac_cyt_veg, 0.25), quantile(bdm_neu5Ac_veg, 0.25), quantile(bdm_cyt_veg, 0.25)))
veg_q75 <- as.numeric(c(quantile(bdm_neu5Ac_cyt_veg, 0.75), quantile(bdm_neu5Ac_veg, 0.75), quantile(bdm_cyt_veg, 0.75)))
spore_med <- c(median(bdm_neu5Ac_cyt_spore), median(bdm_neu5Ac_spore), median(bdm_cyt_spore))
spore_q25 <- as.numeric(c(quantile(bdm_neu5Ac_cyt_spore, 0.25), quantile(bdm_neu5Ac_spore, 0.25), quantile(bdm_cyt_spore, 0.25)))
spore_q75 <- as.numeric(c(quantile(bdm_neu5Ac_cyt_spore, 0.75), quantile(bdm_neu5Ac_spore, 0.75), quantile(bdm_cyt_spore, 0.75)))

# Transform
veg_med <- log10(veg_med+ 1)
spore_med <- log10(spore_med + 1)
veg_q25 <- log10(veg_q25 + 1)
spore_q25 <- log10(spore_q25 + 1)
veg_q75 <- log10(veg_q75 + 1)
spore_q75 <- log10(spore_q75 + 1)
all_med <- rbind(veg_med, spore_med)

# Fix LOD for visualization
veg_q25[veg_q25 == 0] <- 2
spore_q25[spore_q25 == 0] <- 2

#----------------------------------------------------------------------------------------------#

# Generate figures
hispore_col <- 'darkorchid2'
lospore_col <- 'aquamarine2'
veg_col <- 'deepskyblue'
spore_col <- 'tomato'
if (!require(plotrix)) install.packages('plotrix')
library(plotrix)
if (!require(scales)) install.packages('scales')
library(scales)

pdf(file='~/Desktop/repos/Jenior_Cdifficile_2019/results/figures/Figure_S4.pdf', width=8, height=6)
layout(matrix(c(1,2,2,
                3,4,5), nrow=2, ncol=3, byrow=TRUE))

# Biomass
par(mar=c(4,5,1,1), xpd=FALSE, las=1, mgp=c(2,0.7,0), lwd=2)
boxplot(highspore_biomass, at=0.5, xlim=c(0,2), ylab='Simulated Biomass Flux', 
        ylim=c(80,120), yaxt='n', col=hispore_col, outline=FALSE, cex.lab=1.2,
        boxlwd=2, medlwd=2, staplelwd=2, whisklwd=2, whisklty=1, width=0.5)
boxplot(lowspore_biomass, at=1.5, add=TRUE, yaxt='n', col=lospore_col,
        boxlwd=2, medlwd=2, staplelwd=2, whisklwd=2, whisklty=1, width=0.5, outline=FALSE)
axis(side=2, at=seq(80,120,10), labels=c(0,seq(90,120,10)), cex.axis=0.9, lwd=1.7)
axis.break(2, 83)
par(xpd=TRUE)
text(x=c(0.5,1.5), y=73, labels=c('High Spores','Low Spores'), cex=1.2)
text(x=-0.6, y=120.5, 'A', cex=1.5, font=2)
par(xpd=FALSE)

# Random forest results
par(mar=c(2.5, 2.5, 0.5, 0.5), mgp=c(1.4, 0.5, 0), xpd=FALSE, lwd=1.7)
dotchart(rf_mda$mda,  xlab='Mean Decrease Accuracy (%)', xlim=c(0,20),  
         pch=16, lwd=1.7, xaxs='i', pt.cex=0.1, cex=0.8)
text(x=-0.025, y=seq(1.3,10.3,1), labels=rf_mda$name, cex=1, pos=4)
points(x=rf_mda$mda, y=c(1:10), pch=16, cex=1.5)
par(xpd=TRUE)
text(x=-1.5, y=10.5, 'B', cex=1.5, font=2)
par(xpd=FALSE)

# Metabolomics
select_and_plot('5-aminovalerate', title='5-Aminovalerate', best_ylim=5)
par(xpd=TRUE)
text(x=0.1, y=5, 'C', cex=1.5, font=2)
par(xpd=FALSE)

select_and_plot('N-acetylneuraminate', title='N-Acetylneuraminic acid', best_ylim=2.5)
par(xpd=TRUE)
text(x=0.1, y=2.5, 'D', cex=1.5, font=2)
par(xpd=FALSE)

select_and_plot('cytidine', title='Cytidine', best_ylim=1.5)
par(xpd=TRUE)
text(x=0.1, y=1.5, 'E', cex=1.5, font=2)
par(xpd=FALSE)

dev.off()



