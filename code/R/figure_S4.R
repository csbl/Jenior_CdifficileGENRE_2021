
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
biomass_pval <- wilcox.test(highspore_biomass, lowspore_biomass, exact=FALSE)$p.value

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
rf_obj <- randomForest(condition ~ ., data=spore_fluxes, importance=TRUE, err.rate=TRUE, ntree=1500, mtry=15)
rf_obj <- importance(rf_obj, type=1, scale=TRUE)
rf_mda <- as.data.frame(subset(rf_obj, rf_obj > (abs(min(rf_obj)))))
rm(condition)

# Subset to most informative features
rf_mda <- as.data.frame(subset(rf_obj, rf_obj > (abs(min(rf_obj)))))
rf_mda$reaction <- rownames(rf_mda)
rf_mda <- rf_mda[order(-rf_mda$MeanDecreaseAccuracy),]
rf_mda <- rf_mda[c(1:20),]
rm(rf_obj)
write.table(rf_mda, file='~/Desktop/repos/Jenior_Cdifficile_2019/data/spore_flux_mda.tsv', 
            quote=FALSE, sep='\t', row.names=FALSE, col.names=TRUE)
rf_mda <- read.delim('~/Desktop/repos/Jenior_Cdifficile_2019/data/spore_flux_mda.tsv', sep='\t', header=TRUE)
rf_mda <- rf_mda[order(rf_mda$mda),]
rf_mda$name <- gsub('_', ' ', rf_mda$name)

# Generate figures
hispore_col <- 'darkorchid2'
lospore_col <- 'aquamarine2'
library(plotrix)

pdf(file='~/Desktop/repos/Jenior_Cdifficile_2019/results/figures/Figure_S4.pdf', width=6, height=6)
layout(matrix(c(1,2,
                1,2,
                3,3), nrow=3, ncol=2, byrow=TRUE))

par(mar=c(4,3,0.5,0.5), xpd=FALSE, las=1, mgp=c(2,0.7,0), lwd=2)
boxplot(highspore_biomass, at=0.5, xlim=c(0,2), ylab='Simulated Biomass Flux', 
        ylim=c(80,120), yaxt='n', col=hispore_col, outline=FALSE, cex.lab=1.2,
        boxlwd=2, medlwd=2, staplelwd=2, whisklwd=2, whisklty=1, width=0.5)
boxplot(lowspore_biomass, at=1.5, add=TRUE, yaxt='n', col=lospore_col,
        boxlwd=2, medlwd=2, staplelwd=2, whisklwd=2, whisklty=1, width=0.5, outline=FALSE)
axis(side=2, at=seq(80,120,10), labels=c(0,seq(90,120,10)), cex.axis=0.9, lwd=1.7)
axis.break(2, 83)
par(xpd=TRUE)
text(x=c(0.5,1.5), y=74, labels=c('High\nSpores','Low\nSpores'), cex=1.3)
text(x=-0.35, y=120.5, 'A', cex=1.5, font=2)
par(xpd=FALSE)


par(mar=c(2.5, 2.5, 0.5, 0.5), mgp=c(1.4, 0.5, 0), xpd=FALSE, lwd=1.7)
dotchart(rf_mda$mda,  xlab='Mean Decrease Accuracy', xlim=c(0,25),  
         pch=16, lwd=1.7, xaxs='i', pt.cex=0.1, cex=0.8)
text(x=-0.025, y=seq(1.3,20.3,1), labels=rf_mda$name, cex=0.8, pos=4)
points(x=rf_mda$mda, y=c(1:20), pch=16, cex=1.2)
par(xpd=TRUE)
text(x=-3, y=20.5, 'B', cex=1.5, font=2)
par(xpd=FALSE)




# C. extra in vitro results....



dev.off()
