
# Read in data
rough_fluxes <- read.delim('~/Desktop/repos/Jenior_Cdifficile_2019/data/transcript/tamayo_etal/phase_variation/riptide_rough_maxfit/flux_samples.tsv', sep='\t', header=TRUE, row.names=1)
smooth_fluxes <- read.delim('~/Desktop/repos/Jenior_Cdifficile_2019/data/transcript/tamayo_etal/phase_variation/riptide_smooth_maxfit/flux_samples.tsv', sep='\t', header=TRUE, row.names=1)

# Format data
shared_rxns <- intersect(colnames(rough_fluxes), colnames(smooth_fluxes))
rough_fluxes <- rough_fluxes[,shared_rxns]
smooth_fluxes <- smooth_fluxes[,shared_rxns]
rm(shared_rxns)

# Remove Biomass components
rough_fluxes[,c('biomass','dna_rxn','rna_rxn','peptidoglycan_rxn','protein_rxn','lipid_rxn','cofactor_rxn','teichoicacid_rxn','cellwall_rxn')] <- NULL
smooth_fluxes[,c('biomass','dna_rxn','rna_rxn','peptidoglycan_rxn','protein_rxn','lipid_rxn','cofactor_rxn','teichoicacid_rxn','cellwall_rxn')] <- NULL

# Subsample data
sample_size <- min(c(nrow(rough_fluxes), nrow(smooth_fluxes), 250))
sub_sample <- sample(1:min(c(nrow(rough_fluxes), nrow(smooth_fluxes))), sample_size, replace=FALSE)
rough_fluxes <- rough_fluxes[sub_sample,]
smooth_fluxes <- smooth_fluxes[sub_sample,]
rm(sample_size, sub_sample)

# Add column and merge
phase_fluxes <- as.data.frame(rbind(rough_fluxes, smooth_fluxes))
rownames(phase_fluxes) <- c(paste0('rough_', c(1:nrow(rough_fluxes))), paste0('smooth_', c(1:nrow(smooth_fluxes))))
rough_names <- paste0('rough_', c(1:nrow(rough_fluxes)))
smooth_names <- paste0('smooth_', c(1:nrow(smooth_fluxes)))
rownames(phase_fluxes) <- c(rough_names, smooth_names)
phase_fluxes <- phase_fluxes + abs(min(phase_fluxes))

# Create metadata
rough_metadata <- cbind(rough_names, rep('rough', length(rough_names)))
smooth_metadata <- cbind(smooth_names, rep('smooth', length(smooth_names)))
phase_metadata <- rbind(rough_metadata, smooth_metadata)
colnames(phase_metadata) <- c('label', 'phase')
phase_metadata <- as.data.frame(phase_metadata)
rm(rough_metadata, smooth_metadata)

# Calculate dissimilarity (Bray-Curtis)
library(vegan)
flux_dist <- vegdist(phase_fluxes, method='bray')
flux_nmds <- as.data.frame(metaMDS(flux_dist, k=2, trymax=25)$points)

# Center points
flux_x <- (abs(max(flux_nmds$MDS1)) - abs(min(flux_nmds$MDS1))) / 2
flux_y <- (abs(max(flux_nmds$MDS2)) - abs(min(flux_nmds$MDS2))) / 2
flux_nmds$MDS1 <- flux_nmds$MDS1 - flux_x
flux_nmds$MDS2 <- flux_nmds$MDS2 - flux_y
flux_xlim <- max(abs(max(flux_nmds$MDS1)), abs(min(flux_nmds$MDS1))) + 0.01
flux_ylim <- max(abs(max(flux_nmds$MDS2)), abs(min(flux_nmds$MDS2))) + 0.01
#write.table(flux_nmds, file='~/Desktop/repos/Jenior_Cdifficile_2019/data/phase_flux_nmds.tsv', 
#            quote=FALSE, sep='\t', row.names=TRUE, col.names=TRUE)
flux_nmds <- read.delim('~/Desktop/repos/Jenior_Cdifficile_2019/data/phase_flux_nmds.tsv', sep='\t', 
                        header=TRUE, row.names=1)

# Subset axes
rough_names <- paste0('rough_', c(1:nrow(rough_fluxes)))
smooth_names <- paste0('smooth_', c(1:nrow(smooth_fluxes)))
rough_nmds_points <- subset(flux_nmds, rownames(flux_nmds) %in% rough_names)
smooth_nmds_points <- subset(flux_nmds, rownames(flux_nmds) %in% smooth_names)
rm(rough_names, smooth_names)

# Statistical testing (permANOVA)
test <- merge(x=phase_metadata, y=phase_fluxes, by.x='label', by.y='row.names')
rownames(test) <- test$label
test$label <- NULL
phase_pval <- adonis(flux_dist ~ phase, data=test, perm=999, method='bray')
phase_pval <- phase_pval$aov.tab[[6]][1]
rm(test)

# Supervised learning
library(randomForest)
flux_groups <- as.factor(c(rep('rough', nrow(rough_fluxes)), rep('smooth', nrow(smooth_fluxes))))
rf_obj <- randomForest(flux_groups ~ ., data=phase_fluxes, importance=TRUE, err.rate=TRUE, ntree=1500, mtry=15)
rf_obj <- importance(rf_obj, type=1, scale=TRUE)
rf_mda <- as.data.frame(subset(rf_obj, rf_obj > (abs(min(rf_obj)))))
rf_mda$reaction <- rownames(rf_mda)
rm(rough_fluxes, smooth_fluxes, flux_groups)

# Subset for plotting
rf_mda <- rf_mda[order(-rf_mda$MeanDecreaseAccuracy),]
rf_mda <- rf_mda[c(1:20),]
rf_mda <- rf_mda[order(rf_mda$MeanDecreaseAccuracy),]
write.table(rf_mda, file='~/Desktop/repos/Jenior_Cdifficile_2019/data/phase_flux_mda.tsv', 
            quote=FALSE, sep='\t', row.names=TRUE, col.names=TRUE)
rf_mda <- read.delim('~/Desktop/repos/Jenior_Cdifficile_2019/data/phase_flux_mda.tsv', sep='\t', header=TRUE)
rf_mda$name <- gsub('_', ' ', rf_mda$name)
rf_mda <- rf_mda[order(rf_mda$mda),]

#--------------------------------------------------------------------------------------------------#

# Generate figures
rough_col <- 'sienna2'
smooth_col <- 'lightblue1'
library(scales)

pdf(file='~/Desktop/repos/Jenior_Cdifficile_2019/results/figures/Figure_S2.pdf', width=9, height=4.5)
layout(matrix(c(1,2), nrow=1, ncol=2, byrow=TRUE))

par(mar=c(3.5,3.5,1,1), las=1, mgp=c(2.2,0.7,0), lwd=2)
plot(x=flux_nmds$MDS1, y=flux_nmds$MDS2, xlim=c(-0.035,0.035), ylim=c(-0.025, 0.025), lwd.tick=2,
     xlab='NMDS Axis 1', ylab='NMDS Axis 2', pch=19, cex.lab=1.1, cex=0, cex.axis=0.8, xaxs='i', yaxs='i')
points(x=rough_nmds_points$MDS1, y=rough_nmds_points$MDS2, bg=alpha(rough_col,0.8), pch=21, cex=1.7)
points(x=smooth_nmds_points$MDS1, y=smooth_nmds_points$MDS2, bg=alpha(smooth_col,0.8), pch=21, cex=1.7)
legend('topleft', legend=c('Rough','Smooth'), bg='white',
       pt.bg=c(rough_col, smooth_col), pch=21, pt.cex=1.5, cex=1, box.lwd=2)
legend('bottomright', as.expression(bquote(paste(italic('p'),'-value = 0.001 ***'))), bty='n', pt.cex=0, cex=0.9)
par(xpd=TRUE)
text(x=-0.04, y=0.025, 'A', cex=1.2, font=2)
par(xpd=FALSE)

par(mar=c(2.5, 2.5, 0.5, 0.5), mgp=c(1.4, 0.5, 0), xpd=FALSE, lwd=1.7)
dotchart(rf_mda$mda,  xlab='Mean Decrease Accuracy', xlim=c(0,25),  
         pch=16, lwd=1.7, xaxs='i', pt.cex=0.1, cex=0.8)
text(x=-0.025, y=seq(1.3,20.3,1), labels=rf_mda$name, cex=0.75, pos=4)
points(x=rf_mda$mda, y=c(1:20), pch=16, cex=1.2)
par(xpd=TRUE)
text(x=-3, y=20.5, 'B', cex=1.2, font=2)
par(xpd=FALSE)

dev.off()
