
# Read in data
rough_fluxes <- read.delim('~/Desktop/repos/Jenior_Cdifficile_2019/data/transcript/tamayo_etal/phase_variation/riptide_rough_maxfit/flux_samples.tsv', sep='\t', header=TRUE, row.names=1)
smooth_fluxes <- read.delim('~/Desktop/repos/Jenior_Cdifficile_2019/data/transcript/tamayo_etal/phase_variation/riptide_smooth_maxfit/flux_samples.tsv', sep='\t', header=TRUE, row.names=1)

# Identify context-specific reactions
rough_only <- setdiff(colnames(rough_fluxes), colnames(smooth_fluxes))
smooth_only <- setdiff(colnames(smooth_fluxes), colnames(rough_fluxes))

# Get median absolute flux for each
rough_med_flux <- c()
for (x in rough_only) {rough_med_flux <- c(rough_med_flux, abs(median(rough_fluxes[,x])))}
rough_only_flux <- as.data.frame(cbind(rough_only, rough_med_flux, rep('rough',length(rough_only))))
colnames(rough_only_flux) <- c('reaction','abs_med_flux','context')
smooth_med_flux <- c()
for (x in smooth_only) {smooth_med_flux <- c(smooth_med_flux, abs(median(smooth_fluxes[,x])))}
smooth_only_flux <- as.data.frame(cbind(smooth_only, smooth_med_flux, rep('smooth',length(smooth_only))))
colnames(smooth_only_flux) <- c('reaction','abs_med_flux','context')
unique_fluxes <- as.data.frame(rbind(rough_only_flux, smooth_only_flux))
rm(rough_only, rough_med_flux, rough_only_flux, smooth_only, smooth_med_flux, smooth_only_flux, x)
write.table(unique_fluxes, file='~/Desktop/repos/Jenior_Cdifficile_2019/data/phase_unique_flux.tsv', 
            quote=FALSE, sep='\t', row.names=FALSE, col.names=TRUE)

# Biomass flux
smooth_biomass <- as.vector(smooth_fluxes[,'biomass'])
rough_biomass <- as.vector(rough_fluxes[,'biomass'])

# Format data for intersection analysis
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
rf_mda <- rf_mda[c(1:15),]
rf_mda <- rf_mda[order(rf_mda$MeanDecreaseAccuracy),]
write.table(rf_mda, file='~/Desktop/repos/Jenior_Cdifficile_2019/data/phase_flux_mda.tsv', 
            quote=FALSE, sep='\t', row.names=TRUE, col.names=TRUE)

#--------------------------------------------------------------------------------------------------#

# Generate figures
smooth_col <- 'lightsteelblue2'
rough_col <- 'red3'
library(scales)

pdf(file='~/Desktop/repos/Jenior_Cdifficile_2019/results/figures/Figure_S2.pdf', width=5, height=3.5)
layout(matrix(c(1,2,2), nrow=1, ncol=3, byrow=TRUE))

par(mar=c(5,3,1,1), xpd=FALSE, las=1, mgp=c(1.6,0.7,0), lwd=2)
boxplot(smooth_biomass, at=0.5, xlim=c(0,2), ylab='Simulated Biomass Flux', 
        ylim=c(0,60), yaxt='n', col=smooth_col, outline=FALSE, cex.lab=1.2,
        boxlwd=2, medlwd=2, staplelwd=2, whisklwd=2, whisklty=1, width=0.5)
boxplot(rough_biomass, at=1.5, add=TRUE, yaxt='n', col=rough_col,
        boxlwd=2, medlwd=2, staplelwd=2, whisklwd=2, whisklty=1, width=0.5, outline=FALSE)
axis(side=2, at=seq(0,60,10), cex.axis=0.9, lwd=1.7)
segments(x0=0.5, y0=56, x1=1.5, lwd=2)
text(x=1, y=59, '***', cex=1.5, font=2)
par(xpd=TRUE)
text(x=c(0.5,1.5), y=-10, labels=c('Smooth','Rough'), cex=1.4, srt=55)
text(x=-0.7, y=63, 'A', cex=1.4, font=2)
par(xpd=FALSE)

par(mar=c(3.5,3.5,1,1), las=1, mgp=c(2.2,0.7,0), lwd=2)
plot(x=flux_nmds$MDS1, y=flux_nmds$MDS2, xlim=c(-0.035,0.035), ylim=c(-0.025, 0.025), lwd.tick=2,
     xlab='NMDS Axis 1', ylab='NMDS Axis 2', pch=19, cex.lab=1.1, cex=0, cex.axis=0.8, xaxs='i', yaxs='i')
points(x=rough_nmds_points$MDS1, y=rough_nmds_points$MDS2, bg=alpha(rough_col,0.8), pch=21, cex=1.7)
points(x=smooth_nmds_points$MDS1, y=smooth_nmds_points$MDS2, bg=alpha(smooth_col,0.8), pch=21, cex=1.7)
legend('topleft', legend=c('Rough','Smooth'), bg='white',
       pt.bg=c(rough_col, smooth_col), pch=21, pt.cex=1.5, cex=1, box.lwd=2)
legend('bottomright', as.expression(bquote(paste(italic('p'),'-value = 0.001 ***'))), bty='n', pt.cex=0, cex=0.9)
par(xpd=TRUE)
text(x=-0.04, y=0.025, 'B', cex=1.4, font=2)
par(xpd=FALSE)

dev.off()
