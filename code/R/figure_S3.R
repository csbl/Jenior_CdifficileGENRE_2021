
# Doubling times
mutant_doubling <- as.data.frame(t(read.delim('~/Desktop/active_projects/cmrRST/cmr_growthrate.tsv', sep='\t', header=F, row.names=1)))

# Flux sampling files
rough_samples <- '~/Desktop/repos/Jenior_Cdifficile_2019/data/transcript/tamayo_etal/phase_variation/riptide_rough/flux_samples.tsv'
smooth_samples <- '~/Desktop/repos/Jenior_Cdifficile_2019/data/transcript/tamayo_etal/phase_variation/riptide_smooth/flux_samples.tsv'

# Read in data
rough_samples <- read.delim(rough_samples, sep='\t', header=TRUE)
rough_samples$X <- NULL
smooth_samples <- read.delim(smooth_samples, sep='\t', header=TRUE)
smooth_samples$X <- NULL

# Format data
overlap <- intersect(colnames(smooth_samples), colnames(rough_samples))
rough_samples <- rough_samples[, overlap]
smooth_samples <- smooth_samples[, overlap]
rm(overlap)

# Subsample data
sample_size <- min(c(nrow(rough_samples), nrow(smooth_samples), 250))
sub_sample <- sample(1:min(c(nrow(rough_samples), nrow(smooth_samples))), sample_size, replace=FALSE)
rough_samples <- rough_samples[sub_sample,]
smooth_samples <- smooth_samples[sub_sample,]
rm(sample_size, sub_sample)

# Format row names
rough_names <- paste('rough_', 1:nrow(rough_samples), sep='')
rownames(rough_samples) <- rough_names
smooth_names <- paste('smooth_', 1:nrow(smooth_samples), sep='')
rownames(smooth_samples) <- smooth_names

# Create metadata
rough_metadata <- cbind(rough_names, rep('rough', length(rough_names)))
smooth_metadata <- cbind(smooth_names, rep('smooth', length(smooth_names)))
metadata <- rbind(rough_metadata, smooth_metadata)
colnames(metadata) <- c('label', 'phase')
metadata <- as.data.frame(metadata)
rm(rough_metadata, smooth_metadata)

# Merge data and prep for unsupervised learning
all_samples <- rbind(rough_samples, smooth_samples)
all_samples <- all_samples + abs(min(all_samples))

# Calculate dissimilarity (Bray-Curtis)
library(vegan)
flux_dist <- vegdist(all_samples, method='bray') # Bray-Curtis

# Mean wwithin-group dissimilarity
flux_groups <- as.factor(c(rep('rough',nrow(rough_samples)), rep('smooth',nrow(smooth_samples))))
meandist(flux_dist, grouping=flux_groups)
#             rough     smooth
# rough  0.04235277 0.04734441
# smooth 0.04734441 0.03592208

# Unsupervised learning (NMDS)
flux_nmds <- as.data.frame(metaMDS(flux_dist, k=2, trymax=25)$points)

# Center points
flux_x <- (abs(max(flux_nmds$MDS1)) - abs(min(flux_nmds$MDS1))) / 2
flux_y <- (abs(max(flux_nmds$MDS2)) - abs(min(flux_nmds$MDS2))) / 2
flux_nmds$MDS1 <- flux_nmds$MDS1 - flux_x
flux_nmds$MDS2 <- flux_nmds$MDS2 - flux_y
rm(flux_x, flux_y)

# Subset axes
rough_nmds_points <- subset(flux_nmds, rownames(flux_nmds) %in% rough_names)
smooth_nmds_points <- subset(flux_nmds, rownames(flux_nmds) %in% smooth_names)
rm(rough_names, smooth_names)

# Statistical testing (permANOVA)
test <- merge(x=metadata, y=all_samples, by.x='label', by.y='row.names')
rownames(test) <- test$label
test$label <- NULL
test$treatment <- NULL
pval <- adonis(flux_dist ~ phase, data=test, perm=999, method='bray')
pval <- pval$aov.tab[[6]][1]
pval <- as.character(round(pval, 3))
rm(all_samples, test, flux_dist, metadata)

# Flux sampling files
cmrON_samples <- '~/Desktop/active_projects/cmrRST/cmr_ON/flux_samples.tsv'
cmrOFF_samples <- '~/Desktop/active_projects/cmrRST/cmr_OFF/flux_samples.tsv'

# Read in data
cmrON_samples <- read.delim(cmrON_samples, sep='\t', header=TRUE)
cmrON_samples$X <- NULL
cmrOFF_samples <- read.delim(cmrOFF_samples, sep='\t', header=TRUE)
cmrOFF_samples$X <- NULL

# Format data
overlap <- intersect(colnames(cmrOFF_samples), colnames(cmrON_samples))
cmrON_samples <- cmrON_samples[, overlap]
cmrOFF_samples <- cmrOFF_samples[, overlap]
rm(overlap)

# Subsample data
sample_size <- min(c(nrow(cmrON_samples), nrow(cmrOFF_samples), 250))
sub_sample <- sample(1:min(c(nrow(cmrON_samples), nrow(cmrOFF_samples))), sample_size, replace=FALSE)
cmrON_samples <- cmrON_samples[sub_sample,]
cmrOFF_samples <- cmrOFF_samples[sub_sample,]
rm(sample_size, sub_sample)

# Format row names
cmrON_names <- paste('cmrON_', 1:nrow(cmrON_samples), sep='')
rownames(cmrON_samples) <- cmrON_names
cmrOFF_names <- paste('cmrOFF_', 1:nrow(cmrOFF_samples), sep='')
rownames(cmrOFF_samples) <- cmrOFF_names

# Create metadata
cmrON_metadata <- cbind(cmrON_names, rep('cmrON', length(cmrON_names)))
cmrOFF_metadata <- cbind(cmrOFF_names, rep('cmrOFF', length(cmrOFF_names)))
cmr_metadata <- rbind(cmrON_metadata, cmrOFF_metadata)
colnames(cmr_metadata) <- c('label', 'phase')
cmr_metadata <- as.data.frame(cmr_metadata)
rm(cmrON_metadata, cmrOFF_metadata)

# Merge data and prep for unsupervised learning
cmr_all_samples <- rbind(cmrON_samples, cmrOFF_samples)
cmr_all_samples <- cmr_all_samples + abs(min(cmr_all_samples))

# Calculate dissimilarity (Bray-Curtis)
cmr_flux_dist <- vegdist(cmr_all_samples, method='bray') # Bray-Curtis

# Mean wwithin-group dissimilarity
cmr_flux_groups <- as.factor(c(rep('cmrON',nrow(cmrON_samples)), rep('cmrOFF',nrow(cmrOFF_samples))))
meandist(cmr_flux_dist, grouping=cmr_flux_groups)
#             cmrON     cmrOFF
# cmrON  0.03823604 0.05805795
# cmrOFF 0.05805795 0.05438550

# Unsupervised learning (NMDS)
cmr_flux_nmds <- as.data.frame(metaMDS(cmr_flux_dist, k=2, trymax=25)$points)

# Center points
flux_x <- (abs(max(cmr_flux_nmds$MDS1)) - abs(min(cmr_flux_nmds$MDS1))) / 2
flux_y <- (abs(max(cmr_flux_nmds$MDS2)) - abs(min(cmr_flux_nmds$MDS2))) / 2
cmr_flux_nmds$MDS1 <- cmr_flux_nmds$MDS1 - flux_x
cmr_flux_nmds$MDS2 <- cmr_flux_nmds$MDS2 - flux_y
rm(flux_x, flux_y)

# Subset axes
cmrON_nmds_points <- subset(cmr_flux_nmds, rownames(cmr_flux_nmds) %in% cmrON_names)
cmrOFF_nmds_points <- subset(cmr_flux_nmds, rownames(cmr_flux_nmds) %in% cmrOFF_names)
rm(cmrON_names, cmrOFF_names)

# Statistical testing (permANOVA)
test <- merge(x=cmr_metadata, y=cmr_all_samples, by.x='label', by.y='row.names')
rownames(test) <- test$label
test$label <- NULL
test$treatment <- NULL
cmr_pval <- adonis(cmr_flux_dist ~ phase, data=test, perm=999, method='bray')
cmr_pval <- cmr_pval$aov.tab[[6]][1]
cmr_pval <- as.character(round(cmr_pval, 3))
rm(cmr_all_samples, test, cmr_flux_dist, cmr_metadata)

# Generate figures
library(scales)
library(plotrix)
smooth_col <- 'white'
rough_col <- 'cornflowerblue'

# Flux sample dissimilarity
png(filename='~/Desktop/repos/Jenior_Cdifficile_2019/results/figures/Figure_S3A.png', 
    units='in', width=5, height=5, res=300)
par(mar=c(3.5,3.5,0.5,0.5), las=1, mgp=c(2.2,0.7,0), lwd=1.7)
plot(x=flux_nmds$MDS1, y=flux_nmds$MDS2, xlim=c(-0.045,0.045), ylim=c(-0.045, 0.045),
     xlab='NMDS Axis 1', ylab='NMDS Axis 2', pch=19, cex.lab=1.1, cex=0, cex.axis=0.9)
points(x=smooth_nmds_points$MDS1, y=smooth_nmds_points$MDS2, bg=alpha(smooth_col,0.8), pch=21, cex=1.7)
points(x=rough_nmds_points$MDS1, y=rough_nmds_points$MDS2, bg=alpha(rough_col,0.8), pch=21, cex=1.7)
legend('topleft', legend=c('Smooth variant','Rough variant'), 
       pt.bg=c(alpha(smooth_col,0.8),alpha(rough_col,0.8)), pch=21, pt.cex=1.7, box.lwd=2)
legend('bottomright', legend='Shared Core Metabolism', pt.cex=0, bty='n')
legend('topright', legend='p-value < 0.001', pt.cex=0, bty='n', cex=0.9)
box(lwd=2)
dev.off()

# Biomass flux samples
png(filename='~/Desktop/repos/Jenior_Cdifficile_2019/results/figures/Figure_S3B.png', 
    units='in', width=2, height=4, res=300)
par(mar=c(2,3,0.5,0.5), xpd=FALSE, las=1, mgp=c(2,0.8,0), lwd=1.7)
boxplot(mutant_doubling$OFF, at=0.5, xlim=c(0,2), ylab='Sampled Doubling Time (min)', 
        ylim=c(40,60), yaxt='n', col=smooth_col, outline=FALSE,
        boxlwd=2, medlwd=2, staplelwd=2, whisklwd=2, whisklty=1, width=0.9)
boxplot(mutant_doubling$ON, at=1.5, xlim=c(0,2), add=TRUE, yaxt='n', col=rough_col,
        boxlwd=2, medlwd=2, staplelwd=2, whisklwd=2, whisklty=1, width=0.9, outline=FALSE)
axis(side=2, at=seq(40,60,2), labels=c(0,seq(42,60,2)), cex.axis=0.9, lwd=1.7)
axis.break(2, 41, style='slash')
segments(x0=0.5, y0=55, x1=1.5, lwd=2)
text(x=1, y=56, 'n.s.', cex=0.9)
mtext(c('cmr OFF','cmr ON'), side=1, padj=0.5, adj=c(0.1,0.9), cex=0.8)
dev.off()

# Mutant flux sample dissimilarity
png(filename='~/Desktop/repos/Jenior_Cdifficile_2019/results/figures/Figure_S3C.png', 
    units='in', width=5, height=5, res=300)
par(mar=c(3.5,3.5,0.5,0.5), las=1, mgp=c(2.2,0.7,0), lwd=1.7)
plot(x=cmr_flux_nmds$MDS1, y=cmr_flux_nmds$MDS2, xlim=c(-0.05,0.05), ylim=c(-0.05, 0.05),
     xlab='NMDS Axis 1', ylab='NMDS Axis 2', pch=19, cex.lab=1.1, cex=0, cex.axis=0.9)
points(x=cmrOFF_nmds_points$MDS1, y=cmrOFF_nmds_points$MDS2, bg=alpha(smooth_col,0.8), pch=21, cex=1.7)
points(x=cmrON_nmds_points$MDS1, y=cmrON_nmds_points$MDS2, bg=alpha(rough_col,0.8), pch=21, cex=1.7)
legend('bottomleft', legend=c('cmr ON-locked mutant','cmr OFF-locked mutant'), 
       pt.bg=c(alpha(smooth_col,0.8),alpha(rough_col,0.8)), pch=21, pt.cex=1.7, box.lwd=2)
legend('bottomright', legend='Shared Core Metabolism', pt.cex=0, bty='n')
legend('topright', legend='p-value = 0.001', pt.cex=0, bty='n', cex=0.9)
box(lwd=2)
dev.off()

