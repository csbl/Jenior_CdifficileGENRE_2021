
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
#rough  0.04464993 0.04671530
#smooth 0.04671530 0.04322881

# Unsupervised learning (NMDS)
flux_nmds <- as.data.frame(metaMDS(flux_dist, k=2, trymax=25)$points)

# Center points
flux_x <- (abs(max(flux_nmds$MDS1)) - abs(min(flux_nmds$MDS1))) / 2
flux_y <- (abs(max(flux_nmds$MDS2)) - abs(min(flux_nmds$MDS2))) / 2
flux_nmds$MDS1 <- flux_nmds$MDS1 - flux_x
flux_nmds$MDS2 <- flux_nmds$MDS2 - flux_y
flux_x <- max(abs(max(flux_nmds$MDS1)), abs(min(flux_nmds$MDS1))) + 0.01
flux_y <- max(abs(max(flux_nmds$MDS2)), abs(min(flux_nmds$MDS2))) + 0.01

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

# Ordination of shared reaction flux samples
library(scales)
png(filename='~/Desktop/repos/Jenior_Cdifficile_2019/results/figures/figure_S3.png', units='in', width=4.5, height=4.5, res=300)
par(mar=c(3.5,3.5,0.5,0.5), las=1, mgp=c(2.2,0.7,0), lwd=1.7)
plot(x=flux_nmds$MDS1, y=flux_nmds$MDS2, xlim=c(-0.045,0.045), ylim=c(-0.045, 0.045),
     xlab='NMDS Axis 1', ylab='NMDS Axis 2', pch=19, cex.lab=1.1, cex=0, cex.axis=0.9)
points(x=smooth_nmds_points$MDS1, y=smooth_nmds_points$MDS2, bg=alpha('chocolate2',0.8), pch=21, cex=1.7)
points(x=rough_nmds_points$MDS1, y=rough_nmds_points$MDS2, bg=alpha('cyan3',0.8), pch=21, cex=1.7)
legend('topleft', legend=c('Smooth','Rough'), 
       pt.bg=c(alpha('chocolate2',0.8),alpha('cyan3',0.8)), pch=21, pt.cex=1.7, box.lwd=2)
#legend('bottomleft', legend=substitute(paste(italic('p'), '-value = 0.001**')), pt.cex=0, bty='n')
legend('bottomright', legend='Shared Core Metabolism', pt.cex=0, bty='n')
box(lwd=2)
dev.off()

