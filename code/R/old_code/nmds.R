
library(vegan)

# Metabolomes
metabolome <- '~/Desktop/repos/Jenior_Metatranscriptomics_mSphere_2018/data/metabolome/scaled_intensities.log10.tsv'

# 16S
shared_otu <- '~/Desktop/repos/Jenior_Metatranscriptomics_mSphere_2018/data/16S_analysis/all_treatments.0.03.unique_list.conventional.shared'
otu_tax <- '~/Desktop/repos/Jenior_Metatranscriptomics_mSphere_2018/data/16S_analysis/formatted.all_treatments.0.03.cons.taxonomy'

# Metadata
metadata <- '~/Desktop/repos/Jenior_Metatranscriptomics_mSphere_2018/data/metadata.tsv'

#----------------#

# Read in data

# Metabolomes
metabolome <- read.delim(metabolome, sep='\t', header=TRUE)

# 16S
shared_otu <- read.delim(shared_otu, sep='\t', header=T, row.names=2)
shared_otu$numOtus <- NULL
shared_otu$label <- NULL
shared_otu <- shared_otu[!rownames(shared_otu) %in% c('CefC5M2'), ]  # Remove possible contaminated sample
shared_otu <- shared_otu[,!(names(shared_otu) %in% c('Otu0004','Otu0308'))] # Remove residual C. difficile OTUs
otu_tax <- read.delim(otu_tax, sep='\t', header=T, row.names=1)

# Metadata
metadata <- read.delim(metadata, sep='\t', header=T, row.names=1)

#-------------------------------------------------------------------------------------------------------------------------#

# Format data

# Metadata
metadata$type <- NULL
metadata$cage <- NULL
metadata$mouse <- NULL
metadata$gender <- NULL

# Metabolomes
metabolome$BIOCHEMICAL <- gsub('_', ' ', metabolome$BIOCHEMICAL)
rownames(metabolome) <- metabolome$BIOCHEMICAL
metabolome$BIOCHEMICAL <- NULL
metabolome$PUBCHEM <- NULL
metabolome$KEGG <- NULL
metabolome$SUB_PATHWAY <- NULL
metabolome$SUPER_PATHWAY <- NULL
metabolome <- as.data.frame(t(metabolome))
metabolome <- merge(metadata, metabolome, by='row.names')
rownames(metabolome) <- metabolome$Row.names
metabolome$Row.names <- NULL
metabolome <- subset(metabolome, abx != 'germfree')
metabolome <- subset(metabolome, abx != 'none')
metabolome <- subset(metabolome, infection == 'mock')
metabolome$susceptibility <- NULL
metabolome$infection <- NULL

# 16S
subSize <- min(rowSums(shared_otu))
shared_otu <- t(shared_otu)
for (x in 1:ncol(shared_otu)) {
  shared_otu[,x] <- as.vector(rrarefy(shared_otu[,x], sample=subSize))
}
rm(subSize)
shared_otu <- as.data.frame(t(shared_otu))
otu_tax$genus <- gsub('_', ' ', otu_tax$genus)
otu_tax$genus <- gsub('Ruminococcus2', 'Ruminococcus', otu_tax$genus)
otu_tax$taxon <- paste(otu_tax$genus, otu_tax$OTU.1, sep='_')
otu_tax$phylum <- NULL
otu_tax$genus <- NULL
otu_tax$OTU.1 <- NULL
shared_otu <- t(shared_otu)
shared_otu <- merge(otu_tax, shared_otu, by='row.names')
rownames(shared_otu) <- shared_otu$Row.names
shared_otu$Row.names <- NULL
rownames(shared_otu) <- shared_otu$taxon
shared_otu$taxon <- NULL
shared_otu <- t(shared_otu)
shared_otu <- merge(metadata, shared_otu, by='row.names')
rownames(shared_otu) <- shared_otu$Row.names
shared_otu$Row.names <- NULL
shared_otu <- subset(shared_otu, abx != 'germfree')
shared_otu <- subset(shared_otu, abx != 'none')
shared_otu <- subset(shared_otu, infection == 'mock')
shared_otu$susceptibility <- NULL
shared_otu$infection <- NULL
rm(otu_tax)

#-------------------------------------------------------------------------------------------------------------------------#

# Ordination analysis

# Metabolome - 728 metabolites
metabolome_dist <- vegdist(metabolome[,2:ncol(metabolome)], method='bray') # Bray-Curtis
metabolome_nmds <- as.data.frame(metaMDS(metabolome_dist, k=2, trymax=100)$points)
# Subset NMDS axes to color points
rownames(metabolome_nmds) <- rownames(metabolome)
metabolome_nmds <- merge(metadata, metabolome_nmds, by='row.names')
rownames(metabolome_nmds) <- metabolome_nmds$Row.names
metabolome_nmds$Row.names <- NULL
metabolome_nmds_strep <- subset(metabolome_nmds, abx == 'streptomycin')
metabolome_nmds_cef <- subset(metabolome_nmds, abx == 'cefoperazone')
metabolome_nmds_clinda <- subset(metabolome_nmds, abx == 'clindamycin')
# permANOVA
metabolome_permANOVA_pval <- adonis(metabolome_dist ~ metabolome$abx, metabolome, perm=999)$aov.tab[[6]][1]
metabolome_permANOVA_pval <- as.character(round(metabolome_permANOVA_pval, 3))
rm(metabolome_dist)

# Find maxima and minima of axes
metabolome_x_axes <- c(min(metabolome_nmds$MDS1), max(metabolome_nmds$MDS1))
metabolome_y_axes <- c(min(metabolome_nmds$MDS2), max(metabolome_nmds$MDS2))

# 16S - 810 OTUs
otu_dist <- vegdist(shared_otu[,2:ncol(shared_otu)], method='bray') # Bray-Curtis
otu_nmds <- as.data.frame(metaMDS(otu_dist, k=2, trymax=100)$points)
otu_nmds$MDS1 <- otu_nmds$MDS1 - 0.1
# Subset NMDS axes to color points
rownames(otu_nmds) <- rownames(shared_otu)
otu_nmds <- merge(metadata, otu_nmds, by='row.names')
rownames(otu_nmds) <- otu_nmds$Row.names
otu_nmds$Row.names <- NULL
otu_nmds_strep <- subset(otu_nmds, abx == 'streptomycin')
otu_nmds_cef <- subset(otu_nmds, abx == 'cefoperazone')
otu_nmds_clinda <- subset(otu_nmds, abx == 'clindamycin')
# permANOVA
otu_permANOVA_pval <- adonis(otu_dist ~ shared_otu$abx, shared_otu, perm=999)$aov.tab[[6]][1]
otu_permANOVA_pval <- as.character(round(otu_permANOVA_pval, 3))
rm(otu_dist)
rm(metadata)

# Find maxima and minima of axes
otu_x_axes <- c(min(otu_nmds$MDS1), max(otu_nmds$MDS1))
otu_y_axes <- c(min(otu_nmds$MDS2), max(otu_nmds$MDS2))

#-------------------------------------------------------------------------------------------------------------------------#

# Plot the figure
png(filename='~/Desktop/nmds.png', units='in', width=11, height=5.5, res=300)
layout(matrix(c(1,2), nrow=1, ncol=2, byrow=TRUE))
par(mar=c(4,4,1,3), las=1, mgp=c(2.8,0.75,0))

# 16S
plot(x=otu_nmds$MDS1, y=otu_nmds$MDS2, xlim=c(-0.6,0.6), ylim=c(-0.4,0.4),
     xlab='NMDS axis 1', ylab='NMDS axis 2', pch=19, cex.axis=1.2, cex.lab=1.2)
points(x=otu_nmds_strep$MDS1, y=otu_nmds_strep$MDS2, bg='gray', pch=21, cex=2, lwd=1.2)
points(x=otu_nmds_cef$MDS1, y=otu_nmds_cef$MDS2, bg='dodgerblue3', pch=21, cex=2, lwd=1.2)
points(x=otu_nmds_clinda$MDS1, y=otu_nmds_clinda$MDS2, bg='firebrick3', pch=21, cex=2, lwd=1.2)
legend('bottomleft', legend=c('Between Antibiotics', as.expression(bquote(paste(italic('p'),' < 0.001')))), 
       pch=1, cex=1.2, pt.cex=0, bty='n')
legend('bottomright', legend=c('Streptomycin','Cefoperazone','Clindamycin'), 
       pt.bg=c('gray', 'dodgerblue3', 'firebrick3'), pch=21, cex=1.1, pt.cex=2)
legend('topleft', legend='16S amplicon sequencing', bty='n', cex=1.2, pt.cex=0)
box(lwd=2)

# Metabolome
plot(x=metabolome_nmds$MDS1, y=metabolome_nmds$MDS2, xlim=c(-0.4,0.4), ylim=c(-0.25,0.25),
     xlab='NMDS axis 1', ylab='NMDS axis 2', pch=19, cex.axis=1.2, cex.lab=1.2)
points(x=metabolome_nmds_strep$MDS1, y=metabolome_nmds_strep$MDS2, bg='gray', pch=21, cex=2, lwd=1.2)
points(x=metabolome_nmds_cef$MDS1, y=metabolome_nmds_cef$MDS2, bg='dodgerblue3', pch=21, cex=2, lwd=1.2)
points(x=metabolome_nmds_clinda$MDS1, y=metabolome_nmds_clinda$MDS2, bg='firebrick3', pch=21, cex=2, lwd=1.2)
legend('bottomleft', legend=c('Between Antibiotics', as.expression(bquote(paste(italic('p'),' < 0.001')))), 
       pch=1, cex=1.2, pt.cex=0, bty='n')
legend('bottomright', legend=c('Streptomycin','Cefoperazone','Clindamycin'), 
       pt.bg=c('gray', 'dodgerblue3', 'firebrick3'), pch=21, cex=1.1, pt.cex=2)
legend('topleft', legend='Metabolomics', bty='n', cex=1.2, pt.cex=0)
box(lwd=2)

dev.off()


