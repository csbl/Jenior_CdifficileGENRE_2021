
# Set up environment
rm(list=ls())
gc()

# Load in functions
starting_dir <- getwd()
source('~/Desktop/repositories/Cdiff_modeling/src/init.R')

# Output plot name
ordination_plot <- 'results/figures/nmds.pdf'
features_plot <- 'results/figures/aucrf.pdf'
family_plot <- 'results/figures/barplot.pdf'
legend_plot <- 'results/figures/legend.pdf'

# Data
shared <- 'data/16S/cdf_models.opti_mcc.0.03.shared'
tax <- 'data/16S/cdf_models.opti_mcc.0.03.cons.format.taxonomy'
shared_family <- 'data/16S/all_treatments.family.subsample.shared'
taxonomy_family <- 'data/16S/all_treatments.family.cons.format.taxonomy'
metadata <- 'data/metadata.tsv'

#----------------#

# Taxonomy
tax <- read.delim(tax, sep='\t', header=T, row.names=1)
taxonomy_family <- read.delim(taxonomy_family, sep='\t', header=T)

# Shared file
shared <- read.delim(shared, sep='\t', header=T, row.names=2)
shared_family <- read.delim(shared_family, sep='\t', header=T, row.names=2)

# Metadata
metadata <- read.delim(metadata, sep='\t', header=T, row.names=1)

#----------------#

# Reformatting

# Metadata
metadata$cage <- NULL
metadata$mouse <- NULL
metadata$gender <- NULL
metadata$clearance <- NULL
metadata$infection <- NULL

# Shared
shared$label <- NULL
shared$numOtus <- NULL
shared <- rarefyOTU(shared)
shared <- filtOTU(shared, 3)
shared <- mergeByRow(metadata, shared)
shared <- subset(shared, type == 'conventional')
shared$type <- NULL
metadata <- subset(metadata, type == 'conventional')
metadata$type <- NULL
shared_family <- shared_family[!rownames(shared_family) %in% c('CefC5M2'), ]  # Remove contaminated sample
shared_family <- shared_family[ ,!(names(shared_family) == 'Otu008')] # Remove residual C. difficile OTU
shared_family$numOtus <- NULL
shared_family$label <- NULL

# Remove C. difficile OTUs
cdf <- c('Otu0004','Otu0117','Otu0171','Otu0321','Otu0335','Otu0371',
         'Otu0403','Otu0537','Otu0549','Otu0592','Otu1179','Otu1296')
shared <- shared[,which(!colnames(shared) %in% cdf)]
tax <- tax[which(!rownames(tax) %in% cdf),]
rm(cdf)

#-------------------------------------------------------------------------------------------------------------------------#

# Distance analysis
otu_dist <- designdist(shared[,3:ncol(shared)], method='1-(J/(A+B-J))', terms='quadratic', abcd=FALSE) # Theta-YC
# permANOVA
susceptibility_pval <- adonis(otu_dist ~ shared$susceptibility, shared[,3:ncol(shared)], perm=999)$aov.tab[[6]][1]
susceptibility_pval <- as.character(round(susceptibility_pval, 3))
abx_pval <- adonis(otu_dist ~ shared$abx, shared[,3:ncol(shared)], perm=999)$aov.tab[[6]][1]
abx_pval <- as.character(round(abx_pval, 3))
shared$abx <- NULL
# Ordination
otu_nmds <- as.data.frame(metaMDS(otu_dist, k=2, trymax=100)$points)
otu_nmds$MDS1 <- otu_nmds$MDS1 - 0.07 # Center points
otu_nmds$MDS2 <- otu_nmds$MDS2 - 0.04 # Center points
rm(otu_dist)
# Subset NMDS axes to color points
rownames(otu_nmds) <- rownames(shared)
otu_nmds <- mergeByRow(metadata, otu_nmds)
otu_nmds_strep <- subset(otu_nmds, abx == 'streptomycin')
otu_nmds_cef <- subset(otu_nmds, abx == 'cefoperazone')
otu_nmds_clinda <- subset(otu_nmds, abx == 'clindamycin')
otu_nmds_resistant <- subset(otu_nmds, susceptibility == 'resistant')
otu_nmds_susceptible <- subset(otu_nmds, susceptibility == 'susceptible')
# Calculate centroids
otu_resistant_centroids <- aggregate(cbind(otu_nmds_resistant$MDS1, otu_nmds_resistant$MDS2)~otu_nmds_resistant$susceptibility, data=otu_nmds_resistant, mean)
otu_susceptible_centroids <- aggregate(cbind(otu_nmds_susceptible$MDS1, otu_nmds_susceptible$MDS2)~otu_nmds_susceptible$susceptibility, data=otu_nmds_susceptible, mean)

#-------------------------------------------------------------------------------------------------------------------------#

# Feature selection and subset shared by significant features
aucrf_feat <- aucrfSus(shared) # 0% OOB
sig_shared <- as.data.frame(t(subset(t(shared), colnames(shared) %in% rownames(aucrf_feat))))

# Break into separate shared files for statistical comparison
sig_shared$susceptibility <- shared$susceptibility
sus_shared <- subset(sig_shared, susceptibility == 'susceptible')
sus_shared$susceptibility <- NULL
for (i in 1:ncol(sus_shared)) {
  sus_shared[,i] <- as.numeric(as.character(sus_shared[,i]))
}
res_shared <- subset(sig_shared, susceptibility == 'resistant')
res_shared$susceptibility <- NULL
for (i in 1:ncol(res_shared)) {
  res_shared[,i] <- as.numeric(as.character(res_shared[,i]))
}
rm(sig_shared)

# Find significant differences
abund_pval <- c()
for (i in 1:ncol(sus_shared)){abund_pval[i] <- wilcox.test(sus_shared[,i], res_shared[,i], exact=FALSE)$p.value}
abund_pval <- round(p.adjust(abund_pval, method='BH'), 5)
# Log transform values
sus_shared <- log10(sus_shared + 1)
res_shared <- log10(res_shared + 1)

#-------------#

# Add names from tax
sus_shared <- merge(t(sus_shared), tax, by='row.names')
sus_shared$Row.names <- NULL
sus_shared$genus <- c('Porphyromonadaceae unclassified(100) 1', 'Porphyromonadaceae_unclassified(100) 2', 
                      'Lachnospiraceae unclassified(100) 1', 'Clostridium XlVa(100)', 
                      'Lachnospiraceae unclassified(100) 2')
rownames(sus_shared) <- sus_shared$genus
sus_shared$phylum <- NULL
sus_shared$class <- NULL
sus_shared$order <- NULL
sus_shared$family <- NULL
sus_shared$genus <- NULL
sus_shared <- as.data.frame(t(sus_shared))
res_shared <- merge(t(res_shared), tax, by='row.names')
res_shared$Row.names <- NULL
res_shared$genus <- c('Porphyromonadaceae unclassified(100) 1', 'Porphyromonadaceae_unclassified(100) 2', 
                      'Lachnospiraceae unclassified(100) 1', 'Clostridium XlVa(100)', 
                      'Lachnospiraceae unclassified(100) 2')
rownames(res_shared) <- res_shared$genus
res_shared$phylum <- NULL
res_shared$class <- NULL
res_shared$order <- NULL
res_shared$family <- NULL
res_shared$genus <- NULL
res_shared <- as.data.frame(t(res_shared))


# Family-level shared file
shared_family <- mergeByRow(metadata, shared_family)
shared_family$susceptibility <- NULL
shared_family <- aggregate(. ~ abx, data=shared_family, FUN=median)
rownames(shared_family) <- shared_family$abx
shared_family$abx <- NULL

# Convert to relative abundance
relabund_family <- (shared_family / rowSums(shared_family)) * 100
rm(shared_family)

# Bin lowly abundant OTUs into an 'Other' category
relabund_family[relabund_family < 1] <- 0
relabund_family <- relabund_family[, which(colSums(relabund_family) > 0)]
top_otus <- colnames(relabund_family)

# Subset family-level taxonomy
taxonomy_family <- subset(taxonomy_family, taxonomy_family$otu %in% top_otus)
taxonomy_family$phylum <- gsub('_', ' ', taxonomy_family$phylum)
taxonomy_family$family <- gsub('_', ' ', taxonomy_family$family)
taxonomy_family <- taxonomy_family[order(taxonomy_family$family),]
rm(top_otus)

# Calculate left over abundances
relabund_family$Other <- 100 - rowSums(relabund_family)

# Define group colors
taxonomy_family[] <- lapply(taxonomy_family, as.character)
taxonomy_family <- rbind(taxonomy_family, c('Other', 'Other', 'Other (<1% each)'))
taxonomy_family <- taxonomy_family[order(taxonomy_family$phylum), ]
family_colors <- c('mediumblue','dodgerblue2','powderblue',
                   'darkred','firebrick3','red2','tomato2','salmon',
                   'navajowhite4',
                   'darkgoldenrod1')
taxonomy_family$color <- family_colors
other <- taxonomy_family[which(taxonomy_family$family == 'Other (<1% each)'),]
taxonomy_family <- subset(taxonomy_family, family != 'Other (<1% each)')
taxonomy_family <- rbind(taxonomy_family, other)
rm(other)

# Add empty columns for plotting and sort table
relabund_family <- relabund_family[, taxonomy_family$otu] # reorder shared according to taxonomy
relabund_family$abx <- c("cefoperazone","clindamycin","none","streptomycin")
relabund_family <- relabund_family[order(match(relabund_family$abx, c('none','streptomycin','cefoperazone','clindamycin'))),] # sort by treatment
relabund_family$abx <- NULL

#-------------------------------------------------------------------------------------------------------------------------#

# Generate figures

# NMDS
pdf(file=ordination_plot, width=7, height=6)
par(mar=c(3,3,1,1), las=1, mgp=c(2,0.7,0))
plot(x=otu_nmds$MDS1, y=otu_nmds$MDS2, xlim=c(-0.6,0.6), ylim=c(-0.6,0.6),
     xlab='NMDS axis 1', ylab='NMDS axis 2', pch=19, cex.axis=1.2, cex.lab=1.2, cex=0.2)
segments(x0=otu_nmds_resistant$MDS1, y0=otu_nmds_resistant$MDS2, x1=otu_resistant_centroids[1,2], y1=otu_resistant_centroids[1,3], col='grey25')
points(x=otu_nmds_resistant$MDS1, y=otu_nmds_resistant$MDS2, bg=resistant_col, pch=21, cex=2.1, lwd=1.2)
segments(x0=otu_nmds_susceptible$MDS1, y0=otu_nmds_susceptible$MDS2, x1=otu_susceptible_centroids[1,2], y1=otu_susceptible_centroids[1,3], col='gray25')
points(x=otu_nmds_strep$MDS1, y=otu_nmds_strep$MDS2, bg=susceptible_col, pch=24, cex=1.8, lwd=1.2)
points(x=otu_nmds_cef$MDS1, y=otu_nmds_cef$MDS2, bg=susceptible_col, pch=23, cex=2.1, lwd=1.2)
points(x=otu_nmds_clinda$MDS1, y=otu_nmds_clinda$MDS2, bg=susceptible_col, pch=22, cex=1.8, lwd=1.2)
legend('bottomright', legend=c(as.expression(bquote(paste(italic('p'),' = 0.001')))), 
       pch=1, cex=1.6, pt.cex=0, bty='n')
legend('topleft', legend=c('Susceptible','Resistant'), 
       pt.bg=c(susceptible_col, resistant_col), pch=22, cex=1.3, pt.cex=2)
legend('bottomleft', legend=c('Streptomycin','Cefoperazone','Clindamycin'), 
       pch=c(17,15,18), cex=1.3, pt.cex=c(2,2,3))
box(lwd=2)
dev.off()

# Feature selection results
pdf(file=features_plot, width=4, height=6)
multiStripchart(sus_shared, res_shared,
                'Susceptible','Resistant',
                susceptible_col, resistant_col,
                abund_pval, 'Relative Abundance (Log10)')
dev.off()

# Family-level phylotype bar chart
pdf(file=family_plot, width=6, height=4)
par(mar=c(2,4,1,1), mgp=c(2.5, 0.25, 0), new=FALSE, xpd=FALSE)
barplot(t(rev(relabund_family)), col=rev(taxonomy_family$color), yaxt='n', xaxt='n', cex.lab=1.3,
        ylim=c(0,100), ylab='Relative Abundance', cex.names=1.2, space=0)
box(lwd=1.5)
axis(side=2, at=seq(0,100,20), labels=c('0%','20%','40%','60%','80%','100%'), tick=FALSE, las=1)
abline(h=c(20,40,60,80), lty=2)
mtext(c('No Antibiotics','Streptomycin','Cefoperazone','Clindamycin'), font=2,
      side=1, at=c(0.5,1.5,2.5,3.5), adj=0.5, padj=1, cex=0.9, col='black')
dev.off()

# Create a figure legend in empty plot
pdf(file=legend_plot, width=3.5, height=4)
par(mar=c(0,0,0,1))
plot(0, type='n', ylim=c(-10,10), xlim=c(5,5), ylab='', xlab='', xaxt='n', yaxt='n', axes=FALSE)
legend('right', legend=taxonomy_family$family, pt.bg=taxonomy_family$color, 
       pch=22, pt.cex=2, cex=0.9, bty='n')
# Add in phylum classifications
segments(x0=rep(4.45,3), x1=rep(4.45,3), 
         y0=c(4.9,1.7,-2.75), y1=c(2.2,-2.25,-3.55), 
         lwd=3) # vertical
text(x=rep(3.7,3), y=c(3.5,-0.3,-3.15), cex=0.9,
     labels=c('Bacteroidetes', 'Firmicutes', 'Proteobacteria'))
text(x=c(3.75,5.5), y=-6, labels=c('Phylum','Family'), cex=1.1, font=2, xpd=TRUE)
dev.off()

#-------------------------------------------------------------------------------------------------------------------------#

# Clean up
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
setwd(starting_dir)
rm(list=ls())
gc()
