
set.seed(9861)

# Flux sampling files
highspore_samples <- '~/Desktop/repos/Jenior_Cdifficile_2019/data/contextualized_iCdG709/riptide_highspore_maxfit/flux_samples.tsv'
lowspore_samples <- '~/Desktop/repos/Jenior_Cdifficile_2019/data/contextualized_iCdG709/riptide_lowspore_maxfit/flux_samples.tsv'

# Read in data
highspore_samples <- read.delim(highspore_samples, sep='\t', header=TRUE)
highspore_samples$X <- NULL
lowspore_samples <- read.delim(lowspore_samples, sep='\t', header=TRUE)
lowspore_samples$X <- NULL

# Limit to overlapping reactions
overlap <- intersect(colnames(lowspore_samples), colnames(highspore_samples))
highspore_samples <- highspore_samples[, overlap]
lowspore_samples <- lowspore_samples[, overlap]
rm(overlap)

# Subsample data
sample_size <- min(c(nrow(highspore_samples), nrow(lowspore_samples), 250))
sub_sample <- sample(1:min(c(nrow(highspore_samples), nrow(lowspore_samples))), sample_size, replace=FALSE)
highspore_samples <- highspore_samples[sub_sample,]
lowspore_samples <- lowspore_samples[sub_sample,]
rm(sample_size, sub_sample)

# Format row names
highspore_names <- paste('highspore_', 1:nrow(highspore_samples), sep='')
rownames(highspore_samples) <- highspore_names
lowspore_names <- paste('lowspore_', 1:nrow(lowspore_samples), sep='')
rownames(lowspore_samples) <- lowspore_names

# Create metadata
highspore_metadata <- cbind(highspore_names, rep('highspore', length(highspore_names)), rep('cleared', length(highspore_names)))
lowspore_metadata <- cbind(lowspore_names, rep('lowspore', length(lowspore_names)), rep('colonized', length(lowspore_names)))
spore_metadata <- rbind(highspore_metadata, lowspore_metadata)
colnames(spore_metadata) <- c('label', 'sporulation','clearance')
spore_metadata <- as.data.frame(spore_metadata)
rm(highspore_metadata, lowspore_metadata)

# Merge data and prep for unsupervised learning
all_samples <- rbind(highspore_samples, lowspore_samples)
flux_groups <- as.factor(c(rep('highspore', nrow(highspore_samples)), rep('lowspore', nrow(lowspore_samples))))
all_samples_adj <- all_samples + abs(min(all_samples))

library(randomForest)
rf_obj <- randomForest(flux_groups ~ ., data=all_samples, importance=TRUE, err.rate=TRUE, ntree=1500, mtry=15)
rf_obj <- importance(rf_obj, type=1, scale=TRUE)
rf_mda <- as.data.frame(subset(rf_obj, rf_obj > (abs(min(rf_obj)))))
all_samples_adj <- all_samples_adj[,rownames(rf_mda)]

# Calculate dissimilarity (Bray-Curtis)
library(vegan)
flux_dist <- vegdist(all_samples_adj, method='bray') # Bray-Curtis

# Mean wwithin-group dissimilarity
meandist(flux_dist, grouping=flux_groups)
#           highspore   lowspore
# highspore 0.0212113 0.02562170
# lowspore  0.0256217 0.02425617

# Unsupervised learning (NMDS)
flux_nmds <- as.data.frame(metaMDS(flux_dist, k=2, trymax=25)$points)

# Center points
flux_x <- (abs(max(flux_nmds$MDS1)) - abs(min(flux_nmds$MDS1))) / 2
flux_y <- (abs(max(flux_nmds$MDS2)) - abs(min(flux_nmds$MDS2))) / 2
flux_nmds$MDS1 <- flux_nmds$MDS1 - flux_x
flux_nmds$MDS2 <- flux_nmds$MDS2 - flux_y
flux_xlim <- max(abs(max(flux_nmds$MDS1)), abs(min(flux_nmds$MDS1))) + 0.01
flux_ylim <- max(abs(max(flux_nmds$MDS2)), abs(min(flux_nmds$MDS2))) + 0.01
write.table(flux_nmds, file='~/Desktop/repos/Jenior_Cdifficile_2019/data/spore_flux_nmds.tsv', 
            quote=FALSE, sep='\t', row.names=TRUE, col.names=TRUE)

flux_nmds <- read.delim('~/Desktop/repos/Jenior_Cdifficile_2019/data/spore_flux_nmds.tsv', sep='\t', 
                        header=TRUE, row.names=1)

# Subset axes
highspore_nmds_points <- subset(flux_nmds, rownames(flux_nmds) %in% highspore_names)
lowspore_nmds_points <- subset(flux_nmds, rownames(flux_nmds) %in% lowspore_names)
rm(highspore_names, lowspore_names)

# Statistical testing (permANOVA)
test <- merge(x=spore_metadata, y=all_samples_adj, by.x='label', by.y='row.names')
rownames(test) <- test$label
test$label <- NULL
test$clearance <- NULL
spore_pval <- adonis(flux_dist ~ sporulation, data=test, perm=999, method='bray')
spore_pval <- spore_pval$aov.tab[[6]][1]

#----------------------------------------------------------------------#

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

check_diff <- function(metabolite_name) {
  metabolite <- metabolome[, c(1,2,which(colnames(metabolome) %in% c(metabolite_name)))]
  metabolite <- subset(metabolite, abx %in% c('cefoperazone','clindamycin'))
  colnames(metabolite) <- c('abx','infection','intensity')
  
  hi_mock <- subset(subset(metabolite, abx=='clindamycin'), infection=='mock')[,3]
  hi_630 <- subset(subset(metabolite, abx=='clindamycin'), infection=='630')[,3]
  lo_mock <- subset(subset(metabolite, abx=='cefoperazone'), infection=='mock')[,3]
  lo_630 <- subset(subset(metabolite, abx=='cefoperazone'), infection=='630')[,3]
  
  hi_mean_diff <- mean(hi_630) - mean(hi_mock)
  lo_mean_diff <- mean(lo_630) - mean(lo_mock)
  hi_pval <- wilcox.test(hi_mock, hi_630, exact=FALSE)$p.value
  lo_pval <- wilcox.test(lo_mock, lo_630, exact=FALSE)$p.value
  
  print(paste('high mean diff:', hi_mean_diff))
  print(paste('high p-value:', hi_pval))
  print(paste('low mean diff:', lo_mean_diff))
  print(paste('low p-value:', lo_pval))
}

check_diff('cytidine')
check_diff('ornithine')
check_diff('serine')
check_diff('proline')
check_diff('N-acetylneuraminate')
check_diff('5-aminovalerate')
check_diff('2,3-dihydroxyisovalerate')

#----------------------------------------------------------------------#

# Biomass contribution importance

# Read and format data
contributions <- read.delim('/home/mjenior/Desktop/repos/Jenior_Cdifficile_2019/data/invivo_biomass_contributions.tsv', sep='\t', header=TRUE)
high_spore <- subset(contributions, spores =='high')
high_spore_uptake <- subset(high_spore, direction =='consume')
high_spore_uptake <- high_spore_uptake[order(high_spore_uptake$effect),]
high_spore_efflux <- subset(high_spore, direction =='produce')
high_spore_efflux <- high_spore_efflux[order(high_spore_efflux$effect),]
high_spore <- as.data.frame(rbind(high_spore_efflux, high_spore_uptake))
low_spore <- subset(contributions, spores =='low')
low_spore_uptake <- subset(low_spore, direction =='consume')
low_spore_uptake <- low_spore_uptake[order(low_spore_uptake$effect),]
low_spore_efflux <- subset(low_spore, direction =='produce')
low_spore_efflux <- low_spore_efflux[order(low_spore_efflux$effect),]
low_spore <- as.data.frame(rbind(low_spore_efflux, low_spore_uptake))
rm(contributions)

#----------------------------------------------------------------------#

# Generate figure
library(scales)
consume_col <- 'white'
produce_col <- 'darkgray'
hispore_col <- 'darkorchid2'
lospore_col <- 'aquamarine2'

pdf(file='~/Desktop/repos/Jenior_Cdifficile_2019/results/figures/Figure_4.pdf', width=9, height=6)
layout(matrix(c(1,2,3,
                1,4,5), nrow=2, ncol=3, byrow=TRUE))

par(mar=c(15.5,3.5,11,0.5), las=1, mgp=c(2.2,0.7,0), lwd=2)
plot(x=flux_nmds$MDS1, y=flux_nmds$MDS2, xlim=c(-0.014,0.014), ylim=c(-0.014, 0.014), lwd.tick=2,
     xlab='NMDS Axis 1', ylab='NMDS Axis 2', pch=19, cex.lab=1.1, cex=0, cex.axis=0.8, xaxs='i', yaxs='i')
points(x=lowspore_nmds_points$MDS1, y=lowspore_nmds_points$MDS2, bg=alpha(lospore_col,0.8), pch=21, cex=1.7)
points(x=highspore_nmds_points$MDS1, y=highspore_nmds_points$MDS2, bg=alpha(hispore_col,0.8), pch=21, cex=1.7)
legend('topleft', legend=c('High Sporulation','Low Sporulation'), bg='white',
       pt.bg=c(hispore_col, lospore_col), pch=21, pt.cex=1.5, cex=1, box.lwd=2)
text(x=0, y=-0.013, as.expression(bquote(paste(italic('p'),'-value = 0.001 ***'))), cex=0.9, pos=4)
box()
par(xpd=TRUE)
text(x=-0.017, y=0.014, 'A', cex=1.5, font=2)
par(xpd=FALSE)

par(mar=c(3,6,1.5,2.7), xpd=FALSE, las=1, mgp=c(1.5,0.5,0), lwd=2)
barplot(high_spore$effect, horiz=TRUE, xlab='Impact on Biomass (%)', xlim=c(0,6), lwd.tick=2, main='High Sporulation',
        names.arg=high_spore$metabolite, cex.names=0.9, cex.axis=0.8, cex.main=0.9,
        col=c(rep(produce_col,nrow(high_spore_efflux)),rep(consume_col,nrow(high_spore_uptake))))
text(x=5.6, y=4.5, as.expression(bquote(paste(italic('in silico'),' Predictions'))), srt=-90, cex=0.9)
legend('bottomright', legend=c('Consumption','Production'), bg='white', pt.lwd=1.5,
       pt.bg=c(consume_col,produce_col), pch=22, pt.cex=1.5, cex=0.9, bty='n')
box(lwd=2.5)
par(xpd=TRUE)
#rect(xleft=6.3, xright=7.0, ytop=8.4, ybottom=7.6, col='white')
#text(x=6.65, y=8, 'ND', cex=0.6)
text(x=6.65, y=8, 'Not\nMeasured', cex=0.6)
#rect(xleft=6.3, xright=7.0, ytop=7.2, ybottom=6.4, col='white')
#text(x=6.65, y=6.8, 'ND', cex=0.6)
text(x=6.65, y=6.8, 'Not\nMeasured', cex=0.6)
rect(xleft=6.3, xright=7.0, ytop=5.95, ybottom=5.15, col='deepskyblue2')
rect(xleft=6.3, xright=7.0, ytop=4.75, ybottom=3.95, col='deepskyblue2')
text(x=6.65, y=4.35, '*', cex=1.1, font=2)
rect(xleft=6.3, xright=7.0, ytop=3.5, ybottom=2.7, col='deepskyblue2')
rect(xleft=6.3, xright=7.0, ytop=2.3, ybottom=1.5, col='firebrick2')
text(x=6.65, y=1.9, '*', cex=1.1, font=2)
rect(xleft=6.3, xright=7.0, ytop=1.05, ybottom=0.25, col='firebrick2')
text(x=6.65, y=0.65, '*', cex=1.1, font=2)
text(x=-1.2, y=8.75, 'B', cex=1.5, font=2)
par(xpd=FALSE)

par(mar=c(0,0,0,0), lwd=1.5)
plot(0, type='n', ylim=c(0,50), xlim=c(0,5), ylab='', xlab='', xaxt='n', yaxt='n', axes=FALSE)
legend('left', legend=c('Decreased during CDI', 'Increased during CDI'), cex=0.9,
       pt.bg=c('deepskyblue2','firebrick2'), pch=22, pt.cex=1.9, bty='n', 
       title=as.expression(bquote(paste(italic('in vivo'),' Metabolomics'))))

par(mar=c(3,6,1.5,2.7), xpd=FALSE, las=1, mgp=c(1.5,0.5,0), lwd=2)
barplot(low_spore$effect, horiz=TRUE, xlab='Impact on Biomass (%)', xlim=c(0,10), lwd.tick=2, main='Low Sporulation',
        names.arg=low_spore$metabolite, cex.names=0.9, cex.axis=0.8, cex.main=0.9,
        col=c(rep(produce_col,nrow(low_spore_efflux)),rep(consume_col,nrow(low_spore_uptake))))
text(x=9.4, y=5, as.expression(bquote(paste(italic('in silico'),' Predictions'))), srt=-90, cex=0.9)
legend('bottomright', legend=c('Consumption','Production'), bg='white', pt.lwd=1.5,
       pt.bg=c(consume_col,produce_col), pch=22, pt.cex=1.5, cex=0.9, bty='n')
box(lwd=2.5)
par(xpd=TRUE)
rect(xleft=10.5, xright=11.5, ytop=9.6, ybottom=8.8, col='deepskyblue2')
#rect(xleft=10.5, xright=11.5, ytop=8.4, ybottom=7.6, col='white')
#text(x=11, y=8, 'ND', cex=0.6)
text(x=11, y=8, 'Not\nMeasured', cex=0.6)
rect(xleft=10.5, xright=11.5, ytop=7.2, ybottom=6.4, col='deepskyblue2')
text(x=11, y=6.8, '*', cex=1.1, font=2)
rect(xleft=10.5, xright=11.5, ytop=5.95, ybottom=5.15, col='deepskyblue2')
text(x=11, y=5.5, '*', cex=1.1, font=2)
rect(xleft=10.5, xright=11.5, ytop=4.75, ybottom=3.95, col='deepskyblue2')
text(x=11, y=4.35, '*', cex=1.1, font=2)
#rect(xleft=10.5, xright=11.5, ytop=3.5, ybottom=2.7, col='white')
#text(x=11, y=3.1, 'ND', cex=0.6)
text(x=11, y=3.1, 'Not\nMeasured', cex=0.6)
rect(xleft=10.5, xright=11.5, ytop=2.3, ybottom=1.5, col='deepskyblue2')
rect(xleft=10.5, xright=11.5, ytop=1.05, ybottom=0.25, col='firebrick2')
text(x=11, y=0.65, '*', cex=1.1, font=2)
text(x=-2, y=10, 'C', cex=1.5, font=2)
par(xpd=FALSE)

par(mar=c(0,0,0,0), lwd=1.5)
plot(0, type='n', ylim=c(0,50), xlim=c(0,5), ylab='', xlab='', xaxt='n', yaxt='n', axes=FALSE)
legend('left', legend=c('Decreased during CDI', 'Increased during CDI'), cex=0.9,
       pt.bg=c('deepskyblue2','firebrick2'), pch=22, pt.cex=1.9, bty='n', 
       title=as.expression(bquote(paste(italic('in vivo'),' Metabolomics'))))

dev.off()


