
# Start with clean environment
rm(list=ls())
gc()

# Flux sampling files
rough1 <- '~/Desktop/tamayo_analysis/riptide_rough1/flux_samples.tsv'
rough2 <- '~/Desktop/tamayo_analysis/riptide_rough2/flux_samples.tsv'
smooth2 <- '~/Desktop/tamayo_analysis/riptide_smooth2/flux_samples.tsv'
smooth3 <- '~/Desktop/tamayo_analysis/riptide_smooth3/flux_samples.tsv'

# Read in data
rough1 <- read.delim(rough1, sep='\t', header=TRUE)
rough1$X <- NULL
rough2 <- read.delim(rough2, sep='\t', header=TRUE)
rough2$X <- NULL
smooth2 <- read.delim(smooth2, sep='\t', header=TRUE)
smooth2$X <- NULL
smooth3 <- read.delim(smooth3, sep='\t', header=TRUE)
smooth3$X <- NULL

# Combine datasets
rough_cols <- intersect(colnames(rough1), colnames(rough2))
rough1 <- rough1[,rough_cols]
rough2 <- rough2[,rough_cols]
rough <- rbind(rough1, rough2)
rm(rough_cols)
smooth_cols <- intersect(colnames(smooth2), colnames(smooth3))
smooth2 <- smooth2[,smooth_cols]
smooth3 <- smooth3[,smooth_cols]
smooth <- rbind(smooth2, smooth3)
rm(smooth_cols)

# Format row names
rough1_names <- paste('rough1_', 1:nrow(rough1), sep='')
rownames(rough1) <- rough1_names
rough2_names <- paste('rough2_', 1:nrow(rough2), sep='')
rownames(rough2) <- rough2_names
smooth2_names <- paste('smooth2_', 1:nrow(smooth2), sep='')
rownames(smooth2) <- smooth2_names
smooth3_names <- paste('smooth3_', 1:nrow(smooth3), sep='')
rownames(smooth3) <- smooth3_names
rough_names <- paste('rough_', 1:nrow(rough), sep='')
rownames(rough) <- rough_names
smooth_names <- paste('smooth_', 1:nrow(smooth), sep='')
rownames(smooth) <- smooth_names

# Create metadata
rough_metadata <- cbind(rough_names, rep('rough', length(rough_names)))
smooth_metadata <- cbind(smooth_names, rep('smooth', length(smooth_names)))
metadata <- rbind(rough_metadata, smooth_metadata)
colnames(metadata) <- c('label', 'group')
metadata <- as.data.frame(metadata)
metadata$group <- as.factor(metadata$group)
rm(rough_metadata, smooth_metadata)
rough1_metadata <- cbind(rough1_names, rep('rough1', length(rough1_names)))
rough2_metadata <- cbind(rough2_names, rep('rough2', length(rough2_names)))
rough_metadata <- rbind(rough1_metadata, rough2_metadata)
colnames(rough_metadata) <- c('label', 'group')
rough_metadata <- as.data.frame(rough_metadata)
rough_metadata$group <- as.factor(rough_metadata$group)
rm(rough1_metadata, rough2_metadata)
smooth2_metadata <- cbind(smooth2_names, rep('smooth2', length(smooth2_names)))
smooth3_metadata <- cbind(smooth3_names, rep('smooth3', length(smooth3_names)))
smooth_metadata <- rbind(smooth2_metadata, smooth3_metadata)
colnames(smooth_metadata) <- c('label', 'group')
smooth_metadata <- as.data.frame(smooth_metadata)
smooth_metadata$group <- as.factor(smooth_metadata$group)
rm(smooth2_metadata, smooth3_metadata)

# Find overlapping reactions
shared_rxns <- intersect(colnames(rough), colnames(smooth))
rough <- rough[,shared_rxns]
smooth <- smooth[,shared_rxns]
shared_rxns <- intersect(colnames(rough1), colnames(rough2))
rough1 <- rough1[,shared_rxns]
rough1 <- rough1[,shared_rxns]
shared_rxns <- intersect(colnames(smooth2), colnames(smooth3))
smooth2 <- smooth2[,shared_rxns]
smooth3 <- smooth3[,shared_rxns]
rm(shared_rxns)

# Subsample data
sub_sample <- 150
sub_sample <- sample(1:500, sub_sample, replace=FALSE)
rough <- rough[sub_sample,]
smooth <- smooth[sub_sample,]
rough1 <- rough1[sub_sample,]
rough2 <- rough2[sub_sample,]
smooth2 <- smooth2[sub_sample,]
smooth3 <- smooth3[sub_sample,]
rm(sub_sample)

# Merge data
flux_samples <- rbind(rough, smooth)
rm(rough, smooth)
rough_flux_samples <- rbind(rough1, rough2)
rm(rough1, rough2)
smooth_flux_samples <- rbind(smooth2, smooth3)
rm(smooth2, smooth3)

# Ordination analysis
library(vegan)
library(ape)
#flux_samples <- flux_samples + abs(min(flux_samples))
#flux_dist <- vegdist(flux_samples, method='bray') # Bray-Curtis
flux_dist <- designdist(flux_samples, method='1-(J/(A+B-J))', terms='quadratic', abcd=FALSE) # Theta-YC
ord_points <- as.data.frame(metaMDS(flux_dist, k=2, trymax=100)$points)
#rough_flux_samples <- rough_flux_samples + abs(min(rough_flux_samples))
#rough_flux_dist <- vegdist(rough_flux_samples, method='bray') # Bray-Curtis
rough_flux_dist <- designdist(rough_flux_samples, method='1-(J/(A+B-J))', terms='quadratic', abcd=FALSE) # Theta-YC
rough_ord_points <- as.data.frame(metaMDS(rough_flux_dist, k=2, trymax=100)$points)
#smooth_flux_samples <- smooth_flux_samples + abs(min(smooth_flux_samples))
#smooth_flux_dist <- vegdist(smooth_flux_samples, method='bray') # Bray-Curtis
smooth_flux_dist <- designdist(smooth_flux_samples, method='1-(J/(A+B-J))', terms='quadratic', abcd=FALSE) # Theta-YC
smooth_ord_points <- as.data.frame(metaMDS(smooth_flux_dist, k=2, trymax=100)$points)

# Center points
flux_x <- (abs(max(ord_points[,1])) - abs(min(ord_points[,1]))) / 2
flux_y <- (abs(max(ord_points[,2])) - abs(min(ord_points[,2]))) / 2
ord_points[,1] <- ord_points[,1] - flux_x
ord_points[,2] <- ord_points[,2] - flux_y
x_lim <- c(min(ord_points[,1])-0.004, max(ord_points[,1])+0.004)
x_lim <- round(x_lim, digits=2)
y_lim <- c(min(ord_points[,2])-0.004, max(ord_points[,2])+0.004)
y_lim <- round(y_lim, digits=2)
rough_flux_x <- (abs(max(rough_ord_points[,1])) - abs(min(rough_ord_points[,1]))) / 2
rough_flux_y <- (abs(max(rough_ord_points[,2])) - abs(min(rough_ord_points[,2]))) / 2
rough_ord_points[,1] <- rough_ord_points[,1] - rough_flux_x
rough_ord_points[,2] <- rough_ord_points[,2] - rough_flux_y
rough_x_lim <- c(min(rough_ord_points[,1])-0.004, max(rough_ord_points[,1])+0.004)
rough_x_lim <- round(rough_x_lim, digits=2)
rough_y_lim <- c(min(rough_ord_points[,2])-0.004, max(rough_ord_points[,2])+0.004)
rough_y_lim <- round(rough_y_lim, digits=2)
smooth_flux_x <- (abs(max(smooth_ord_points[,1])) - abs(min(smooth_ord_points[,1]))) / 2
smooth_flux_y <- (abs(max(smooth_ord_points[,2])) - abs(min(smooth_ord_points[,2]))) / 2
smooth_ord_points[,1] <- smooth_ord_points[,1] - smooth_flux_x
smooth_ord_points[,2] <- smooth_ord_points[,2] - smooth_flux_y
smooth_x_lim <- c(min(smooth_ord_points[,1])-0.004, max(smooth_ord_points[,1])+0.004)
smooth_x_lim <- round(smooth_x_lim, digits=2)
smooth_y_lim <- c(min(smooth_ord_points[,2])-0.004, max(smooth_ord_points[,2])+0.004)
smooth_y_lim <- round(smooth_y_lim, digits=2)

# Subset axes
rough1_ord_points <- subset(rough_ord_points, rownames(rough_ord_points) %in% rough1_names)
rough2_ord_points <- subset(rough_ord_points, rownames(rough_ord_points) %in% rough2_names)
smooth2_ord_points <- subset(smooth_ord_points, rownames(smooth_ord_points) %in% smooth2_names)
smooth3_ord_points <- subset(smooth_ord_points, rownames(smooth_ord_points) %in% smooth3_names)
rough_ord_points <- subset(ord_points, rownames(ord_points) %in% rough_names)
smooth_ord_points <- subset(ord_points, rownames(ord_points) %in% smooth_names)
rm(smooth_names, smooth2_names, smooth3_names, 
   rough_names, rough1_names, rough2_names)

# Statistical testing (permANOVA)
test <- merge(x=metadata, y=flux_samples, by.x='label', by.y='row.names')
test$label <- NULL
all_pval <- adonis2(flux_dist ~ group, data=test)
all_pval <- all_pval[[5]][1]
all_pval <- paste('p =', as.character(round(all_pval, 3)))
test <- merge(x=rough_metadata, y=rough_flux_samples, by.x='label', by.y='row.names')
test$label <- NULL
rough_pval <- adonis2(rough_flux_dist ~ group, data=test)
rough_pval <- rough_pval[[5]][1]
rough_pval <- paste('p =', as.character(round(rough_pval, 3)))
test <- merge(x=smooth_metadata, y=smooth_flux_samples, by.x='label', by.y='row.names')
test$label <- NULL
smooth_pval <- adonis2(smooth_flux_dist ~ group, data=test)
smooth_pval <- smooth_pval[[5]][1]
smooth_pval <- paste('p =', as.character(round(smooth_pval, 3)))
rm(test)
# Use previously calculated p-values
all_pval <- 'p = 0.001'
rough_pval <- 'p = 0.191'
smooth_pval <- 'p = 0.001'

# Generate figure
library(scales)
png(filename='~/Desktop/tamayo_analysis/fluxsample_ord.png', units='in', width=6, height=6, res=300)
par(mar=c(3.5,3.5,0.5,0.5), las=1, mgp=c(2.15,0.5,0))
plot(x=ord_points[,1], y=ord_points[,2], xlim=x_lim, ylim=y_lim,
     xlab='NMDS 1', ylab='NMDS 2', pch=19, cex.lab=1.4, cex=0)
points(x=smooth_ord_points[,1], y=smooth_ord_points[,2], bg=alpha('chartreuse1',0.75), pch=21, cex=2.4, lwd=1.5)
points(x=rough_ord_points[,1], y=rough_ord_points[,2], bg=alpha('blue3',0.75), pch=21, cex=2.4, lwd=1.5)
legend('topleft', legend=c('Smooth','Rough'), 
       pt.bg=c('chartreuse1','blue3'), pch=21, pt.cex=2, pt.lwd=1.5, cex=1.1, bty='n')
legend('topright', legend=all_pval, cex=0.9, pt.cex=0, bty='n')
legend('bottomright', legend='Shared Reaction Flux Dissimilarity', cex=0.9, pt.cex=0, bty='n')
box(lwd=2)
dev.off()

png(filename='~/Desktop/tamayo_analysis/fluxsample_ord_groups.png', units='in', width=12, height=6, res=300)
layout(matrix(c(1,2), nrow=1, ncol=2, byrow=TRUE))
par(mar=c(3.5,3.5,0.5,0.5), las=1, mgp=c(2.15,0.5,0), lwd=2)
plot(x=rough_ord_points[,1], y=rough_ord_points[,2], xlim=rough_x_lim, ylim=rough_y_lim,
     xlab='NMDS 1', ylab='NMDS 2', pch=19, cex.lab=1.4, cex=0)
points(x=rough1_ord_points[,1], y=rough1_ord_points[,2], bg=alpha('firebrick',0.75), pch=21, cex=2.4, lwd=1.5)
points(x=rough2_ord_points[,1], y=rough2_ord_points[,2], bg=alpha('blueviolet',0.75), pch=21, cex=2.4, lwd=1.5)
legend('topleft', legend=c('Rough 1','Rough 2'), 
       pt.bg=c('firebrick','blueviolet'), pch=21, pt.cex=2, pt.lwd=1.5, cex=1.1, bty='n')
legend('topright', legend=rough_pval, cex=0.9, pt.cex=0, bty='n')
plot(x=smooth_ord_points[,1], y=smooth_ord_points[,2], xlim=smooth_x_lim, ylim=smooth_y_lim,
     xlab='NMDS 1', ylab='NMDS 2', pch=19, cex.lab=1.4, cex=0)
points(x=smooth2_ord_points[,1], y=smooth2_ord_points[,2], bg=alpha('dodgerblue3',0.75), pch=21, cex=2.4, lwd=1.5)
points(x=smooth3_ord_points[,1], y=smooth3_ord_points[,2], bg=alpha('goldenrod2',0.75), pch=21, cex=2.4, lwd=1.5)
legend('topleft', legend=c('Smooth 2','Smooth 3'), 
       pt.bg=c('dodgerblue3','goldenrod2'), pch=21, pt.cex=2, pt.lwd=1.5, cex=1.1, bty='n')
legend('topright', legend=smooth_pval, cex=0.9, pt.cex=0, bty='n')
dev.off()

