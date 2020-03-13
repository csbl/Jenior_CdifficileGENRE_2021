
# Set up environment
rm(list=ls())
gc()

library(randomForest)
library(vegan)
library(ape)

# Pick which samples will be taken
small_sample <- sample(c(1:10000), 100, replace=FALSE)
large_sample <- sample(c(1:10000), 300, replace=FALSE)

# Load in  and format data
lb_aerobic_samples <- read.delim(file='/home/mjenior/Desktop/repos/Cdiff_modeling/data/flux_samples/LB_aerobic.flux_samples.format.tsv', sep='\t', header=TRUE)
colnames(lb_aerobic_samples) <- make.names(colnames(lb_aerobic_samples))
lb_aerobic_samples <- lb_aerobic_samples[small_sample,]
lb_aerobic_samples$BIOMASS_Ec_iJO1366_WT_53p95M <- NULL
lb_aerobic_samples$condition <- 'in vitro'

m9_aerobic_samples <- read.delim(file='/home/mjenior/Desktop/repos/Cdiff_modeling/data/flux_samples/M9_aerobic.flux_samples.format.tsv', sep='\t', header=TRUE)
colnames(m9_aerobic_samples) <- make.names(colnames(m9_aerobic_samples))
m9_aerobic_samples <- m9_aerobic_samples[small_sample,]
m9_aerobic_samples$BIOMASS_Ec_iJO1366_WT_53p95M <- NULL
m9_aerobic_samples$condition <- 'in vitro'

m9_anaerobic_samples <- read.delim(file='/home/mjenior/Desktop/repos/Cdiff_modeling/data/flux_samples/M9_anaerobic.flux_samples.format.tsv', sep='\t', header=TRUE)
colnames(m9_anaerobic_samples) <- make.names(colnames(m9_anaerobic_samples))
m9_anaerobic_samples <- m9_anaerobic_samples[small_sample,]
m9_anaerobic_samples$BIOMASS_Ec_iJO1366_WT_53p95M <- NULL
m9_anaerobic_samples$condition <- 'in vitro'
rm(small_sample)

# Save metadata
lb_metadata <- lb_aerobic_samples[, c('sample', 'condition')]
lb_metadata$media <- 'LB'
rownames(lb_metadata) <- lb_metadata$sample
lb_metadata$sample <- NULL
rownames(lb_aerobic_samples) <- lb_aerobic_samples$sample
lb_aerobic_samples$sample <- NULL
lb_aerobic_samples$condition <- NULL
lb_aerobic_samples <- as.data.frame(t(apply(lb_aerobic_samples, 2, as.numeric)))
m9_aerobic_metadata <- m9_aerobic_samples[, c('sample', 'condition')]
m9_aerobic_metadata$media <- 'm9_aerobic'
rownames(m9_aerobic_metadata) <- m9_aerobic_metadata$sample
m9_aerobic_metadata$sample <- NULL
rownames(m9_aerobic_samples) <- m9_aerobic_samples$sample
m9_aerobic_samples$sample <- NULL
m9_aerobic_samples$condition <- NULL
m9_aerobic_samples <- as.data.frame(t(apply(m9_aerobic_samples, 2, as.numeric)))
m9_anaerobic_metadata <- m9_anaerobic_samples[, c('sample', 'condition')]
m9_anaerobic_metadata$media <- 'm9_anaerobic'
rownames(m9_anaerobic_metadata) <- m9_anaerobic_metadata$sample
m9_anaerobic_metadata$sample <- NULL
rownames(m9_anaerobic_samples) <- m9_anaerobic_samples$sample
m9_anaerobic_samples$sample <- NULL
m9_anaerobic_samples$condition <- NULL
m9_anaerobic_samples <- as.data.frame(t(apply(m9_anaerobic_samples, 2, as.numeric)))
invitro_metadata <- as.data.frame(rbind(lb_metadata, m9_aerobic_metadata, m9_anaerobic_metadata))
rm(lb_metadata, m9_aerobic_metadata, m9_anaerobic_metadata)

# Merge in vitro samples
invitro_samples <- merge(lb_aerobic_samples, m9_aerobic_samples, by='row.names')
rownames(invitro_samples) <- invitro_samples$Row.names
invitro_samples$Row.names <- NULL
invitro_samples <- merge(invitro_samples, m9_anaerobic_samples, by='row.names')
rownames(invitro_samples) <- invitro_samples$Row.names
invitro_samples$Row.names <- NULL
invitro_samples <- as.data.frame(t(invitro_samples))
rownames(invitro_samples) <- rownames(invitro_metadata)
rm(lb_aerobic_samples, m9_aerobic_samples, m9_anaerobic_samples)

# Load and format in vivo samples
invivo_samples <- read.delim(file='/home/mjenior/Desktop/repos/Cdiff_modeling/data/flux_samples/clinda_k12.flux_samples.format.tsv', sep='\t', header=TRUE)
colnames(invivo_samples) <- make.names(colnames(invivo_samples))
invivo_samples <- invivo_samples[large_sample,]
invivo_samples$BIOMASS_Ec_iJO1366_WT_53p95M <- NULL
invivo_samples$condition <- 'in vivo'
invivo_metadata <- invivo_samples[, c('sample', 'condition')]
invivo_metadata$media <- 'none'
rownames(invivo_metadata) <- invivo_metadata$sample
invivo_metadata$sample <- NULL
invivo_samples$sample <- NULL
invivo_samples$condition <- NULL
rownames(invivo_samples) <- rownames(invivo_metadata)

# Combine datasets
invitro_samples <- as.data.frame(t(invitro_samples))
invivo_samples <- as.data.frame(t(invivo_samples))
all_samples <- merge(invitro_samples, invivo_samples, by='row.names')
rownames(all_samples) <- all_samples$Row.names
all_samples$Row.names <- NULL
all_samples <- as.data.frame(t(all_samples))
all_metadata <- as.data.frame(rbind(invivo_metadata, invitro_metadata))
rm(invivo_metadata, invitro_metadata)

# Transform data due to negative values and calculate axes
all_samples <- all_samples + max(all_samples)
samples_dist <- vegdist(all_samples, method='bray')
samples_pcoa_obj <- pcoa(samples_dist, correction='none', rn=NULL)
samples_pcoa_points <- as.data.frame(samples_pcoa_obj$vectors[,c(1,2)])
samples_pcoa_values <- samples_pcoa_obj$values
samples_pcoa_values <- samples_pcoa_values$Relative_eig[c(1,2)] * 100.0
colnames(samples_pcoa_points) <- round(samples_pcoa_values, digits=2)
rm(samples_pcoa_obj, samples_pcoa_values)

# Combine with metadata
samples_pcoa_points <- merge(samples_pcoa_points, all_metadata, by='row.names')
rownames(samples_pcoa_points) <- samples_pcoa_points$Row.names
samples_pcoa_points$Row.names <- NULL

# Calculate differences
all_samples <- merge(all_metadata, all_samples, by='row.names')
rownames(all_samples) <- all_samples$Row.names
all_samples$Row.names <- NULL
all_samples$media <- NULL
condition_permANOVA_pval <- adonis(samples_dist ~ condition, all_samples, perm=999)$aov.tab[[6]][1]
all_samples$condition <- NULL
all_samples <- merge(all_metadata, all_samples, by='row.names')
rownames(all_samples) <- all_samples$Row.names
all_samples$Row.names <- NULL
all_samples$condition <- NULL
media_permANOVA_pval <- adonis(samples_dist ~ media, all_samples, perm=999)$aov.tab[[6]][1]
rm(all_samples, samples_dist)

# Prep for plotting
samples_pcoa_points[,2] <- samples_pcoa_points[,2] + 0.001 # center points
x_axis_lab <- paste ('PC1 (', as.character(colnames(samples_pcoa_points)[1]), '%)', sep='')
y_axis_lab <- paste ('PC2 (', as.character(colnames(samples_pcoa_points)[2]), '%)', sep='')

# Subset points
vivo_points <- subset(samples_pcoa_points, condition == 'in vivo')
vitro_points <- subset(samples_pcoa_points, condition == 'in vitro')
lb_points <- subset(samples_pcoa_points, media == 'LB')
m9_aer_points <- subset(samples_pcoa_points, media == 'm9_aerobic')
m9_an_points <- subset(samples_pcoa_points, media == 'm9_anaerobic')

# Generate figures
png(filename='~/Desktop/repos/Cdiff_modeling/results/figure_5b.png', units='in', width=4, height=5, res=300)
par(mar=c(3,3,0.5,0.5), las=1, mgp=c(2,0.75,0))
plot(x=samples_pcoa_points[,2], y=samples_pcoa_points[,1], xlim=c(-0.008,0.008), ylim=c(-0.012,0.012),
     xlab=x_axis_lab, ylab=y_axis_lab, pch=19, xaxt='n', yaxt='n', cex.lab=1.2)
axis(1, at=c(-0.006,0,0.006), labels=c('-0.006','0.0','0.006'), cex.axis=0.9, cex.lab=1.2) 
axis(2, at=c(-0.01,0,0.01), labels=c('-0.01','0.0','0.01'), cex.axis=0.9, cex.lab=1.2) 
points(x=vivo_points[,2], y=vivo_points[,1], bg='orange', pch=21, cex=1.5, lwd=1.2)
points(x=vitro_points[,2], y=vitro_points[,1], bg='blue2', pch=21, cex=1.5, lwd=1.2)
legend('bottomright', legend=c(as.expression(bquote(paste(italic('p'),' = ',as.character(condition_permANOVA_pval))))), 
       pch=1, cex=1.2, pt.cex=0, bty='n')
legend('topleft', legend=c(as.expression(bquote(italic('In vitro'))), as.expression(bquote(italic('In vivo')))), 
       pt.bg=c('blue2', 'orange'), pch=21, pt.cex=1.5, cex=0.9)
box(lwd=2)
dev.off()

png(filename='~/Desktop/repos/Cdiff_modeling/results/supp_ordination.png', units='in', width=4, height=5, res=300)
par(mar=c(3,3,0.5,0.5), las=1, mgp=c(2,0.75,0))
plot(x=samples_pcoa_points[,2], y=samples_pcoa_points[,1], xlim=c(-0.008,0.008), ylim=c(-0.012,0.012),
     xlab=x_axis_lab, ylab=y_axis_lab, pch=19, xaxt='n', yaxt='n', cex.lab=1.2)
axis(1, at=c(-0.006,0,0.006), labels=c('-0.006','0.0','0.006'), cex.axis=0.9, cex.lab=1.2) 
axis(2, at=c(-0.01,0,0.01), labels=c('-0.01','0.0','0.01'), cex.axis=0.9, cex.lab=1.2)
points(x=vivo_points[,2], y=vivo_points[,1], bg='orange', pch=21, cex=1.5, lwd=1.2)
points(x=lb_points[,2], y=lb_points[,1], bg='chartreuse3', pch=21, cex=1.5, lwd=1.2)
points(x=m9_aer_points[,2], y=m9_aer_points[,1], bg='cyan3', pch=21, cex=1.5, lwd=1.2)
points(x=m9_an_points[,2], y=m9_an_points[,1], bg='darkorchid3', pch=21, cex=1.5, lwd=1.2)
legend('bottomright', legend=c(as.expression(bquote(paste(italic('p'),' = 0.001')))), 
       pch=1, cex=1.2, pt.cex=0, bty='n')
legend('topleft', legend=c(as.expression(bquote(italic('In vivo'))),'LB (Aerobic)','M9 (Aerobic)','M9 (Anaerobic)'), 
       pt.bg=c('orange','chartreuse3','cyan3','darkorchid3'), pch=21, pt.cex=1.5, cex=0.9)
box(lwd=2)
dev.off()

rm(samples_pcoa_points, vivo_points, vitro_points, lb_points, m9_aer_points, m9_an_points, none_points, all_metadata,
   media_permANOVA_pval, condition_permANOVA_pval, x_axis_lab, y_axis_lab, invitro_samples, invivo_samples)

#-----------------------------------------------------------------------------------------#

# Load and format samples
pfba_samples <- read.delim(file='/home/mjenior/Desktop/repos/Cdiff_modeling/data/flux_samples/pFBA.flux_samples.format.tsv', sep='\t', header=TRUE, row.names=1)
colnames(pfba_samples) <- make.names(colnames(pfba_samples))
pfba_samples <- pfba_samples[large_sample,]
pfba_samples$BIOMASS_Ec_iJO1366_WT_53p95M <- NULL
invivo_samples <- read.delim(file='/home/mjenior/Desktop/repos/Cdiff_modeling/data/flux_samples/clinda_k12.flux_samples.format.tsv', sep='\t', header=TRUE, row.names=1)
colnames(invivo_samples) <- make.names(colnames(invivo_samples))
invivo_samples <- invivo_samples[large_sample,]
invivo_samples$BIOMASS_Ec_iJO1366_WT_53p95M <- NULL

# Calculate relative flux within datasets (to sum of flux)
pfba_samples <- sweep(pfba_samples, 1, rowSums(pfba_samples), FUN='/')
invivo_samples <- sweep(invivo_samples, 1, rowSums(invivo_samples), FUN='/')

# Combine datasets
pfba_samples$condition <- 'pfba'
invivo_samples$condition <- 'in vivo'
pfba_samples <- as.data.frame(t(pfba_samples))
invivo_samples <- as.data.frame(t(invivo_samples))
all_samples <- merge(pfba_samples, invivo_samples, by='row.names')
rownames(all_samples) <- all_samples$Row.names
all_samples$Row.names <- NULL
all_samples <- as.data.frame(t(all_samples))

# Prep for machine learning
condition <- as.factor(as.character(all_samples$condition))
all_samples$condition <- NULL
all_samples <- as.data.frame(apply(all_samples, 2, as.numeric))

# Run AUCRF and retreive data of interest



aucrf_obj <- AUCRF(condition ~ ., data=all_samples, pdel=0.05, k0=5, ranking='MDA')
aucrf_oob <- aucrf_obj$RFopt
aucrf_oob <- aucrf_oob$err.rate
aucrf_oob <- as.character(round(median(aucrf_oob[,1]) * 100, 3))
aucrf_obj <- as.character(OptimalSet(aucrf_obj)$Name)
pfba_samples <- subset(all_samples, condition == 0)[, aucrf_obj]
pfba_samples <- as.data.frame(apply(pfba_samples, 2, as.numeric))
invivo_samples <- subset(all_samples, condition == 1)[, aucrf_obj]
invivo_samples <- as.data.frame(apply(invivo_samples, 2, as.numeric))
rm(all_samples, aucrf_obj)

# Find significant differences
pvals <- c()
pfba_samples <- as.data.frame(apply(pfba_samples, 2, as.numeric))
invivo_samples <- as.data.frame(apply(invivo_samples, 2, as.numeric))
for (i in 1:ncol(pfba_samples)){pvals[i] <- wilcox.test(pfba_samples[,i], invivo_samples[,i], exact=FALSE)$p.value}
pvals <- round(p.adjust(pvals, method='BH'), 3)
rm(i)

# Generates plot for significant differences in flux values
multiStripchart <- function(plot_file, fluxes1, fluxes2, treatments, treatment_col1, treatment_col2){
  
  pdf(file=plot_file, width=10, height=ncol(fluxes1)*1.5)
  layout(matrix(c(1:(ncol(fluxes1)+2)), nrow=(ncol(fluxes1)+2), ncol=1, byrow = TRUE))
  
  par(mar=c(0.2, 0, 0, 1), mgp=c(2.3, 0.75, 0), xpd=FALSE)
  plot(0, type='n', axes=FALSE, xlab='', ylab='', xlim=c(-10,10), ylim=c(-5,5))
  
  legend('bottomright', legend=treatments, bty='n',
         pt.bg=c(treatment_col1, treatment_col2), pch=21, cex=1.2, pt.cex=2, ncol=2)
  
  par(mar=c(0.2, 2, 0.2, 1), mgp=c(2.3, 0.75, 0), xpd=FALSE, yaxs='i')
  for(i in c(1:(ncol(fluxes1)))){
    xmax <- round(max(abs(c(fluxes1[,i], fluxes2[,i]))),3) # Adaptive x-axis limit
    plot(0, type='n', xlab='', ylab='', xaxt='n', yaxt='n', xlim=c(-xmax,xmax), ylim=c(0.3,1.8))
    abline(v=0, lty=5, lwd=1.5, col='gray35')
    stripchart(at=1.2, jitter(fluxes1[,i], amount=1e-5), 
               pch=21, bg=treatment_col1, method='jitter', jitter=0.12, cex=2, add=TRUE)
    stripchart(at=0.66, jitter(fluxes2[,i], amount=1e-5), 
               pch=21, bg=treatment_col2, method='jitter', jitter=0.12, cex=2, add=TRUE)
    box(lwd=2)
    legend('topright', legend=colnames(fluxes1)[i], pch=1, cex=1.3, pt.cex=0, bty='n')
    axis(1, at=c(-xmax,xmax), NA, cex.axis=0.8, tck=0.015)
    text(x=c(-xmax,xmax), y=0.4, labels=c(as.character(-xmax), as.character(xmax)), cex=0.9)

    segments(median(fluxes1[,i]), 1.03, median(fluxes1[,i]), 1.37, lwd=2.5)
    segments(median(fluxes2[,i]), 0.49, median(fluxes2[,i]), 0.83, lwd=2.5)
  }
  
  par(mar=c(0, 0, 0, 0))
  plot(0, type='n', axes=FALSE, xlab='', ylab='', xlim=c(-10,10), ylim=c(-5,5))
  text(x=0, y=4, labels='Sampled Flux Value', cex=1.4)
  dev.off()
}

# Generate figure
multiStripchart('~/Desktop/flux_plot.pdf', invivo_samples, pfba_samples, c('In vivo','pFBA'), 'orange', '#d817ff')

