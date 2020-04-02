
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

# Subsample data
sub_sample <- sample(1:500, 100, replace=FALSE)
rough1 <- rough1[sub_sample,]
rough2 <- rough2[sub_sample,]
smooth2 <- smooth2[sub_sample,]
smooth3 <- smooth3[sub_sample,]
rm(sub_sample)

# Combine datasets
rough_cols <- intersect(colnames(rough1), colnames(rough2))
rough1 <- rough1[,rough_cols]
rough2 <- rough2[,rough_cols]
rough <- rbind(rough1, rough2)
rm(rough_cols, rough1, rough2)
smooth_cols <- intersect(colnames(smooth2), colnames(smooth3))
smooth2 <- smooth2[,smooth_cols]
smooth3 <- smooth3[,smooth_cols]
smooth <- rbind(smooth2, smooth3)
rm(smooth_cols, smooth2, smooth3)

# Format row names
rownames(rough) <- paste('rough_', 1:nrow(rough), sep='')
rownames(smooth) <- paste('smooth_', 1:nrow(smooth), sep='')

# Find overlapping reactions
shared_rxns <- intersect(colnames(rough), colnames(smooth))
rough <- rough[,shared_rxns]
smooth <- smooth[,shared_rxns]
rm(shared_rxns)

# Merge and format data
rough$condition <- 'rough'
smooth$condition <- 'smooth'
all_samples <- rbind(rough, smooth)

# Format class data
conditions <- all_samples$condition
conditions <- as.factor(conditions)
conditions <- droplevels(conditions)
all_samples$condition <- NULL

# Calculate optimal parameters
mTries <- round(sqrt(ncol(all_samples) - 1))
nTrees <- (ncol(all_samples) - 1) * 9

# Run random forest and get MDA values
set.seed(906801)
library(randomForest)
modelRF <- randomForest(conditions ~ ., data=all_samples, 
                          importance=TRUE, replace=FALSE, err.rate=TRUE, mtry=mTries, ntree=nTrees)
importanceRF <- as.data.frame(importance(modelRF))
rm(modelRF, mTries, nTrees, conditions, all_samples)

# Filter to significant features (Strobl 2002) and sort
sig_importanceRF <- as.data.frame(subset(importanceRF, importanceRF$MeanDecreaseAccuracy >= abs(min(importanceRF$MeanDecreaseAccuracy))))
sig_importanceRF <- sig_importanceRF[order(-sig_importanceRF$MeanDecreaseAccuracy), ]
rm(importanceRF)

# Compare flux distributions
select_sig_importanceRF <- sig_importanceRF[1:6,] # Large jump in MDA and Gini for top 6
sig_reactions <- rownames(select_sig_importanceRF)
rough <- rough[, sig_reactions]
smooth <- smooth[, sig_reactions]
rxn_names <- cbind(c('rxn07124_c', 'ID008_c', 'EX_cpd03170_e', 'EX_cpd00339_e', 'rxn12566_c', 'rxn20606_c'),
                   c('p-Hydroxyphenylacetate decarboxylase','2,3-Dihydroxyphenylpropanoate diffusion',
                     '4-Hydroxymandelate exchange','5-Aminopentanoate exchange',
                     '5-Aminopentanoate transport via proton symport','D-proline reductase'))
rxn_names <- as.data.frame(rxn_names)
colnames(rxn_names) <- c('id', 'name')
select_sig_importanceRF <- merge(x=select_sig_importanceRF, y=rxn_names, by.x='row.names', by.y='id')
select_sig_importanceRF <- select_sig_importanceRF[order(-select_sig_importanceRF$MeanDecreaseAccuracy), ]
rm(rxn_names)

# Test for significant differences
phase_pval <- c()
for (x in 1:ncol(rough)){phase_pval[x] <- wilcox.test(rough[,x], smooth[,x], exact=FALSE)$p.value}
phase_pval <- p.adjust(phase_pval, method='BH')



# Generate figures
png(filename='~/Desktop/tamayo_analysis/mda.png', units='in', width=6, height=10, res=300)
par(mar=c(3, 1, 1, 1), mgp=c(1.8, 0.5, 0), xpd=FALSE, yaxs='i', lwd=2)
dotchart(rev(sig_importanceRF$MeanDecreaseAccuracy), labels=rev(rownames(sig_importanceRF)),
         pch=21, bg='chartreuse3', pt.cex=1.8, xlim=c(0,30),xlab='MDA')
dev.off()

png(filename='~/Desktop/tamayo_analysis/select_mda.png', units='in', width=6, height=5, res=300)
par(mar=c(3, 1, 1, 1), mgp=c(1.8, 0.5, 0), xpd=FALSE, yaxs='i', lwd=2)
dotchart(rev(select_sig_importanceRF$MeanDecreaseAccuracy), labels=rev(select_sig_importanceRF$name),
         pch=21, bg='chartreuse3', pt.cex=1.8, xlim=c(20,30),xlab='MDA')
dev.off()

png(filename='~/Desktop/tamayo_analysis/fluxes.png', units='in', width=8, height=6, res=300)
par(mar=c(3, 1, 1, 1), mgp=c(1.8, 0.5, 0), xpd=FALSE, yaxs='i', lwd=2)
plot(1, type='n', ylim=c(0.5,18.5), xlim=c(-1000,1000), 
     ylab='', xlab='Sampled Flux', yaxt='n', cex.lab=0.9)
abline(v=0, lwd=3, lty=3, col='gray')
box()
index <- 2
for(i in c(1:ncol(rough))){
  stripchart(at=index+0.4, rough[,i], 
             pch=21, bg='white', method='jitter', jitter=0.12, cex=1.5, add=TRUE)
  stripchart(at=index-0.8, smooth[,i], 
             pch=21, bg='mediumorchid4', method='jitter', jitter=0.12, cex=1.5, add=TRUE)
  if (i != ncol(rough)){
    abline(h=index+1.5, lty=2)
  }
  segments(median(rough[,i]), index+0.8, median(rough[,i]), index, lwd=3)
  segments(median(smooth[,i]), index-1.2, median(smooth[,i]), index-0.4, lwd=3)
  index <- index + 3
}
text(x=-550, y=c(3,6,9,12,15,18), labels=rev(select_sig_importanceRF$name))
mtext('OOB Error = 0%', side=1, cex=0.75, padj=3.4, adj=1)
mtext(c('***','***','***','***','***','***'), side=4, at=c(2,5,8,11,14,17), cex=1.5, font=2, padj=0.25)
dev.off()

