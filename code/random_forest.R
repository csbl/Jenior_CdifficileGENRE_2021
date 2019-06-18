
# Set up environment
rm(list=ls())
gc()

# Load in functions
starting_dir <- getwd()
source('~/Desktop/Repositories/Jenior_Metatranscriptomics_mSphere_2018/code/R/functions.R')

# Output plot name
plot_ab <- 'results/figures/figure_2ac.pdf'
metabolome_plot <- 'results/figures/figure_2d.pdf'
otu_plot <- 'results/figures/figure_2b.pdf'

# Metabolomes
metabolome <- 'data/metabolome/scaled_intensities.log10.tsv'

# 16S
shared_otu <- 'data/16S_analysis/all_treatments.0.03.unique_list.conventional.shared'
otu_tax <- 'data/16S_analysis/formatted.all_treatments.0.03.cons.taxonomy'

# Metadata
metadata <- 'data/metadata.tsv'

# Merge by row name function
row_merge <- function(data_1, data_2){
  
  clean_merged <- merge(data_1, data_2, by='row.names', all.y=TRUE)
  rownames(clean_merged) <- clean_merged$Row.names
  clean_merged$Row.names <- NULL
  
  return(clean_merged)
}

#----------------#

# Read in data

# Metabolomes
metabolome <- read.delim(metabolome, sep='\t', header=TRUE)

# 16S
shared_otu <- read.delim(shared_otu, sep='\t', header=T, row.names=2)
shared_otu$numOtus <- NULL
shared_otu$label <- NULL
shared_otu <- shared_otu[!rownames(shared_otu) %in% c('CefC5M2'), ]  # Remove possible contaminated sample
cdf_otu <- shared_otu[,(names(shared_otu) %in% c('Otu0004','Otu0308'))] 
shared_otu <- shared_otu[,!(names(shared_otu) %in% c('Otu0004','Otu0308'))] # Remove residual C. difficile OTUs
shared_noncdf <- subset(shared_otu, rowSums(cdf_otu) != 0)
cdf_otu <- subset(cdf_otu, rowSums(cdf_otu) != 0)
cdf_percent <- as.character(round(mean(rowSums(cdf_otu) / rowSums(shared_noncdf)) * 100, 3))
rm(cdf_otu, shared_noncdf)
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
metabolome <- row_merge(metadata, metabolome)
metabolome <- subset(metabolome, abx != 'germfree')
metabolome <- subset(metabolome, infection == 'mock')
metabolome$abx <- NULL
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
shared_otu <- row_merge(otu_tax, shared_otu)
rownames(shared_otu) <- shared_otu$taxon
shared_otu$taxon <- NULL
shared_otu <- t(shared_otu)
shared_otu <- row_merge(metadata, shared_otu)
shared_otu <- subset(shared_otu, abx != 'germfree')
shared_otu <- subset(shared_otu, infection == 'mock')
shared_otu$abx <- NULL
shared_otu$infection <- NULL
rm(otu_tax)

#-------------------------------------------------------------------------------------------------------------------------#

# Feature selection

rf_feature_select <- function(training_data, column){
  
  # Calculate optimal parameters
  mTries <- round(sqrt(ncol(training_data) - 1))
  nTrees <- (ncol(training_data) - 1) * 9
  
  # Format class data
  classes <- training_data[, column]
  classes <- as.factor(classes)
  classes <- droplevels(classes)
  training_data[, column] <- NULL
  
  # Run random forest and get MDA values
  set.seed(906801)
  modelRF <- randomForest(classes ~ ., data=training_data, 
                          importance=TRUE, replace=FALSE, err.rate=TRUE, mtry=mTries, ntree=nTrees)
  featRF <- importance(modelRF, type=1)
  featRF <- as.data.frame(featRF)
  featRF$features <- rownames(featRF)

  # Filter to significant features (Strobl 2002, p <= ~0.01) and sort
  sigFeatRF <- as.data.frame(subset(featRF, featRF$MeanDecreaseAccuracy >= (abs(min(featRF$MeanDecreaseAccuracy))) * 5.0))
  sigFeatRF <- sigFeatRF[order(-sigFeatRF$MeanDecreaseAccuracy), ]
  finalFeat <- training_data[, which(rownames(sigFeatRF) %in% colnames(training_data))]
  
  return(sigFeatRF)
}

check_prediction <- function(training_data, column){

  # Get class info
  feature <- training_data[, column]
  training_data[, column] <- NULL
  feature <- as.factor(feature)
  feature <- droplevels(feature)
  
  # Run random forest
  set.seed(906801)
  randomForest(feature ~ ., data=training_data, replace=FALSE)
}

# Metabolome
colnames(metabolome) <- make.names(colnames(metabolome))
metabolome_rf <- rf_feature_select(metabolome, 'susceptibility')
# Get top features with MDA breakpoint
breakpoint <- round(median(as.numeric(format(diff(metabolome_rf$MeanDecreaseAccuracy), scientific=FALSE))),4)
breakpoint <- min(metabolome_rf$MeanDecreaseAccuracy[which(round(diff(metabolome_rf$MeanDecreaseAccuracy),4) == breakpoint)])
metabolome_rf_top <- subset(metabolome_rf, MeanDecreaseAccuracy >= breakpoint)
# Check OOB top features
save_columns <- c('susceptibility', rownames(metabolome_rf_top))
metabolome <- metabolome[, save_columns]
check_prediction(metabolome, 'susceptibility') # = 0%
# Find significant differences
susceptible_metabolome <- subset(metabolome, susceptibility == 'susceptible')
susceptible_metabolome$susceptibility <- NULL
resistant_metabolome <- subset(metabolome, susceptibility == 'resistant')
resistant_metabolome$susceptibility <- NULL
metabolome_pval <- c()
for (i in 1:ncol(susceptible_metabolome)){metabolome_pval[i] <- wilcox.test(susceptible_metabolome[,i], resistant_metabolome[,i], exact=FALSE)$p.value}
round(p.adjust(metabolome_pval, method='BH'), 4)
rm(metabolome_pval, metabolome, save_columns, breakpoint, metabolome_rf)

# 16S
colnames(shared_otu) <- make.names(colnames(shared_otu))
otu_rf <- rf_feature_select(shared_otu, 'susceptibility')
# Get top features with MDA breakpoint
breakpoint <- round(median(as.numeric(format(diff(otu_rf$MeanDecreaseAccuracy), scientific=FALSE))),4)
breakpoint <- min(otu_rf$MeanDecreaseAccuracy[which(round(diff(otu_rf$MeanDecreaseAccuracy),4) == breakpoint)])
otu_rf_top <- subset(otu_rf, MeanDecreaseAccuracy >= breakpoint)
# Check OOB for top features
save_columns <- c('susceptibility', rownames(otu_rf_top))
shared_otu <- shared_otu[, save_columns]
check_prediction(shared_otu, 'susceptibility') # = 0%
# Find significant differences
susceptible_otu <- subset(shared_otu, susceptibility == 'susceptible')
susceptible_otu$susceptibility <- NULL
resistant_otu <- subset(shared_otu, susceptibility == 'resistant')
resistant_otu$susceptibility <- NULL
otu_pval <- c()
for (i in 1:ncol(susceptible_otu)){otu_pval[i] <- wilcox.test(susceptible_otu[,i], resistant_otu[,i], exact=FALSE)$p.value}
round(p.adjust(otu_pval, method='BH'), 4)
rm(otu_pval, shared_otu, save_columns, breakpoint, otu_rf)

#-------------#

# Reformat names to be more human readable
# Metabolome
colnames(susceptible_metabolome) <- c("N-Acetylglutamine","Allantoic acid","Lactobionate","4-Hydroxyhippuric acid","Glutaconic acid",
                                      "Oxindole","5-Aminovaleric acid","Uracil","Methylmalonic acid","Undecanedioic acid","Pinitol","N-Acetylproline",
                                      "Inosine","Glutaminylleucine","Guanidinopropanoic acid","Dimethylglycine")
colnames(resistant_metabolome) <- c("N-Acetylglutamine","Allantoic acid","Lactobionic acid","4-Hydroxyhippuric acid","Glutaconic acid",
                                    "Oxindole","5-Aminovaleric acid","Uracil","Methylmalonic acid","Undecanedioic acid","Pinitol","N-Acetylproline",
                                    "Inosine","Glutaminylleucine","Guanidinopropanoic acid","Dimethylglycine")

# 16S
formatted_names <- gsub('_\\.', ' ', colnames(susceptible_otu))
genera <- sapply(strsplit(formatted_names, ' '), `[`, 1)
genera <- gsub('\\.Shigella', '', genera)
genera <- gsub('\\.unclassified', ' unclassified', genera)
genera <- gsub('\\.XlVa', ' XIVa', genera)
otu <- sapply(strsplit(formatted_names, ' '), `[`, 2)
otu <- gsub('\\.', ')', otu)
otu <- gsub('O', '(O', otu)
formatted_names <- as.vector(lapply(1:length(genera), function(i) paste(genera[i], ' ', otu[i], sep='')))
rm(genera, otu)

#-------------------------------------------------------------------------------------------------------------------------#


# Group by broader taxonomic category or metabolic pathway for plotting




# Feature selection results
# 16S
pdf(file=otu_plot, width=4, height=6)
par(mar=c(3, 1, 1.25, 1), mgp=c(1.8, 0.5, 0), xpd=FALSE, yaxs='i')
plot(1, type='n', ylim=c(0.5,nrow(susceptible_otu)+6.5), xlim=c(0,100), 
     ylab='', xlab='Relative Abundance (%)', xaxt='n', yaxt='n', cex.lab=0.9)
index <- 2
for(i in c(1:ncol(susceptible_otu))){
  stripchart(at=index+0.4, susceptible_otu[,i], 
             pch=21, bg='white', method='jitter', jitter=0.12, cex=1.5, add=TRUE)
  stripchart(at=index-0.8, resistant_otu[,i], 
             pch=21, bg='mediumorchid4', method='jitter', jitter=0.12, cex=1.5, add=TRUE)
  if (i != ncol(susceptible_otu)){
    abline(h=index+1.5, lty=2)
  }
  # Medians
  segments(median(susceptible_otu[,i]), index+0.8, median(susceptible_otu[,i]), index, lwd=2)
  segments(median(resistant_otu[,i]), index-1.2, median(resistant_otu[,i]), index-0.4, lwd=2)
  index <- index + 3
}
axis(1, at=c(0,20,40,60,80,100), labels=c(0,20,40,60,80,100), cex.axis=0.8) 
box()
text(x=c(82,82,80,81,81), y=c(3,6,9,12,15), labels=do.call(expression, formatted_names), cex=0.75)
mtext('OOB Error = 0%', side=1, cex=0.75, padj=3.4, adj=1)
mtext(c('*','*','n.s.','n.s.','n.s.'), side=4, at=c(2,5,8,11,14), # Significance
      cex=c(1.5,1.5,0.8,0.8,0.8), font=c(2,2,1,1,1), padj=c(0.25,0.25,-0.5,-0.5,-0.5))
dev.off()








# Metabolome
multiStripchart(metabolome_plot, cleared_metabolome, colonized_metabolome, metabolome_pval, metabolome_aucrf_oob, 
                'Cleared', 'Colonized', 'white', 'darkorchid4', '', 'black', 
                as.list(colnames(cleared_metabolome)), expression(paste('Scaled Intensity (',log[10],')')))









