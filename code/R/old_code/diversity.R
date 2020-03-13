
# Load dependencies
deps <- c('vegan')
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  } 
  library(dep, verbose=FALSE, character.only=TRUE)
}
rm(dep)

# 16S
shared_otu <- read.delim('~/Desktop/repos/Jenior_Metatranscriptomics_mSphere_2018/data/16S_analysis/all_treatments.0.03.unique_list.conventional.shared', 
                         sep='\t', header=T, row.names=2)
shared_otu$numOtus <- NULL
shared_otu$label <- NULL
shared_otu <- shared_otu[!rownames(shared_otu) %in% c('CefC5M2'), ]  # Remove possible contaminated sample
cdf_otu <- shared_otu[,(names(shared_otu) %in% c('Otu0004','Otu0308'))] 
shared_otu <- shared_otu[,!(names(shared_otu) %in% c('Otu0004','Otu0308'))] # Remove residual C. difficile OTUs
shared_noncdf <- subset(shared_otu, rowSums(cdf_otu) != 0)
cdf_otu <- subset(cdf_otu, rowSums(cdf_otu) != 0)
cdf_percent <- as.character(round(mean(rowSums(cdf_otu) / rowSums(shared_noncdf)) * 100, 3))
rm(cdf_otu, shared_noncdf)
otu_tax <- read.delim('~/Desktop/repos/Jenior_Metatranscriptomics_mSphere_2018/data/16S_analysis/formatted.all_treatments.0.03.cons.taxonomy', 
                      sep='\t', header=T, row.names=1)

# Metadata
metadata <- read.delim('~/Desktop/repos/Jenior_Metatranscriptomics_mSphere_2018/data/metadata.tsv', 
                       sep='\t', header=T, row.names=1)


# Format data

# Metadata
metadata$type <- NULL
metadata$cage <- NULL
metadata$mouse <- NULL
metadata$gender <- NULL
metadata$clearance <- c(rep('colonized',18), rep('cleared',18), rep('colonized',18),
                        rep('mock',12), rep('colonized',18))

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
shared_otu$susceptibility <- NULL
shared_otu$infection <- NULL
shared_otu$clearance <- NULL
rm(metadata, otu_tax)

# Subset data
noabx <- subset(shared_otu, abx == 'none')
noabx$abx <- NULL
cef <- subset(shared_otu, abx == 'cefoperazone')
cef$abx <- NULL
clinda <- subset(shared_otu, abx == 'clindamycin')
clinda$abx <- NULL
strep <- subset(shared_otu, abx == 'streptomycin')
strep$abx <- NULL
rm(shared_otu)


# Calculate inverse-Simpson diversity
noabx <- as.vector(diversity(noabx, index='invsimpson'))
cef <- as.vector(diversity(cef, index='invsimpson'))
clinda <- as.vector(diversity(clinda, index='invsimpson'))
strep <- as.vector(diversity(strep, index='invsimpson'))


# Test differences
pvals <- c(round(wilcox.test(noabx, cef, exact=FALSE)$p.value, 3),
           round(wilcox.test(noabx, clinda, exact=FALSE)$p.value, 3),
           round(wilcox.test(noabx, strep, exact=FALSE)$p.value, 3),
           round(wilcox.test(cef, clinda, exact=FALSE)$p.value, 3),
           round(wilcox.test(cef, strep, exact=FALSE)$p.value, 3),
           round(wilcox.test(clinda, strep, exact=FALSE)$p.value, 3))
pvals <- p.adjust(pvals, method='BH')           


#Generate figure
pdf(file='~/Desktop/diversity.pdf', width=5, height=4)
par(mar=c(5,4,1,1), mgp=c(2.2, 0.75, 0), las=1)
stripchart(noabx, bg='black', xlim=c(0.25,2.25), ylim=c(0,25), pch=21,
           vertical=TRUE, at=0.5, xaxt='n', ylab='Inv. Simpson Diversity', lwd=2,
           cex=1.7, method='jitter', jitter=0.12)
stripchart(strep, bg='red', xlim=c(0.25,2.25), ylim=c(0,25), pch=21,
           vertical=TRUE, at=1, xaxt='n', yaxt='n', ylab='', lwd=2,
           cex=1.7, method='jitter', jitter=0.12, add=TRUE)
stripchart(cef, bg='blue', xlim=c(0.25,2.25), ylim=c(0,25), pch=21,
           vertical=TRUE, at=1.5, xaxt='n', yaxt='n', ylab='', lwd=2,
           cex=1.7, method='jitter', jitter=0.12, add=TRUE)
stripchart(clinda, bg='green', xlim=c(0.25,2.25), ylim=c(0,25), pch=21,
           vertical=TRUE, at=2, xaxt='n', yaxt='n', ylab='', lwd=2,
           cex=1.7, method='jitter', jitter=0.12, add=TRUE)
box(lwd=2)
text(c(0.2,0.7,1.2,1.7), -9.3, adj = 0, srt=45, xpd = TRUE, font=2,
     labels=c('No Antibiotics','Streptomycin','Cefoperazone','Clindamycin'),
     col='black')
dev.off()
