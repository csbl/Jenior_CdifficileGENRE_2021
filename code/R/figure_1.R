
# Set up environment
source('~/Desktop/repositories/Cdiff_modeling/src/init.R')

#----------------#

# Read in data

# Metadata
metadata <- read.delim(metadata, sep='\t', header=T, row.names=1)
metadata$type <- NULL
metadata$cage <- NULL
metadata$mouse <- NULL
metadata$gender <- NULL
metadata$clearance <- NULL

# 16S
shared <- read.delim(shared, sep='\t', header=T, row.names=2)
shared$numOtus <- NULL
shared$label <- NULL
tax <- read.delim(tax, sep='\t', header=T, row.names=1)

#----------------#

# Format data
shared <- rarefyOTU(shared)
shared <- filtOTU(shared, 3)
shared <- mergeByRow(metadata, shared)
shared <- subset(shared, abx != 'germfree')
shared <- subset(shared, infection == 'mock')
shared$infection <- NULL
strep_shared <- subset(shared, abx %in% c('streptomycin','none'))
strep_shared$susceptibility <- NULL
strep_shared <- filtOTU(strep_shared, 3)
cef_shared <- subset(shared, abx %in% c('cefoperazone','none'))
cef_shared$susceptibility <- NULL
cef_shared <- filtOTU(cef_shared, 3)
clinda_shared <- subset(shared, abx %in% c('clindamycin','none'))
clinda_shared$susceptibility <- NULL
clinda_shared <- filtOTU(clinda_shared, 3)

#-------------------------------------------------------------------------------------------------------------------------#

# Theta-YC distance analysis
otu_dist <- designdist(shared[,3:ncol(shared)], method='1-(J/(A+B-J))', terms='quadratic', abcd=FALSE)

# permANOVA
susceptibility_pval <- adonis(otu_dist ~ shared$susceptibility, shared[,3:ncol(shared)], perm=999)$aov.tab[[6]][1]
susceptibility_pval <- as.character(round(susceptibility_pval, 3))
abx_pval <- adonis(otu_dist ~ shared$abx, shared[,3:ncol(shared)], perm=999)$aov.tab[[6]][1]
abx_pval <- as.character(round(abx_pval, 3))

# NMDS ordination
otu_nmds <- as.data.frame(metaMDS(otu_dist, k=2, trymax=100)$points)
otu_nmds$MDS1 <- otu_nmds$MDS1 - 0.06

# Subset points
otu_nmds <- mergeByRow(metadata, otu_nmds)
otu_nmds$infection <- NULL
susceptible_nmds <- subset(otu_nmds, susceptibility == 'susceptible')
susceptible_nmds$abx <- NULL
resistant_nmds <- subset(otu_nmds, susceptibility == 'resistant')
resistant_nmds$abx <- NULL
streptomycin_nmds <- subset(otu_nmds, abx == 'streptomycin')
streptomycin_nmds$susceptibility <- NULL
cefoperazone_nmds <- subset(otu_nmds, abx == 'cefoperazone')
cefoperazone_nmds$susceptibility <- NULL
clindamycin_nmds <- subset(otu_nmds, abx == 'clindamycin')
clindamycin_nmds$susceptibility <- NULL
noabx_nmds <- subset(otu_nmds, abx == 'none')
noabx_nmds$susceptibility <- NULL
rm(metadata)

# Calculate centroids
susceptible_centroid <- centroidMDS(susceptible_nmds)
resistant_centroid <- centroidMDS(resistant_nmds)
streptomycin_centroid <- centroidMDS(streptomycin_nmds)
cefoperazone_centroid <- centroidMDS(cefoperazone_nmds)
clindamycin_centroid <- centroidMDS(clindamycin_nmds)
noabx_centroid <- centroidMDS(noabx_nmds)

#-------------------------------------------------------------------------------------------------------------------------#

# Feature selection

# 16S
# RF feature selection
shared_rf <- rfSus(shared) # 0% OOB

# Find significant differences
res_shared <- subset(shared, susceptibility == 'resistant')
colnames(res_shared) <- make.names(colnames(res_shared))
res_shared <- res_shared[,which(colnames(res_shared) %in% rownames(shared_rf))]
sus_shared <- subset(shared, susceptibility == 'susceptible')
colnames(sus_shared) <- make.names(colnames(sus_shared))
sus_shared <- sus_shared[,which(colnames(sus_shared) %in% rownames(shared_rf))]
rm(shared)
shared_pval <- c()
for (i in 1:ncol(res_shared)){shared_pval[i] <- wilcox.test(res_shared[,i], sus_shared[,i], exact=FALSE)$p.value}
shared_pval <- as.numeric(p.adjust(shared_pval, method='BH'))

# Subset to significant values
res_shared <- res_shared[,which(shared_pval <= 0.05)]
sus_shared <- sus_shared[,which(shared_pval <= 0.05)]

# Subset to OTUs with greater median abundance in resistant communities
keep <- c()
medians <- c()
y <- 1
for (x in 1:ncol(res_shared)) {
  if (median(res_shared[,x]) > median(sus_shared[,x])) {
    medians[y] <- median(res_shared[,x]) - median(sus_shared[,x])
    keep[y] <- x
    y <- y + 1
  }
}
res_shared <- res_shared[,keep]
sus_shared <- sus_shared[,keep]
rm(x, y, keep)

# Create importance + abundance change table
change <- as.data.frame(cbind(colnames(res_shared), medians))
colnames(change) <- c('otu','diff')
imp_diff <- merge(shared_rf, change, by.x='row.names', by.y='otu')
rownames(imp_diff) <- imp_diff$Row.names
imp_diff$Row.names <- NULL
imp_diff$MeanDecreaseAccuracy <- as.numeric(as.character(imp_diff$MeanDecreaseAccuracy))
imp_diff$diff <- as.numeric(as.character(imp_diff$diff))
imp_diff$diff <- log10(imp_diff$diff)
rm(shared_rf, change)

# Log transform values
#res_shared <- log10(res_shared + 1)
#sus_shared <- log10(sus_shared + 1)

# Reformat taxonomy
tax <- tax[which(rownames(tax) %in% colnames(res_shared)), ]
tax$genus_otu <- paste(tax$genus, rownames(tax), sep='_')
tax$phylum <- NULL
tax$class <- NULL
tax$order <- NULL
tax$family <- NULL
tax$genus <- NULL

# Replace row names with genera + OTUs
res_shared <- t(res_shared)
res_shared <- mergeByRow(tax, res_shared)
rownames(res_shared) <- res_shared$genus_otu
res_shared$genus_otu <- NULL
res_shared <- as.data.frame(t(res_shared))
sus_shared <- t(sus_shared)
sus_shared <- mergeByRow(tax, sus_shared)
rownames(sus_shared) <- sus_shared$genus_otu
sus_shared$genus_otu <- NULL
sus_shared <- as.data.frame(t(sus_shared))

#---------#

# Individual pretreatment groups
temp <- rfAbx(strep_shared)
strep_rf <- temp$sigFeatRF
strep_shared <- temp$subset_data
strep <- subset(strep_shared, abx == 'streptomycin')
strep$abx <- NULL
strep_ctrl <- subset(strep_shared, abx == 'none')
strep_ctrl$abx <- NULL
rm(strep_shared)
strep_diff <- sigDiff(strep_ctrl, strep)
rm(strep_ctrl, strep)
temp <- rfAbx(cef_shared)
cef_rf <- temp$sigFeatRF
cef_shared <- temp$subset_data
cef <- subset(cef_shared, abx == 'cefoperazone')
cef$abx <- NULL
cef_ctrl <- subset(cef_shared, abx == 'none')
cef_ctrl$abx <- NULL
rm(cef_shared)
cef_diff <- sigDiff(cef_ctrl, cef)
rm(cef_ctrl, cef)
temp <- rfAbx(clinda_shared)
clinda_rf <- temp$sigFeatRF
clinda_shared <- temp$subset_data
clinda <- subset(clinda_shared, abx == 'clindamycin')
clinda$abx <- NULL
clinda_ctrl <- subset(clinda_shared, abx == 'none')
clinda_ctrl$abx <- NULL
rm(clinda_shared)
clinda_diff <- sigDiff(clinda_ctrl, clinda)
rm(clinda_ctrl, clinda)
rm(temp)



strep_diff <- merge(strep_diff, strep_rf, by.x='names', by.y='features')
rownames(strep_diff) <- strep_diff$names
strep_diff$names <- NULL
strep_diff$diff <- log10(as.numeric(as.character(strep_diff$diff)))
rm(strep_rf)


cef_diff <- merge(cef_diff, cef_rf, by.x='names', by.y='features')
rownames(cef_diff) <- cef_diff$names
cef_diff$names <- NULL
cef_diff$diff <- log10(as.numeric(as.character(cef_diff$diff)))
rm(cef_rf)



clinda_diff <- merge(clinda_diff, clinda_rf, by.x='names', by.y='features')
rownames(clinda_diff) <- clinda_diff$names
clinda_diff$names <- NULL
clinda_diff$diff <- log10(as.numeric(as.character(clinda_diff$diff)))
rm(clinda_rf)



# Streptomycin vs Resistant
summary(lm(diff ~ MeanDecreaseAccuracy, data=strep_diff))
# Cefoperazone vs Resistant
summary(lm(diff ~ MeanDecreaseAccuracy, data=cef_diff))
# Clindamycin vs Resistant
summary(lm(diff ~ MeanDecreaseAccuracy, data=clinda_diff))
# All susceptible vs Resistant
summary(lm(diff ~ MeanDecreaseAccuracy, data=imp_diff))





#-------------------------------------------------------------------------------------------------------------------------#

# Plot
pdf(file=fig1Plot, width=10, height=5)
layout(matrix(c(1,2), nrow=1, ncol=2, byrow = TRUE))
par(mar=c(4,4,1,1), las=1, mgp=c(2.5,0.75,0), las=1)

plot(x=otu_nmds$MDS1, y=otu_nmds$MDS2, xlim=c(-0.5,0.5), ylim=c(-0.5,0.5),
     xlab='NMDS axis 1', ylab='NMDS axis 2', pch=19, cex.axis=1.2, cex.lab=1.2, cex=0.2)
segments(x0=resistant_nmds$MDS1, y0=resistant_nmds$MDS2, 
         x1=resistant_centroid[1,2], y1=resistant_centroid[1,3], col='grey25', lwd=1.2)
points(x=resistant_nmds$MDS1, y=resistant_nmds$MDS2, bg=resistant_col, pch=21, cex=1.8, lwd=1.5)
segments(x0=susceptible_nmds$MDS1, y0=susceptible_nmds$MDS2, 
         x1=susceptible_centroid[1,2], y1=susceptible_centroid[1,3], col='gray25', lwd=1.2)
points(x=susceptible_nmds$MDS1, y=susceptible_nmds$MDS2, bg=susceptible_col, pch=21, cex=1.8, lwd=1.5)
mtext('A', side=2, line=2, las=2, adj=1, padj=-9.5, cex=1.7, font=2)
legend('bottomright', legend=c('Resistant','Susceptible'), 
       pt.bg=c(resistant_col, susceptible_col), pch=21, pt.cex=2, bty='n')
legend('bottomleft', legend=as.expression(bquote(paste(italic('p'),' = ', .(susceptibility_pval), '***'))), 
       pch=1, cex=1.2, pt.cex=0, bty='n')
box()

plot(x=otu_nmds$MDS1, y=otu_nmds$MDS2, xlim=c(-0.5,0.5), ylim=c(-0.5,0.5),
     xlab='NMDS axis 1', ylab='NMDS axis 2', pch=19, cex.axis=1.2, cex.lab=1.2, cex=0.2)
segments(x0=streptomycin_nmds$MDS1, y0=streptomycin_nmds$MDS2, 
         x1=streptomycin_centroid[1,2], y1=streptomycin_centroid[1,3], col='grey25', lwd=1.2)
points(x=streptomycin_nmds$MDS1, y=streptomycin_nmds$MDS2, bg=streptomycin_col, pch=21, cex=1.8, lwd=1.5)
segments(x0=cefoperazone_nmds$MDS1, y0=cefoperazone_nmds$MDS2, 
         x1=cefoperazone_centroid[1,2], y1=cefoperazone_centroid[1,3], col='grey25', lwd=1.2)
points(x=cefoperazone_nmds$MDS1, y=cefoperazone_nmds$MDS2, bg=cefoperazone_col, pch=21, cex=1.8, lwd=1.5)
segments(x0=clindamycin_nmds$MDS1, y0=clindamycin_nmds$MDS2, 
         x1=clindamycin_centroid[1,2], y1=clindamycin_centroid[1,3], col='grey25', lwd=1.2)
points(x=clindamycin_nmds$MDS1, y=clindamycin_nmds$MDS2, bg=clindamycin_col, pch=21, cex=1.8, lwd=1.5)
segments(x0=noabx_nmds$MDS1, y0=noabx_nmds$MDS2, 
         x1=noabx_centroid[1,2], y1=noabx_centroid[1,3], col='grey25', lwd=1.2)
points(x=noabx_nmds$MDS1, y=noabx_nmds$MDS2, bg=noabx_col, pch=21, cex=1.8, lwd=1.5)
mtext('B', side=2, line=2, las=2, adj=1, padj=-9.5, cex=1.7, font=2)
legend('bottomright', legend=c('Untreated','Cefoperazone','Streptomycin','Clindamycin'), 
       pt.bg=c(noabx_col, cefoperazone_col, streptomycin_col, clindamycin_col), pch=21, pt.cex=2, bty='n')
legend('bottomleft', legend=as.expression(bquote(paste(italic('p'),' = ', .(abx_pval), '***'))), 
       pch=1, cex=1.2, pt.cex=0, bty='n')
box()

dev.off()

#----------------#



par(mar=c(4,4,1,1), las=1, mgp=c(2.5,0.75,0), las=1)

plot(x=imp_diff$MeanDecreaseAccuracy, y=imp_diff$diff, xlim=c(0,4), ylim=c(0,4), cex=0.8,
     xlab='Mean Decrease Accuracy', ylab='Log10 Abundance Decrease', pch=19, cex.axis=1.2, cex.lab=1.2)


abline(lm(imp_diff$diff ~ imp_diff$MeanDecreaseAccuracy), col="red", lwd=2)


box()


plot(fit)


# Porphyromonadaceae_unclassified(100)_Otu0005 -  Muribaculum intestinale strain YL27 16S ribosomal RNA, partial sequence 416 	416 	100% 	2e-116	96%
# Porphyromonadaceae_unclassified(100)_Otu0009 -  Muribaculum intestinale strain YL27 16S ribosomal RNA, partial sequence 350 	350 	100% 	3e-96	92%
# Alistipes(100)_Otu0012 - Alistipes onderdonkii strain JCM 16771 16S ribosomal RNA gene, partial sequence 396 	396 	100% 	3e-110	95%
# Lachnospiraceae_unclassified(100)_Otu0014 -  Clostridium aerotolerans strain DSM 5434 16S ribosomal RNA gene, partial sequence 361 	361 	100% 	1e-99	93%
# Porphyromonadaceae_unclassified(100)_Otu0017 -  Muribaculum intestinale strain YL27 16S ribosomal RNA, partial sequence 466 	466 	100% 	2e-131	100%
# Porphyromonadaceae_unclassified(100)_Otu0013 -  Muribaculum intestinale strain YL27 16S ribosomal RNA, partial sequence 377 	377 	100% 	1e-104	94%




#-------------------------------------------------------------------------------------------------------------------------#

# Clean up
#for (dep in deps){
#  pkg <- paste('package:', dep, sep='')
#  detach(pkg, character.only = TRUE)
#}
#setwd(starting_dir)
#rm(list=ls())
#gc()
