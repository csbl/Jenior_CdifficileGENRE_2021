
# Flux sampling files
clinda_samples <- '~/Desktop/repos/Jenior_Cdifficile_2019/data/contextualized_iCdG698/riptide_high_spore/flux_samples.tsv'
strep_samples <- '~/Desktop/repos/Jenior_Cdifficile_2019/data/contextualized_iCdG698/riptide_low_spore/flux_samples.tsv'

# Read in data
clinda_samples <- read.delim(clinda_samples, sep='\t', header=TRUE)
clinda_samples$X <- NULL
strep_samples <- read.delim(strep_samples, sep='\t', header=TRUE)
strep_samples$X <- NULL

# Format data
overlap <- intersect(colnames(strep_samples), colnames(clinda_samples))
clinda_samples <- clinda_samples[, overlap]
strep_samples <- strep_samples[, overlap]
rm(overlap)

# Subsample data
sample_size <- min(c(nrow(clinda_samples), nrow(strep_samples), 250))
sub_sample <- sample(1:min(c(nrow(clinda_samples), nrow(strep_samples))), sample_size, replace=FALSE)
clinda_samples <- clinda_samples[sub_sample,]
strep_samples <- strep_samples[sub_sample,]
rm(sample_size, sub_sample)

# Format row names
clinda_names <- paste('clinda_', 1:nrow(clinda_samples), sep='')
rownames(clinda_samples) <- clinda_names
strep_names <- paste('strep_', 1:nrow(strep_samples), sep='')
rownames(strep_samples) <- strep_names

# Create metadata
clinda_metadata <- cbind(clinda_names, rep('clinda', length(clinda_names)), rep('cleared', length(clinda_names)))
strep_metadata <- cbind(strep_names, rep('strep', length(strep_names)), rep('colonized', length(strep_names)))
metadata <- rbind(clinda_metadata, strep_metadata)
colnames(metadata) <- c('label', 'treatment','clearance')
metadata <- as.data.frame(metadata)
rm(clinda_metadata, strep_metadata)

# Merge data and prep for unsupervised learning
all_samples <- rbind(clinda_samples, strep_samples)
all_samples <- all_samples + abs(min(all_samples))
    
# Calculate dissimilarity (Bray-Curtis)
library(vegan)
flux_dist <- vegdist(all_samples, method='bray') # Bray-Curtis

# Mean wwithin-group dissimilarity
flux_groups <- as.factor(c(rep('clinda',nrow(clinda_samples)), rep('strep',nrow(strep_samples))))
meandist(flux_dist, grouping=flux_groups)
#        clinda      strep
# clinda 0.01981804 0.03645153
# strep  0.03645153 0.02875940

# Unsupervised learning (NMDS)
flux_nmds <- as.data.frame(metaMDS(flux_dist, k=2, trymax=25)$points)

# Center points
flux_x <- (abs(max(flux_nmds$MDS1)) - abs(min(flux_nmds$MDS1))) / 2
flux_y <- (abs(max(flux_nmds$MDS2)) - abs(min(flux_nmds$MDS2))) / 2
flux_nmds$MDS1 <- flux_nmds$MDS1 - flux_x
flux_nmds$MDS2 <- flux_nmds$MDS2 - flux_y
flux_xlim <- max(abs(max(flux_nmds$MDS1)), abs(min(flux_nmds$MDS1))) + 0.01
flux_ylim <- max(abs(max(flux_nmds$MDS2)), abs(min(flux_nmds$MDS2))) + 0.01

# Subset axes
clinda_nmds_points <- subset(flux_nmds, rownames(flux_nmds) %in% clinda_names)
strep_nmds_points <- subset(flux_nmds, rownames(flux_nmds) %in% strep_names)
rm(clinda_names, strep_names)

# Statistical testing (permANOVA)
test <- merge(x=metadata, y=all_samples, by.x='label', by.y='row.names')
rownames(test) <- test$label
test$label <- NULL
test$treatment <- NULL
spore_pval <- adonis(flux_dist ~ clearance, data=test, perm=99, method='bray')
spore_pval <- spore_pval$aov.tab[[6]][1]
spore_pval <- as.character(round(spore_pval, 3))
rm(all_samples, test, flux_dist, metadata)

#-------------------------------------#

# Limit to Biomass flux
clinda_samples <- clinda_samples[,'biomass']
strep_biomass <- strep_samples[,'biomass']

# Test differences
biomass_pval <- round(wilcox.test(clinda_samples, strep_biomass, exact=FALSE)$p.value, 3)

# Convert to doubling time
clinda_doubling <- (1 / clinda_samples) * 3600
strep_doubling <- (1 / strep_biomass) * 3600

# Center on plot area
clinda_doubling <- clinda_doubling - 35
strep_doubling <- strep_doubling - 35

#-------------------------------------#


# Flux sampling files
clinda_samples <- '~/Desktop/repos/Jenior_Cdifficile_2019/data/contextualized_iCdG698/riptide_high_spore/flux_samples.tsv'
strep_samples <- '~/Desktop/repos/Jenior_Cdifficile_2019/data/contextualized_iCdG698/riptide_low_spore/flux_samples.tsv'

# Read in data
clinda_samples <- read.delim(clinda_samples, sep='\t', header=TRUE)
clinda_samples$X <- NULL
strep_samples <- read.delim(strep_samples, sep='\t', header=TRUE)
strep_samples$X <- NULL

# Subset to shared reactions
overlap <- intersect(colnames(strep_samples), colnames(clinda_samples))
clinda_samples <- clinda_samples[, overlap]
strep_samples <- strep_samples[, overlap]
rm(overlap)

# Merge data for supervised learning
clinda_samples$condition <- 1
strep_samples$condition <- 0
all_samples <- rbind(clinda_samples, strep_samples)
all_samples$condition <- as.factor(all_samples$condition)
rm(clinda_samples, strep_samples)

# Run AUCRF and obtain feature lists
library(AUCRF)
set.seed(906801)
all_aucrf <- AUCRF(condition ~ ., data=all_samples, pdel=0, k0=10)
print(all_aucrf)
#rm(all_samples)

# Assemble feature table
top_rxns_importance <- all_aucrf$ranking[1:all_aucrf$Kopt]
rf_rxns <- as.data.frame(cbind(labels(top_rxns_importance), as.vector(top_rxns_importance)))
colnames(rf_rxns) <- c('id','mda')
rf_rxns$mda <- as.numeric(as.character(rf_rxns$mda))
write.table(rf_rxns, file='~/Desktop/repos/Jenior_Cdifficile_2019/data/aucrf_invivo.tsv', 
            quote=FALSE, sep='\t', row.names=FALSE, col.names=TRUE)
#rm(all_aucrf, top_rxns_importance)

# Import pre-run data
rf_rxns <- read.delim('/home/mjenior/Desktop/repos/Jenior_Cdifficile_2019/data/aucrf_invivo.tsv', sep='\t', header=TRUE)
rf_rxns$name <- as.character(rf_rxns$name)
rf_rxns$name <- gsub('_', ' ', rf_rxns$name)
rf_rxns$id <- as.character(rf_rxns$id)
rf_rxns$label <- paste0(rf_rxns$name , ' (', rf_rxns$id,')')
rf_rxns$mda <- round(as.numeric(as.character(rf_rxns$mda)), 2)

# Rank and subset data for plotting
rf_rxns <- rf_rxns[order(rf_rxns$mda),]
rf_clinda_samples <- clinda_samples[, rf_rxns$id]
rf_strep_samples <- strep_samples[, rf_rxns$id]

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

select_and_plot <- function(metabolite_name, best_ylim=0, correction=FALSE, title=NA, across=FALSE, title_cex=0.75) {
    metabolite <- metabolome[, c(1,2,which(colnames(metabolome) %in% c(metabolite_name)))]
    metabolite <- subset(metabolite, abx %in% c('cefoperazone','clindamycin'))
    colnames(metabolite) <- c('abx','infection','intensity')
    
    if (correction == TRUE) {metabolite$intensity <- metabolite$intensity - min(metabolite$intensity)}
    if (is.na(title)) {title <- metabolite_name}
    
    a_mock <- subset(subset(metabolite, abx=='clindamycin'), infection=='mock')[,3]
    a_630 <- subset(subset(metabolite, abx=='clindamycin'), infection=='630')[,3]
    b_mock <- subset(subset(metabolite, abx=='cefoperazone'), infection=='mock')[,3]
    b_630 <- subset(subset(metabolite, abx=='cefoperazone'), infection=='630')[,3]
    colnames(metabolite) <- c('abx','infection','intensity')

    if (best_ylim == 0) {best_ylim <- ceiling(as.numeric(quantile(metabolite[,3], 0.95)))}
    if (across == TRUE) {best_ylim <- best_ylim + (best_ylim * 0.5)}
    
    par(mar=c(2.6,2.8,1.5,0.5), xpd=FALSE, mgp=c(1.6,0.7,0), lwd=1.5, xaxt='n', las=1)
    boxplot(a_mock, at=1, xlim=c(0.6,3.4), ylim=c(0,best_ylim), col='firebrick2', cex.axis=0.7, cex.lab=0.8, 
            xlab="", ylab="Scaled Intensity (Log10)", outcex=0, whisklty=1, medlwd=2,  main=title, cex.main=title_cex) 
    boxplot(a_630, at=1.5, xlim=c(0.6,3.4), ylim=c(0,best_ylim), col='firebrick2', 
            xlab="", ylab="", outcex=0, whisklty=1, medlwd=2, yaxt='n', add=TRUE) 
    if (across == FALSE) {abline(v=2)}
    boxplot(b_mock, at=2.5, xlim=c(0.6,3.4), ylim=c(0,best_ylim), col='blue4', 
            xlab="", ylab="", outcex=0, whisklty=1, medlwd=2, yaxt='n', add=TRUE) 
    boxplot(b_630, at=3, xlim=c(0.6,3.4), ylim=c(0,best_ylim), col='blue4', 
            xlab="", ylab="", outcex=0, whisklty=1, medlwd=2, yaxt='n', add=TRUE) 
    mtext(c('High\nSpores','Low\nSpores'), side=1, at=c(1.25,2.75), padj=1.3, cex=0.65)
    mtext(c('CDI:','-','+','-','+'), side=1, at=c(0.5,1,1.5,2.5,3), padj=0.3, cex=c(0.5,0.8,0.8,0.8,0.8))

    a_pval <- round(wilcox.test(a_mock, a_630, exact=FALSE)$p.value,3)
    b_pval <- round(wilcox.test(b_mock, b_630, exact=FALSE)$p.value,3)
    ab_pval <- round(wilcox.test(a_mock, b_mock, exact=FALSE)$p.value,3)
    ba_pval <- round(wilcox.test(a_630, b_630, exact=FALSE)$p.value,3)
    if (across == TRUE) {fctr <- 0.7} else {fctr <- 0.96}
    if (a_pval <= 0.05) {
        segments(x0=1, y0=best_ylim*fctr, x1=1.5, y1=best_ylim*fctr, lwd=1.5)
        text(x=1.25, y=best_ylim*(fctr+0.04), '*', font=2, cex=1.2)}
    if (b_pval <= 0.05) {
        segments(x0=2.5, y0=best_ylim*fctr, x1=3, y1=best_ylim*fctr, lwd=1.5)
        text(x=2.75, y=best_ylim*(fctr+0.04), '*', font=2, cex=1.2)}
    if (across == TRUE) {
        if (ab_pval <= 0.05) {
            segments(x0=1, y0=best_ylim*0.96, x1=2.5, y1=best_ylim*0.96, lwd=1.5)
            text(x=1.75, y=best_ylim, '*', font=2, cex=1.2)}
        if (ba_pval <= 0.05) {
            segments(x0=1.5, y0=best_ylim*0.86, x1=3, y1=best_ylim*0.86, lwd=1.5)
            text(x=2.25, y=best_ylim*0.9, '*', font=2, cex=1.2)}
    }
}

#---------------------------------------------------------------------#

# Generate figure panels

# Generate figure
library(vioplot)
library(plotrix)

# Doubling time
png(filename='~/Desktop/repos/Jenior_Cdifficile_2019/results/figures/Figure_4A.png', 
    units='in', width=2, height=4, res=300)
par(mar=c(3,3,0.5,0.5), xpd=FALSE, las=1, mgp=c(2,0.8,0), lwd=2.5)
boxplot(clinda_doubling, at=0.5, xlim=c(0,2), ylab='Predicted Doubling Time (min)', 
        ylim=c(0,70), yaxt='n', col='firebrick2', 
        boxlwd=2, medlwd=2, staplelwd=2, whisklwd=2, whisklty=1, width=0.9)
boxplot(strep_doubling, at=1.5, xlim=c(0,2), add=TRUE, yaxt='n', col='blue3',
        boxlwd=2, medlwd=2, staplelwd=2, whisklwd=2, whisklty=1, width=0.9)
axis(side=2, at=seq(0,70,10), labels=c(0,seq(20,80,10)), cex.axis=0.9)
axis.break(2, 5, style='slash')
segments(x0=0.5, y0=65, x1=1.5)
text(x=1, y=68, '***', font=2, cex=1.5)
par(xpd=TRUE)
text(x=c(0.5,1.5), y=-10, labels=c('High\nSpores','Low\nSpores'), cex=1)
par(xpd=FALSE)
dev.off()

# Ordination of shared reaction flux samples
png(filename='~/Desktop/repos/Jenior_Cdifficile_2019/results/figures/Figure_4B.png', 
    units='in', width=4.5, height=4.5, res=300)
library(scales)
par(mar=c(3.5,3.5,0.5,0.5), las=1, mgp=c(2.2,0.7,0), lwd=2)
plot(x=flux_nmds$MDS1, y=flux_nmds$MDS2, xlim=c(-0.035,0.035), ylim=c(-0.035, 0.035),
     xlab='NMDS Axis 1', ylab='NMDS Axis 2', pch=19, cex.lab=1.1, cex=0, cex.axis=0.9)
points(x=strep_nmds_points$MDS1, y=strep_nmds_points$MDS2, bg=alpha('blue3',0.8), pch=21, cex=1.7)
points(x=clinda_nmds_points$MDS1, y=clinda_nmds_points$MDS2, bg=alpha('firebrick2',0.8), pch=21, cex=1.7)
legend('topright', legend=c('Low Sporulation','High Sporulation'), 
       pt.bg=c(alpha('blue3',0.8),alpha('firebrick2',0.8)), pch=21, pt.cex=1.7, box.lwd=2)
legend('bottomleft', legend=substitute(paste(italic('p'), '-value = 0.001***')), pt.cex=0, bty='n')
legend('bottomright', legend='Core Metabolism', pt.cex=0, bty='n')
box(lwd=2.5)
dev.off()

# Flux samples from AUCRF features
library(plotrix)
png(filename='~/Desktop/repos/Jenior_Cdifficile_2019/results/figures/Figure_4C.png', 
    units='in', width=2.5, height=4, res=300)
par(mar=c(2.5, 0.5, 0.5, 0.5), mgp=c(1.4, 0.5, 0), xpd=FALSE, lwd=1.7, xaxt='n')
dotchart(rf_rxns$mda, color='darkorchid4', xlab='Mean Decrease Accuracy (%)', xlim=c(15,30),  
         pch=16, lwd=1.7, xaxs='i', pt.cex=1.5, cex=0.9)
par(xaxt='s')
axis(1, at=seq(15,30,5), labels=c(0,seq(20,30,5)), cex.axis=0.7) 
axis.break(1, 17, style='slash')
text(x=15, y=seq(1.3,10.3,1), labels=rf_rxns$name, cex=0.79, pos=4)
dev.off()

# Main body
png(filename='~/Desktop/repos/Jenior_Cdifficile_2019/results/figures/Figure_4EF.png', 
    units='in', width=3.5, height=3.5, res=300)
layout(matrix(c(1,2), nrow=1, ncol=2, byrow=TRUE))
select_and_plot('2,3-dihydroxyisovalerate', title='Isovalerate', across=TRUE, best_ylim=2.25)
select_and_plot('N-acetylneuraminate', title='N-Acetylneuraminate', across=TRUE)
dev.off()


select_and_plot('imidazole_propionate', across=TRUE, best_ylim=1.5)


# Supplement
hi_glcnac <- as.vector(clinda_samples[,'EX_cpd00122_e'])
lo_glcnac <- as.vector(strep_samples[,'EX_cpd00122_e'])
lo_mannac <- as.vector(strep_samples[,'EX_cpd00492_e'])
lo_mannac <- subset(lo_mannac, lo_mannac >= 0)

pvals <- c()
for (x in c(1:1000)) {
    test_1 <- sample(hi_glcnac, size=10)
    test_2 <- sample(lo_glcnac, size=10)
    pvals[x] <- round(wilcox.test(test_1, test_2, exact=FALSE)$p.value, 3)
}
glcnac_pval <- median(pvals)

library(vioplot)
png(filename='~/Desktop/repos/Jenior_Cdifficile_2019/results/figures/Figure_S4AB.png', 
    units='in', width=3.5, height=3.5, res=300)
layout(matrix(c(1,2), nrow=1, ncol=2, byrow=TRUE))
par(mar=c(3, 2, 1, 1), mgp=c(1.4, 0.5, 0), xpd=FALSE, lwd=1.7, xaxt='n')

vioplot(hi_glcnac, lo_glcnac, col=c('firebrick2','blue4'), 
        ylim=c(0,1000), lwd=2, drawRect=FALSE, yaxt='n')
axis(side=2, at=seq(0,1000,200), cex.axis=0.5, lwd=2)
box(lwd=2)
par(xpd=TRUE)
text(x=c(1,2), y=-150, labels=c('High\nSpores','Low\nSpores'), cex=0.8)
text(x=1.5, y=1080, labels='N-Acetylglucosamine', cex=0.8, font=2)
text(x=-0.2, y=500, labels='Predicted Secretion Flux', cex=0.8, srt=90)
par(xpd=FALSE)

vioplot(0, lo_mannac, col=c('firebrick2','blue4'), 
        ylim=c(0,1000), lwd=2, drawRect=FALSE, yaxt='n')
axis(side=2, at=seq(0,1000,200), cex.axis=0.5, lwd=2)
box(lwd=2)
text(x=1, y=50, labels='inactive', cex=0.7)
par(xpd=TRUE)
text(x=c(1,2), y=-150, labels=c('High\nSpores','Low\nSpores'), cex=0.8)
text(x=1.5, y=1080, labels='N-Acetylmannosamine', cex=0.8, font=2)
text(x=-0.2, y=500, labels='Predicted Secretion Flux', cex=0.8, srt=90)
par(xpd=FALSE)

dev.off()


png(filename='~/Desktop/repos/Jenior_Cdifficile_2019/results/figures/Figure_S4C.png', 
    units='in', width=3, height=3, res=300)
select_and_plot('glucose', title='D-Glucose', best_ylim=6, across=TRUE)
dev.off()

