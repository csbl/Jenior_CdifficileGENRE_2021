
# Flux sampling files
clinda_samples <- '~/Desktop/repos/Jenior_Cdifficile_2019/data/contextualized_iCdG693/riptide_med_tox/flux_samples.tsv'
strep_samples <- '~/Desktop/repos/Jenior_Cdifficile_2019/data/contextualized_iCdG693/riptide_low_tox/flux_samples.tsv'

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
# clinda = 0.0333
# strep = 0.0425

# Unsupervised learning (NMDS)
flux_nmds <- as.data.frame(metaMDS(flux_dist, k=2, trymax=25)$points)

# Center points
flux_x <- (abs(max(flux_nmds$MDS1)) - abs(min(flux_nmds$MDS1))) / 2
flux_y <- (abs(max(flux_nmds$MDS2)) - abs(min(flux_nmds$MDS2))) / 2
flux_nmds$MDS1 <- flux_nmds$MDS1 - flux_x
flux_nmds$MDS2 <- flux_nmds$MDS2 - flux_y
flux_x <- max(abs(max(flux_nmds$MDS1)), abs(min(flux_nmds$MDS1))) + 0.01
flux_y <- max(abs(max(flux_nmds$MDS2)), abs(min(flux_nmds$MDS2))) + 0.01

# Subset axes
clinda_nmds_points <- subset(flux_nmds, rownames(flux_nmds) %in% clinda_names)
strep_nmds_points <- subset(flux_nmds, rownames(flux_nmds) %in% strep_names)
rm(clinda_names, strep_names)

# Statistical testing (permANOVA)
test <- merge(x=metadata, y=all_samples, by.x='label', by.y='row.names')
rownames(test) <- test$label
test$label <- NULL
test$treatment <- NULL
clearance_pval <- adonis(flux_dist ~ clearance, data=test, perm=99, method='bray')
clearance_pval <- clearance_pval$aov.tab[[6]][1]
clearance_pval <- as.character(round(clearance_pval, 3))
rm(all_samples, test, flux_dist, metadata)

#-------------------------------------#

# Flux sampling files
clinda_samples <- '~/Desktop/repos/Jenior_Cdifficile_2019/data/contextualized_iCdG693/riptide_med_tox/flux_samples.tsv'
strep_samples <- '~/Desktop/repos/Jenior_Cdifficile_2019/data/contextualized_iCdG693/riptide_low_tox/flux_samples.tsv'

# Read in data
clinda_samples <- read.delim(clinda_samples, sep='\t', header=TRUE)
clinda_samples$X <- NULL
strep_samples <- read.delim(strep_samples, sep='\t', header=TRUE)
strep_samples$X <- NULL
rxn_names <- read.delim(rxn_names, sep='\t', header=TRUE)

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
#rm(clinda_samples, strep_samples)

# Run AUCRF and obtain feature lists
#library(AUCRF)
#set.seed(906801)
#all_aucrf <- AUCRF(condition ~ ., data=all_samples, pdel=0, k0=10)
#print(all_aucrf)
rm(all_samples)

# Assemble feature table
#top_rxns_importance <- all_aucrf$ranking[1:all_aucrf$Kopt]
#rf_rxns <- as.data.frame(cbind(labels(top_rxns_importance), as.vector(top_rxns_importance)))
#colnames(rf_rxns) <- c('id','mda')
#rf_rxns$mda <- as.numeric(as.character(rf_rxns$mda))
#rm(all_aucrf, top_rxns_importance)

# Import pre-run data
rf_rxns <- read.delim('/home/mjenior/Desktop/repos/Jenior_Cdifficile_2019/data/aucrf_invivo.tsv', sep='\t', header=TRUE)
rf_rxns$name <- as.character(rf_rxns$name)
rf_rxns$name <- gsub('_', ' ', rf_rxns$name)
rf_rxns$id <- as.character(rf_rxns$id)
rf_rxns$mda <- round(as.numeric(as.character(rf_rxns$mda)), 2)
rf_rxns$label <- paste0(rf_rxns$name , ' (MDA: ', rf_rxns$mda,')')
rf_rxns <- rf_rxns[order(rf_rxns$mda),]

# Subset to most informative
top_rf_rxns <- subset(rf_rxns, mda >= 30)
top_rf_rxns <- top_rf_rxns[order(top_rf_rxns$mda),]

# Subset data for plotting
rf_clinda_samples <- clinda_samples[, top_rf_rxns$id]
rf_strep_samples <- strep_samples[, top_rf_rxns$id]

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

select_and_plot <- function(metabolite_name, best_ylim=0, correction=FALSE, title=NA, panel=FALSE, across=FALSE) {
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
    
    par(mar=c(3.3,3,1.5,0.5), xpd=FALSE, mgp=c(1.9,0.7,0), lwd=1.7, xaxt='n', las=1)
    boxplot(a_mock, at=1, xlim=c(0.6,3.4), ylim=c(0,best_ylim), col='firebrick2', cex.axis=0.8,
            xlab="", ylab="Scaled Intesity", outcex=0, whisklty=1, medlwd=2,  main=title, cex.main=0.8) 
    boxplot(a_630, at=1.5, xlim=c(0.6,3.4), ylim=c(0,best_ylim), col='firebrick2', 
            xlab="", ylab="Scaled Intesity", outcex=0, whisklty=1, medlwd=2, yaxt='n', add=TRUE) 
    if (across == FALSE) {abline(v=2)}
    boxplot(b_mock, at=2.5, xlim=c(0.6,3.4), ylim=c(0,best_ylim), col='blue3', 
            xlab="", ylab="Scaled Intesity", outcex=0, whisklty=1, medlwd=2, yaxt='n',add=TRUE) 
    boxplot(b_630, at=3, xlim=c(0.6,3.4), ylim=c(0,best_ylim), col='blue3', 
            xlab="", ylab="Scaled Intesity", outcex=0, whisklty=1, medlwd=2, yaxt='n',add=TRUE) 
    mtext(c('High\nSpores','Low\nSpores'), side=1, at=c(1.25,2.75), padj=1.5, cex=0.9)
    mtext(c('CDI:','-','+','-','+'), side=1, at=c(0.5,1,1.5,2.5,3), padj=0.5, cex=c(0.9,1.2,1.2,1.2,1.2))
    if (panel != FALSE) {mtext(text=panel, side=2, font=2, cex=1.1, padj=-8, adj=3.5)}
    
    a_pval <- round(wilcox.test(a_mock, a_630, exact=FALSE)$p.value,3)
    b_pval <- round(wilcox.test(b_mock, b_630, exact=FALSE)$p.value,3)
    ab_pval <- round(wilcox.test(a_mock, b_mock, exact=FALSE)$p.value,3)
    ba_pval <- round(wilcox.test(a_630, b_630, exact=FALSE)$p.value,3)
    if (across == TRUE) {fctr <- 0.7} else {fctr <- 0.96}
    if (a_pval <= 0.05) {
        segments(x0=1, y0=best_ylim*fctr, x1=1.5, y1=best_ylim*fctr, lwd=1.7)
        text(x=1.25, y=best_ylim, '*', font=2, cex=1.5)}
    if (b_pval <= 0.05) {
        segments(x0=2.5, y0=best_ylim*fctr, x1=3, y1=best_ylim*fctr, lwd=1.7)
        text(x=2.75, y=best_ylim*(fctr+0.04), '*', font=2, cex=1.5)}
    if (across == TRUE) {
        if (ab_pval <= 0.05) {
            segments(x0=1, y0=best_ylim*0.96, x1=2.5, y1=best_ylim*0.96, lwd=1.7)
            text(x=1.75, y=best_ylim, '*', font=2, cex=1.5)}
        if (ba_pval <= 0.05) {
            segments(x0=1.5, y0=best_ylim*0.86, x1=3, y1=best_ylim*0.86, lwd=1.7)
            text(x=2.25, y=best_ylim*0.9, '*', font=2, cex=1.5)}
    }
}

#---------------------------------------------------------------------#

# Generate figure panels

# Ordination of shared reaction flux samples
png(filename='~/Desktop/repos/Jenior_Cdifficile_2019/results/figures/figure_4A.png', 
    units='in', width=4.5, height=4.5, res=300)
library(scales)
par(mar=c(3.5,3.5,0.5,0.5), las=1, mgp=c(2.2,0.7,0), lwd=1.7)
plot(x=flux_nmds$MDS1, y=flux_nmds$MDS2, xlim=c(-0.045,0.045), ylim=c(-0.045, 0.045),
     xlab='NMDS Axis 1', ylab='NMDS Axis 2', pch=19, cex.lab=1.1, cex=0, cex.axis=0.9)
points(x=strep_nmds_points$MDS1, y=strep_nmds_points$MDS2, bg=alpha('blue3',0.8), pch=21, cex=1.7)
points(x=clinda_nmds_points$MDS1, y=clinda_nmds_points$MDS2, bg=alpha('firebrick2',0.8), pch=21, cex=1.7)
legend('topleft', legend=c('Low Sporulation','High Sporulation'), 
       pt.bg=c(alpha('blue3',0.8),alpha('firebrick2',0.8)), pch=21, pt.cex=1.7, box.lwd=2)
legend('bottomleft', legend=substitute(paste(italic('p'), '-value = 0.01*')), pt.cex=0, bty='n')
legend('bottomright', legend='Shared Core Metabolism', pt.cex=0, bty='n')
box(lwd=2)
dev.off()

# Flux samples from AUCRF features
library(plotrix)
png(filename='~/Desktop/repos/Jenior_Cdifficile_2019/results/figures/figure_4B.png', 
    units='in', width=2.5, height=4, res=300)
par(mar=c(2.5, 0.5, 0.5, 0.5), mgp=c(1.4, 0.5, 0), xpd=FALSE, lwd=1.7, xaxt='n')
dotchart(rf_rxns$mda, color='darkorchid4', xlab='Mean Decrease Accuracy (%)', xlim=c(10,40),  
         pch=16, lwd=1.7, xaxs='i', pt.cex=1.5, cex=0.9)
par(xaxt='s')
axis(1, at=seq(10,40,5), labels=c(0,seq(15,40,5)), cex.axis=0.7) 
axis.break(1, 12.5, style='slash')
text(x=9, y=seq(1.3,10.3,1), labels=rf_rxns$name, cex=0.75, pos=4)
dev.off()

# Metabolite concentrations
png(filename='~/Desktop/repos/Jenior_Cdifficile_2019/results/figures/figure_4D.png', 
    units='in', width=4, height=3, res=300)
layout(matrix(c(1,2), nrow=1, ncol=2, byrow=TRUE))
select_and_plot('N-acetylglucosamine/N-acetylgalactosamine', 0.3, correction=TRUE, title='N-Acetylglucosamine', panel='D')
select_and_plot('p-cresol_sulfate', 3.0, correction=TRUE, title='p-Cresol', panel='E')
dev.off()

png(filename='~/Desktop/repos/Jenior_Cdifficile_2019/results/figures/Figure_S4.png', 
    units='in', width=4, height=3, res=300)
layout(matrix(c(1,2), nrow=1, ncol=2, byrow=TRUE))
select_and_plot('methionine', title='Methionine', panel='A', across=TRUE)
select_and_plot('glutamate', title='Glutamate', panel='B', across=TRUE)
dev.off()

