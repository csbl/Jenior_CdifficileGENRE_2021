
# Set up environment
rm(list=ls())
gc()

# Load in functions
starting_dir <- getwd()
source('~/Desktop/repos/Jenior_Metatranscriptomics_mSphere_2018/code/R/functions.R')
library(packcircles)
library(ggplot2)
library(viridis)
library(gridExtra)

# Output plot name
plot_16S <- '~/Desktop/plot_16S.png'
plot_metabolome <- '~/Desktop/plot_metabolome.png'
plot_metagenome <- '~/Desktop/plot_metagenome.png'

# Metabolomes
metabolome <- 'data/metabolome/scaled_intensities.log10.tsv'

# 16S
shared_otu <- 'data/16S_analysis/all_treatments.0.03.unique_list.conventional.shared'
otu_tax <- 'data/16S_analysis/formatted.all_treatments.0.03.cons.taxonomy'

# Metadata
metadata <- 'data/metadata.tsv'

# Colors
strep_col <- 'red3'
cef_col <- 'royalblue3'
clinda_col <- 'springgreen4'
noabx_col <- 'gray20'

#----------------#

# Read in data

# Metabolomes
metabolome <- read.delim(metabolome, sep='\t', header=TRUE)

# 16S
shared_otu <- read.delim(shared_otu, sep='\t', header=T, row.names=2)
shared_otu <- shared_otu[!rownames(shared_otu) %in% c('CefC5M2'), ]  # Remove possible contaminated sample
shared_otu <- shared_otu[,!(names(shared_otu) %in% c('Otu0004','Otu0308'))] # Remove residual C. difficile OTUs
shared_otu$numOtus <- NULL
shared_otu$label <- NULL
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
metabolome <- clean_merge(metadata, metabolome)
metabolome <- subset(metabolome, abx != 'germfree')
metabolome <- subset(metabolome, infection == 'mock')
metabolome$abx <- NULL
metabolome$infection <- NULL

# 16S
otu_tax$genus <- gsub('_', ' ', otu_tax$genus)
otu_tax$genus <- gsub('Ruminococcus2', 'Ruminococcus', otu_tax$genus)
otu_tax$taxon <- paste(otu_tax$genus, otu_tax$OTU.1, sep='_')
otu_tax$phylum <- NULL
otu_tax$genus <- NULL
otu_tax$OTU.1 <- NULL
shared_otu <- t(shared_otu)
shared_otu <- clean_merge(otu_tax, shared_otu)
rownames(shared_otu) <- shared_otu$taxon
shared_otu$taxon <- NULL
shared_otu <- t(shared_otu)
shared_otu <- clean_merge(metadata, shared_otu)
shared_otu <- subset(shared_otu, abx != 'germfree')
shared_otu <- subset(shared_otu, infection == 'mock')
shared_otu$infection <- NULL
rm(otu_tax)

# Metagenomes
cefoperazone <- c(1350,1055,620,1354,923,4266)
streptomycin <- c(365,759,1243,1663,2533,274,811,278,371,295,1485,1959,1660,1890)
clindamycin <- c(1889,699,1533,1917)
no_abx <- c(2054,2311,233,1331,1647,2233,1255,2018,8815,256,1500,1846,2304,957,395,1882,257,786,2331,488,2731,465,994,398,1092,1400,255,1026,1677,3084,1851,1773,1202,393,379,3230,1580,877,1196,995,3407,905,580,490,494,607,213,2219,315,8323,1023,585,1525,284,3175,727,1018,316,5214,490,959,761,1278,845,1101,826,1466,520,359,196,3115,9663,1658,2367,1209,1503,2193,1074,1650,1208,2923,261,206,3440,3583,2471,233)
max(cefoperazone)
max(streptomycin)
max(clindamycin)
max(no_abx)
max_genes <- 9663
cefoperazone <- c(max_genes, sort(cefoperazone, decreasing=TRUE))
streptomycin <- c(max_genes, sort(streptomycin, decreasing=TRUE))
clindamycin <- c(max_genes, sort(clindamycin, decreasing=TRUE))
no_abx <- sort(no_abx, decreasing=TRUE)

#-------------------------------------------------------------------------------------------------------------------------#

# Generate color palettes
family_colors <- c('chartreuse3', # Coriobacteriaceae
                   'mediumblue','royalblue3','dodgerblue2','deepskyblue','powderblue', # Bacteroidetes
                   'darkred','firebrick3','red2','brown3','tomato2','coral3','salmon', # Firmicutes
                   'darkgoldenrod1','#CCCC00', # Proteobacteria
                   'darkmagenta') # Verrucomicrobiaceae

clinda_colors <- c('chartreuse', 'darkgoldenrod1','coral3','deepskyblue','mediumblue')
cef_colors <- c('chartreuse', 'coral3','mediumblue','salmon','brown3','deepskyblue','tomato2')
strep_colors <- c('chartreuse', 'mediumblue','dodgerblue2','brown3','tomato2','coral3',
                  'darkgoldenrod1','darkmagenta','chartreuse3','salmon',
                  'firebrick3','#CCCC00','darkred','white','white')
noabx_colors <- c('dodgerblue2','firebrick3','powderblue','mediumblue','salmon',
                  'coral3','brown3','red2','tomato2',sample(family_colors,68,replace=T),
                  rep('white',10))

#-------------------------------------------------------------------------------------------------------------------------#

# Ordination analysis

# 16S - 810 OTUs
otu_dist <- vegdist(shared_otu[,2:ncol(shared_otu)], method='bray') # Bray-Curtis
otu_nmds <- as.data.frame(metaMDS(otu_dist, k=2, trymax=100)$points)
otu_nmds$MDS1 <- otu_nmds$MDS1 * -1
otu_nmds$MDS1 <- otu_nmds$MDS1 + 0.03
otu_nmds$MDS2 <- otu_nmds$MDS2 * -1
otu_nmds$MDS2 <- otu_nmds$MDS2 + 0.1
# Subset NMDS axes to color points
rownames(otu_nmds) <- rownames(shared_otu)
otu_nmds <- clean_merge(metadata, otu_nmds)
otu_nmds_susceptible <- subset(otu_nmds, abx != 'none')
otu_nmds_strep <- subset(otu_nmds, abx == 'streptomycin')
otu_nmds_cef <- subset(otu_nmds, abx == 'cefoperazone')
otu_nmds_clinda <- subset(otu_nmds, abx == 'clindamycin')
otu_nmds_resistant <- subset(otu_nmds, abx == 'none')
otu_res_centroids <- aggregate(cbind(otu_nmds_resistant$MDS1,otu_nmds_resistant$MDS2)~otu_nmds_resistant$infection, data=otu_nmds_resistant, mean)
otu_sus_centroids <- aggregate(cbind(otu_nmds_susceptible$MDS1,otu_nmds_susceptible$MDS2)~otu_nmds_susceptible$infection, data=otu_nmds_susceptible, mean)
rm(otu_dist, metadata)

#-------------------------------------------------------------------------------------------------------------------------#

# Feature selection

# Metabolome
# AUCRF feature selection/reduction (to 0% OOB)
colnames(metabolome) <- make.names(colnames(metabolome))
metabolome_aucrf <- aucrfSusceptibility(metabolome)
# Get OOB
metabolome_aucrf_oob <- metabolome_aucrf$RFopt
metabolome_aucrf_oob <- metabolome_aucrf_oob$err.rate
metabolome_aucrf_oob <- as.character(round(median(metabolome_aucrf_oob[,1]) * 100, 3))
# Get features
metabolome_aucrf <- as.character(OptimalSet(metabolome_aucrf)$Name)
res_metabolome <- subset(metabolome, susceptibility == 'resistant')[, metabolome_aucrf]
res_metabolome$susceptibility <- NULL
sus_metabolome <- subset(metabolome, susceptibility == 'susceptible')[, metabolome_aucrf]
sus_metabolome$susceptibility <- NULL
rm(metabolome, metabolome_aucrf)
# Find significant differences
metabolome_pval <- c()
for (i in 1:ncol(res_metabolome)){metabolome_pval[i] <- wilcox.test(res_metabolome[,i], sus_metabolome[,i], exact=FALSE)$p.value}
metabolome_pval <- round(p.adjust(metabolome_pval, method='BH'), 4)

# Reformat names to be more human readable
colnames(res_metabolome) <- c('Nudifloramide','N-Acetylproline','Sebacate / Decanedioate','Hyodeoxycholate','Murideoxycholate')
colnames(sus_metabolome) <- c('Nudifloramide','N-Acetylproline','Sebacate / Decanedioate','Hyodeoxycholate','Murideoxycholate')

#-------------------------------------------------------------------------------------------------------------------------#

# Function to process data and generate bubble plots
bubble_plot <- function(sizes, current=1, spacing=1, scale_text=FALSE, pal) {
  
  end <- length(sizes) + current - 1
  print(end+1)
  data <- data.frame(bin=paste('OGU\n', current:end), genes=sizes)
  
  packing <- circleProgressiveLayout(sizes, sizetype='area')
  packing$radius <- spacing * packing$radius
  data <- cbind(data, packing)
  data.gg <- circleLayoutVertices(packing, npoints=50)
  
  if (scale_text == FALSE) {
    txt_size=1
  } else {txt_size=sizes}
  
  bubbles <- ggplot() + 
    geom_polygon(data=data.gg, aes(x, y, group=id, fill=as.factor(id)), color='black') +
    scale_fill_manual(values=pal) +
    geom_text(data=data, aes(x, y, size=txt_size, label=bin), color='black') +
    theme_void() +
    theme(legend.position='none') +
    coord_equal()
  
  return(bubbles)
}

# Create figures
no_abx_panel <- bubble_plot(no_abx, scale_text=TRUE, pal=noabx_colors)
strep_panel <- bubble_plot(streptomycin, current=88, spacing=0.95, scale_text=TRUE, pal=strep_colors)
cef_panel <- bubble_plot(cefoperazone, current=103, spacing=0.95, pal=cef_colors)
clinda_panel <- bubble_plot(clindamycin, current=110, spacing=0.95, pal=clinda_colors)

#-------------------------------------------------------------------------------------------------------------------------#

# Generate figures

# 16S
png(filename=plot_16S, width=5, height=5, units='in', res=300)
par(mar=c(4,4,1,1), las=1, mgp=c(2.8,0.75,0))
plot(x=otu_nmds$MDS1, y=otu_nmds$MDS2, xlim=c(-0.4,0.4), ylim=c(-0.5,0.5),
     xlab='NMDS axis 1', ylab='NMDS axis 2', pch=19, cex.axis=1.2, cex.lab=1.2)
segments(x0=otu_nmds_susceptible$MDS1, y0=otu_nmds_susceptible$MDS2, x1=otu_sus_centroids[1,2], y1=otu_sus_centroids[1,3], col='gray30')
points(x=otu_nmds_strep$MDS1, y=otu_nmds_strep$MDS2, col=strep_col, pch=16, cex=2, lwd=1.2)
points(x=otu_nmds_cef$MDS1, y=otu_nmds_cef$MDS2, col=cef_col, pch=16, cex=2, lwd=1.2)
points(x=otu_nmds_clinda$MDS1, y=otu_nmds_clinda$MDS2, col=clinda_col, pch=16, cex=2, lwd=1.2)
segments(x0=otu_nmds_resistant$MDS1, y0=otu_nmds_resistant$MDS2, x1=otu_res_centroids[1,2], y1=otu_res_centroids[1,3], col='gray30')
points(x=otu_nmds_resistant$MDS1, y=otu_nmds_resistant$MDS2, col=noabx_col, pch=16, cex=2, lwd=1.2)
legend('bottomleft', legend=c('Resistant vs Susceptible', as.expression(bquote(paste(italic('p'),' << 0.001')))), 
       pch=1, pt.cex=0, bty='n')
legend('bottomright', legend=c('Untreated','Streptomycin','Cefoperazone','Clindamycin'), 
       col=c(noabx_col, strep_col, cef_col, clinda_col), pch=16, pt.cex=2)
legend('topleft', legend='16S rRNA gene region sequencing', bty='n', pt.cex=0)
box(lwd=2)
dev.off()

#---------------#

# Metabolomes
multiStripchart_png(plot_metabolome, res_metabolome, sus_metabolome, metabolome_pval, metabolome_aucrf_oob, 
                'Untreated', 'All Antibiotic Groups', noabx_col, 'maroon4', '', 'black', 
                as.list(colnames(res_metabolome)), expression(paste('Scaled Intensity (',log[10],')')))

#---------------#

# Metagenomes
png(filename=plot_metagenome, width=10, height=10, units='in', res=300)
grid.arrange(no_abx_panel, strep_panel, cef_panel, clinda_panel, nrow=2)
dev.off()

#-------------------------------------------------------------------------------------------------------------------------#

# Clean up
#for (dep in deps){
#  pkg <- paste('package:', dep, sep='')
#  detach(pkg, character.only = TRUE)
#}
#setwd(starting_dir)
#rm(list=ls())
#gc()
