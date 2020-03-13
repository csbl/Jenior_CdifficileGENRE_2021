
# Set up environment
rm(list=ls())
gc()

# Load in functions
starting_dir <- getwd()
source('~/Desktop/repositories/Cdiff_modeling/src/init.R')

# Output plot name
plot_file <- 'results/figures/figure_3.pdf'

# Input raw data
# Metabolomes
metabolome <- 'data/metabolome/scaled_intensities.log10.tsv'
# Metadata
metadata <- 'data/metadata.tsv'

#----------------#

# Read in data

# Metabolomes
metabolome <- read.delim(metabolome, sep='\t', header=TRUE)
# Metadata
metadata <- read.delim(metadata, sep='\t', header=T, row.names=1)

#-------------------------------------------------------------------------------------------------------------------------#

# Format data

# Metadata
metadata$cage <- NULL
metadata$mouse <- NULL
metadata$gender <- NULL
metadata$susceptibility <- NULL
metadata$clearance <- NULL

# Metabolomes
metabolome$BIOCHEMICAL <- gsub('_', ' ', metabolome$BIOCHEMICAL)
rownames(metabolome) <- metabolome$BIOCHEMICAL
metabolome$BIOCHEMICAL <- NULL
metabolome$PUBCHEM <- NULL
metabolome$KEGG <- NULL
metabolome$SUB_PATHWAY <- NULL
metabolome$SUPER_PATHWAY <- NULL
metabolome <- as.data.frame(t(metabolome))

#-------------------------------------------------------------------------------------------------------------------------#

# Combine with metadata 
metabolome <- mergeByRow(metadata, metabolome)
rm(metadata)

# Subset to treatment groups
metabolome <- subset(metabolome, type == 'conventional')
metabolome$type <- NULL

cef <- subset(metabolome, abx == 'cefoperazone')
cef$abx <- NULL
strep <- subset(metabolome, abx == 'streptomycin')
strep$abx <- NULL
clinda <- subset(metabolome, abx == 'clindamycin')
clinda$abx <- NULL
noabx_mock <- subset(metabolome, abx == 'none')
noabx_mock$abx <- NULL
noabx_mock$infection <- NULL
rm(metabolome)

cef_630 <- subset(cef, infection == '630')
cef_630$infection <- NULL
cef_mock <- subset(cef, infection == 'mock')
cef_mock$infection <- NULL
rm(cef)
strep_630 <- subset(strep, infection == '630')
strep_630$infection <- NULL
strep_mock <- subset(strep, infection == 'mock')
strep_mock$infection <- NULL
rm(strep)
clinda_630 <- subset(clinda, infection == '630')
clinda_630$infection <- NULL
clinda_mock <- subset(clinda, infection == 'mock')
clinda_mock$infection <- NULL
rm(clinda)

#-------------------------------------------------------------------------------------------------------------------------#

# Subset to metabolites of interest
#shared_metabolites <- c("xanthine","glycerophosphoethanolamine","arabonate/xylonate",
#                        "fumarate","arabinose","mevalonate","malate","N-acetylneuraminate",
#                        "maltose","diacetylchitobiose","alpha-ketoglutarate",
#                        "2-keto-3-deoxy-gluconate","sucrose")

shared_metabolites <- c("succinate","arabonate/xylonate","N-acetylneuraminate","raffinose",
                        "mevalonate","fructose")


cef_strep_metabolites <- c("N-acetylglucosamine/N-acetylgalactosamine","gluconate","succinate")
cef_clinda_metabolites <- c("mannose")
strep_clinda_metabolites <- c("lactate","trimethylamine N-oxide","galactitol (dulcitol)")

cef_metabolites <- c("glucose","nicotinate","betaine")
strep_metabolites <- c("alanine")
clinda_metabolites <- c("fructose","raffinose","cytosine")

#-------------#

# Shared across all groups
cef_630_shared <- cef_630[,which(colnames(cef_630) %in% shared_metabolites)]
cef_mock_shared <- cef_mock[,which(colnames(cef_mock) %in% shared_metabolites)]
strep_630_shared <- strep_630[,which(colnames(strep_630) %in% shared_metabolites)]
strep_mock_shared <- strep_mock[,which(colnames(strep_mock) %in% shared_metabolites)]
clinda_630_shared <- clinda_630[,which(colnames(clinda_630) %in% shared_metabolites)]
clinda_mock_shared <- clinda_mock[,which(colnames(clinda_mock) %in% shared_metabolites)]
noabx_mock_shared <- noabx_mock[,which(colnames(noabx_mock) %in% shared_metabolites)]

# Most interesting combination
cef_strep_630 <- cef_630[,which(colnames(cef_630) %in% cef_strep_metabolites)]
cef_strep_mock <- cef_mock[,which(colnames(cef_mock) %in% cef_strep_metabolites)]
strep_cef_630 <- cef_630[,which(colnames(cef_630) %in% cef_strep_metabolites)]
strep_cef_mock <- strep_mock[,which(colnames(strep_mock) %in% cef_strep_metabolites)]
cef_strep_noabx <- noabx_mock[,which(colnames(noabx_mock) %in% cef_strep_metabolites)]

# Most interesting solo group
clinda_630_alone <- clinda_630[,which(colnames(clinda_630) %in% clinda_metabolites)]
clinda_mock_alone <- clinda_mock[,which(colnames(clinda_mock) %in% clinda_metabolites)]
clinda_mock_noabx <- noabx_mock[,which(colnames(noabx_mock) %in% clinda_metabolites)]

rm(cef_630, cef_mock, clinda_630, clinda_mock, strep_630, strep_mock, noabx_mock, 
   shared_metabolites, cef_strep_metabolites, cef_clinda_metabolites, strep_clinda_metabolites, 
   strep_metabolites, clinda_metabolites, cef_metabolites)

#-------------#



#-------------------------------------------------------------------------------------------------------------------------#

# Plotting
pdf(file=plot_file, width=14, height=7)
layout(matrix(c(1,2,3,
                4,5,6),
              nrow=2, ncol=3, byrow=TRUE))

for (index in 1:ncol(noabx_mock_shared)) {
  metabolitePlot(index, noabx_mock_shared, 
                 strep_mock_shared, strep_630_shared,
                 cef_mock_shared, cef_630_shared,
                 clinda_mock_shared, clinda_630_shared,
                 '')
}
dev.off()

# arabonate
# glycerophosphoethanolamine
# maltose
# mevalonate
# N-acetylneuraminate
# sucrose
# xanthine




#-------------------------------------------------------------------------------------------------------------------------#

# Clean up
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
setwd(starting_dir)
rm(list=ls())
gc()
