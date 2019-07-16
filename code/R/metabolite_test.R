
# Set up environment
rm(list=ls())
gc()

# Load in functions
source('~/Desktop/old_repos/Jenior_Metatranscriptomics_mSphere_2018/code/R/functions.R')

# Metabolomes
metabolome <- 'data/metabolome/scaled_intensities.log10.tsv'

# Metadata
metadata <- 'data/metadata.tsv'

# Metabolite of interest
metabolite <- 'glutamine'


# Test differences function
testDiff <- function(data){
  infected <- subset(data, infection == '630')
  infected <- as.vector(infected[,2])
  mock <- subset(data, infection == 'mock')
  mock <- as.vector(mock[,2])
  
  med_infected <- median(infected)
  med_mock <- median(mock)
  pval <- wilcox.test(infected, mock, exact=FALSE)$p.value
  
  print('Infected median:')
  print(med_infected)
  print('Mock median:')
  print(med_mock)
  print('p-value:')
  print(pval)
}

#----------------#

# Read in data

# Metabolomes
metabolome <- read.delim(metabolome, sep='\t', header=TRUE)

# Metadata
metadata <- read.delim(metadata, sep='\t', header=T, row.names=1)

#-------------------------------------------------------------------------------------------------------------------------#

# Format data

# Metadata
metadata$type <- NULL
metadata$cage <- NULL
metadata$mouse <- NULL
metadata$gender <- NULL
metadata$clearance <- c(rep('colonized',18), rep('cleared',18), rep('colonized',18),
                        rep('mock',12), rep('colonized',18))

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
untreated <- subset(metabolome, abx == 'none')[,metabolite]
metabolome <- subset(metabolome, abx != 'none')
metabolome <- metabolome[,c('abx','infection', metabolite)]
rm(metadata)

# Subset
cef <- subset(metabolome, abx == 'cefoperazone')
cef$abx <- NULL
clinda <- subset(metabolome, abx == 'clindamycin')
clinda$abx <- NULL
strep <- subset(metabolome, abx == 'streptomycin')
strep$abx <- NULL
gf <- subset(metabolome, abx == 'germfree')
gf$abx <- NULL
rm(metabolome)

# Test
testDiff(cef)
testDiff(clinda)
testDiff(strep)
testDiff(gf)
print('Untreated median:')
print(median(untreated))

