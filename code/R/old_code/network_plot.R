
# Load dependencies
deps <- c('vegan', 'igraph', 'ggplot2', 'shape');
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  }
  library(dep, verbose=FALSE, character.only=TRUE)
}
set.seed(42)

# Colors specific nodes in network
color_points <- function(graph, reactions, metabolites, new_color){
  for (x in reactions) {
    if (x %in% V(graph)$name) {
      V(graph)[x]$color <- adjustcolor(new_color, alpha.f=0.8)
    }
  }
  for (y in metabolites) {
    if (y %in% V(graph)$name) {
      V(graph)[y]$color <- new_color
    }
  }
  return(graph)
}

#-------------------------------------------------------------------------------------------------------------------------------------#

# Define file variables for network plot
network <- '~/Desktop/repositories/Cdiff_modeling/data/icdj794_net_rxn.tsv'

# Read in metabolic network data
network <- read.delim(network, header=FALSE, sep='\t')

# Format directed graph
network <- graph.data.frame(network, directed=TRUE)

# Simplify graph
network <- simplify(network)

# Get largest component
#network <- decompose.graph(network)[[1]]

# Prune "core" metabolites - >25 edges
core <- c('cpd00067_c','cpd00001_c','cpd00002_c','cpd00003_c','cpd00004_c','cpd00009_c',
          'cpd00008_c','cpd00010_c','cpd00012_c','cpd00011_c','cpd00020_c','cpd00013_c',
          'cpd00006_c','cpd00005_c','cpd00080_c','cpd00022_c','cpd00023_c','cpd11493_c',
          'cpd00014_c','cpd00046_c','cpd00026_c','cpd00018_c','cpd00067_e','cpd00061_c',
          'cpd00052_c','cpd00072_c','cpd02090_c','cpd00041_c','cpd11492_c','cpd00908_c',
          'cpd00102_c')
network <- delete_vertices(network, core)

#------------------------------------------------------------------#

# Prep for plotting

# Color and resize graph
V(network)$color <- ifelse(grepl('rxn', V(network)$name), adjustcolor('white', alpha.f=0.65), adjustcolor('gray25', alpha.f=0.85))
E(network)$color <- 'gray40' # Color edges
V(network)$size <- ifelse(grepl('rxn', V(network)$name), 8, 4)

# Set subsystem colors
# Anaerobic metabolism & auxotrophies - darkred
anaer_rxn <- c('rxn00799_c','rxn08168_c','rxn30075_c','rxn27328_c','EX_cpd00107_e','rxn05243_c','EX_cpd00069_e','rxn18579_c','rxn18584_c','rxn10144_c','rxn10883_c','rxn34357_c','EX_cpd00041_e','rxn08686_c','','rxn09031_c','rxn05467_c','rxn08186_c','rxn08952_c','rxn10161_c','rxn05175_c','rxn08304_c','rxn05175_c','rxn30210_c','EX_cpd00084_e','rxn05237_c','R03018_c','EX_cpd00644_e','EX_cpd00263_e','COG3601_c','rxn05255_c','rxn05466_c','rxn10826_c','EX_cpd00133_e','rxn11334_c','EX_cpd00443_e','rxn10421_c')
anaer_cpd <- c('cpd00107_e','cpd00069_e','cpd00066_e','cpd00132_e','cpd00054_e','cpd00041_e','cpd00084_e','cpd00644_e','cpd00263_e','cpd00013_e','cpd00242_e','cpd00133_e','cpd00443_e')
network <- color_points(network, anaer_rxn, anaer_cpd, 'darkred')

# Stickland fermentation - royalblue3
stick_rxn <- c('rxn05165_c','ENOG4108HXH_c','EX_cpd29317_e','rxn05165_2_c','rxn10862_c','rxn32054_c','rxn05164_c','EX_cpd00119_e','rxn05169_c','rxn05155_c','rxn05151_c','rxn14152_c','rxn00804_c','R01174_1_c','ID001_c','EX_cpd05178_e','R01214_c','rxn10618_c','R01174_2_c','ID002_c','EX_cpd01711_e','rxn01575_c','R01174_3_c','ID003_c','rxn10563_c','R01174_4_c','ID004_c','EX_cpd00489_e','rxn05061_c','ID005_c','EX_cpd01042_e','rxn10562_c','R01174_5_c','rxn11336_c','EX_cpd00430_e','rxn20606_c','rxn12566_c','rxn11676_c','rxn13022_c','K20025_c','R11462_c','R11076_c','ID006_c','R06782_c','rxn01610_c','ID007_c','rxn07124_c','ID008_c','EX_cpd02501_e')
stick_cpd <- c('cpd00129_e','cpd00567_e','cpd29317_e','cpd00033_e','cpd00119_e','cpd00119_c','cpd00161_e','cpd00053_e','cpd05178_c','cpd05178_e','cpd00156_c','cpd00123_c','cpd00481_c','cpd01711_e','cpd00322_c','cpd00760_c','cpd19585_e','cpd03165_c','cpd00489_c','cpd00489_e','cpd01042_c','cpd01042_e','cpd00452_c','cpd00430_c','cpd00430_e','cpd00339_e','cpd00200_c','C21400_c','C21399_c','C21090_c','C21399_e','cpd07326_c','cpd03343_c','cpd03343_e','cpd03170_c','cpd03170_e')
network <- color_points(network, stick_rxn, stick_cpd, 'royalblue3')

# Carbohydrate acquisition - forestgreen
carb_rxn <- c('EX_cpd05161_e','rxn12619_c','rxn18524_c','R01103_c','ID010_c','rxn00629_c','rxn00634_c')
carb_cpd <- c('cpd00082_c','cpd00082_e','cpd00027_e','cpd00027_c','cpd00314_c','cpd00314_e','cpd00588_e','cpd00588_c','cpd00105_e','cpd00105_c','cpd01030_e','cpd01030_c','cpd00794_e','cpd00794_c','cpd05161_e','cpd05161_c','cpd20885_c','cpd00076_c','cpd00314_c','cpd00588_c')
network <- color_points(network, carb_rxn, carb_cpd, 'forestgreen')

# Short-chain fatty acid production - red
scfa_rxn <- c('EX_cpd00036_e','rxn09271_c','rxn25164_c','rxn16140_c','rxn02167_c','rxn16586_c','rxn13713_c','rxn00871_c','ID009_c','EX_cpd00029_e','rxn08061_c','rxn30338_c','rxn05602_c','rxn10171_c','rxn00602_c','rxn00669_c','rxn09172_c','EX_cpd00141_e')
scfa_cpd <- c('cpd00036_e','cpd00199_c','cpd07946_c','cpd00120_c','cpd00728_c','cpd00211_e','cpd00029_e','cpd00159_e','cpd00221_e','cpd00141_e')
network <- color_points(network, scfa_rxn, scfa_cpd, 'red')

# Host-derived glycan catabolism - purple1
nac_rxn <- c('EX_cpd00232_e','rxn29571_c','EX_cpd00492_e','rxn05486_c','rxn01484_c','rxn02262_c','rxn01620_c','rxn04937_c')
nac_cpd <- c('cpd00232_e','cpd00492_e','cpd00751_e','cpd00334_c','cpd00232_c','cpd00492_c','cpd00293_c','cpd00288_c','cpd01912_c','cpd00428_c')
network <- color_points(network, nac_rxn, nac_cpd, 'purple1')

# Original Biomass objective - goldenrod3
bio_rxn <- c('rxn03141_c','dna_rxn','rna_rxn','protein_rxn','teichoicacid_rxn','peptidoglycan_rxn','cellwall_rxn','lipid_rxn','cofactor_rxn','biomass_rxn')
bio_cpd <- c('cpd17042_c','cpd17043_c','cpd17041_c','cpd12894_c','cpd16661_c','cellwall_c','cpd11852_c','cofactor_c','cpd11416_c')
network <- color_points(network, bio_rxn, bio_cpd, 'goldenrod3')

# Gapfilled
gap_rxn <- c('rxn00010_c','rxn00285_c','rxn00341_c','rxn00530_c','rxn00553_c','rxn00579_c','rxn00605_c',
             'rxn00620_c','rxn00777_c','rxn00958_c','rxn00970_c','rxn01014_c','rxn01105_c','rxn01182_c',
             'rxn01256_c','rxn01359_c','rxn01520_c','rxn01538_c','rxn01539_c','rxn01644_c','rxn01669_c',
             'rxn01953_c','rxn02011_c','rxn02269_c','rxn02286_c','rxn02476_c','rxn02928_c','rxn03164_c',
             'rxn03167_c','rxn03852_c','rxn05039_c','rxn05227_c','rxn05291_c','rxn05293_c','rxn05744_c',
             'rxn07466_c','rxn08040_c','rxn08131_c','rxn09063_c','rxn09065_c','rxn09119_c','rxn09124_c',
             'rxn09126_c','rxn09142_c','rxn09145_c','rxn09147_c','rxn09163_c','rxn09341_c','rxn09348_c',
             'rxn09429_c','rxn10030_c','rxn10036_c','rxn10053_c','rxn10175_c','rxn10484_c','rxn10663_c',
             'rxn10840_c','rxn12239_c','rxn12676_c','rxn13012_c','rxn13139_c','rxn13140_c','rxn13251_c','rxn13914_c')
gap_cpd <- c('')
network <- color_points(network, gap_rxn, gap_cpd, 'blue')

# Calculate layout
layout4 <- layout_on_sphere(network)

#------------------------------------------------------------------#

# Set up plotting environment
pdf(file='~/Desktop/repositories/Cdiff_modeling/results/icdj794.pdf', width=6, height=6)
par(mar=c(0,0,0,0), xpd=TRUE)
plot(network, vertex.label=NA, layout=layout4, edge.arrow.size=0.4, edge.arrow.width=0.6)
dev.off()

#------------------------------------------------------------------#

# Clean up
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
rm(list=ls())
gc()

