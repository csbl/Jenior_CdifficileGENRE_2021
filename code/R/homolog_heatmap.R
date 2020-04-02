
# Read in data
genes <- read.delim('~/Desktop/repos/Jenior_Cdifficile_2019/data/all_gene_alignment.tsv', header=TRUE)

# Subset to genes included in at least one GENRE
genes$genre <- genes$iCdR758_included + genes$iCdG791_included
genes <- subset(genes, genre > 0)
genes <- subset(genes, genre < 2)
genes$genre <- NULL

# Remove unnecessary columns
genes$cdR20291_homolog <- NULL
genes$cd630_homolog <- NULL
genes$cdR20291_gene <- NULL
genes$cd630_gene <- NULL

# Reformat for heatmap
rownames(genes) <- make.unique(as.character(genes$gene_name))
genes$gene_name <- NULL
genes <- as.matrix(genes)

# Generate figure
heatmap(genes, col=c('orange','purple'))


