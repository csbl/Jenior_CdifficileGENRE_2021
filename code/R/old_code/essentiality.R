
library(gplots)

# Load in data
essential <- as.data.frame(t(read.delim(file='~/Desktop/repos/Cdiff_modeling/data/essentiality.tsv', sep='\t', header=TRUE, row.names=1)))

# Assess differential groups of essentiality
test <- as.data.frame(t(essential))
# Group 1
m9_anaerobic_only <- subset(test, m9_anaerobic == 2)
m9_anaerobic_only <- subset(m9_anaerobic_only, lb_aerobic == 0)
m9_anaerobic_only <- subset(m9_anaerobic_only, m9_aerobic == 0)
m9_anaerobic_only <- subset(m9_anaerobic_only, pfba == 0)
m9_anaerobic_only <- rownames(m9_anaerobic_only)
print(m9_anaerobic_only)
rm(m9_anaerobic_only)
# Group 2
lb_aerobic_only <- subset(test, lb_aerobic == 2)
lb_aerobic_only <- subset(lb_aerobic_only, m9_aerobic == 0)
lb_aerobic_only <- subset(lb_aerobic_only, m9_anaerobic == 0)
lb_aerobic_only <- subset(lb_aerobic_only, pfba == 0)
lb_aerobic_only <- rownames(lb_aerobic_only)
print(lb_aerobic_only)
rm(lb_aerobic_only)
# Group 3
m9_aerobic_only <- subset(test, m9_aerobic == 2)
m9_aerobic_only <- subset(m9_aerobic_only, lb_aerobic == 0)
m9_aerobic_only <- subset(m9_aerobic_only, m9_anaerobic == 0)
m9_aerobic_only <- subset(m9_aerobic_only, pfba == 0)
m9_aerobic_only <- rownames(m9_aerobic_only)
print(m9_aerobic_only)
rm(m9_aerobic_only)
# Group 4
m9_anaerobic_non <- subset(test, m9_anaerobic == 0)
m9_anaerobic_non <- subset(m9_anaerobic_non, lb_aerobic == 2)
m9_anaerobic_non <- subset(m9_anaerobic_non, m9_aerobic == 2)
m9_anaerobic_non <- subset(m9_anaerobic_non, pfba == 2)
m9_anaerobic_non <- rownames(m9_anaerobic_non)
print(m9_anaerobic_non)
rm(m9_anaerobic_non)
rm(test)

# Format data for plotting
genres <- rownames(essential)
essential <- as.matrix(sapply(essential, as.numeric))
rownames(essential) <- genres
rm(genres)

# Generate figure
pdf(file='~/Desktop/repos/Cdiff_modeling/results/essentiality.pdf', width=6, height=4)
heatmap.2(essential, col=c('black','white','forestgreen'), dendrogram='none', density.info='none', 
          trace='none', key=FALSE, margins=c(1,1), Rowv=FALSE, sepwidth=c(0.05,0.01), 
          sepcolor='black', colsep=1:ncol(essential), rowsep=1:nrow(essential), labCol=FALSE)
dev.off()

# Clean up
rm(list=ls())
gc()
