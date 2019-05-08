
# Set up environment
rm(list=ls())
gc()

# Load dependencies
deps <- c('vegan', 'igraph', 'ggplot2', 'shape');
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  }
  library(dep, verbose=FALSE, character.only=TRUE)
}
set.seed(42)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Define network object
network <- '~/Desktop/repositories/Cdiff_modeling/cdf_met_pairs.tsv'
network <- read.delim(network, header=FALSE, sep='\t')

# Define energy and small molecule metabolites to omit
omit <- c('cpd00001_c','cpd00002_c','cpd00003_c','cpd00004_c','cpd00005_c','cpd00006_c',
          'cpd00008_c','cpd00009_c','cpd00038_c','cpd00031_c','cpd00067_c','cpd00012_c',
          'cpd00018_c','cpd00052_c','cpd00096_c','cpd00357_c','cpd00297_c','cpd00011_c',
          'cpd00010_c','cpd00020_c')

# Format directed graph
network <- graph.data.frame(network, directed=TRUE)
network <- simplify(network) # Simplify graph
network <- decompose.graph(network) # Decompose graph

# Separate the largest component from small islands
large <- which.max(sapply(network, vcount))
largest_component <- network[[large]]
network[[large]] <- NULL # Leave only small components
rm(large)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Analysis

# Large Component
# Find degree centrality
graph_degree <- as.data.frame(degree(largest_component, v=V(largest_component), mode='all'))
# Calculate betweenness centrality
graph_betweenness <- as.data.frame(betweenness(largest_component))
# Calculate closeness centrality
graph_closeness_in <- as.data.frame(closeness(largest_component, vids=V(largest_component), mode='in'))
graph_closeness_out <- as.data.frame(closeness(largest_component, vids=V(largest_component), mode='out'))
graph_closeness_total <- as.data.frame(closeness(largest_component, vids=V(largest_component), mode='total'))
# Merge characteristic tables
graph_topology <- merge(graph_degree, graph_betweenness, by='row.names')
rownames(graph_topology) <- graph_topology$Row.names
graph_topology$Row.names <- NULL
graph_topology <- merge(graph_topology, graph_closeness_total, by='row.names')
rownames(graph_topology) <- graph_topology$Row.names
graph_topology$Row.names <- NULL
colnames(graph_topology) <- c('degree_centrality','betweenness_centrality','closeness_centrality')
rm(graph_degree, graph_betweenness, graph_closeness_in, graph_closeness_out, graph_closeness_total)
# Remove uninformative metabolites
graph_topology <- subset(graph_topology, !rownames(graph_topology) %in% omit)
graph_topology <- subset(graph_topology, graph_topology$betweenness_centrality > 0)

# Correlate centralities
r_value <- cor(x=graph_topology$degree_centrality, y=graph_topology$betweenness_centrality, method='spearman')
p_value <- cor.test(x=graph_topology$degree_centrality, y=graph_topology$betweenness_centrality, method='spearman', exact=FALSE)$p.value

plot(x=graph_topology$degree_centrality, y=graph_topology$betweenness_centrality)
abline(lm(graph_topology$betweenness_centrality ~ graph_topology$degree_centrality), col="red")






#------------#

# Examine sub-components
for (x in 1:length(network)) {
  print(V(network[[x]])$name)
}





#-------------------------------------------------------------------------------------------------------------------------------------#





# Write tables to files, ranked by betweenness
table_file <- '~/Desktop/'
write.table(enzyme_topology, file=table_file, quote=FALSE, sep='\t', row.names=FALSE)
table_file <- '~/Desktop/'
write.table(substrate_topology, file=table_file, quote=FALSE, sep='\t', row.names=FALSE)
rm(table_file, graph_topology, enzyme_topology, substrate_topology)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Transform largest component graph for plotting

# Get nodes for visualization
ko_simp <- as.vector(grep('K', V(largest_whole_graph)$name, value=TRUE))
substrate_simp <- as.vector(grep('C', V(largest_whole_graph)$name, value=TRUE))
nodes <- c(ko_simp, substrate_simp)
summary(largest_whole_graph)
print(length(ko_simp))
print(length(substrate_simp))

# Remove zeroes so transformation doesn't return negative infinity
ko[,2][ko[,2] == 0] <- 1
ko[,2] <- log10(ko[,2])
ko[,2][ko[,2] < 0.4] <- 0.4

# Scale points by number of transcripts mapped
ko <- as.data.frame(subset(ko, ko[,1] %in% ko_simp))
ko <- ko[match(ko_simp, ko$KO_code),]
mappings <- c(as.vector(ko[,2] * 5), rep(2, length(substrate_simp)))
node_size <- cbind.data.frame(nodes, mappings)
node_size <- node_size[match(V(largest_whole_graph)$name, node_size$nodes),]
node_size <- setNames(as.numeric(node_size[,2]), as.character(node_size[,1]))
V(largest_whole_graph)$size <- as.matrix(node_size)
rm(raw_graph, ko, ko_simp, substrate_simp, node_size, mappings, nodes, simple_graph)

# Color graph
V(largest_whole_graph)$color <- ifelse(grepl('K', V(largest_whole_graph)$name), adjustcolor('firebrick3', alpha.f=0.6), 'blue3') # Color nodes
E(largest_whole_graph)$color <- 'gray15' # Color edges