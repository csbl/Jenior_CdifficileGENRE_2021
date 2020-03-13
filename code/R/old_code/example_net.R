#install.packages('igraph')
library(igraph)

# Define file variables for network plot
network <- '~/Desktop/iJO1366_net.tsv'

# Read in metabolic network data
network <- read.delim(network, header=FALSE, sep='\t')

# Format directed graph
raw_graph <- graph.data.frame(network, directed=TRUE)
rm(network)

# Simplify graph
simple_graph <- simplify(raw_graph)
rm(raw_graph)

# Decompose graph
decomp_whole_graph <- decompose.graph(simple_graph)
rm(simple_graph)

# Get largest component
largest_component <- which.max(sapply(decomp_whole_graph, vcount))
largest_whole_graph <- decomp_whole_graph[[largest_component]]
rm(decomp_whole_graph, largest_component)

# Find strongly-connected components
largest_scc <- rownames(as.data.frame(clusters(largest_whole_graph, mode='strong')[1]))
scc_graph <- subgraph(largest_whole_graph, largest_scc)
rm(largest_whole_graph, largest_scc)

# Remove core metabolites (most connected)
graph_degree <- as.data.frame(degree(scc_graph, v=V(scc_graph), mode='all'))
colnames(graph_degree) <- 'degree'
scc_sub <- rownames(subset(graph_degree, degree < 8))
sub_graph <- subgraph(scc_graph, scc_sub)
graph_degree <- as.data.frame(degree(sub_graph, v=V(sub_graph), mode='all'))
colnames(graph_degree) <- 'degree'
scc_sub <- rownames(subset(graph_degree, degree > 1))
final_graph <- subgraph(sub_graph, scc_sub)
rm(graph_degree, scc_graph, scc_sub)


largest_scc <- rownames(as.data.frame(clusters(final_graph, mode='strong')[1]))
scc_graph <- subgraph(final_graph, largest_scc)
rm(largest_whole_graph, largest_scc)

# Set visualization parameters
V(final_graph)$color <- 'dodgerblue4'
E(final_graph)$color <- 'gray10' 
V(final_graph)$size <- 5
optlayout <- layout_with_graphopt(final_graph)

# Create figure
pdf(file='~/Desktop/ecoli_scc.pdf', width=5, height=5)
par(mar=c(0,0,0,0), xpd=TRUE)
plot(final_graph, vertex.label=NA, layout=optlayout, edge.arrow.size=0.25, edge.width=1.5)
dev.off()
