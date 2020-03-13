
library(ggraph)
library(igraph)
library(tidyverse)
library(viridis)

# We need a data frame giving a hierarchical structure. Let's consider the flare dataset:
edges=flare$edges
vertices=flare$vertices
mygraph <- graph_from_data_frame( edges, vertices=vertices )

# Hide the first level (right)
ggraph(mygraph, layout='circlepack', weight='size') + 
  geom_node_circle(aes(fill=as.factor(depth), color=as.factor(depth) )) +
  scale_fill_manual(values=c('0'='white', '1'=viridis(4)[1], '2'=viridis(4)[2], '3'=viridis(4)[3], '4'=viridis(4)[4])) +
  scale_color_manual( values=c('0'='white', '1'='black', '2'='black', '3'='black', '4'='black') ) +
  theme_void() + 
  theme(legend.position='FALSE') 

# Second one: add 2 first levels
ggraph(mygraph, layout='circlepack', weight='size') + 
  geom_node_circle(aes(fill=as.factor(depth), color=as.factor(depth) )) +
  scale_fill_manual(values=c('0'='white', '1'='white', '2'=magma(4)[2], '3'=magma(4)[3], '4'=magma(4)[4])) +
  scale_color_manual( values=c('0'='white', '1'='white', '2'='black', '3'='black', '4'='black')) +
  theme_void() + 
  theme(legend.position='FALSE')





# Add the data.tree library
library(data.tree)

# Rebuild the data
edges=flare$edges
vertices=flare$vertices

# Transform it in a 'tree' format
tree <- FromDataFrameNetwork(edges)

# Then I can easily get the level of each node, and add it to the initial data frame:
mylevels=data.frame(name=tree$Get('name'), level=tree$Get('level') )
vertices=vertices %>% left_join(., mylevels, by=c('name'='name'))

# Now we can add label for level1 and 2 only for example:
vertices=vertices %>% mutate(new_label=ifelse(level==2, shortName, NA))
mygraph <- graph_from_data_frame( edges, vertices=vertices )

# Make the graph
ggraph(mygraph, layout='circlepack', weight='size') + 
  geom_node_circle(aes(fill=as.factor(depth), color=as.factor(depth) )) +
  scale_fill_manual(values=c('0'='white', '1'=viridis(4)[1], '2'=viridis(4)[2], '3'=viridis(4)[3], '4'=viridis(4)[4])) +
  scale_color_manual(values=c('0'='white', '1'='black', '2'='black', '3'='black', '4'='black')) +
  geom_node_label(aes(label=new_label), size=4) +
  theme_void() + 
  theme(legend.position='FALSE', plot.margin=unit(rep(0,4), 'cm'))





