
# libraries
library(packcircles)
library(ggplot2)
library(viridis)
library(gridExtra)

# Raw data
cefoperazone <- c(1350,1055,620,1354,923,4266)
streptomycin <- c(365,759,1243,1663,2533,274,811,278,371,295,1485,1959,1660,1890)
clindamycin <- c(1889,699,1533,1917)
no_abx <- c(2054,2311,233,1331,1647,2233,1255,2018,8815,256,1500,1846,2304,957,395,1882,257,786,2331,488,2731,465,994,398,1092,1400,255,1026,1677,3084,1851,1773,1202,393,379,3230,1580,877,1196,995,3407,905,580,490,494,607,213,2219,315,8323,1023,585,1525,284,3175,727,1018,316,5214,490,959,761,1278,845,1101,826,1466,520,359,196,3115,9663,1658,2367,1209,1503,2193,1074,1650,1208,2923,261,206,3440,3583,2471,233)

# Format data
max_genes <- max(c(cefoperazone, streptomycin, clindamycin, no_abx))
min_genes <- min(c(cefoperazone, streptomycin, clindamycin, no_abx))
cefoperazone <- c(sort(cefoperazone, decreasing=TRUE),max_genes)
streptomycin <- c(sort(streptomycin, decreasing=TRUE),max_genes)
clindamycin <- c(sort(clindamycin, decreasing=TRUE),max_genes)
no_abx <- c(sort(no_abx, decreasing=TRUE),max_genes)

# Function to process data and generate bubble plots
bubble_plot <- function(sizes, spacing=1, scale_text=FALSE) {
  
  data <- data.frame(bin=paste('MAG\n', 1:length(sizes)), genes=sizes)
  
  packing <- circleProgressiveLayout(sizes, sizetype='area')
  packing$radius <- spacing*packing$radius
  data <- cbind(data, packing)
  data.gg <- circleLayoutVertices(packing, npoints=50)
  
  if (scale_text == FALSE) {
    txt_size=1
  } else {txt_size=sizes}
  
  bubbles <- ggplot() + 
    geom_polygon(data=data.gg, aes(x, y, group=id, fill=id), color='black', alpha=0.6) +
    scale_fill_viridis() +
    geom_text(data=data, aes(x, y, size=txt_size, label=bin), color='black') +
    theme_void() +
    theme(legend.position='none') +
    coord_equal()
  
  return(bubbles)
}

# Create figure
cef_panel <- bubble_plot(cefoperazone, spacing=0.95)
strep_panel <- bubble_plot(streptomycin, spacing=0.95)
clinda_panel <- bubble_plot(clindamycin, spacing=0.95)
no_abx_panel <- bubble_plot(no_abx, scale_text=TRUE)



pdf(file='~/Desktop/metaG_bubbles.pdf', width=10, height=10)

grid.arrange(no_abx_panel, cef_panel, 
             strep_panel, clinda_panel, nrow = 2)

dev.off()


#============================================#

# Grouped bubble testing

library(ggraph)
library(igraph)
library(tidyverse)
data(flare)

# We need a data frame giving a hierarchical structure. Let's consider the flare dataset:
edges=flare$edges

# Usually we associate another dataset that give information about each node of the dataset:
vertices=flare$vertices

# Then we have to make a 'graph' object using the igraph library:
mygraph <- graph_from_data_frame(edges, vertices=vertices)

# Make the plot
ggraph(mygraph, layout='circlepack') + 
  geom_node_circle() +
  theme_void()



# Hide the first level (right)
ggraph(mygraph, layout='circlepack', weight='size') + 
  geom_node_circle(aes(fill=as.factor(depth), color=as.factor(depth) )) +
  scale_fill_manual(values=c('0'='white', '1'=viridis(4)[1], '2'=viridis(4)[2], '3'=viridis(4)[3], '4'=viridis(4)[4])) +
  scale_color_manual( values=c('0'='white', '1'='black', '2'='black', '3'='black', '4'='black') ) +
  theme_void() + 
  theme(legend.position='FALSE') 



