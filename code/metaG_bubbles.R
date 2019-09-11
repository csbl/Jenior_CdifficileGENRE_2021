
# libraries
library(packcircles)
library(ggplot2)
library(viridis)
library(gridExtra)

# Raw data
cefoperazone <- c(1350,1055,620,1354,923,4266)
streptomycin <- c(365,759,1243,1663,2533,274,811,278,371,295,1485,1959,1660,1890)
clindamycin <- c(1889,699,1533,1917)
no_abx <- c(2054,2311,233,1331,1647,2233,1255,2018,8815,256,1500,1846,2304,957,395,1882,
            257,786,2331,488,2731,465,994,398,1092,1400,255,1026,1677,3084,1851,1773,1202,
            393,379,3230,1580,877,1196,995,3407,905,580,490,494,607,213,2219,315,8323,1023,
            585,1525,284,3175,727,1018,316,5214,490,959,761,1278,845,1101,826,1466,520,359,
            196,3115,9663,1658,2367,1209,1503,2193,1074,1650,1208,2923,261,206,3440,3583,2471,233)

# Define taxon colors
cef_colors <- c('red','blue','green','orange')
strep_colors <- c('red','blue','green','orange')
clinda_colors <- c('red','blue','green','orange')
noabx_colors <- c('red','blue','green','orange')


# Find min/max bin sizes
max_genes <- max(c(cefoperazone, streptomycin, clindamycin, no_abx))
min_genes <- min(c(cefoperazone, streptomycin, clindamycin, no_abx))

# Format data
cefoperazone <- sort(cefoperazone, decreasing=TRUE)
streptomycin <- sort(streptomycin, decreasing=TRUE)
clindamycin <- sort(clindamycin, decreasing=TRUE)
no_abx <- sort(no_abx, decreasing=TRUE)



# Function to process data and generate bubble plots
bubble_plot <- function(sizes, max_size, tax_colors=c('none')) {
  
  if (tax_colors[1] == 'none') {tax_colors = viridis(n=length(sizes))}
  tax_colors <- c(tax_colors, 'white')
  circle_colors <- c(rep('black', length(sizes)), 'white')
  sizes <- c(sizes, max_size)
  
  data <- data.frame(bin=1:length(sizes), genes=sizes)
  packing <- circleProgressiveLayout(sizes, sizetype='area')
  packing$radius <- 0.9 * packing$radius
  packing$tax_colors <- tax_colors
  data <- cbind(data, packing)
  data.gg <- circleLayoutVertices(packing, npoints=50)
  
  if (length(sizes) > 14) {txt_size=sizes} else {txt_size=4}
  if (length(txt_size) > 1) {
    bubbles <- ggplot() +
      geom_polygon(data=data.gg, aes(x, y, group=id, col=factor(id), fill=factor(id)), color='white') +
      scale_color_manual(values=tax_colors) +
      scale_fill_manual(values=tax_colors) +
      geom_text(data=data, aes(x, y, size=txt_size, label=round(genes), fontface=2), color='white') +
      theme_void() +
      theme(legend.position='none') +
      coord_equal()
  } else {
    bubbles <- ggplot() +
      geom_polygon(data=data.gg, aes(x, y, group=id, col=factor(id), fill=factor(id)), color='white') +
      scale_color_manual(values=tax_colors) +
      scale_fill_manual(values=tax_colors) +
      geom_text(data=data, aes(x, y, label=round(genes), fontface=2), color='white', size=txt_size) +
      theme_void() +
      theme(legend.position='none') +
      coord_equal()
  }
  
  
  return(bubbles)
}


# Create figure

pdf(file='~/Desktop/metaG_bubbles.pdf', width=9, height=9)

cef_panel <- bubble_plot(cefoperazone, max_size=max_genes)
strep_panel <- bubble_plot(streptomycin, max_size=max_genes)
clinda_panel <- bubble_plot(clindamycin, max_size=max_genes)
no_abx_panel <- bubble_plot(no_abx, max_size=max_genes)

grid.arrange(no_abx_panel, strep_panel, cef_panel, clinda_panel, 
             layout_matrix = rbind(c(1,1,1,2), 
                                   c(1,1,1,3), 
                                   c(1,1,1,4)))

dev.off()


