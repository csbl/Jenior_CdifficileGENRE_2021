


# Campos M et al. (2014) A constant size extension drives bacterial cell size homeostasis. 
# Cell 159:1433–46. doi: 10.1016/j.cell.2014.11.022. p.1439 right column top paragraph
# bionumbers.hms.harvard.edu/bionumber.aspx?id=111767

library(vioplot)
library(scales)
library(vegan)
library(plotrix)

#------------------------------------------------------------------------------

# Load in data and format
growth_rates <- as.data.frame(t(read.delim(file='~/Desktop/cdf_growth_rates.tsv', sep='\t', header=FALSE, row.names=1)))
pfba <- as.numeric(growth_rates$base_pfba)
cef <- as.numeric(growth_rates$cef)
clinda <- as.numeric(growth_rates$clinda)
strep <- as.numeric(growth_rates$strep)
rm(growth_rates)

# Subsample data
pfba <- sample(pfba, 500, replace=FALSE)
cef <- sample(cef, 500, replace=FALSE)
clinda <- sample(clinda, 500, replace=FALSE)
strep <- sample(strep, 500, replace=FALSE)

# Find ceiling for data
max_rate <- max(c(max(pfba),max(cef),max(clinda),max(strep)))

# Calculate significant differences
pvals <- p.adjust(c(round(wilcox.test(pfba, cef, exact=FALSE)$p.value,5),
                    round(wilcox.test(pfba, strep, exact=FALSE)$p.value,5),
                    round(wilcox.test(pfba, clinda, exact=FALSE)$p.value,5),
                    round(wilcox.test(strep, clinda, exact=FALSE)$p.value,5),
                    round(wilcox.test(strep, cef, exact=FALSE)$p.value,5),
                    round(wilcox.test(cef, clinda, exact=FALSE)$p.value,5)), method='BH')

# Generate figure
png(filename='~/Desktop/cdf_growth.png', units='in', width=6, height=4, res=300)

par(mar=c(3,3.3,1,1), xpd=FALSE, las=1, mgp=c(2,0.75,0), lwd=2)
plot(0,0,type="n",xlim=c(0.5,4.5), ylim=c(0,0.8),  xaxt = 'n', xlab='', 
     ylab=expression(paste('Calculated Growth Rate (hr'^'-1',')')))
text(x=3.5, y=0.01, labels='RIPTiDe-contextualized', cex=1.2)
vioplot(pfba, at=1, col='white', lwd=2, drawRect=FALSE, add=TRUE)
vioplot(cef, at=2, col='dodgerblue3', lwd=2, drawRect=FALSE, add=TRUE)
vioplot(clinda, at=3, col='firebrick3', lwd=2, drawRect=FALSE, add=TRUE)
vioplot(strep, at=4, col='gray', lwd=2, drawRect=FALSE, add=TRUE)
mtext(c('Unweighted\npFBA','Cefoperazone','Clindamycin','Streptomycin'), 
      side=1, at=c(1,2,3,4), padj=1)
dev.off()



