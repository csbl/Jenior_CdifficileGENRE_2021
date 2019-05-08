
# Base model

# Load data
base <- read.delim('~/Desktop/repos/Cdiff_modeling/data/base_shadow.tsv', header=TRUE, sep='\t', row.names=1)
# Generate figure
png(filename='~/Desktop/repos/Cdiff_modeling/results/base_shadow.png', width=850, height=850, res=300)
par(mar=c(4,6,1,1), las=1, xpd=FALSE, mgp=c(2.5,1,0), xaxs='r')
plot(0, type='n', xaxt='n', yaxt='n', xlim=c(-10,10), ylim=c(0,8), pch=15, xlab='pFBA Shadow Price', ylab='')
axis(1, at=c(-10,0,10))
abline(v=0, lwd=1.5, lty=5)
box(lwd=1.5)
# Add data
axis(2, at=seq(0.5,7.5,1), labels=rev(rownames(base)), tck=0, cex.axis=1.2)
points(x=rev(base$shadow_price), y=seq(0.5,7.5,1), pch=19, cex=1.5, col='black')
segments(x0=0, y0=seq(0.5,7.5,1), x1=rev(base$shadow_price), y1=seq(0.5,7.5,1), lwd=3, col='black')
dev.off()

#--------------------------------------#

# Glucose minimal media

# Load data
glucose <- read.delim('~/Desktop/repos/Cdiff_modeling/data/glucose_shadow.tsv', header=TRUE, sep='\t', row.names=1)
# Generate figure
png(filename='~/Desktop/repos/Cdiff_modeling/results/glucose_shadow.png', width=850, height=850, res=300)
par(mar=c(4,6,1,1), las=1, xpd=FALSE, mgp=c(2.5,1,0), xaxs='r')
plot(0, type='n', xaxt='n', yaxt='n', xlim=c(-10,10), ylim=c(0,5), pch=15, xlab='pFBA Shadow Price', ylab='')
axis(1, at=c(-10,0,10))
abline(v=0, lwd=1.5, lty=5)
box(lwd=1.5)
# Add data
axis(2, at=seq(0.5,4.5,1), labels=rev(rownames(glucose)), tck=0, cex.axis=1.2)
points(x=rev(glucose$shadow_price), y=seq(0.5,4.5,1), pch=19, cex=1.5, col='#8e7cc3')
segments(x0=0, y0=seq(0.5,4.5,1), x1=rev(glucose$shadow_price), y1=seq(0.5,4.5,1), lwd=3, col='#8e7cc3')
dev.off()

#--------------------------------------#

# Peptide minimal media

# Load data
peptide <- read.delim('~/Desktop/repos/Cdiff_modeling/data/peptide_shadow.tsv', header=TRUE, sep='\t', row.names=1)
# Generate figure
png(filename='~/Desktop/repos/Cdiff_modeling/results/peptide_shadow.png', width=850, height=850, res=300)
par(mar=c(4,6,1,1), las=1, xpd=FALSE, mgp=c(2.5,1,0), xaxs='r')
plot(0, type='n', xaxt='n', yaxt='n', xlim=c(-15,15), ylim=c(0,7), pch=15, xlab='pFBA Shadow Price', ylab='')
axis(1, at=c(-15,0,15))
abline(v=0, lwd=1.5, lty=5)
box(lwd=1.5)
# Add data
axis(2, at=seq(0.5,6.5,1), labels=rev(rownames(peptide)), tck=0, cex.axis=1.2)
points(x=rev(peptide$shadow_price), y=seq(0.5,6.5,1), pch=19, cex=1.5, col='#e69138')
segments(x0=0, y0=seq(0.5,6.5,1), x1=rev(peptide$shadow_price), y1=seq(0.5,6.5,1), lwd=3, col='#e69138')
dev.off()
