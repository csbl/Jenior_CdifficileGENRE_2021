

models <- rev(c('iMLTC806cdf','icdf834','iHD992','iCN900','cd630_PATRIC'))
doubling <- rev(c(3.6,24.57,3.6,576.0,327.6))
freeMass <- rev(c(17,2,693,0,0))
blocked <- rev(c(192,66,71,78,0))
massImbal <- rev(c(34,144,0,132,0))




png(filename='~/Desktop/quality.png', units='in', width=8, height=8, res=300)
layout(matrix(c(1,2,
                3,4), nrow=2, ncol=2, byrow=TRUE))

par(mar=c(4,7,1,1), mgp=c(2.5, 1, 0), las=1, lwd=2)
barplot(doubling, horiz=TRUE, names.arg=models, xlab='Doubling Time (min)', xlim=c(0,600),cex.axis=1.2, col='blue3')
box()

barplot(freeMass, horiz=TRUE, names.arg=models, xlab='Free Cytosolic Metabolites', xlim=c(0,700), cex.axis=1.2, col='chartreuse3')
box()

barplot(blocked, horiz=TRUE, names.arg=models, xlab='No Gene Blocked Reactions', xlim=c(0,200), cex.axis=1.2, col='darkorange2')
box()

barplot(massImbal, horiz=TRUE, names.arg=models, xlab='Mass Imbalanced Reactions', xlim=c(0,150),cex.axis=1.2, col='firebrick3')
box()

dev.off()



