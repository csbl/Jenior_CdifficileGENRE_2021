


pfba_isovalerate <- c(59.72, 103.73, 151.69)
pfba_ornithine <- c(-2.01, -1.72, -1.42)
strep_d_mannose <- c(-102.288, -20.546, 56.17)
cef_d_lactate <- c(-241.79, -133.57, -28.566)

# Log2
pfba_isovalerate <- c(5.90, 6.69, 7.24)
pfba_ornithine <- c(-1.01, -0.78, -0.51)
strep_d_mannose <- c(-6.68, -4.36, 5.81)
cef_d_lactate <- c(-7.92, -7.06, -4.84)



# Generate figure
png(filename='~/Desktop/cdf_exchange.png', units='in', width=3, height=4, res=300)
par(mar=c(3,7,1,1), xpd=FALSE, las=1, mgp=c(2,0.75,0), lwd=2)
plot(0, 0, type='n', xlim=c(-8,8), ylim=c(0,4),  yaxt='n', xlab='Sampled Exchange Flux (Log10)', ylab='', cex.lab=0.7)
rect(xleft=5.9, ybottom=3.25, xright=7.24, ytop=3.75, col='white', border='black')
rect(xleft=-1.01, ybottom=2.25, xright=-0.51, ytop=2.75, col='white', border='black')
rect(xleft=-6.68, ybottom=1.25, xright=5.81, ytop=1.75, col='gray', border='black')
rect(xleft=-7.92, ybottom=0.25, xright=-4.84, ytop=0.75, col='dodgerblue3', border='black')
abline(v=0, lty=5)
mtext(rev(c('Isovalerate\npFBA','Ornithine\npFBA','D-Mannose\nCefoperazone','D-Lactate\nStreptomycin')), 
      side=2, at=c(0.5,1.5,2.5,3.5), adj=1.2)
dev.off()

