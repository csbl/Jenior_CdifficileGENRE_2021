

quality <- read.delim('~/Desktop/repositories/Cdiff_modeling/data/quality.tsv', sep='\t', header=T, row.names=1)



pdf(file='~/Desktop/repositories/Cdiff_modeling/results/topology.pdf', width=7, height=5)
par(mar=c(3,4,1,1), las=1)
barplot(as.matrix(quality[1:3,]), ylab="Total", col=c("darkblue","red",'chartreuse3'), beside=TRUE, ylim=c(0,2000))
box(lwd=1.5)
abline(v=16.5,lwd=1.5,lty=5)
text(x=10.5, y=1250, '*', cex=3)
legend('top', legend=c('Metabolites','Reactions','Genes'), 
       pt.bg=c("darkblue","red",'chartreuse3'), pch=22, cex=1.2, pt.cex=2.5)
dev.off()


gpr <- as.vector((quality[4,] / quality[2,]) * 100)
mass <- as.vector((quality[6,] / quality[2,]) * 100)
charge <- as.vector((quality[7,] / quality[2,]) * 100)
blocked <- as.vector((quality[8,] / quality[2,]) * 100)
free <- as.vector((quality[5,] / quality[1,]) * 100)
pdf(file='~/Desktop/repositories/Cdiff_modeling/results/quality_percent.pdf', width=7, height=5)
par(mar=c(3,4,1,1), las=1)
barplot(as.matrix(rbind(gpr,mass,charge,blocked,free)), xlim=c(1,29.5),
        ylab="Percent (%)", col=c('darkorchid3','darkorange2','royalblue2','olivedrab2','gray70'), beside=TRUE, ylim=c(0,100))
box(lwd=1.5)
abline(v=c(6.5,12.5,18.5,24.5), lwd=1.5, lty=2)
text(x=c(14.5,15.5), y=3, 'NA', cex=0.7)
legend('topright', legend=c('Rxns /w no GPR','Mass imbalanced rxns','Charge imbalanced rxns','Blocked rxns', 'Free metabolites'), 
       pt.bg=c('darkorchid3','darkorange2','royalblue2','olivedrab2','gray70'), bg='white', pch=22, cex=0.9, pt.cex=1.5)
dev.off()

gpr <- as.vector(quality[4,])
mass <- as.vector(quality[6,])
charge <- as.vector(quality[7,])
blocked <- as.vector(quality[8,])
free <- as.vector(quality[5,])
pdf(file='~/Desktop/repositories/Cdiff_modeling/results/quality_total.pdf', width=7, height=5)
par(mar=c(3,4,1,1), las=1)
barplot(as.matrix(rbind(gpr,mass,charge,blocked,free)), xlim=c(1,29.5),
        ylab="Total", col=c('darkorchid3','darkorange2','royalblue2','olivedrab2','gray70'), beside=TRUE, ylim=c(0,1250))
box(lwd=1.5)
abline(v=c(6.5,12.5,18.5,24.5), lwd=1.5, lty=2)
text(x=c(14.5,15.5), y=30, 'NA', cex=0.7)
legend('topright', legend=c('Rxns /w no GPR','Mass imbalanced rxns','Charge imbalanced rxns','Blocked rxns', 'Free metabolites'), 
       pt.bg=c('darkorchid3','darkorange2','royalblue2','olivedrab2','gray70'), bg='white', pch=22, cex=0.9, pt.cex=1.5)
dev.off()



doubling <- as.matrix(quality[9,])
pdf(file='~/Desktop/repositories/Cdiff_modeling/results/doubling_complete.pdf', width=7, height=5)
par(mar=c(3,4,1,1), las=1)
barplot(doubling, col='forestgreen', ylim=c(0,110), ylab='Minutes')
box(lwd=1.5)
text(x=1,y=103,'Complete media', cex=1.5)
abline(h=38, lwd=2, col='red')
dev.off()


doubling <- as.matrix(quality[10,])
pdf(file='~/Desktop/repositories/Cdiff_modeling/results/doubling_mdm.pdf', width=7, height=5)
par(mar=c(3,4,1,1), las=1)
barplot(doubling, col='forestgreen', ylim=c(0,160), ylab='Minutes')
box(lwd=1.5)
text(x=1,y=150,'Minimal media', cex=1.5)
abline(h=64, lwd=2, col='red')
dev.off()

