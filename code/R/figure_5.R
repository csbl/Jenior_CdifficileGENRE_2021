
# Growth curve data
growth <- read.delim('~/Desktop/repos/Jenior_Cdifficile_2019/data/630_BDM_growthcurve.tsv', sep='\t', header=TRUE)

collateNorm <- function(col_num, data) {
  blank <- paste0('A', as.character(col_num))
  wells <- c(paste0('B', as.character(col_num)), paste0('C', as.character(col_num)), 
             paste0('D', as.character(col_num)), paste0('E', as.character(col_num)))
  blank <- as.matrix(data[,blank])
  summ_growth <- as.matrix(data[,wells])
  blank <- apply(blank, 1, quantile, probs=c(0.25,0.5,0.75))
  summ_growth <- apply(summ_growth, 1, quantile, probs=c(0.25,0.5,0.75))
  summ_growth <- abs(summ_growth - blank)
  summ_growth <- t(summ_growth)
  
  return(summ_growth)}

# Collate and normalize results
bhi <- collateNorm(1, growth)
bdm <- collateNorm(2, growth)
bdm_gluc <- collateNorm(3, growth)
bdm_neu_cyt <- collateNorm(4, growth)
bdm_neu <- collateNorm(5, growth)
bdm_cyt <- collateNorm(6, growth)

# Format for statistical testing
collateDTW <- function(group1_col, group2_col, data) {
  blank <- paste0('A', as.character(group1_col))
  wells <- c(paste0('B', as.character(group1_col)), paste0('C', as.character(group1_col)), 
               paste0('D', as.character(group1_col)), paste0('E', as.character(group1_col)))
  blank <- as.vector(data[,blank])
  growth1 <- as.data.frame(data[,wells])
  for (x in 1:ncol(growth1)) {
    growth1[,x] <- abs(growth1[,x] - blank)}
  growth1 <- t(growth1)
  growth1 <- as.data.frame(cbind(rep('group1', nrow(growth1)), growth1))
  
  blank <- paste0('A', as.character(group2_col))
  wells <- c(paste0('B', as.character(group2_col)), paste0('C', as.character(group2_col)), 
             paste0('D', as.character(group2_col)), paste0('E', as.character(group2_col)))
  blank <- as.vector(data[,blank])
  growth2 <- as.data.frame(data[,wells])
  for (x in 1:ncol(growth2)) {
    growth2[,x] <- abs(growth2[,x] - blank)}
  growth2 <- t(growth2)
  growth2 <- as.data.frame(cbind(rep('group2', nrow(growth2)), growth2))
  
  comb_growth <- as.data.frame(rbind(growth1, growth2))
  colnames(comb_growth)[1] <- 'condition'
  rownames(comb_growth) <- NULL
  comb_growth$condition <- as.factor(comb_growth$condition)
  comb_growth <- droplevels(comb_growth)
  write.table(comb_growth, file='~/Desktop/temp_od.tsv', 
              quote=FALSE, sep='\t', row.names=FALSE, col.names=TRUE)
  comb_growth <- read.delim('~/Desktop/temp_od.tsv', sep='\t', header=TRUE)
  
  return(comb_growth)}

# Calculate significant differences
library(dtw)
library(vegan)
gluc_nocarb <- collateDTW(1, 2, growth)
gluc_nocarb_dist <- dist(gluc_nocarb[,2:ncol(gluc_nocarb)], method='DTW')
gluc_nocarb_pval <- adonis(gluc_nocarb_dist ~ gluc_nocarb$condition, gluc_nocarb[,2:ncol(gluc_nocarb)], perm=999)$aov.tab[[6]][1]
rm(gluc_nocarb, gluc_nocarb_dist)

gluc_neucyt <- collateDTW(1, 3, growth)
gluc_neucyt_dist <- dist(gluc_neucyt[,2:ncol(gluc_neucyt)], method='DTW')
gluc_neucyt_pval <- adonis(gluc_neucyt_dist ~ gluc_neucyt$condition, gluc_neucyt[,2:ncol(gluc_neucyt)], perm=999)$aov.tab[[6]][1]
rm(gluc_neucyt, gluc_neucyt_dist)

gluc_neu <- collateDTW(1, 4, growth)
gluc_neu_dist <- dist(gluc_neu[,2:ncol(gluc_neu)], method='DTW')
gluc_neu_pval <- adonis(gluc_neu_dist ~ gluc_neu$condition, gluc_neu[,2:ncol(gluc_neu)], perm=999)$aov.tab[[6]][1]
rm(gluc_neu, gluc_neu_dist)

gluc_cyt <- collateDTW(1, 5, growth)
gluc_cyt_dist <- dist(gluc_cyt[,2:ncol(gluc_cyt)], method='DTW')
gluc_cyt_pval <- adonis(gluc_cyt_dist ~ gluc_cyt$condition, gluc_cyt[,2:ncol(gluc_cyt)], perm=999)$aov.tab[[6]][1]
rm(gluc_cyt, gluc_cyt_dist)

neu_cyt <- collateDTW(4, 5, growth)
neu_cyt_dist <- dist(neu_cyt[,2:ncol(neu_cyt)], method='DTW')
neu_cyt_pval <- adonis(neu_cyt_dist ~ neu_cyt$condition, neu_cyt[,2:ncol(neu_cyt)], perm=999)$aov.tab[[6]][1]
rm(neu_cyt, neu_cyt_dist)

neucyt_neu <- collateDTW(3, 4, growth)
neucyt_neu_dist <- dist(neucyt_neu[,2:ncol(neucyt_neu)], method='DTW')
neucyt_neu_pval <- adonis(neucyt_neu_dist ~ neucyt_neu$condition, neucyt_neu[,2:ncol(neucyt_neu)], perm=999)$aov.tab[[6]][1]
rm(neucyt_neu, neucyt_neu_dist)

neucyt_cyt <- collateDTW(3, 5, growth)
neucyt_cyt_dist <- dist(neucyt_cyt[,2:ncol(neucyt_cyt)], method='DTW')
neucyt_cyt_pval <- adonis(neucyt_cyt_dist ~ neucyt_cyt$condition, neucyt_cyt[,2:ncol(neucyt_cyt)], perm=999)$aov.tab[[6]][1]
rm(neucyt_cyt, neucyt_cyt_dist)
rm(growth)

#-----------------------------#

# Sporulation assay
sporulation <- read.delim('~/Desktop/repos/Jenior_Cdifficile_2019/data/veg_spore_cfu.tsv', sep='\t', header=TRUE)

# Subset to experimental groups
bdm_veg <- as.vector(sporulation[,'bdm_veg'])
bdm_gluc_veg <- as.vector(sporulation[,'bdm_gluc_veg'])
bdm_neu5Ac_cyt_veg <- as.vector(sporulation[,'bdm_neu5Ac_cyt_veg'])
bdm_neu5Ac_veg <- as.vector(sporulation[,'bdm_neu5Ac_veg'])
bdm_cyt_veg <- as.vector(sporulation[,'bdm_cyt_veg'])
bdm_spore <- as.vector(sporulation[,'bdm_spore'])
bdm_gluc_spore <- as.vector(sporulation[,'bdm_gluc_spore'])
bdm_neu5Ac_cyt_spore <- as.vector(sporulation[,'bdm_neu5Ac_cyt_spore'])
bdm_neu5Ac_spore <- as.vector(sporulation[,'bdm_neu5Ac_spore'])
bdm_cyt_spore <- as.vector(sporulation[,'bdm_cyt_spore'])
rm(sporulation)

# Calculate log ratio
control <- log10(bdm_spore) / log10(bdm_veg)
control[which(!is.finite(control))] <- 0
gluc <- log10(bdm_gluc_spore) / log10(bdm_gluc_veg)
gluc[which(!is.finite(gluc))] <- 0
neu5Ac_cyt <- log10(bdm_neu5Ac_cyt_spore) / log10(bdm_neu5Ac_cyt_veg)
neu5Ac_cyt[which(!is.finite(neu5Ac_cyt))] <- 0
neu5Ac <- log10(bdm_neu5Ac_spore) / log10(bdm_neu5Ac_veg)
neu5Ac[which(!is.finite(neu5Ac))] <- 0
cyt <- log10(bdm_cyt_spore) / log10(bdm_cyt_veg)
cyt[which(!is.finite(cyt))] <- 0
rm(bdm_veg, bdm_gluc_veg, bdm_neu5Ac_cyt_veg, bdm_neu5Ac_veg, bdm_cyt_veg, 
   bdm_spore, bdm_gluc_spore, bdm_neu5Ac_cyt_spore, bdm_neu5Ac_spore, bdm_cyt_spore)

# Test relevant differences
pvals <- p.adjust(c(wilcox.test(control, gluc, exact=FALSE)$p.value,
                    wilcox.test(control, neu5Ac_cyt, exact=FALSE)$p.value,
                    wilcox.test(control, neu5Ac, exact=FALSE)$p.value,
                    wilcox.test(control, cyt, exact=FALSE)$p.value), method='BH')

# Assemble summary statistics
meds <- c(median(control), median(gluc), median(neu5Ac_cyt), median(neu5Ac), median(cyt))
q25s <- as.numeric(c(quantile(control, 0.25), quantile(gluc, 0.25), quantile(neu5Ac_cyt, 0.25), quantile(neu5Ac, 0.25), quantile(cyt, 0.25)))
q75s <- as.numeric(c(quantile(control, 0.75), quantile(gluc, 0.75), quantile(neu5Ac_cyt, 0.75), quantile(neu5Ac, 0.75), quantile(cyt, 0.75)))
rm(control, gluc, neu5Ac_cyt, neu5Ac, cyt)

#----------------------------------------------------------------------------------------------#

# Generate figures
veg_col <- 'gray90'
spore_col <- 'lightslateblue'
library(scales)

pdf(file='~/Desktop/repos/Jenior_Cdifficile_2019/results/figures/Figure_5.pdf', width=4, height=6)
layout(matrix(c(1,
                2), nrow=2, ncol=1, byrow=TRUE))

# Growth curve
par(mar=c(3,3,0.5,2.3), mgp=c(1.8, 0.7, 0), new=FALSE, xpd=FALSE, lwd=2, las=1, xaxs='i', yaxs='i')
plot(bdm_gluc[,2], type='l', ylab='OD600', xaxt='n', yaxt='n', xlab='Hours', ylim=c(0,0.25), xlim=c(1,176.5))
axis(side=1, at=seq(1,195,19.5), labels=seq(0,18,2), cex.axis=0.7, lwd=2)
axis(side=2, at=seq(0,0.25,0.05), cex.axis=0.7, lwd=2)
abline(h=c(0.05,0.1,0.15,0.2), col='gray', lty=5, lwd=2)
polygon(c(1:195, rev(1:195)), c(bdm[,3], rev(bdm[,1])), col=alpha('black', 0.25), border=NA)
lines(bdm[,2], lwd=3, col='black')
polygon(c(1:195, rev(1:195)), c(bdm_gluc[,3], rev(bdm_gluc[,1])), col=alpha('dodgerblue3', 0.25), border=NA)
lines(bdm_gluc[,2], lwd=3, col='dodgerblue3')
polygon(c(1:195, rev(1:195)), c(bdm_neu_cyt[,3], rev(bdm_neu_cyt[,1])), col=alpha('forestgreen', 0.25), border=NA)
lines(bdm_neu_cyt[,2], lwd=3, col='forestgreen')
polygon(c(1:195, rev(1:195)), c(bdm_neu[,3], rev(bdm_neu[,1])), col=alpha('darkorange2', 0.25), border=NA)
lines(bdm_neu[,2], lwd=3, col='darkorange2')
polygon(c(1:195, rev(1:195)), c(bdm_cyt[,3], rev(bdm_cyt[,1])), col=alpha('firebrick3', 0.25), border=NA)
lines(bdm_cyt[,2], lwd=3, col='firebrick3')
legend('topleft', legend=c('BDM (no additional substrate)',
                           'BDM + D-Glucose',
                           'BDM + N-Acetylneuraminate + Cytidine',
                           'BDM + N-Acetylneuraminate',
                           'BDM + Cytidine'), 
       col=c('black','dodgerblue3','forestgreen','darkorange2','firebrick3'), 
       lwd=2.5, cex=0.5, bg='white')
box()
par(xpd=TRUE)
segments(x0=c(179,179,179,179,183), x1=c(183,183,183,188,188), y0=c(0.115,0.104,0.09,0.048,0.1065))
segments(x0=c(183,188), y0=c(0.115,0.1065), y1=c(0.09,0.048))
text(x=192, y=0.07725, '*', cex=1.3, font=2)
text(x=-28, y=0.25, 'A', cex=1.2, font=2)
par(xpd=FALSE)

# Sporulation
par(mar=c(3,3,1,1), xpd=FALSE, las=1, mgp=c(1.8,0.8,0), lwd=2)
barplot(meds, xlim=c(0,6.2), ylim=c(0,0.7), ylab='Log10 Ratio Spores / Vegetative cells', 
        col='gray72', beside=TRUE, cex.axis=0.7, yaxt='n', cex.lab=0.7)
axis(side=2, at=seq(0,0.6,0.2), cex.axis=0.7, lwd=2)
box()
segments(x0=c(0.7,1.9,3.1,4.3,5.5), y0=q25s, y1=q75s)
segments(x0=c(0.9,2.1,3.3,4.5,5.7), y0=q25s, x1=c(0.5,1.7,2.9,4.1,5.3))
segments(x0=c(0.9,2.1,3.3,4.5,5.7), y0=q75s, x1=c(0.5,1.7,2.9,4.1,5.3))
segments(x0=0.7, x1=c(1.9,3.1,4.3,5.5), y0=c(0.65,0.61,0.57,0.53))
text(x=c(1.3,1.9,2.5,3.1), y=c(0.67,0.63,0.59,0.55), '*', font=2, cex=1)
par(xpd=TRUE)
text(x=c(0.7,1.9,3.1,4.3,5.5), cex=0.7, y=-0.08,
     labels=c('BDM\nno additives','BDM\n+ Glucose',
              'BDM\n+ Neu5Ac\n+ Cytidine',
              'BDM\n+ Neu5Ac','BDM\n+ Cytidine'))
text(x=-0.9, y=0.7, 'B', cex=1.2, font=2)
par(xpd=FALSE)

dev.off()
