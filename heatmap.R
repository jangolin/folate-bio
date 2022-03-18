library(ComplexHeatmap)
library(circlize)


fa <- read.csv('C:/Users/Joshua Zhang/OneDrive - The University of Melbourne/201802 BIOM30003/_qPCR data/R-heatmapdata.csv', stringsAsFactors = F)
row.names(fa) <- fa$ï..Gene
fa <- fa[,2:31]
fa1 <- fa[,1:9]
fa2 <- fa[,13:30]

fa <- cbind(fa1, fa2)
names(fa) <- c(rep('HDCM', 3), rep('LDCM', 3), rep('UCM', 3), rep('Glc-', 3), rep('FA', 3), rep('THF', 3), rep('Low pH', 3), rep('DMSO', 3), rep('NaOH', 3))

fa_norm <- t(scale(t(fa), center = F))

png("C:/Users/Joshua Zhang/OneDrive - The University of Melbourne/201802 BIOM30003/_qPCR data/heatmap-final-clustered.png", units="in", width=10, height=4.5, res=300)
Heatmap(fa_norm, col = colorRamp2(c(min(fa), 0, max(fa)), c('green', 'black', 'red')), cluster_rows = T, cluster_columns = T)
dev.off()
png("C:/Users/Joshua Zhang/OneDrive - The University of Melbourne/201802 BIOM30003/_qPCR data/heatmap-final.png", units="in", width=10, height=4, res=300)
Heatmap(fa_norm, col = colorRamp2(c(min(fa), 0, max(fa)), c('green', 'black', 'red')), cluster_rows = F, cluster_columns = F)
dev.off()