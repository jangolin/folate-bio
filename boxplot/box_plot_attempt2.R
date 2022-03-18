library(ggplot2)
library(gridExtra)

t <- read.csv(file="C:/Users/Joshua Zhang/OneDrive - The University of Melbourne/201802 BIOM30003/_qPCR data/box_plotdata/data2.csv", sep=",")
row.names(t) <- t$ï..
t <- t[,2:91]
names(t) <- c(rep('HD', 9), rep('LD', 9), rep('U', 9), rep('FM', 9), rep('Glc-', 9), rep('FA', 9), rep('THF', 9), rep('low pH', 9), rep('DMSO', 9), rep('NaOH', 9))

n <- unique(colnames(t)) 
t_avg <- c() 
for(i in n){ 
  t_avg <- cbind(t_avg, rowMeans(t[, names(t) == i])) 
  } 
colnames(t_avg) <- n 


genes <- rownames(t_avg) 

#t_sd <- read.csv(file="C:/Users/Joshua Zhang/OneDrive - The University of Melbourne/201802 BIOM30003/_qPCR data/box_plotdata/data_sd.csv", sep=",")
#row.names(t_sd) <- t_sd$gene
#t_sd <- t_sd[,2:11]

t_ymin <- read.csv(file="C:/Users/Joshua Zhang/OneDrive - The University of Melbourne/201802 BIOM30003/_qPCR data/box_plotdata/data_ymin.csv", sep=",")
row.names(t_ymin) <- t_ymin$ï..gene
t_ymin <- t_ymin[,2:11]
colnames(t_ymin) <- n 

t_ymax <- read.csv(file="C:/Users/Joshua Zhang/OneDrive - The University of Melbourne/201802 BIOM30003/_qPCR data/box_plotdata/data_ymax.csv", sep=",")
row.names(t_ymax) <- t_ymax$ï..gene
t_ymax <- t_ymax[,2:11]
colnames(t_ymax) <- n 

for(i in 1:length(genes)){
p[[i]] <- ggplot() +
    geom_bar(aes(x = factor(n, levels=unique(n)), y = t_avg[i, ]), stat = 'identity', colour = "black", fill = "skyblue") +
    ggtitle(genes[i]) +
    labs(x = "", y = "log2(FC)") +
    theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
          axis.text = element_text(colour = "black"))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    theme(plot.title = element_text(hjust = 0.5))+
    geom_hline(yintercept=0)
    
print(p[i])
}

