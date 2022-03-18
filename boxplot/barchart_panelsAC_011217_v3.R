#HD Paper New RTPCR:

#   Load the necessary libraries.
library(plyr)
library(dplyr)
library(ggplot2)
library(cowplot)

#  Load the necessary data files.
rtpcr_hdld <- read.csv("~/Dropbox/~~ONDECK/November_NewQPCR/datafiles/barchartv4.csv")
#rtpcr_hdldui <- read.csv("~/Desktop/qpcr/dataforR_hdldui.csv")

colvec <-c("lightcoral", "lightgray", "gray30", "light blue", "mediumorchid")
colfill <-c("firebrick3", "lightcoral", "navy", "mediumpurple3")
colfill2 <-c("lightcoral", "mediumpurple3")
shapvec <- c(15, 16, 17, 18, 4, 3)

rtpcr_hdld$Name <- as.character(rtpcr_hdld$Name)
rtpcr_hdld$Name <- factor(rtpcr_hdld$Name, levels=unique(rtpcr_hdld$Name))
rtpcr_hdld$Assay <- as.character(rtpcr_hdld$Assay)
rtpcr_hdld$Assay <- factor(rtpcr_hdld$Assay, levels=unique(rtpcr_hdld$Assay))

barchart <- ggplot(data = rtpcr_hdld, 
                        aes(x = Name, y = FC, group=Assay, fill=Assay, colour=Assay)) + 
  geom_bar(data = rtpcr_hdld, aes(x = Name, y = FC, group=Assay,fill = Assay, color = Assay),
           size = 0.5, stat = "identity", position = position_dodge(.91), width = 0.8) + 
  geom_errorbar(data = rtpcr_hdld, aes(ymin = FC-Sterr, ymax = FC+Sterr), 
                width = 0.3, color = "black", size = 0.4, position = position_dodge(.91)) +
#  geom_rect(aes(xmin = 0, xmax = 1, ymin = -Inf, ymax = Inf),
#            fill = "pink", alpha = 0.03)+
  scale_fill_manual(name="",breaks=c("HDG", "LDL","HDLD"), 
                     labels=c("HDCM+/LDCM","HDCM/LDCM+","HDCM/LDCM"),values = colvec) +
  xlab("") +
  ylab("Fold change (log2)") +
  scale_color_manual(name="",breaks=c("HDG", "LDL","HDLD"), 
                     labels=c("HDCM+/LDCM","HDCM/LDCM+","HDCM/LDCM"),values = colvec) +
  #  scale_x_discrete(
  #    breaks=c("Cysteine", "Glutamic acid", "Glutamine", "Isoleucine", "Methionine", "Proline", "Tyrosine"),
  #    labels=c("Cys", "Glu", "Gln", "Ile", "Met", "Pro", "Tyr")) +
  theme(legend.position = c(.75,.9),
        legend.title = element_text(size=8),
        legend.text = element_text(size=8),
        axis.text.x  = element_text(size=8, colour="black", angle = 90, hjust = 1),  
        axis.title.x = element_text(size = 8, vjust=0.1),
        axis.text.y = element_text(size=9, colour="black"),
        axis.title.y = element_text(vjust=0.2, size = 9), 
        legend.key.size =  unit(0.3, "cm")
  )
barchart

#LDLD <- ggplot(data = rtpcr_hdld, 
#                aes(x = Name, y = LDL.LD)) + 
#  geom_bar(data = rtpcr_hdld, aes(x = Name, y = LDL.LD),
#           size = 0.5, stat = "identity", position = position_dodge(.91), width = 0.8,
#           fill = "gray", color = "gray") + 
#  geom_errorbar(data = rtpcr_hdld, aes(ymin = LDL.LD-LDL.Ldpval, ymax = LDL.LD+LDL.Ldpval), 
#                width = 0.3, color = "black", size = 0.4, position = position_dodge(.91)) +
  #  geom_rect(aes(xmin = 0, xmax = 1, ymin = -Inf, ymax = Inf),
  #            fill = "pink", alpha = 0.03)+
  #  scale_fill_manual(name="", values = c('black','gray')) +
  #  scale_color_manual(name="", values = c('black','gray')) +
#  xlab("") +
#  ylab("Fold change (log2)") +
  #  scale_x_discrete(
  #    breaks=c("Cysteine", "Glutamic acid", "Glutamine", "Isoleucine", "Methionine", "Proline", "Tyrosine"),
  #    labels=c("Cys", "Glu", "Gln", "Ile", "Met", "Pro", "Tyr")) +
#  theme(legend.position = c(.75,.9),
#        legend.title = element_text(size=8),
#        legend.text = element_text(size=8),
#        axis.text.x  = element_text(size=8, colour="black", angle = 90, hjust = 1),  
#        axis.title.x = element_text(size = 8, vjust=0.1),
#        axis.text.y = element_text(size=9, colour="black"),
#        axis.title.y = element_text(vjust=0.2, size = 9), 
#        legend.key.size =  unit(0.3, "cm")
#  )
#LDLD


pear <- read.csv("~/Dropbox/~~ONDECK/November_NewQPCR/datafiles/pearsonv2.csv")
cor.test(pear$M.logFC,pear$Q.logFC, method="pearson")

##Highlight genes that have an absolute fold change > 2 and a p-value < Bonferroni cut-off
pear$threshold = as.factor(pear$Same=="Y")
#bonferroni or not
#pear$threshold = as.factor(abs(pear$Q.FC) > 2 & pear$Q.pval < 0.05/54)

##Construct the plot object
g = ggplot(data=pear, aes(x=Q.logFC, y=-log10(Q.pval), colour=threshold)) +
  geom_point(alpha=0.4, size=1.75) +
  theme(legend.position = "none",
        axis.text.x  = element_text(size=8, colour="black", angle = 90, hjust = 1),  
        axis.title.x = element_text(size = 8, vjust=0.1),
        axis.text.y = element_text(size=9, colour="black"),
        axis.title.y = element_text(vjust=0.2, size = 9), 
        legend.key.size =  unit(0.3, "cm")
        ) +
  xlim(c(-1.5, 1.5)) + ylim(c(0, 6)) +
  xlab("log2 fold change") + ylab("-log10 p-value")
g

gtext <- g + geom_text(aes(x=pear$Q.logFC[1], y=-log10(pear$Q.pval[1]),
                  label=pear$Name[1]), colour="black",size=2.5)
gtext <- gtext + geom_text(aes(x=pear$Q.logFC[3], y=-log10(pear$Q.pval[3]),
                  label=pear$Name[3]), colour="black",size=2.5)
gtext <- gtext + geom_text(aes(x=pear$Q.logFC[2], y=-log10(pear$Q.pval[2]),
                               label=pear$Name[2]), colour="black",size=2.5)
gtext <- gtext + geom_text(aes(x=pear$Q.logFC[14], y=-log10(pear$Q.pval[14]),
                               label=pear$Name[14]), colour="black",size=2.5)
gtext <- gtext + geom_text(aes(x=pear$Q.logFC[23], y=-log10(pear$Q.pval[23]),
                               label=pear$Name[23]), colour="black",size=2.5)
gtext
gtext_ann <- gtext + geom_hline(yintercept=c(10^0.05), linetype="dotted")
gtext_ann

pdf("Figure8_FINAL.pdf",width=4,height=3)
gtext_ann
dev.off()



multig <- plot_grid(barchart, gtext, labels=c("A", "B"), 
                       ncol=2, rel_widths=c(2,1), label_size=10)
multig

pdf("newFigure8_011217.pdf", width=12, height=3)
multig
dev.off()


#pdf("Figure9_301117.pdf",width=8,height=4,paper='special') 
#gcompare_HDLD
#dev.off()







#p-value = 0.026
#Correlation coefficient by pearson correlation = 0.36


