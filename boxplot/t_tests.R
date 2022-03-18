
z <- read.csv("C:/Users/Joshua Zhang/OneDrive - The University of Melbourne/201802 BIOM30003/_qPCR data/ddCt_ttest.csv")
z <- z[,2:10]
b <- c()

for (w in c(1, 10, 19, 28, 37, 46, 55, 64, 73, 82, 91)) {
  for (q in w:(w+8)) {
    p <- t.test(as.numeric(as.vector(z[q, ])), as.numeric(as.vector(z[w, ])), paired = FALSE)$p.value
    b <- cbind(b, p)
  }
}

write.csv(b, "C:/Users/Joshua Zhang/OneDrive - The University of Melbourne/201802 BIOM30003/_qPCR data/ddCtttestval.csv")

z2 <- read.csv("C:/Users/Joshua Zhang/OneDrive - The University of Melbourne/201802 BIOM30003/_qPCR data/box_plotdata/dCt_ttest.csv")
z2 <- z2[,2:10]
b2 <- c()

for (w in c(4, 14, 24, 34, 44, 54, 64, 74, 84, 94, 104)) {
  for (q in (w-3):(w+6)) {
    p2 <- t.test(as.numeric(as.vector(z2[q, ])), as.numeric(as.vector(z2[w, ])), paired = FALSE)$p.value
    b2 <- cbind(b2, p2)
  }
}

write.csv(b2, "C:/Users/Joshua Zhang/OneDrive - The University of Melbourne/201802 BIOM30003/_qPCR data/dCtttestval.csv")

a <- read.csv("C:/Users/Joshua Zhang/OneDrive - The University of Melbourne/201802 BIOM30003/_qPCR data/box_plotdata/DHFR-TS.csv")
t.test(as.numeric(as.vector(a$THF)), as.numeric(as.vector(a$DMSO)), var.equal = FALSE, paired = FALSE)$p.value
t.test(as.numeric(as.vector(a$FA)), as.numeric(as.vector(a$NaOH)), var.equal = FALSE, paired = FALSE)$p.value
t.test(as.numeric(as.vector(a$Glc.)), as.numeric(as.vector(a$HD)), var.equal = FALSE, paired = FALSE)$p.value

a_ <- read.csv("C:/Users/Joshua Zhang/OneDrive - The University of Melbourne/201802 BIOM30003/_qPCR data/box_plotdata/data2_transposed.csv")
a_ <- a_[,2:12]
colnames(a_) <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11)

v1 <- c()
v2 <- c()
v3 <- c()
for (g in c(1:11)) {
  p1_ <- t.test(as.numeric(as.vector(a_[1:9, g])), as.numeric(as.vector(a_[37:45, g])), var.equal = FALSE, paired = FALSE)$p.value
  v1 <- cbind(v1, p1_)
  p2_ <- t.test(as.numeric(as.vector(a_[55:63, g])), as.numeric(as.vector(a_[73:81, g])), var.equal = FALSE, paired = FALSE)$p.value
  v2 <- cbind (v2, p2_)
  p3_ <- t.test(as.numeric(as.vector(a_[46:54, g])), as.numeric(as.vector(a_[82:90, g])), var.equal = FALSE, paired = FALSE)$p.value
  v3 <- cbind (v3, p3_)
}
v <- rbind(v1, v2, v3)
write.csv(v, "C:/Users/Joshua Zhang/OneDrive - The University of Melbourne/201802 BIOM30003/_qPCR data/ttestpairs.csv")

#ANOVA & Tukey

#DHFRTS
anova <- read.csv(file="C:/Users/Joshua Zhang/OneDrive - The University of Melbourne/201802 BIOM30003/_qPCR data/box_plotdata/data2_transposed.csv", header = TRUE, sep=",")
anova <- rbind(anova[1:27,], anova[37:90,])
colnames(anova) <- c("condition", "DHFRTS", "PTPS", "PFAM", "DHPSPPPK", "PBAS", "GCH1", "HAD3", "SHMT", "MTFMT", "DHFSFPGS", "ADCL")

aovv <- aov(cbind(DHFRTS, PTPS, PFAM, DHPSPPPK, PBAS, GCH1, HAD3, SHMT, MTFMT, DHFSFPGS, ADCL) ~ condition, data = anova)


dhfrts_aov <- aov(DHFRTS ~ condition, data = anova)
dhfrts_aovsum <- summary(dhfrts_aov)
dhfrts_hsd <- TukeyHSD(dhfrts_aov)
#dhfrts_aovsum
#dhfrts_hsd
#plot(dhfrts_hsd, xlim=c(0,2), las=1)
dhfrts_ <- dhfrts_hsd$condition
dhfrts_ <- cbind(dhfrts_[, 1], dhfrts_[, 4])
colnames(dhfrts_) <- c("diff", "p")


#PTPS
ptps_aov <- aov(PTPS ~ condition, data = anova)
ptps_aovsum <- summary(ptps_aov)
ptps_hsd <- TukeyHSD(ptps_aov)
#ptps_aovsum
#ptps_hsd
#plot(ptps_hsd, xlim=c(0,2), las=1)
ptps_ <- ptps_hsd$condition
ptps_ <- cbind(ptps_[, 1], ptps_[, 4])
colnames(ptps_) <- c("diff", "p")


#PFAM
pfam_aov <- aov(PFAM ~ condition, data = anova)
pfam_aovsum <- summary(pfam_aov)
pfam_hsd <- TukeyHSD(pfam_aov)
#pfam_aovsum
#pfam_hsd
#plot(pfam_hsd, xlim=c(0,2), las=1)
pfam_ <- pfam_hsd$condition
pfam_ <- cbind(pfam_[, 1], pfam_[, 4])
colnames(pfam_) <- c("diff", "p")


#DHPSPPPK
dhpspppk_aov <- aov(DHPSPPPK ~ condition, data = anova)
dhpspppk_aovsum <- summary(dhpspppk_aov)
dhpspppk_hsd <- TukeyHSD(dhpspppk_aov)
#dhpspppk_aovsum
#dhpspppk_hsd
#plot(dhpspppk_hsd, xlim=c(0,2), las=1)
dhpspppk_ <- dhpspppk_hsd$condition
dhpspppk_ <- cbind(dhpspppk_[, 1], dhpspppk_[, 4])
colnames(dhpspppk_) <- c("diff", "p")


#PBAS
pbas_aov <- aov(PBAS ~ condition, data = anova)
pbas_aovsum <- summary(pbas_aov)
pbas_hsd <- TukeyHSD(pbas_aov)
#pbas_aovsum
#pbas_hsd
#plot(pbas_hsd, xlim=c(0,2), las=1)
pbas_ <- pbas_hsd$condition
pbas_ <- cbind(pbas_[, 1], pbas_[, 4])
colnames(pbas_) <- c("diff", "p")


#GCH1
gch1_aov <- aov(GCH1 ~ condition, data = anova)
gch1_aovsum <- summary(gch1_aov)
gch1_hsd <- TukeyHSD(gch1_aov)
#gch1_aovsum
#gch1_hsd
#plot(gch1_hsd, xlim=c(0,2), las=1)
gch1_ <- gch1_hsd$condition
gch1_ <- cbind(gch1_[, 1], gch1_[, 4])
colnames(gch1_) <- c("diff", "p")


#HAD3
had3_aov <- aov(HAD3 ~ condition, data = anova)
had3_aovsum <- summary(had3_aov)
had3_hsd <- TukeyHSD(had3_aov)
#had3_aovsum
#had3_hsd
#plot(had3_hsd, xlim=c(0,2), las=1)
had3_ <- had3_hsd$condition
had3_ <- cbind(had3_[, 1], had3_[, 4])
colnames(had3_) <- c("diff", "p")


#SHMT
shmt_aov <- aov(SHMT ~ condition, data = anova)
shmt_aovsum <- summary(shmt_aov)
shmt_hsd <- TukeyHSD(shmt_aov)
#shmt_aovsum
#shmt_hsd
#plot(shmt_hsd, xlim=c(0,2), las=1)
shmt_ <- shmt_hsd$condition
shmt_ <- cbind(shmt_[, 1], shmt_[, 4])
colnames(shmt_) <- c("diff", "p")


#MTFMT
mtfmt_aov <- aov(MTFMT ~ condition, data = anova)
mtfmt_aovsum <- summary(mtfmt_aov)
mtfmt_hsd <- TukeyHSD(mtfmt_aov)
#mtfmt_aovsum
#mtfmt_hsd
#plot(mtfmt_hsd, xlim=c(0,2), las=1)
mtfmt_ <- mtfmt_hsd$condition
mtfmt_ <- cbind(mtfmt_[, 1], mtfmt_[, 4])
colnames(mtfmt_) <- c("diff", "p")


#DHFSFPGS
dhfsfpgs_aov <- aov(DHFSFPGS ~ condition, data = anova)
dhfsfpgs_aovsum <- summary(dhfsfpgs_aov)
dhfsfpgs_hsd <- TukeyHSD(dhfsfpgs_aov)
#dhfsfpgs_aovsum
#dhfsfpgs_hsd
#plot(dhfsfpgs_hsd, xlim=c(0,2), las=1)
dhfsfpgs_ <- dhfsfpgs_hsd$condition
dhfsfpgs_ <- cbind(dhfsfpgs_[, 1], dhfsfpgs_[, 4])
colnames(dhfsfpgs_) <- c("diff", "p")


#ADCL
adcl_aov <- aov(ADCL ~ condition, data = anova)
adcl_aovsum <- summary(adcl_aov)
adcl_hsd <- TukeyHSD(adcl_aov)
#adcl_aovsum
#adcl_hsd
#plot(adcl_hsd, xlim=c(0,2), las=1)
adcl_ <- adcl_hsd$condition
adcl_ <- cbind(adcl_[, 1], adcl_[, 4])
colnames(adcl_) <- c("diff", "p")



library(ggplot2)

png("C:/Users/Joshua Zhang/OneDrive - The University of Melbourne/201802 BIOM30003/_qPCR data/box_plotdata/figures_boxplots/volcano.png", units="in", width=4, height=4, res=400)
ggplot()+
  geom_point(data = as.data.frame(dhfrts_), aes(x=diff, y=-log10(p)), alpha=0.4, size=1.75, colour="black")+
  geom_point(data = as.data.frame(ptps_), aes(x=diff, y=-log10(p)), alpha=0.4, size=1.75, colour="black")+
  geom_point(data = as.data.frame(pfam_), aes(x=diff, y=-log10(p)), alpha=0.4, size=1.75, colour="black")+
  geom_point(data = as.data.frame(dhpspppk_), aes(x=diff, y=-log10(p)), alpha=0.4, size=1.75, colour="black")+
  geom_point(data = as.data.frame(pbas_), aes(x=diff, y=-log10(p)), alpha=0.4, size=1.75, colour="black")+
  geom_point(data = as.data.frame(gch1_), aes(x=diff, y=-log10(p)), alpha=0.4, size=1.75, colour="black")+
  geom_point(data = as.data.frame(had3_), aes(x=diff, y=-log10(p)), alpha=0.4, size=1.75, colour="black")+
  geom_point(data = as.data.frame(shmt_), aes(x=diff, y=-log10(p)), alpha=0.4, size=1.75, colour="black")+
  geom_point(data = as.data.frame(mtfmt_), aes(x=diff, y=-log10(p)), alpha=0.4, size=1.75, colour="black")+
  geom_point(data = as.data.frame(dhfsfpgs_), aes(x=diff, y=-log10(p)), alpha=0.4, size=1.75, colour="black")+
  geom_point(data = as.data.frame(adcl_), aes(x=diff, y=-log10(p)), alpha=0.4, size=1.75, colour="black")+
  theme_bw()+
  theme(legend.position = "none",
        axis.text.x  = element_text(size=8, colour="black", angle = 90, hjust = 1),  
        axis.title.x = element_text(size = 8, vjust=0.1),
        axis.text.y = element_text(size=9, colour="black"),
        axis.title.y = element_text(vjust=0.2, size = 9), 
        legend.key.size =  unit(0.3, "cm")
  ) +
  xlim(c(-1.5, 1.5)) + ylim(c(0, 15)) +
  geom_text(data = as.data.frame(adcl_), aes(x=diff[16], y=-log10(p[16])), label="ADCL:HD/Glc-", vjust=-1, colour="black",size=2.5)+
  geom_point(data = as.data.frame(adcl_), aes(x=diff[16], y=-log10(p[16])), alpha=0.4, size=1.75, colour="red")+
  geom_text(data = as.data.frame(shmt_), aes(x=diff[16], y=-log10(p[16])), label="SHMT:HD/Glc-", vjust=-1, colour="black",size=2.5)+
  geom_point(data = as.data.frame(shmt_), aes(x=diff[16], y=-log10(p[16])), alpha=0.4, size=1.75, colour="red")+
  geom_hline(yintercept=-log10(0.05), linetype="dashed")+
  labs(x = expression('log'[2]*'(FC ratio)'), y = expression('-log'[10]*'(P-value)'))
dev.off()

#extra t test
scat <- read.csv(file="C:/Users/Joshua Zhang/OneDrive - The University of Melbourne/201802 BIOM30003/_qPCR data/box_plotdata/scatter_cat.csv", sep = "\t")

g1h1 <- t.test(as.numeric(as.vector(scat$g1)), as.numeric(as.vector(scat$h1)), paired = FALSE)$p.value
g2h2 <- t.test(as.numeric(as.vector(scat$g2)), as.numeric(as.vector(scat$h2)), paired = FALSE)$p.value
  
#cleanedupvolcano
library(ggplot2)

png("C:/Users/Joshua Zhang/OneDrive - The University of Melbourne/201802 BIOM30003/_qPCR data/box_plotdata/figures_boxplots/volcano_clean.png", units="in", width=4, height=4, res=400)
ggplot()+
  geom_point(data = as.data.frame(dhfrts_[c(7,13,16), ]), aes(x=diff, y=-log10(p)), alpha=0.4, size=1.75, colour="black")+
  geom_point(data = as.data.frame(ptps_[c(7,13,16), ]), aes(x=diff, y=-log10(p)), alpha=0.4, size=1.75, colour="black")+
  geom_point(data = as.data.frame(pfam_[c(7,13,16), ]), aes(x=diff, y=-log10(p)), alpha=0.4, size=1.75, colour="black")+
  geom_point(data = as.data.frame(dhpspppk_[c(7,13,16), ]), aes(x=diff, y=-log10(p)), alpha=0.4, size=1.75, colour="black")+
  geom_point(data = as.data.frame(pbas_[c(7,13,16), ]), aes(x=diff, y=-log10(p)), alpha=0.4, size=1.75, colour="black")+
  geom_point(data = as.data.frame(gch1_[c(7,13,16), ]), aes(x=diff, y=-log10(p)), alpha=0.4, size=1.75, colour="black")+
  geom_point(data = as.data.frame(had3_[c(7,13,16), ]), aes(x=diff, y=-log10(p)), alpha=0.4, size=1.75, colour="black")+
  geom_point(data = as.data.frame(shmt_[c(7,13,16), ]), aes(x=diff, y=-log10(p)), alpha=0.4, size=1.75, colour="black")+
  geom_point(data = as.data.frame(mtfmt_[c(7,13,16), ]), aes(x=diff, y=-log10(p)), alpha=0.4, size=1.75, colour="black")+
  geom_point(data = as.data.frame(dhfsfpgs_[c(7,13,16), ]), aes(x=diff, y=-log10(p)), alpha=0.4, size=1.75, colour="black")+
  geom_point(data = as.data.frame(adcl_[c(7,13,16), ]), aes(x=diff, y=-log10(p)), alpha=0.4, size=1.75, colour="black")+
  theme_bw()+
  theme(legend.position = "none",
        axis.text.x  = element_text(size=8, colour="black", angle = 90, hjust = 1),  
        axis.title.x = element_text(size = 8, vjust=0.1),
        axis.text.y = element_text(size=9, colour="black"),
        axis.title.y = element_text(vjust=0.2, size = 9), 
        legend.key.size =  unit(0.3, "cm")
  ) +
  xlim(c(-1.5, 1.5)) + ylim(c(0, 10)) +
  geom_text(data = as.data.frame(adcl_), aes(x=diff[16], y=-log10(p[16])), label="ADCL:HDCM/Glc-", vjust=-1, colour="black",size=2.5)+
  geom_point(data = as.data.frame(adcl_), aes(x=diff[16], y=-log10(p[16])), alpha=0.4, size=1.75, colour="red")+
  geom_text(data = as.data.frame(shmt_), aes(x=diff[16], y=-log10(p[16])), label="SHMT:HDCM/Glc-", vjust=-1, colour="black",size=2.5)+
  geom_point(data = as.data.frame(shmt_), aes(x=diff[16], y=-log10(p[16])), alpha=0.4, size=1.75, colour="red")+
  labs(x = expression('log'[2]*'(FC ratio)'), y = expression('-log'[10]*'(P-value)'))
dev.off()