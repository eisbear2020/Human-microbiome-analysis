## Plots for Pathways ##
rm(list=ls())
getwd()
setwd("C:/Users/pradler/Desktop/Courses/IST Core Course/Pathway analysis")
options(scipen=999)


UC_stats=read.csv("UC_stats.csv", header=T, sep=",", na= "0")
CD_stats=read.csv("CD_stats.csv", header=T, sep=",", na= "0")

UC_stats=read.csv("UC_stats.csv", header=T, sep=",", na= "0")
CD_stats=read.csv("CD_stats.csv", header=T, sep=",", na= "0")

CD_stats = CD_stats[,1:7]

CD_sign = read.csv("CD_stats_sign.csv", header = T, sep=",")
UC_sign = read.csv("UC_stats_sign.csv", header = T, sep=",")

CD_sign = CD_sign[10:13]
CD_sign = CD_sign[1:24,]
UC_sign = UC_sign[10:12]
rownames(CD_sign) = CD_sign$Pathway
rownames(UC_sign) = UC_sign$Pathway



CD_sign = CD_sign[,-1]
UC_sign = UC_sign[,-1]
dev.off()
heatmap.2(as.matrix(CD_sign), Rowv = F, Colv = F,  trace="none", margins=c(13,30), cexCol = 1.5, cexRow = 0.7, keysize = 1.4, na.color = "grey", col=greenred(100)) ## Heatmap of the pValues of the comparison between 
heatmap.2(as.matrix(UC_sign), Rowv = F, Colv = F,  trace="none", margins=c(16,30), cexCol = 1.5, cexRow = 0.8, keysize = 1.4, na.color = "grey") ## Heatmap of the pValues of the comparison between 
