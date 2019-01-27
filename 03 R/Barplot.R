# Stacked Bar Plots for Phylum Proportion # 
rm(list=ls())
getwd()
setwd("C:/Users/pradler/Desktop/Courses/IST Core Course")
options(scipen=999)

Species = read.csv("Phylum Proportion.csv", header = T, sep = ",")
rownames(Species)=Species$Phylum
Species$Count = (Species$Count/982) * 100 ## Normalize the data by total species abundance in percent
Species$Species = "Species"



ggplot(Species, aes(x=Species, y = Count, fill = Phylum)) + 
  geom_bar(stat = "identity") +
  xlab("") +
  ylab("Percentage") +
  theme_bw()
