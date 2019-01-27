### OTUs and Taxons
rm(list=ls())
getwd()
setwd("C:/Users/pradler/Desktop/Courses/IST Core Course/Pathway analysis")
options(scipen=999)


OTU_Identifier = read.csv("OTUs and Taxons.csv", header = T, na="")
sign_parv = read.csv("Sign_Pvalue_Fold_Parvathy.csv", header = T, na="")
all_parv = read.csv("All_pvalue_fold_parvathy.csv", header = T, na="")

names(OTU_Identifier) = c("OTU", "Taxon")

OTU_UC = sign_parv[1:5]
OTU_CD = sign_parv[7:11]

names(OTU_CD)[1] = "OTU"
names(OTU_UC)[1] = "OTU"

OTU_CD=merge(OTU_CD, OTU_Identifier, by = "OTU")
OTU_UC=merge(OTU_UC, OTU_Identifier, by = "OTU")

shared=merge(OTU_CD, OTU_UC, by= "OTU")

all_UC = all_parv[1:5]
all_CD = all_parv[7:11]
names(all_UC)[1] = "OTU"
names(all_CD)[1] = "OTU"

all_UC = merge(all_UC, OTU_Identifier, by="OTU")
all_CD = merge(all_CD, OTU_Identifier, by="OTU")

names(all_UC) = c("OTU", "pvalue", "fold", "log2(Fold)", "-log10(pvalue", "Taxon")
names(all_CD) = c("OTU", "pvalue", "fold", "log2(Fold)", "-log10(pvalue", "Taxon")

x = paste(all_CD$`log2(Fold)`)
y = paste(all_UC$`log2(Fold)`)
x = as.numeric(x)
y=as.numeric(y)
all_UC$`log2(Fold)` = y
all_CD$`log2(Fold)` = x


fold_UC = all_UC$`log2(Fold)`
pvalue_UC= all_UC$pvalue

fold_CD = all_CD$`log2(Fold)`
pvalue_CD= all_CD$pvalue

filter_by_fold_CD = abs(fold_CD) >= 5
filter_by_fold_UC = abs(fold_UC) >= 5

# P-value filter for "statistical" significance
filter_by_pvalue_CD = pvalue_CD <= 0.05
filter_by_pvalue_UC = pvalue_UC <= 0.05


# Combined filter (both biological and statistical)
filter_combined_CD = filter_by_fold_CD & filter_by_pvalue_CD
filter_combined_UC = filter_by_fold_UC & filter_by_pvalue_UC

fold_cutoff = 5
pvalue_cutoff = 0.05

par(mfrow=c(1,2))
## Volcano Plot for UC

plot(fold_UC, -log10(pvalue_UC), main = "UC vs. control (Species)",  xlab= "log2(Fold)", ylab="-log10(P-Value)", xlim=c(-10,10), ylim=c(0,5))
points (fold_UC[filter_combined_UC & fold_UC > 0], 
        -log10(pvalue_UC)[filter_combined_UC & fold_UC > 0],
        pch = 16, col = "green")
points (fold_UC[filter_combined_UC & fold_UC < 0], 
        -log10(pvalue_UC)[filter_combined_UC & fold_UC < 0],
        pch = 16, col = "red")
abline(v = fold_cutoff, col = "green", lwd = 3)
abline(v = -fold_cutoff, col = "red", lwd = 3)
abline(h = -log10(pvalue_cutoff), col = "gray", lwd = 3)

fold_CD = all_CD$`log2(Fold)`
pvalue_CD= all_CD$pvalue

## Volcano Plot for CD
plot(fold_CD, -log10(pvalue_CD), main = "CD vs. Control (Species) ",  xlab= "log2(Fold)", ylab="-log10(P-Value)", xlim=c(-10,10), ylim=c(0,5))
points (fold_CD[filter_combined_CD & fold_CD > 0], 
        -log10(pvalue_CD)[filter_combined_CD & fold_CD > 0],
        pch = 16, col = "green")
points (fold_CD[filter_combined_CD & fold_CD < 0], 
        -log10(pvalue_CD)[filter_combined_CD & fold_CD < 0],
        pch = 16, col = "red")
abline(v = fold_cutoff, col = "green", lwd = 3)
abline(v = -fold_cutoff, col = "red", lwd = 3)
abline(h = -log10(pvalue_cutoff), col = "gray", lwd = 3)

filtered_CD = na.omit(all_CD[filter_combined_CD,])
filtered_UC = na.omit(all_UC[filter_combined_UC,])

filtered_CD=merge(filtered_CD, OTU_Identifier, by = "OTU")
filtered_UC=merge(filtered_UC, OTU_Identifier, by = "OTU") # Unfortunately no species is present in both diseases

write.csv(filtered_CD, file = "Filtered CD species.csv")
write.csv(filtered_UC, file = "Filtered UC species.csv")
