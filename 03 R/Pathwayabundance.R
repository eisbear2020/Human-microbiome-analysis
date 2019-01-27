### Pathway abundance Bacteria Microbiome

rm(list=ls())
getwd()
setwd("C:/Users/pradler/Desktop/Courses/IST Core Course/Pathway analysis")
options(scipen=999)

Pathways <- read_excel("pathabundance_relab.xlsx", na="0")
ID_diagnosis = read.delim("ID_Diagnosis.txt", header = TRUE, sep = "\t", dec = ".") ## ID and diagnosis
ID_pathways = read.delim("Identifier Pathways.txt", header = TRUE, sep = "\t", dec = ".") ## IDs for diagnosis

Pathways$Expression = rowSums(Pathways[3:738], na=T)
Pathways[is.na(Pathways)] <- 0


Pathways = Pathways[rowSums(Pathways[(3:738)], na.rm = T, dims = 1) != 0,] ### Remove pathways that are not expressed at all


## I want to extract the diagnosis information from the ID_diagnosis dataframe for every ID in pathways. Like this, I could split the pathway abundance
## in the control, UC and CD abundances and compare their levels.

result <- merge(ID_pathways, ID_diagnosis, by="ID") ## All Identifiers in the pathways metafile are merged with their respective diagnosis
result = unique(result) ### All duplicates are lost, only 735 Subject IDs with their disease are left


## Extract pathway abundance for ctrl, CD and UC by selecting the Identifier
UC = result[result$diagnosis=="UC",]
CD = result[result$diagnosis=="CD",]
ctrl = result[result$diagnosis=="nonIBD",]

UC = Pathways[ ,c(UC$ID)]
CD = Pathways[ ,c(CD$ID)]
ctrl = Pathways[ ,c(ctrl$ID)] ### Pathway abundance for control, UC and CD

CD=cbind(Pathways$Pathway, Pathways$Species, CD)
ctrl=cbind(Pathways$Pathway, Pathways$Species, ctrl)

colnames(CD)[2] <- "Species"
colnames(CD)[1] <- "Pathway"

colnames(ctrl)[2] <- "Species"
colnames(ctrl)[1] <- "Pathway" ### Somehow the column names got lost. I rename to be sure to have Species, Pathway, ID everywhere

UC$Expression = rowMeans(UC[(3:211)], na.rm=T)
CD$Expression = rowMeans(CD[(3:339)], na.rm=T)
ctrl$Expression = rowMeans(ctrl[(3:189)], na.rm=T) ## Expression levels for each type of disease or control

# Now we have CD, UC or control for each ID that is specified in the pathway abundance file. Next step is to extract the species of interest 
# for control, UC and CD and check if they are differentially abundant.

## Extracting pathways from different species and calculate sum of abundance: UC vs ctrl ##

Mega_UC = UC[UC$Species=="g__Megasphaera.s__Megasphaera_micronuciformis",]
Mega_ctrl = ctrl[ctrl$Species=="g__Megasphaera.s__Megasphaera_micronuciformis",]

Burk_UC = UC[UC$Species=="g__Burkholderiales_noname.s__Burkholderiales_bacterium_1_1_47",]
Burk_ctrl = ctrl[ctrl$Species=="g__Burkholderiales_noname.s__Burkholderiales_bacterium_1_1_47",]

Pseudo_UC = UC[UC$Species=="g__Pseudomonas.s__Pseudomonas_aeruginosa",]
Pseudo_ctrl = ctrl[ctrl$Species=="g__Pseudomonas.s__Pseudomonas_aeruginosa",]

Oscil_UC = UC[grep("Oscillibacter", UC$Species), ]
Oscil_ctrl = ctrl[grep("Oscillibacter", ctrl$Species), ]

Rumino_UC = UC[UC$Species=="g__Ruminococcaceae_noname.s__Ruminococcaceae_bacterium_D16",]
Rumino_ctrl = ctrl[ctrl$Species=="g__Ruminococcaceae_noname.s__Ruminococcaceae_bacterium_D16",]

Prevo_UC = UC[grep("Prevotella", UC$Species), ]
Prevo_ctrl = ctrl[grep("Prevotella", ctrl$Species), ]

UC_species=rbind(Mega_UC, Pseudo_UC, Burk_UC, Oscil_UC, Rumino_UC, Prevo_UC)
ctrl_species=rbind(Mega_ctrl, Pseudo_ctrl, Burk_ctrl, Oscil_ctrl, Rumino_ctrl, Prevo_ctrl)


par(mfrow=c(3,2))# make a plot window
boxplot(Mega_ctrl$Expression, Mega_UC$Expression, main="Mega", names = c("Ctrl", "UC")) ## Abundance of all E.Hallii Pathways. 
boxplot(Burk_ctrl$Expression, Burk_UC$Expression, main="Burk", names = c("Ctrl", "UC")) ## Abundance of all E.Hallii Pathways. 
boxplot(Pseudo_ctrl$Expression, Pseudo_UC$Expression, main="PSeudo", names = c("Ctrl", "UC")) ## Abundance of all E.Hallii Pathways. 
boxplot(Oscil_ctrl$Expression, Oscil_UC$Expression, main="Oscil", names = c("Ctrl", "UC"))
boxplot(Prevo_ctrl$Expression, Prevo_UC$Expression, main="Prevo", names = c("Ctrl", "UC"))
boxplot(Rumino_ctrl$Expression, Rumino_UC$Expression, main="Rumino", names = c("Ctrl", "UC"))


### Get all pathways from 10-20 species. check for any overlapping, expressed pathways with merge and filter function (get rid of all non abundant 
### expressed pathways in the first filtering step, after loading the dataset)


### Check how to plot the pathway abundances in a nice heatmap format, similar to the expression profile. To t-test to compare the means of 
### pathway abundance? Check how they calculated statistics in the paper for this.
UC_test = as.matrix(UC_species[3:211])
ctrl_test = as.matrix(ctrl_species[3:189])


### Check for significant differences in pathway abundances using the Wilcox Test (Ranked Sums, high number, normality can't be assumed)

pvalue = NULL # Empty list for the p-values
Wstat = NULL # Empty list of the W test statistics

for(i in 1 : nrow(UC_test)) { # For each pathway : 
  x = UC_test[i,] # CD pathway number i
  y = ctrl_test[i,] # ctrl pathway number i
  
  if (sum(x)==sum(y)){
    pvalue[i] = 1
    Wstat[i] = 0
  } ### If the sum of expression levels is the same within WT and KO, we assume no significant difference. As R would stop computing, we introduce this if statement
  
  else{
    # Compute wilcox-test between the two conditions: Default test: two sided, assuming non-equal variances
    W = wilcox.test(x, y)
    
    # Put the current p-value in the pvalues list
    pvalue[i] = W$p.value
    # Put the current t-statistic in the tstats list
    Wstat[i] = W$statistic
  }}

head(pvalue)
UC_stat = cbind(as.character(UC_species$Pathway), as.character(UC_species$Species), as.numeric(UC_species$Expression), as.numeric(ctrl_species$Expression),pvalue, Wstat)
head(UC_stat)

colnames(UC_stat) = c("Pathway", "Species", "Expression_UC", "Expression_ctrl" ,"pvalue", "Wstat")
UC_stat$Fold = UC_stat$Expression_ctrl/UC_stat$Expression_UC
write.csv(UC_stat, file = "UC_stats.csv")


### CD ###
Clos_CD = UC[grep("Clostridiales", CD$Species), ]
Clos_ctrl = ctrl[grep("Clostridiales", ctrl$Species), ]

Prevo_CD = UC[grep("Prevotella", CD$Species), ]
Prevo_ctrl = ctrl[grep("Prevotella", ctrl$Species), ]

Ali_CD = UC[grep("Alistipes", CD$Species), ]
Ali_ctrl = ctrl[grep("Alistipes", ctrl$Species), ]

CD_species=rbind(Prevo_CD, Ali_CD, Clos_CD)
ctrl_species=rbind(Prevo_ctrl, Ali_ctrl, Clos_ctrl)


par(mfrow=c(3,2))# make a plot window

boxplot(Prevo_ctrl$Expression, Prevo_CD$Expression, main="Prevo", names = c("Ctrl", "CD"))
boxplot(Ali_ctrl$Expression, Ali_CD$Expression, main="Ali", names = c("Ctrl", "CD"))
boxplot(Clos_ctrl$Expression, Clos_CD$Expression, main="Clos", names = c("Ctrl", "CD"))


### Get all pathways from 10-20 species. check for any overlapping, expressed pathways with merge and filter function (get rid of all non abundant 
### expressed pathways in the first filtering step, after loading the dataset)


### Check how to plot the pathway abundances in a nice heatmap format, similar to the expression profile. To t-test to compare the means of 
### pathway abundance? Check how they calculated statistics in the paper for this.
CD_test = as.matrix(CD_species[3:211])
ctrl_test = as.matrix(ctrl_species[3:189])


### Check for significant differences in pathway abundances using the Wilcox Test (Ranked Sums, high number, normality can't be assumed)

pvalue = NULL # Empty list for the p-values
Wstat = NULL # Empty list of the W test statistics

for(i in 1 : nrow(CD_test)) { # For each pathway : 
  x = CD_test[i,] # CD pathway number i
  y = ctrl_test[i,] # ctrl pathway number i
  
  if (sum(x)==sum(y)){
    pvalue[i] = 1
    Wstat[i] = 0
  } ### If the sum of expression levels is the same within WT and KO, we assume no significant difference. As R would stop computing, we introduce this if statement
  
  else{
    # Compute wilcox-test between the two conditions: Default test: two sided, assuming non-equal variances
    W = wilcox.test(x, y)
    
    # Put the current p-value in the pvalues list
    pvalue[i] = W$p.value
    # Put the current t-statistic in the tstats list
    Wstat[i] = W$statistic
  }}

head(pvalue)
CD_stat = cbind(as.character(CD_species$Pathway), as.character(CD_species$Species), as.numeric(CD_species$Expression), as.numeric(ctrl_species$Expression),pvalue, Wstat)
head(CD_stat)

colnames(CD_stat) = c("Pathway", "Species", "Expression_CD", "Expression_ctrl" ,"pvalue", "Wstat")

write.csv(CD_stat, file = "CD_stats.csv")
