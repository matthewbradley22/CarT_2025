#Load packages
library(Seurat)
library(dplyr)
library(nnet)
library(tidyr)
library(corrplot)

#Make plots work properly
options(bitmapType="cairo") 

#Load in some functions for project
source("mystore/cartdata/scripts/CarT_project_functions.R")

#Load here for t cells with only singlets  (vireo and scdblfinder used), and  no d0 hypoxia groups. 
deseqDat <- LoadSeuratRds("./mystore/cartdata/data/tcell_deseqdat_singlets.rds")

#Assigns cell phase to each cell. Can skip this because deseqDat is now saved with phase included
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
deseqDat <- CellCycleScoring(deseqDat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

day7 <- subset(deseqDat, day == "D7")
day13 <- subset(deseqDat, day == "D13")

CD4 <- subset(deseqDat, CD_pred == "CD4")
CD8 <- subset(deseqDat, CD_pred == "CD8")

#Plot cell cycle
#Saving 9x5 dimensions
totalCountsCD4 = CD4[[]] %>% group_by(day, CAR, hypoxia) %>% summarise(total = n())
barPlotDatCD4 <- CD4[[]] %>% group_by(day, CAR, hypoxia, Phase) %>% summarise(total = n()) 
barPlotDatCD4$overallCounts =rep(totalCountsCD4$total, each = 3)
barPlotDatCD4 %>% mutate(group = paste(day, CAR, hypoxia, sep = "_")) %>% 
  mutate(proportion = total/overallCounts) %>% 
  ggplot(aes(x = factor(group, levels = unique(group)), y = proportion, fill = Phase))+
  geom_bar(stat = 'identity')+
  theme(axis.text.x = element_text(angle =90, vjust = 0.6, size = 13))+
  ggtitle("Cell Phase CD4")

totalCountsCD8 = CD8[[]] %>% group_by(day, CAR, hypoxia) %>% summarise(total = n())
barPlotDatCD8 <- CD8[[]] %>% group_by(day, CAR, hypoxia, Phase) %>% summarise(total = n()) 
barPlotDatCD8$overallCounts =rep(totalCountsCD8$total, each = 3)
barPlotDatCD8 %>% mutate(group = paste(day, CAR, hypoxia, sep = "_")) %>% 
  mutate(proportion = total/overallCounts) %>% 
  ggplot(aes(x = factor(group, levels = unique(group)), y = proportion, fill = Phase))+
  geom_bar(stat = 'identity')+
  theme(axis.text.x = element_text(angle =90, vjust = 0.6, size = 13))+
  ggtitle("Cell Phase CD8")

#Split plots by day instead
totalCounts = day7[[]] %>% group_by(hypoxia, CAR) %>% summarise(total = n())

barPlotDat7 <- day7[[]] %>% group_by(hypoxia, CAR,Phase) %>% summarise(total = n()) 
barPlotDat7$overallCounts =rep(totalCounts$total, each = 3)
barPlotDat7 %>% mutate(group = paste(hypoxia, CAR, sep = "_")) %>% 
  mutate(proportion = total/overallCounts) %>% 
  ggplot(aes(x = group, y = proportion, fill = Phase))+
  geom_bar(stat = 'identity')+
  theme(axis.text.x = element_text(angle =45, vjust = 0.6))+
  ggtitle("Cell Phase Day 7")

totalCounts13 = day13[[]] %>% group_by(hypoxia, CAR) %>% summarise(total = n())
barPlotDat13 <- day13[[]] %>% group_by(hypoxia, CAR,Phase) %>% summarise(total = n()) 
barPlotDat13$overallCounts =rep(totalCounts13$total, each = 3)
barPlotDat13 %>% mutate(group = paste(hypoxia, CAR, sep = "_")) %>% 
  mutate(proportion = total/overallCounts) %>% 
  ggplot(aes(x = group, y = proportion, fill = Phase))+
  geom_bar(stat = 'identity')+
  theme(axis.text.x = element_text(angle =45, vjust = 0.6))+
  ggtitle("Cell Phase Day 13")

#Cell cycle anaysis
test <- multinom(Phase ~ hypoxia + CAR, data = day7[[]])
summary(test)


#Look at variable effects on phase
dds7 <- getPseudoBulkObject(day7, c("hypoxia","CAR", "Phase"))
dds13 <- getPseudoBulkObject(day13, c("hypoxia","CAR","donor_id", "Phase"))

day7TestDat <- data.frame(colData(dds7)) %>% select(hypoxia, CAR, Phase, n) %>% 
  pivot_wider(names_from = Phase, values_from = n) %>% as.data.frame()
rownames(day7TestDat) = paste(day7TestDat$hypoxia, day7TestDat$CAR, sep = "_")
day7Chisq <- day7TestDat %>% select(G1, G2M, S) %>% chisq.test()
residuals(day7Chisq)

corrplot(day7Chisq$observed - day7Chisq$expected, is.corr = F,  cl.ratio=1)

#rerun dds 7 with donor for plotting
dds7 <- getPseudoBulkObject(day7, c("hypoxia", "donor_id","CAR", "Phase"))
#Create PCA plots grouped by various variables
getPCAGrid(dds7, c("hypoxia", "CAR", "donor_id", "Phase"), rowNum = 2)
getPCAGrid(dds13, c("hypoxia", "CAR", "donor_id", "Phase"), rowNum = 2)


#### Examine variables overall effects on model. This was an early plan that we didn't follow up on ####
# Run DESeq2 differential expression analysis
#LRT to compare to reduced formula
#Found it difficult to write a function for these lines because reduced won't take a string
#Example formula "LRT p-value: '~ hypoxia + CAR + day + donor_id' vs '~ CAR + day + donor_id'"

covariateImpact <- function(dat, varSet){
  reducedFormula <- reformulate(termlabels = varSet)
  deseqObj <- DESeq(dat, test="LRT", reduced=reducedFormula)
}

#noHypoxTest <- covariateImpact(dds, c("CAR", "day", "donor_id"))

#Pass vars of interest and call covariateImpact on all combinations of 
#length n-1 where n is number of vars
effectPerVariable <- function(GeneDat, vars){
  #Take all var combos length n-1 
  varCombos = list()
  if(length(vars) > 1){
    for(i in 1:length(vars)){
      combo = vars[-i]
      varCombos[[i]] = combo
    }
    lapply(varCombos, covariateImpact, dat = GeneDat)
  }else{
    stop("Need >1 variables to compare impacts of variables")
  }
}

allVars7 <- effectPerVariable(dds7, c("hypoxia", "CAR", "donor_id"))
allVars13 <- effectPerVariable(dds13, c("hypoxia", "CAR", "donor_id"))

#Check dispersion
plotDispEsts(allVars7[[1]])

GetPVals <- function(ddsModels){
  AllRes <- lapply(ddsModels, FUN = function(x){
    res = results(x)
    res$gene <-  rownames(res)
    res$formula = res@elementMetadata$description[4] #Extract formula to keep track
    res = res[,c("gene","padj","formula")]
  })
}
pValsD7 = GetPVals(allVars7)
pValsD13 = GetPVals(allVars13)

namedPvals <- function(pValDat){
  for(i in 1:length(pValDat)){
    #Surely there's an easier way to automate finding the isolated variable given the formula
    formulaAsString = strsplit(pValDat[[i]]$formula[i], ": |vs|'~|'")[[1]][-c(1,2,4,5)]
    varCounts = table(strsplit(paste(formulaAsString, collapse = ""), " ")[[1]])
    varCounts = varCounts[which(varCounts == 1)]
    isolatedVar = names(varCounts)[2]
    pValDat[[i]]$variable = isolatedVar
  }
  do.call(rbind, pValDat)
}
pValsD7 = data.frame(namedPvals(pValsD7))
pValsD13 = data.frame(namedPvals(pValsD13))

pValsD7Wide <- pValsD7 %>% select(-("formula")) %>% pivot_wider(names_from = variable, values_from = padj)
pValsD13Wide <- pValsD13 %>% select(-("formula")) %>% pivot_wider(names_from = variable, values_from = padj)
#write.csv(pValsD13Wide, "./mystore/cartdata/PValsDay13All.csv")
