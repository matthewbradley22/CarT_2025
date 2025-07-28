#Need to make file of functions to source from
library(tidyr)
library(dplyr)
library(ggplot2)
library(DESeq2)
library(umap)

#Source functions and some packages (Seurat etc...)
source("mystore/cartdata/scripts/CarT_project_functions.R")

#Plot properly
options(bitmapType="cairo") 

#Load in data 
T_cells <-  LoadSeuratRds("./mystore/cartdata/data/tcell_deseqdat_singlets.rds")
T_cells$day = factor(T_cells$day, levels = c("D7", "D13"))
T_cells$CAR_pred = factor(T_cells$CAR_pred)


#Look at top 50 genes
#Load in gene set from paper: 
importantGenes <- read.csv("mystore/cartdata/data/importantGenesWithDetails.csv", header = F)
colnames(importantGenes) <- c("gene", "CD", "type")
importantGenes <- importantGenes[with(importantGenes, order(CD, type)), ]

#Heatmap of important genes
t_cell_bulk <- getPseudoBulkObject(T_cells, c("CAR","hypoxia", "Phase", "donor_id", "CAR_pred", "day", "CD_pred"), return.Seurat = TRUE)
t_cell_bulk[[]] <- t_cell_bulk[[]] %>% mutate("heatmapGroups" = paste(CAR, hypoxia, day, sep = "_"))
t_cell_bulk$heatmapGroups <- gsub('untransduced', 'un', t_cell_bulk$heatmapGroups)
t_cell_bulk$heatmapGroups <- gsub('MBBz', 'MB', t_cell_bulk$heatmapGroups)
t_cell_bulk$heatmapGroups <- gsub('M28z', 'M2', t_cell_bulk$heatmapGroups)
t_cell_bulk$heatmapGroups <- gsub('M1XX', 'M1', t_cell_bulk$heatmapGroups)

cd4_bulk <- subset(t_cell_bulk, CD_pred == "CD4")
cd8_bulk <- subset(t_cell_bulk, CD_pred == "CD8")


#Pseudobulk heatmaps in Fig. 2. Saved dimensions 17 x 20. I just preview in browser and download these, rather than save to hpc2n
DoHeatmap(cd4_bulk, features = importantGenes$gene, group.by = c("heatmapGroups"), size = 10)+
  guides(colour=FALSE) +
  theme(axis.text.y = element_text(size = 25))

DoHeatmap(cd8_bulk, features = importantGenes$gene, group.by = c("heatmapGroups"), size = 10)+
  guides(colour=FALSE) +
  theme(axis.text.y = element_text(size = 25))

avg <- AverageExpression(t_cell_bulk, group.by = "CAR", return.seurat = TRUE)
DoHeatmap(avg, features = importantGenes$gene, draw.lines = F)+
  guides(colour=FALSE)


DoHeatmap(T_cells, features = importantGenes$gene, draw.lines = F, group.by = 'day')+
  guides(colour=FALSE)


#Create violin plot of each group of genes
importantGenes <- importantGenes %>% separate_longer_delim(type, delim = "/")
CD4Genes <- subset(importantGenes, CD == "CD4" | CD == "CD4 and CD8")
CD8Genes <- subset(importantGenes, CD == "CD8" | CD == "CD4 and CD8")

CD4Inflammation <- list(subset(CD4Genes, type == "Inflammation")$gene)
CD4Cytoxicity <- list(subset(CD4Genes, type == "Cytotoxicity")$gene)
CD4Exhaustion <- list(subset(CD4Genes, type == "Exhaustion")$gene)
CD4Inhibition <- list(subset(CD4Genes, type == "Inhibition")$gene)
CD4List = list(CD4Inflammation, CD4Cytoxicity, CD4Exhaustion, CD4Inhibition)
names(CD4List) <- list("CD4Inflammation", "CD4Cytoxicity", "CD4Exhaustion", "CD4Inhibition")

CD8Inflammation <- list(subset(CD8Genes, type == "Inflammation")$gene)
CD8Cytoxicity <- list(subset(CD8Genes, type == "Cytotoxicity")$gene)
CD8Exhaustion <- list(subset(CD8Genes, type == "Exhaustion")$gene)
CD8Inhibition <- list(subset(CD8Genes, type == "Inhibition")$gene)
CD8List = list(CD8Inflammation, CD8Cytoxicity, CD8Exhaustion, CD8Inhibition)
names(CD8List) <- list("CD8Inflammation", "CD8Cytoxicity", "CD8Exhaustion", "CD8Inhibition")

CD4 <- subset(T_cells, CD_pred == "CD4")
CD8 <- subset(T_cells, CD_pred == "CD8")

for(i in 1:length(CD4List)){
  CD4 <- AddModuleScore(CD4, features = CD4List[[i]], name = names(CD4List[i]))
}

for(i in 1:length(CD8List)){
  CD8 <- AddModuleScore(CD8, features = CD8List[[i]], name = names(CD8List[i]))
}

getVlnGrid <- function(dat, datList, rowNum, groupVar){
  plotList = list()
  for(i in 1:length(datList)){
    plotList[[i]] = VlnPlot(dat, features = paste0(names(datList)[i], "1"), group.by = groupVar,
                            pt.size = 0) + ggtitle(names(datList)[i])
  }
  grid.arrange(grobs = plotList, nrow = rowNum)
}

getVlnGrid(CD4, CD4List, 2, "CAR")


#Pseudobulk expression and looking at DEGs by pathway
CD4_deseq <- getPseudoBulkObject(CD4, c("CAR","hypoxia", "Phase", "donor_id", "CAR_pred", "day"))
CD8_deseq <- getPseudoBulkObject(CD8, c("CAR","hypoxia", "Phase", "donor_id", "CAR_pred", "day"))
CD4_deseq <- DESeq(CD4_deseq)
CD8_deseq <- DESeq(CD8_deseq)

CD4Day7 <- subset(CD4, day == "D7")
CD4Day13 <- subset(CD4, day == "D13")
CD8Day7 <- subset(CD8, day == "D7")
CD8Day13 <- subset(CD8, day == "D13")

CD4Day7_deseq <- getPseudoBulkObject(CD4Day7, c("CAR","hypoxia", "Phase", "donor_id"))
CD4Day13_deseq <- getPseudoBulkObject(CD4Day13, c("CAR","hypoxia", "Phase", "donor_id"))
CD4Day7_deseq <- DESeq(CD4Day7_deseq)
CD4Day13_deseq <- DESeq(CD4Day13_deseq)

CD8Day7_deseq <- getPseudoBulkObject(CD8Day7, c("CAR","hypoxia", "Phase", "donor_id"))
CD8Day13_deseq <- getPseudoBulkObject(CD8Day13, c("CAR","hypoxia", "Phase", "donor_id"))
CD8Day7_deseq <- DESeq(CD8Day7_deseq)
CD8Day13_deseq <- DESeq(CD8Day13_deseq)

cd4Day7Genes <- getTopGenes(CD4Day7_deseq, geneList = CD4Genes)
CD4Day13Genes <- getTopGenes(CD4Day13_deseq, geneList = CD4Genes)
cd8Day7Genes <- getTopGenes(CD8Day7_deseq, geneList = CD8Genes)
CD8Day13Genes <- getTopGenes(CD8Day13_deseq, geneList = CD8Genes)

plotTopGenes(cd8Day7Genes) + ggtitle("CD8 Day 7 DEGs")+
  theme(axis.title.y=element_text(size=16))








