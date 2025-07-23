#Need to make file of functions to source from
library(Seurat)
library(tidyr)
library(dplyr)
library(xgboost)
library(ggplot2)
library(DESeq2)
library(umap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(topGO)
library(hgu95av2.db)
library(tidyr)

#Source fucntions
source("mystore/cartdata/scripts/CarT_project_functions.R")

#Make UMAPs plot properly
options(bitmapType="cairo") 

# T_cells <-  LoadSeuratRds("./mystore/cartdata/data/tcell_deseqdat_singlets.rds")
# T_cells$CAR_pred <- as.factor(T_cells$CAR_pred)
# T_cells$day = factor(T_cells$day, levels = c("D7", "D13"))

T_cells_noCC <- LoadSeuratRds('mystore/cartdata/data/T_cells_noCC.rds')
#Can skip ML if loading in ./mystore/cartdata/data/tcell_deseqdat_singlets.rds
#Machine learning approach to find untransduced t cells and distinguish CD4 and 8
#variable features from data either in CAR_classifier.R or CD4CD8_classifier.R 
#depending on which model

#I think default options for saving and loading an xgboost model can change it, 
# so may be best not to use this

# car_classifier <- xgb.load("mystore/cartdata/data/car_classifier.model")
# CAR_genes <- readRDS("mystore/cartdata/data/CAR_classifier_variableGenes.rds")
# CDClassifier_1 <- xgb.load("mystore/cartdata/data/CD_classifier.model")
# CD4CD8_genes <- readRDS("mystore/cartdata/data/CD4CD8Classifier_variableGenes.rds")
# 
# #Quick look at car vs hypoxia'
# T_cells_CAR <- subset(T_cells, CAR_pred == 1 & CAR != 'untransduced')
# table(T_cells_CAR$hypoxia,T_cells_CAR$CAR ) %>% as.data.frame() %>% filter(Var2 != 'untransduced') %>% ggplot(aes(x = Var1, y = Var2, fill = Freq))+
#    geom_tile()+
#   geom_text(aes(label=Freq), col = 'white')+
#   xlab('Hypoxia condition')+
#   ylab('CAR')
# 
# emat <- FetchData(T_cells, vars = CAR_genes) 
# carTMat <- xgb.DMatrix(data = as.matrix(emat))
# 
# pred <- predict(car_classifier, carTMat) #car_classifier from CAR_classifier.R
# pred_outcome <- ifelse(pred < 0.5, 0, 1)
# T_cells$CAR_pred <- pred_outcome
# T_cells[[]] %>% select(day, CAR, CAR_pred) %>% 
#   group_by(day, CAR) %>% summarize(percentExp = sum(CAR_pred)/n()) %>% 
#   ggplot(aes(x = factor(day, levels = c("D7", "D13")), y = percentExp, fill = CAR))+
#   geom_bar(stat = 'identity', position="dodge")+
#   xlab("Day")+
#   ylab("Proportion CAR expression")
# 
# carPreds = T_cells[[]] %>% group_by(CAR, CAR_pred) %>% summarise(amount = n()) %>%
#   pivot_wider(names_from = CAR_pred, values_from = amount)
# colnames(carPreds) <- c("CAR", "Pred no car", "Pred car")
# 
# ematCD_1 <- FetchData(deseqDat, vars = CD4CD8_genes) 
# CDTMat_1 <- xgb.DMatrix(data = as.matrix(ematCD_1))
# pred_1 <- predict(CDClassifier_1, CDTMat_1) 
# preds_outcome_1 <- ifelse(pred_1 < 0.5, 0, 1)
# preds_outcome_1 <- ifelse(preds_outcome == 0, "CD4", "CD8")
# T_cells$CD_pred_1 <- preds_outcome

#Split by day, CD, and CAR, (not phase?)
#Make a bunch of UMAPs by subset. At some point write function to shorten this
#Bystanders
bystanderCells <- subset(T_cells, CAR_pred == 0)

bystander7 <- subset(bystanderCells, day == "D7")
bystander13 <- subset(bystanderCells, day == "D13")

bystander7 <- prepUMAP(bystander7)
bystander13 <- prepUMAP(bystander13)

DimPlot(bystander7) + ggtitle("Day 7 bystanders")
DimPlot(bystander13) + ggtitle("Day 13 bystanders")

#CAR experssers
CARCells <- subset(T_cells, CAR_pred == 1)

CARCells7 <- subset(CARCells, day == "D7")
CARCells13 <- subset(CARCells, day == "D13")

CARCells7 <- prepUMAP(CARCells7)
CARCells13 <- prepUMAP(CARCells13)

DimPlot(CARCells7) + ggtitle("Day 7 CAR cells")
DimPlot(CARCells13) + ggtitle("Day 13 CAR cells")

#CD cells
CD4 <- subset(T_cells_noCC, CD_pred == "CD4")
CD8 <- subset(T_cells, CD_pred == "CD8")

CD4_day7 <- subset(CD4, day == "D7")
CD4_day13 <- subset(CD4, day == "D13")

CD8_day7 <- subset(CD8, day == "D7")
CD8_day13 <- subset(CD8, day == "D13")

CD4_day7 <- prepUMAP(CD4_day7)
CD4_day13 <- prepUMAP(CD4_day13)
CD8_day7 <- prepUMAP(CD8_day7)
CD8_day13 <- prepUMAP(CD8_day13)

DimPlot(CD4_day7) + ggtitle("CD4 day 7")
DimPlot(CD4_day13)+ ggtitle("CD4 day 13")
DimPlot(CD8_day7)+ ggtitle("CD8 day 7")
DimPlot(CD8_day13)+ ggtitle("C8 day 13")


#untransduced analysis split by CD
CD4 <- subset(T_cells, CD_pred == "CD4" & CAR_pred == 0)
CD8 <- subset(T_cells, CD_pred == "CD8" & CAR_pred == 0)



#Function from donor_analysis.R
#Not sure whether or not to include day variable here

#Function to prepare data for topGO analysis
#Get top 5 pcs, get genes and loading value for given PC (which PC)
getTopGODat <- function(seuObj, designVars, whichPC = 1){
  dat <- getPseudoBulkObject(seuObj, designVars)
  dat_loadings <- getPCALoadings(dat, numGenes = 2000, numCom = 5)
  dat_genes <- prepTopGO(dat_loadings, PC = whichPC)
  dat_genes
}

genes_cd4_day7 <- getTopGODat(CD4_day7, c("CAR","hypoxia", "Phase", "donor_id"), whichPC = 3)
genes_cd8_day7 <-  getTopGODat(CD8_day7, c("CAR","hypoxia", "Phase", "donor_id"), whichPC = 1)
genes_cd4_day13 <- getTopGODat(CD4_day13, c("CAR","hypoxia", "Phase", "donor_id"), whichPC = 2)
genes_cd8_day13 <- getTopGODat(CD8_day13, c("CAR","hypoxia", "Phase", "donor_id"), whichPC = 2)

getPCAGrid(CD4_day7, c("hypoxia", "CAR", "Phase", "donor_id"), rowNum = 2)
getPCAGrid(CD8_deseq, c("hypoxia", "CAR", "Phase", "donor_id"), rowNum = 2)
getPCAGrid(CD4_day7_bulk, c("hypoxia", "CAR", "Phase", "donor_id"), rowNum = 2)

#topGO wants a function that gives a score threshold for genes.
#We have genes with loading values and want to take top 100 absolute values

topLoadingGenes <- function (dat) 
{
  topGenes <- (head(dat[order(abs(dat), decreasing = TRUE)], n =100)) #Get top 100 genes by absolute value
  geneScore <- min(abs(topGenes))
  return(abs(dat) >= geneScore) 
}

GoDat <- new("topGOdata",
                     description = "Simple session", ontology = "BP",
                     allGenes = genes_cd4_day7, geneSel = topLoadingGenes,
                     nodeSize = 10,
                     annot = annFUN.org, mapping = "org.Hs.eg.db")

resultFisher <- runTest(GoDat, algorithm = "classic", statistic = "fisher")
allRes <- GenTable(GoDat, classicFisher = resultFisher, topNodes = 50)
showSigOfNodes(GoDat, score(resultFisher), firstSigNodes = 5, useInfo = 'all')


#Look at enrichment of given path
getFisherEnrichment(GoDat, "GO:0009653")

#Look at sig genes that are part of given pathway
path = "GO:0009653"
sigGenesInPath <- intersect(genesInTerm(GoDat, whichGO = path)[[1]],sigGenes(GoDat))
geneNames <- mapIds(org.Hs.eg.db, keys =  sigGenesInPath, 
                 keytype = "ENTREZID", column="SYMBOL")

#Get loading score
data.frame(genes_cd4_day7[which(names(genes_cd4_day7) %in% names(geneNames))]) 

#DESeq comparisons
CD4_deseq <- DESeq(CD4_deseq)
CD8_deseq <- DESeq(CD8_deseq)


resCD4 = deseqOneVsMean(CD4_deseq, "CARuntransduced", c("CARM1XX", "CARM28z", "CARMBBz"))
arrange(data.frame(resCD4), padj) %>% filter(log2FoldChange > 0 & padj < 1E-20) %>% 
  rownames()

resCD8 = deseqOneVsMean(CD8_deseq, "CARMBBz", c("CARM28z", "CARM1XX"))
arrange(data.frame(resCD8), padj)
write.csv(resCD8, "mystore/cartdata/data/CD8CARMBBzvsAll.csv")

#Psuedobulk UMAP

T_cells_UMAP <- getPseudoBulkUMAP(T_cells_noCC, c("CAR","hypoxia", "Phase", "donor_id", "CAR_pred",
                                             "CD_pred", "day"))

ggplot(umap_plot_df,aes(x = X1,y = X2, color = CD_pred)) +
  geom_point() +
  xlab("UMAP 1")+
  ylab("UMAP 2")+
  ggtitle("CD Type (ML prediction)")

#Pseudobulk UMAP by cd
CD4 <- subset(T_cells, CD_pred == "CD4")
CD8 <- subset(T_cells, CD_pred == "CD8")

CD4_deseq <- getPseudoBulkUMAP(CD4, 
                               vars = c("CAR","hypoxia", "Phase", "donor_id", "CAR_pred", "day"))

CD8_deseq <- getPseudoBulkUMAP(CD8, 
                               vars = c("CAR","hypoxia", "Phase", "donor_id", "CAR_pred", "day"))
ggplot(CD4_deseq,aes(x = X1,y = X2, color = CAR)) +
  geom_point() +
  xlab("UMAP 1")+
  ylab("UMAP 2")+
  ggtitle("CD8 CAR")

#Pseudobulk UMAP by day and CD

CD4Day7_umap <- getPseudoBulkUMAP(CD4_day7, c("CAR","hypoxia", "Phase", "donor_id", "CAR_pred"))
CD4Day13_umap <- getPseudoBulkUMAP(CD4_day13, c("CAR","hypoxia", "Phase", "donor_id", "CAR_pred"))
CD8Day7_umap <- getPseudoBulkUMAP(CD8_day7, c("CAR","hypoxia", "Phase", "donor_id", "CAR_pred"))
CD8Day13_umap <- getPseudoBulkUMAP(CD8_day13, c("CAR","hypoxia", "Phase", "donor_id", "CAR_pred"))

ggplot(CD8Day7_umap,aes(x = X1,y = X2, color = CAR)) +
  geom_point() +
  xlab("UMAP 1")+
  ylab("UMAP 2")+
  ggtitle("CD8 Day 13 Phase")

#Explore car trends in cd8 day 7
CD8_day7$CAR<-relevel(CD8_day7$CAR, ref="untransduced")
CD8_day7_bulk <- getPseudoBulkObject(CD8_day7, c("CAR","hypoxia", "Phase", "donor_id", "CAR_pred"))
CD8_day7_bulk <- DESeq(CD8_day7_bulk)


resultsNames(CD8_day7_bulk)
M1vsUntansduced <- lfcShrink(CD8_day7_bulk, coef="CAR_M1XX_vs_untransduced", type="apeglm")
MBvsM1 <- lfcShrink(CD8_day7_bulk, coef="CAR_MBBz_vs_M1XX", type="apeglm")
UntransducedVsM1 <- lfcShrink(CD8_day7_bulk, coef="CAR_untransduced_vs_M1XX", type="apeglm")

printOrderedDEGs <- function(degList, only.pos = TRUE){
  if(only.pos){
    degList[order(degList$padj),]  %>% as.data.frame() %>% filter(log2FoldChange>0) %>% head(n = 100) 
  }else{
    degList[order(degList$padj),]  %>% as.data.frame() %>% head(n = 100) 
  }
  
}
printOrderedDEGs(M1vsUntansduced, only.pos = T)
plotCounts(CD8_day7_bulk, gene="BIRC3", intgroup="CAR")

M2vsM1df <- M2vsM1[order(M2vsM1$padj),]  %>% as.data.frame() %>% filter(log2FoldChange<0) %>% head(n = 100) 
M2vsM1df$gene = rownames(M2vsM1df)
MBvsM1df <- MBvsM1[order(MBvsM1$padj),] %>%  as.data.frame() %>% filter(log2FoldChange<0) %>% head(n = 100)
MBvsM1df$gene = rownames(MBvsM1df)
UntransducedVsM1df <- UntransducedVsM1[order(UntransducedVsM1$padj),]  %>% head(n = 100) %>% as.data.frame()
M1vsM2MB <- intersect(rownames(M2vsM1df), rownames(MBvsM1df))

write.csv(M2vsM1df, "./mystore/cartdata/data/M2vsM1CAR.csv",  row.names=FALSE, quote = FALSE)
plotCounts(CD8Day7, gene="MT-ND5", intgroup="CAR")



getGeneOntResults <- function(geneList){
  geneList <- geneList[order(geneList$padj),]  %>% as.data.frame() %>% filter(log2FoldChange>0) %>% 
    dplyr::select(log2FoldChange, padj)
  
  topGenes <- function () 
  {
    topGenes <- head(geneList, n =100) #Get top 100 genes by absolute value
    return(rownames(topGenes)) 
  }
  
  
  dat <- new("topGOdata",
                  description = "Simple session", ontology = "BP",
                  allGenes = rownames(geneList), geneSel = topGenes,
                  nodeSize = 10,
                  annot = annFUN.org, mapping = "org.Hs.eg.db")
}


#Split donor analysis of each variable
CD4_day7$CAR<-relevel(CD4_day7$CAR, ref="untransduced")
CD4_day7_bulk <- getPseudoBulkObject(CD4_day7, c("CAR","hypoxia", "Phase", "donor_id", "CAR_pred"))
CD4_day7_bulk <- DESeq(CD4_day7_bulk)

getPCAGrid(CD4_day7_bulk, c("CAR","hypoxia", "Phase", "donor_id"), rowNum = 2)
resultsNames(CD4_day7_bulk)
M1xxVsM28z <- lfcShrink(CD4_day7_bulk, coef="CAR_M28z_vs_M1XX", type="apeglm")
M1xxVsM28z[order(M1xxVsM28z$padj),]  %>% as.data.frame() %>% filter(log2FoldChange>0) %>% head(n = 100) 

#Look at viruses
viruses <- c("virus-M28z", "virus-MBBz")
cd4 <- subset(T_cells, CD_pred == 'CD4')
cd8 <- subset(T_cells, CD_pred == 'CD8')

VlnPlot(cd4, features = viruses, assay = 'RNA', slot = 'data', group.by = 'CAR')
VlnPlot(cd8, features = viruses, assay = 'RNA', slot = 'data', group.by = 'CAR') 

cd8_viruses <- cd8[['RNA']]$data[viruses,] %>% as.data.frame() %>% t() %>% as.data.frame()
cd8_viruses$CAR <- unname(cd8$CAR)
CARs <- unique(cd8$CAR)
virusCounts <- lapply(CARs, FUN = function(x){
  virusCounts <- cd8_viruses %>% filter(CAR == x) %>% select(!CAR) %>% colMeans(. != 0)
  virusDat <- data.frame(CAR = x, M28z = virusCounts[1], MBBz = virusCounts[2])
  virusDat
})

virusCounts <- do.call(rbind, virusCounts)
rownames(virusCounts) <- NULL
virusCounts <- virusCounts %>% pivot_longer(!CAR, names_to = 'Virus', values_to = 'PercentExpressing')
ggplot(virusCounts, aes(x = CAR, y = PercentExpressing, fill = Virus))+
  geom_bar(stat = 'identity', position="dodge") + ggtitle('Virus Expression by CAR group (CD8 cells)')+
  ylab('Percent cells expressing virus')

#Separate by hypoxia and day
HHD7 <- subset(T_cells, hypoxia == 'HH' & day == 'D7')
NND7 <- subset(T_cells, hypoxia == 'NN' & day == 'D7')
HHD13 <- subset(T_cells, hypoxia == 'HH' & day == 'D13')
NHD13 <- subset(T_cells, hypoxia == 'NH' & day == 'D13')
NND13 <- subset(T_cells, hypoxia == 'NN' & day == 'D13')

HHD7 <- prepUMAP(HHD7, regressCC = 'CC')
ElbowPlot(HHD7)
HHD7 <- finishUMAP(HHD7, dimNeighbors = 15)
DimPlot(HHD7) + ggtitle('Day 7 HH')
HHD7markers <- FindAllMarkers(HHD7, only.pos = TRUE)
HHD7markers <- filter(HHD7markers, p_val_adj < 0.01)
topClusterGenesHHD7 <- HHD7markers %>% group_by(cluster) %>% slice_max(order_by = p_val_adj, n = 30) %>% .$gene

NND7 <- prepUMAP(NND7, regressCC = 'CC')
ElbowPlot(NND7)
NND7 <- finishUMAP(NND7, dimNeighbors = 15)
DimPlot(NND7) + ggtitle('Day 7 NN')
NND7markers <- FindAllMarkers(NND7, only.pos = TRUE)
NND7markers <- filter(NND7markers, p_val_adj < 0.01)
topClusterGenesNND7 <- NND7markers %>% group_by(cluster) %>% slice_max(order_by = p_val_adj, n = 30) %>% .$gene

HHD13 <- prepUMAP(HHD13, regressCC = 'CC')
ElbowPlot(HHD13)
HHD13 <- finishUMAP(HHD13, dimNeighbors = 15)
DimPlot(HHD13) + ggtitle('Day 13 HH')

NHD13 <- prepUMAP(NHD13, regressCC = 'CC')
ElbowPlot(NHD13)
NHD13 <- finishUMAP(NHD13, dimNeighbors = 15)
DimPlot(NHD13) + ggtitle('Day 13 NH')

NND13 <- prepUMAP(NND13, regressCC = 'CC')
ElbowPlot(NND13)
NND13 <- finishUMAP(NND13, dimNeighbors = 15)
DimPlot(NND13) + ggtitle('Day 13 NN')

