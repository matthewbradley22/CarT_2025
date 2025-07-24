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

T_cells <- LoadSeuratRds("./mystore/cartdata/data/tcell_deseqdat_singlets.rds")

T_cells_noCC <- LoadSeuratRds('mystore/cartdata/data/T_cells_noCC.rds')

#Machine learning approach to distinguish CD4 and 8 cells. Variable features from data 
#either in CAR_classifier.R or CD4CD8_classifier.R depending on which model

#car_classifier <- xgb.load("mystore/cartdata/data/car_classifier.model")
#CAR_genes <- readRDS("mystore/cartdata/data/CAR_classifier_variableGenes.rds")
CDClassifier <- xgb.load("mystore/cartdata/data/CD_classifier.model")
CD4CD8_genes <- readRDS("mystore/cartdata/data/CD4CD8Classifier_variableGenes.rds")

#Use claaifiers to assign CD4/CD8
ematCD <- FetchData(T_cells, vars = CD4CD8_genes)
CDTMat <- xgb.DMatrix(data = as.matrix(ematCD))
pred <- predict(CDClassifier, CDTMat)
preds_outcome <- ifelse(pred < 0.5, 0, 1)
preds_outcome <- ifelse(preds_outcome == 0, "CD4", "CD8")
T_cells$CD_pred <- preds_outcome


#We now label CAR+ cells using CAR expression, rather than machine learning, so commented out ML code

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


####Looking at PCA loadings, running Gene ontology on important genes and plotting results #### 
#Function to prepare data for topGO analysis
#Get top 5 pcs, get genes and loading value for given PC (which PC)
getTopGODat <- function(seuObj, designVars, whichPC = 1){
  dat <- getPseudoBulkObject(seuObj, designVars)
  dat_loadings <- getPCALoadings(dat, numGenes = 2000, numCom = 5)
  dat_genes <- prepTopGO(dat_loadings, PC = whichPC)
  dat_genes
}

CD4_day7 <- subset(T_cells, day == "D7" & CD_pred == 'CD4')

genes_cd4_day7 <- getTopGODat(CD4_day7, c("CAR","hypoxia", "Phase", "donor_id"), whichPC = 3)
genes_cd8_day7 <-  getTopGODat(CD8_day7, c("CAR","hypoxia", "Phase", "donor_id"), whichPC = 1)
genes_cd4_day13 <- getTopGODat(CD4_day13, c("CAR","hypoxia", "Phase", "donor_id"), whichPC = 2)
genes_cd8_day13 <- getTopGODat(CD8_day13, c("CAR","hypoxia", "Phase", "donor_id"), whichPC = 2)

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

####Example comparing one group (untransduced) to mean of other groups####
CD4_bulk <- subset(T_cells, CD_pred == 'CD4')
CD4_bulk <- getPseudoBulkObject(CD4_bulk, designVars = c("CAR","hypoxia", "Phase", "donor_id", "day"))
CD4_deseq <- DESeq(CD4_bulk)

resCD4 = deseqOneVsMean(CD4_deseq, "CARuntransduced", c("CARM1XX", "CARM28z", "CARMBBz"))
arrange(data.frame(resCD4), padj) %>% filter(log2FoldChange > 0 & padj < 1E-20) %>% 
  rownames()




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



