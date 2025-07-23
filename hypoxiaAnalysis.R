#Source functions and packages
source("mystore/cartdata/scripts/CarT_project_functions.R")
library(gprofiler2)
library(tibble)
library(EnhancedVolcano)
library(scGSVA)
library(VennDiagram)
library(UpSetR)
library(SeuratDisk)

#Make UMAPs plot properly
options(bitmapType="cairo")

#Run cellCycleRegressed.R script to get T_cells_noCC or read in:
T_cells_noCC <- LoadSeuratRds('mystore/cartdata/data/T_cells_noCC.rds')
DimPlot(T_cells_noCC, group.by = "Phase")

#Summary look
T_cells_noCC$day = factor(T_cells_noCC$day, levels = c('D7', 'D13'))
T_cells_noCC[[]] %>% group_by(day, CD_pred, hypoxia) %>% summarise(totalHypoxia = n()) %>% 
  mutate(day_cd = factor(paste(day, CD_pred, sep = '_'), levels = c('D7_CD4', 'D7_CD8', 'D13_CD4', 'D13_CD8'))) %>% 
  ggplot(aes (x = day_cd, y = totalHypoxia, fill = hypoxia))+
  geom_bar(stat = 'identity', position="dodge") +
  xlab('Group')+
  ylab('Cell Count')+
  ggtitle('Hypoxia by day and CD')

#Label car t cells by car expression
T_cells_noCC[['RNA']]$counts[c('virus-M28z', 'virus-MBBz'),]
M28zCells <- WhichCells(T_cells_noCC, expression = `virus-M28z` > 0)
MBBzCells <- WhichCells(T_cells_noCC, expression = `virus-MBBz` > 0)
carPosCells <- union(M28zCells, MBBzCells)
T_cells_noCC$carExpression <- ifelse(colnames(T_cells_noCC) %in% carPosCells, "carPos", "carNeg")

#Check that labeling worked and compare to ML approach
sum(colSums(subset(T_cells_noCC, carExpression == 'carNeg')[['RNA']]$data[c('virus-M28z', 'virus-MBBz'),]) == 0)
table(T_cells_noCC$CAR_pred, T_cells_noCC$carExpression)


#Split into 4 groups
day7_CD4 <- subset(T_cells_noCC, day == 'D7' & CD_pred == 'CD4' & CAR_pred == 1 & CAR != 'untransduced')
day7_CD4 <- prepUMAP(day7_CD4)
ElbowPlot(day7_CD4)
day7_CD4 <- finishUMAP(day7_CD4, dimNeighbors = 15)
DimPlot(day7_CD4, group.by = 'hypoxia') + ggtitle('Day 7 CD4 UMAP')
DimHeatmap(day7_CD4, dims = 1:6, cells = 500, balanced = TRUE)

day13_CD4 <- subset(T_cells_noCC, day == 'D13' & CD_pred == 'CD4' & CAR_pred == 1 & CAR != 'untransduced')
day13_CD4 <- prepUMAP(day13_CD4)
ElbowPlot(day13_CD4)
day13_CD4 <- finishUMAP(day13_CD4, dimNeighbors = 20)
DimPlot(day13_CD4, group.by = 'hypoxia') + ggtitle('Day 13 CD4 UMAP')

day7_CD8 <- subset(T_cells_noCC, day == 'D7' & CD_pred == 'CD8' & CAR_pred == 1 & CAR != 'untransduced')
day7_CD8 <- prepUMAP(day7_CD8)
ElbowPlot(day7_CD8)
day7_CD8 <- finishUMAP(day7_CD8, dimNeighbors = 15)
DimPlot(day7_CD8, group.by = 'hypoxia') + ggtitle('Day 7 CD8 UMAP')

day13_CD8 <- subset(T_cells_noCC, day == 'D13' & CD_pred == 'CD8' & CAR_pred == 1 & CAR != 'untransduced')
day13_CD8 <- prepUMAP(day13_CD8)
ElbowPlot(day13_CD8)
day13_CD8 <- finishUMAP(day13_CD8, dimNeighbors = 20)
DimPlot(day13_CD8, group.by = 'hypoxia') + ggtitle('Day 13 CD8 UMAP')

#Pseudobulk objects
day7_CD4_bulk <- getPseudoBulkObject(day7_CD4, c("hypoxia", "donor_id", 'CAR'))
day13_CD4$hypoxia <- relevel(day13_CD4$hypoxia, 'NH') #Only if want to compare to NH
day13_CD4_bulk <- getPseudoBulkObject(day13_CD4, c("hypoxia", "donor_id", 'CAR'))
day7_CD8_bulk <- getPseudoBulkObject(day7_CD8, c("hypoxia", "donor_id", 'CAR'))
day13_CD8$hypoxia <- relevel(day13_CD8$hypoxia, 'NH') #Only if want to compare to NH
day13_CD8_bulk <- getPseudoBulkObject(day13_CD8, c("hypoxia", "donor_id", 'CAR'))

#DeSeq on objects
day7_CD4_bulk <- DESeq(day7_CD4_bulk)
day13_CD4_bulk <- DESeq(day13_CD4_bulk)
day7_CD8_bulk <- DESeq(day7_CD8_bulk)
day13_CD8_bulk <- DESeq(day13_CD8_bulk)

#Make hypoxia comps
res_cd4_7 <- lfcShrink(day7_CD4_bulk, coef="hypoxia_NN_vs_HH", type="apeglm")
cd4_7_NNHH_D <- nrow(res_cd4_7 %>% as.data.frame() %>% filter(padj < 0.01)) #2279
res_cd4_7[order(res_cd4_7$padj),]  %>% as.data.frame() %>% head(n = 1000) %>% 
  write.csv(., './mystore/cartdata/data/HypoxiaDEGLists/CD4_day7_NNvsHH.csv')

EnhancedVolcano(res_cd4_7,
                lab = rownames(res_cd4_7),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'CD4 Day 7 NN vs HH')

res_cd4_13 <- lfcShrink(day13_CD4_bulk, coef="hypoxia_NN_vs_HH", type="apeglm")
cd4_13_NNHH <- nrow(res_cd4_13 %>% as.data.frame() %>% filter(padj < 0.01)) #1321
res_cd4_13[order(res_cd4_13$padj),]  %>% as.data.frame() %>% head(n = 1000) %>% 
  write.csv(., './mystore/cartdata/data/HypoxiaDEGLists/CD4_day13_NNvsHH.csv')

EnhancedVolcano(res_cd4_13,
                lab = rownames(res_cd4_13),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'CD4 Day 13 NN vs HH')

res_cd8_7 <- lfcShrink(day7_CD8_bulk, coef="hypoxia_NN_vs_HH", type="apeglm")
cd8_7_NNHH <- nrow(res_cd8_7 %>% as.data.frame() %>% filter(padj < 0.01)) #2326
res_cd8_7[order(res_cd8_7$padj),]  %>% as.data.frame() %>% head(n = 1000) %>% 
  write.csv(., './mystore/cartdata/data/HypoxiaDEGLists/CD8_day7_NNvsHH.csv')

EnhancedVolcano(res_cd8_7,
                lab = rownames(res_cd8_7),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'CD8 Day 7 NN vs HH')

res_cd8_13 <- lfcShrink(day13_CD8_bulk, coef="hypoxia_NN_vs_HH", type="apeglm")
cd8_13_NNHH <- nrow(res_cd8_13 %>% as.data.frame() %>% filter(padj < 0.01)) #1393
res_cd8_13[order(res_cd8_13$padj),]  %>% as.data.frame() %>% head(n = 1000) %>% 
  write.csv(., './mystore/cartdata/data/HypoxiaDEGLists/CD8_day13_NNvsHH.csv')

EnhancedVolcano(res_cd8_13,
                lab = rownames(res_cd8_13),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'CD8 Day 13 NN vs HH')

#NH Hypoxia comparisons
res_cd4_13_HHNH <- lfcShrink(day13_CD4_bulk, coef="hypoxia_HH_vs_NH", type="apeglm")
cd4_13_HHNH <- nrow(res_cd4_13_HHNH %>% as.data.frame() %>% filter(padj < 0.01)) #353
res_cd4_13_HHNH[order(res_cd4_13_HHNH$padj),]  %>% as.data.frame() %>% head(n = 1000) %>% 
  write.csv(., './mystore/cartdata/data/HypoxiaDEGLists/CD4_day13_NHvsHH.csv')

EnhancedVolcano(res_cd4_13_HHNH,
                lab = rownames(res_cd4_13_HHNH),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'CD4 Day 13 NH vs HH')

res_cd4_13_NNNH <- lfcShrink(day13_CD4_bulk, coef="hypoxia_NN_vs_NH", type="apeglm")
cd4_13_NNNH <- nrow(res_cd4_13_NNNH %>% as.data.frame() %>% filter(padj < 0.01)) #583
res_cd4_13_NNNH[order(res_cd4_13_NNNH$padj),]  %>% as.data.frame() %>% head(n = 1000) %>% 
  write.csv(., './mystore/cartdata/data/HypoxiaDEGLists/CD4_day13_NHvsNN.csv')

EnhancedVolcano(res_cd4_13_NNNH,
                lab = rownames(res_cd4_13_NNNH),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'CD8 Day 13 NH vs HH')

res_cd8_13_HHNH <- lfcShrink(day13_CD8_bulk, coef="hypoxia_HH_vs_NH", type="apeglm")
cd8_13_HHNH <- nrow(res_cd8_13_HHNH %>% as.data.frame() %>% filter(padj < 0.01)) #482
res_cd8_13_HHNH[order(res_cd8_13_HHNH$padj),]  %>% as.data.frame() %>% head(n = 1000) %>% 
  write.csv(., './mystore/cartdata/data/HypoxiaDEGLists/CD8_day13_NHvsHH.csv')

EnhancedVolcano(res_cd8_13_HHNH,
                lab = rownames(res_cd8_13_HHNH),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'CD8 Day 13 NH vs HH')

res_cd8_13_NNNH <- lfcShrink(day13_CD8_bulk, coef="hypoxia_NN_vs_NH", type="apeglm")
cd8_13_NNNH <- nrow(res_cd8_13_NNNH %>% as.data.frame() %>% filter(padj < 0.01)) #625
res_cd8_13_NNNH[order(res_cd8_13_NNNH$padj),]  %>% as.data.frame() %>% head(n = 1000) %>% 
  write.csv(., './mystore/cartdata/data/HypoxiaDEGLists/CD8_day13_NHvsNN.csv')

EnhancedVolcano(res_cd8_13_NNNH,
                lab = rownames(res_cd8_13_NNNH),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'CD4 Day 13 NH vs HH')

#Compare number of DEGs per comparison
degNums = c(cd4_7_NNHH_D, cd8_7_NNHH, cd4_13_NNHH, cd8_13_NNHH, cd4_13_HHNH, cd4_13_NNNH, cd8_13_HHNH, cd8_13_NNNH)
compLevels <- factor(c('cd4_7_NNHH', 'cd8_7_NNHH','cd4_13_NNHH', 'cd8_13_NNHH', 'cd4_13_HHNH', 'cd4_13_NNNH', 'cd8_13_HHNH', 'cd8_13_NNNH'), 
                     levels = c('cd4_7_NNHH', 'cd8_7_NNHH','cd4_13_NNHH', 'cd8_13_NNHH', 'cd4_13_HHNH', 'cd4_13_NNNH', 'cd8_13_HHNH',
                                'cd8_13_NNNH'))
days = factor(c('7', '7', '13', '13', '13', '13', '13', '13'), levels = c('7', '13'))
numDEGS <- data.frame(comparison = compLevels, numDEGS = degNums, day = days)
ggplot(numDEGS, aes(x = comparison, y = numDEGS, fill = days))+
  geom_bar(stat = 'identity')+
  ggtitle('Number of DEGs per group')+
  ylab('Number of DEGs')

#Look at combined genes
res_cd4_7_subset <- res_cd4_7[order(res_cd4_7$padj),]  %>% as.data.frame() %>% head(n = 1000) 
res_cd8_7_subset <- res_cd8_7[order(res_cd8_7$padj),]  %>% as.data.frame() %>% head(n = 1000) 
res_cd4_13_subset <- res_cd4_13[order(res_cd4_13$padj),]  %>% as.data.frame() %>% head(n = 1000) 
res_cd8_13_subset <- res_cd8_13[order(res_cd8_13$padj),]  %>% as.data.frame() %>% head(n = 1000) 

genes <- list(rownames(res_cd4_7_subset), rownames(res_cd8_7_subset), rownames(res_cd4_13_subset), rownames(res_cd8_13_subset))

commonGenes <- Reduce(intersect, genes)
noquote(commonGenes)

#Look at gene scores between hypoxia groups
hsko<-buildAnnot(species="human",keytype="SYMBOL",anntype="KEGG")
res<-scgsva(day7_CD4, hsko)
head(findPathway(res,group = "hypoxia"))
vlnPlot(res,features="Chemokine.signaling.pathway",group_by="hypoxia")
ridgePlot(res,features="Chemokine.signaling.pathway",group_by="hypoxia")
Heatmap(res, group_by="hypoxia")

#Load in gene scores sets from descriptive figures script
scores = list(MemoryScore, proliferationScore, cytotoxicScore, restingScore, activeScore, proInflammatory, lateEffectorDiff)
names(scores) = c('MemoryScore', 'proliferationScore', 'cytotoxicScore', 'restingScore', 'activeScore', 'proInflammatory', 'lateEffectorDiff')
# for(i in 1:length(scores)){
#   day7_CD4 <- AddModuleScore(day7_CD4, scores[[i]], name = names(scores)[i])
#   day13_CD4 <- AddModuleScore(day13_CD4, scores[[i]], name = names(scores)[i])
#   day7_CD8 <- AddModuleScore(day7_CD8, scores[[i]], name = names(scores)[i])
#   day13_CD8 <- AddModuleScore(day13_CD8, scores[[i]], name = names(scores)[i])
# }
# 
# day7_CD4$hypoxia <- factor(day7_CD4$hypoxia, levels = c('NN', 'HH'))
# 
# VlnPlot(day13_CD4, features = 'restingScore1', group.by = 'hypoxia')+ 
#   theme(legend.position = 'none',
#         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#   ggtitle('CD4 Day 7 proliferation score')

DotPlot(object = day13_CD8, features = unique(unlist(scores)), group.by = 'hypoxia')+
   theme(axis.text.x = element_text(angle = 75, vjust = 0.7))+
  ggtitle('CD8 Day 13, notable genes')


#CAR proportions in hypoxia groups
table(day7_CD4$hypoxia, day7_CD4$CAR) %>% as.data.frame() %>% 
  ggplot(aes(x = Var1, y = Var2, fill = Freq))+
  geom_tile()+
  geom_text(aes(label=Freq), col = 'white')+
  ggtitle('Hypoxia CAR table. CD4 Day 7')

table(day13_CD4$hypoxia, day13_CD4$CAR) %>% as.data.frame() %>% 
  ggplot(aes(x = Var1, y = Var2, fill = Freq))+
  geom_tile()+
  geom_text(aes(label=Freq), col = 'white')+
  ggtitle('Hypoxia CAR table. CD4 Day 13')

table(day7_CD8$hypoxia, day7_CD8$CAR) %>% as.data.frame() %>% 
  ggplot(aes(x = Var1, y = Var2, fill = Freq))+
  geom_tile()+
  geom_text(aes(label=Freq), col = 'white')+
  ggtitle('Hypoxia CAR table. CD8 Day 7')

table(day13_CD8$hypoxia, day13_CD8$CAR) %>% as.data.frame() %>% 
  ggplot(aes(x = Var1, y = Var2, fill = Freq))+
  geom_tile()+
  geom_text(aes(label=Freq), col = 'white')+
  ggtitle('Hypoxia CAR table. CD8 Day 13')

#Look at untransduced cells
#Split into 4 groups
day7_CD4_noCAR <- subset(T_cells_noCC, day == 'D7' & CD_pred == 'CD4' & CAR_pred == 0)
day7_CD4_noCAR <- prepUMAP(day7_CD4_noCAR)
ElbowPlot(day7_CD4_noCAR)
day7_CD4_noCAR <- finishUMAP(day7_CD4_noCAR, dimNeighbors = 15)
DimPlot(day7_CD4_noCAR, group.by = 'hypoxia') + ggtitle('Day 7 CD4 UMAP')
DimHeatmap(day7_CD4_noCAR, dims = 1:6, cells = 500, balanced = TRUE)


#Looking at gene sets from GO results
cytoSignalling <- c('CHUK','DHX9','HMGB1','HSPA5','HSPA8','HSP90AA1','KPNA2','MT2A','PSMC5','RANBP2','UBB',
                    'USP14','HNRNPDL','SAMHD1') #Seem highly expressed at day 7 and in NN'
cytoSignallingInImmuneSystem <- c('CSF2RB','S1PR1','HMGB1','OAS2','EIF2AK2','SOS1','YWHAZ','H3C11','IFITM1','XAF1','CCR2','ITGA4',
                                  'H2BC17','H4C3','H2AC12','RIPOR2','KLF2','RNF125')
responseToHypoxia <- c('BNIP3','BNIP3L','HSPG2','LTA','MT3','PGK1','NOL3','NDRG1','FAM162A','HILPDA','ANGPTL4','EGLN3','PLIN2')

aerobicGlycolysis <- c('ALDOA','ENO1','GAPDH','GPI','LDHA','PGK1','TPI1','ALDOC',
                       'ENO2','STAT3','PLIN2','BNIP3','ITGB2','EGLN3','GLUL','SAMD4A','GARS1')


DotPlot(T_cells_noCC, features = cytoSignalling, group.by = 'CAR')

DimPlot(T_cells_noCC, label = TRUE, group.by = 'day')

DoHeatmap(day7_CD4, features = MemoryScore[[1]], slot = 'data') 
DotPlot(day7_CD4, features = cytotoxicScore[[1]], group.by = 'hypoxia') 

DotPlot(day7_CD4, features = 'SELL', group.by = 'hypoxia') 
VlnPlot(day7_CD4, features = 'SELL', group.by = 'hypoxia') 

CD4 <- subset(T_cells_noCC, CD_pred == 'CD4' & CAR != 'untransduced' & CAR_pred == 1)
CD8 <- subset(T_cells_noCC, CD_pred == 'CD8' & CAR != 'untransduced' & CAR_pred == 1)
CarPositives <- subset(T_cells_noCC, CAR != 'untransduced' & CAR_pred == 1)
DotPlot(CD4 ,cytoSignallingInImmuneSystem, group.by = 'hypoxia', assay = 'RNA') +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle('Day 13 CAR+ CD4 Cytokine Signalling genes')

DotPlot(CarPositives, responseToHypoxia, group.by = 'hypoxia', assay = 'RNA') +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle('CAR+ Response to hypoxia genes')

DotPlot(CarPositives, cytoSignallingInImmuneSystem, group.by = 'hypoxia', assay = 'RNA') +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle('CAR+ Immune Cytokine genes')

DotPlot(CarPositives, aerobicGlycolysis, group.by = 'hypoxia', assay = 'RNA') +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle('CAR+ Immune Cytokine genes')
# plinGenes <- rownames(day7_CD4[['RNA']]$data[grep('^PLIN', rownames(day7_CD4[['RNA']]$data)),])
# DotPlot(day7_CD4, plinGenes, group.by = 'hypoxia')

#Split above groups further by CAR
day7_CD4_M1XX <- subset(day7_CD4, CAR == 'M1XX')
day7_CD4_M1XX <- prepUMAP(day7_CD4_M1XX)
ElbowPlot(day7_CD4_M1XX)
day7_CD4_M1XX <- finishUMAP(day7_CD4_M1XX, dimNeighbors = 15)
DimPlot(day7_CD4_M1XX, group.by = 'hypoxia') + ggtitle('Day 7 CD4 M1XX UMAP')


day7_CD4_M28z <- subset(day7_CD4, CAR == 'M28z')
day7_CD4_M28z <- prepUMAP(day7_CD4_M28z)
ElbowPlot(day7_CD4_M28z)
day7_CD4_M28z <- finishUMAP(day7_CD4_M28z, dimNeighbors = 15)
DimPlot(day7_CD4_M28z, group.by = 'hypoxia') + ggtitle('Day 7 M28z CD4 UMAP')

day7_CD4_MBBz <- subset(day7_CD4, CAR == 'MBBz')
day7_CD4_MBBz <- prepUMAP(day7_CD4_MBBz)
ElbowPlot(day7_CD4_MBBz)
day7_CD4_MBBz <- finishUMAP(day7_CD4_MBBz, dimNeighbors = 15)
DimPlot(day7_CD4_MBBz, group.by = 'hypoxia') + ggtitle('Day 7 MBBz CD4 UMAP')

day7_CD4_M1XX_bulk <- getPseudoBulkObject(day7_CD4_M1XX, c("hypoxia", "donor_id"))
day7_CD4_M28z_bulk <- getPseudoBulkObject(day7_CD4_M28z, c("hypoxia", "donor_id"))
day7_CD4_MBBz_bulk <- getPseudoBulkObject(day7_CD4_MBBz, c("hypoxia", "donor_id"))

day7_CD4_M1XX_bulk <- DESeq(day7_CD4_M1XX_bulk)
day7_CD4_M28z_bulk <- DESeq(day7_CD4_M28z_bulk)
day7_CD4_MBBz_bulk <- DESeq(day7_CD4_MBBz_bulk)

res_day7_CD4_M1XX_bulk <- lfcShrink(day7_CD4_M1XX_bulk, coef="hypoxia_NN_vs_HH", type="apeglm")
res_day7_CD4_M1XX_bulk[order(res_day7_CD4_M1XX_bulk$padj),]  %>% as.data.frame() %>% head(n = 1000)# %>% 
 # write.csv(., './mystore/cartdata/data/HypoxiaDEGLists/CD8_day7_NNvsHH.csv')
cd4_7_M1XXD <- rownames(res_day7_CD4_M1XX_bulk %>% as.data.frame() %>% filter(padj < 0.01))

EnhancedVolcano(res_day7_CD4_M1XX_bulk,
                lab = rownames(res_day7_CD4_M1XX_bulk),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'CD4 Day 7 M1XX NN vs HH')

res_day7_CD4_M28z_bulk <- lfcShrink(day7_CD4_M28z_bulk, coef="hypoxia_NN_vs_HH", type="apeglm")
res_day7_CD4_M28z_bulk[order(res_day7_CD4_M28z_bulk$padj),]  %>% as.data.frame() %>% head(n = 1000)# %>% 
# write.csv(., './mystore/cartdata/data/HypoxiaDEGLists/CD8_day7_NNvsHH.csv')
cd4_7_M28zD <- rownames(res_day7_CD4_M28z_bulk %>% as.data.frame() %>% filter(padj < 0.01))
EnhancedVolcano(res_day7_CD4_M28z_bulk,
                lab = rownames(res_day7_CD4_M28z_bulk),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'CD4 Day 7 M28z NN vs HH')

res_day7_CD4_MBBz_bulk <- lfcShrink(day7_CD4_MBBz_bulk, coef="hypoxia_NN_vs_HH", type="apeglm")
res_day7_CD4_MBBz_bulk[order(res_day7_CD4_MBBz_bulk$padj),]  %>% as.data.frame() %>% head(n = 1000)# %>% 
# write.csv(., './mystore/cartdata/data/HypoxiaDEGLists/CD8_day7_NNvsHH.csv')
cd4_7_MBBzD <- rownames(res_day7_CD4_MBBz_bulk %>% as.data.frame() %>% filter(padj < 0.01))
EnhancedVolcano(res_day7_CD4_MBBz_bulk,
                lab = rownames(res_day7_CD4_MBBz_bulk),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'CD4 Day 7 MBBz NN vs HH')

commonGenes <- Reduce(intersect, list(cd4_7_M1XXD, cd4_7_M28zD, cd4_7_MBBzD))

geneInput <- c(M1XX = length(cd4_7_M1XXD), M28z = length(cd4_7_M28zD), MBBz = length(cd4_7_MBBzD), 
               'M1XX&M28z' = 42, 'M1XX&MBBz' = 94, 'M28z&MBBz' = 44, 'M1XX&M28z&MBBz' = 25
               )
upset(fromExpression(geneInput), order.by = "freq")


#Look at responder vs nonresponder genes from Saren single cell rna analysis reveals cell-intrinsic functions...
responder <- c("CCL3", 'LAG3', 'IFNG', 'KLRG1', 'CCL4', 'HIF1A', 'NKG7', 'GNLY', 'TOX', 'CXCL5', 'IL13')
DotPlot(day7_CD8, responder, group.by = 'CAR')
day7_CD8 <- AddModuleScore(day7_CD8, features = list(responder), name = 'responderGenes')
DotPlot(day7_CD8, responder, group.by = 'CAR', scale = FALSE)
DoHeatmap(day7_CD8, features = responder, group.by = 'CAR')
