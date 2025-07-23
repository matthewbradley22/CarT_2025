#Source functions and packages
source("mystore/cartdata/scripts/CarT_project_functions.R") #If this fails can check the IfLoadingSeuratCrashesR script
library(gprofiler2)
library(tibble)

#Make UMAPs plot properly
options(bitmapType="cairo") 

#Load data from beginning (created in donor_analysis script)
T_cells <-  LoadSeuratRds("./mystore/cartdata/data/tcell_deseqdat_singlets.rds")

#Load t cells with cell cycle already regressed out
T_cells_noCC <- LoadSeuratRds('mystore/cartdata/data/T_cells_noCC.rds')

#Start from beginning here
#Cell cycle genes
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

T_cells$day = factor(T_cells$day, levels = c("D7", "D13"))

T_cells <- NormalizeData(T_cells)
T_cells <- FindVariableFeatures(T_cells, selection.method = "vst")
T_cells <- ScaleData(T_cells, features = c(VariableFeatures(T_cells)))
T_cells <- RunPCA(T_cells, features = VariableFeatures(object = T_cells))

T_cells <- CellCycleScoring(T_cells, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

T_cells_noCC <- ScaleData(T_cells, vars.to.regress = c("S.Score", "G2M.Score", "Phase"), features = c(VariableFeatures(T_cells)))
T_cells_noCC <- RunPCA(T_cells_noCC)
T_cells_noCC <- FindNeighbors(T_cells_noCC, dims = 1:25)
T_cells_noCC <- FindClusters(T_cells_noCC, resolution = 0.5)
T_cells_noCC <- RunUMAP(T_cells_noCC, dims = 1:25)

DimPlot(T_cells_noCC, group.by = "Phase")
DimPlot(T_cells, group.by = "Phase")

# T_cells_noCC[[]] %>% filter(CAR_pred == 1 & CAR == 'M1XX') %>% group_by(day, CD_pred) %>% dplyr::summarise(amount = n()) %>% 
#   ggplot(aes(x = day, y = amount, fill = CD_pred))+
#   geom_bar(position="dodge", stat="identity")+
#   ggtitle('M1XX')


T_cells_noCC$CAR_pred <- as.factor(T_cells_noCC$CAR_pred)
bulk_noCC <- getPseudoBulkObject(T_cells_noCC, c("CAR","hypoxia", "Phase", "donor_id", "carExpression", "day", "CD_pred"), return.Seurat = TRUE)
cd4_bulk_noCC = subset(bulk_noCC, CD_pred == "CD4")
cd8_bulk_noCC = subset(bulk_noCC, CD_pred == "CD8")

cd4_bulk_noCC <- FindVariableFeatures(cd4_bulk_noCC)
cd4_bulk_noCC <- ScaleData(cd4_bulk_noCC, vars.to.regress = c("Phase"))
cd4_bulk_noCC <- RunPCA(cd4_bulk_noCC)
ElbowPlot(cd4_bulk_noCC, ndims = 30)
cd4_bulk_noCC <- FindNeighbors(cd4_bulk_noCC, dims = 1:20)
cd4_bulk_noCC <- FindClusters(cd4_bulk_noCC, resolution = 0.5)
cd4_bulk_noCC <- RunUMAP(cd4_bulk_noCC, dims = 1:20)

cd8_bulk_noCC <- FindVariableFeatures(cd8_bulk_noCC)
cd8_bulk_noCC <- ScaleData(cd8_bulk_noCC, vars.to.regress = c("Phase"))
cd8_bulk_noCC <- RunPCA(cd8_bulk_noCC)
ElbowPlot(cd8_bulk_noCC, ndims = 30)
cd8_bulk_noCC <- FindNeighbors(cd8_bulk_noCC, dims = 1:20)
cd8_bulk_noCC <- FindClusters(cd8_bulk_noCC, resolution = 0.5)
cd8_bulk_noCC <- RunUMAP(cd8_bulk_noCC, dims = 1:20)
#saveRDS(cd4_bulk_noCC, "/pfs/stor10/users/home/m/mb223/mystore/cartdata/data/cd4_bulk_noCC.rds")

#saveRDS(cd8_bulk_noCC, "/pfs/stor10/users/home/m/mb223/mystore/cartdata/data/cd8_bulk_noCC.rds")

#Plotting for figure 2
DimPlot(cd8_bulk_noCC, group.by = "day")+ggtitle("CD8 Day UMAP")

# #CD cells
# CD4 <- subset(T_cells, CD_pred == "CD4")
# CD8 <- subset(T_cells, CD_pred == "CD8")
# 
# CD4_day7 <- subset(CD4, day == "D7")
# CD4_day13 <- subset(CD4, day == "D13")
# CD8_day7 <- subset(CD8, day == "D7")
# CD8_day13 <- subset(CD8, day == "D13")
# 
# #CD4_day13 <- getPseudoBulkObject(CD4_day13, c("CAR","hypoxia", "Phase", "donor_id"), return.Seurat = TRUE)
# CD8_day13 <- CellCycleScoring(CD8_day13, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
# CD8_day13 <- prepUMAP(CD8_day13, regressCC = TRUE)
# 
# varList = c("CAR","hypoxia", "Phase", "donor_id")
# plotList = list()
# for(i in 1:length(varList)){
#   plotList[[i]] = DimPlot(CD8_day13, group.by = varList[[i]])
# }
# grid.arrange(grobs = plotList, nrow = 2)
# 
# CD4_day13 <- prepUMAP(CD4_day13)
# CD8_day7 <- prepUMAP(CD8_day7)
# CD4_day13 <- prepUMAP(CD4_day13)
# 
# CD4_day7Bulk <- getPseudoBulkObject(CD4_day7, c("CAR","hypoxia", "Phase", "donor_id"), return.Seurat = TRUE)
# CD4_day7Bulk <- CellCycleScoring(CD4_day7Bulk, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
# CD4_day7Bulk <- prepUMAP(CD4_day7Bulk, regressCC = TRUE)
# DimPlot(CD4_day7Bulk, group.by = "Phase")
# 
# CD4Bulk <- getPseudoBulkObject(CD4, c("CAR","hypoxia", "Phase", "donor_id", "day"), return.Seurat = TRUE)
# CD4Bulk <- CellCycleScoring(CD4Bulk, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
# CD4Bulk <- prepUMAP(CD4Bulk, regressCC = TRUE)
# DimPlot(CD4Bulk, group.by = "Phase")
# 
# CD8Bulk <- getPseudoBulkObject(CD8, c("CAR","hypoxia", "Phase", "donor_id", "day"), return.Seurat = TRUE)
# CD8Bulk <- CellCycleScoring(CD8Bulk, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
# CD8Bulk <- prepUMAP(CD8Bulk, regressCC = TRUE)
# DimPlot(CD8Bulk, group.by = "Phase")
# 
# varList = c("CAR","hypoxia", "Phase", "donor_id", "day")
# plotList = list()
# for(i in 1:length(varList)){
#   plotList[[i]] = DimPlot(CD4Bulk, group.by = varList[[i]])
# }
# 
# grid.arrange(grobs = plotList, nrow = 2)

#Regress out cell cycle from deseq2 dat

#Can start here if loading in cd bulk data directly
cd4_bulk_noCC_CAR = subset(cd4_bulk_noCC, CAR_pred == 1)
cd8_bulk_noCC_CAR = subset(cd8_bulk_noCC, CAR_pred == 1)
cd4_bulk_noCC_CAR[[]]<- cd4_bulk_noCC_CAR[[]] %>% mutate(specificGroups = paste(day, CAR, hypoxia, sep='_'))
cd8_bulk_noCC_CAR[[]]<- cd8_bulk_noCC_CAR[[]] %>% mutate(specificGroups = paste(day, CAR, hypoxia, sep='_'))
T_cells[[]]<- T_cells[[]] %>% mutate(specificGroups = paste(day, CAR, hypoxia, sep='_'))

cd4_bulk_d7 = subset(cd4_bulk, day == 'D7')
cd4_bulk_d13= subset(cd4_bulk, day == 'D13')
cd8_bulk_d7 = subset(cd8_bulk, day == 'D7')
cd8_bulk_d13= subset(cd8_bulk, day == 'D13')
cd4_bulk_d7[[]]<- cd4_bulk_d7[[]] %>% mutate(specificGroups = paste(CAR, hypoxia, sep='_'))
cd4_bulk_d13[[]]<- cd4_bulk_d13[[]] %>% mutate(specificGroups = paste(CAR, hypoxia, sep='_'))
cd8_bulk_d7[[]]<- cd8_bulk_d7[[]] %>% mutate(specificGroups = paste(CAR, hypoxia, sep='_'))
cd8_bulk_d13[[]]<- cd8_bulk_d13[[]] %>% mutate(specificGroups = paste(CAR, hypoxia, sep='_'))
#Gene lists from https://www.science.org/doi/full/10.1126/sciadv.adp4008
#cytoGenes = list(c('GZMB', 'PRF1', 'FASLG'))
activatedGenes = list(c('TNFRSF9', 'TBX21',  'ZBED2'))
inflammatoryGenes = list(c('CRTAM', 'IFNG', 'TNF', 'CSF2', 'XCL1',  'XCL2'))
restingMemGenes <- list(c('TCF7', 'CCR7', 'LEF1', 'SELL'))
MemoryScore <- list(c('CCR7', 'TCF7', 'Il7R', 'AQP3', 'CD27', 'LTB')) #From Bai https://www.nature.com/articles/s41586-024-07762-w
proliferationScore <- list(c('MKI67', 'TYMS', 'TOP2A', 'ASPM')) #From Bai
cytotoxicScore <- list(c( 'GZMA', 'GZMB', 'GZMH', 'GNLY', 'PRF1', 'NKG7')) #Some of these associated with nonresponders in Haradhvala paper
lateEffectorDiff <- list(c('HOPX', 'ENTPD1', 'LAG3', 'HAVCR3',  'GNLY'))

cd4_bulk_noCC_CAR <- AddModuleScore(cd4_bulk_noCC_CAR, inflammatoryGenes, name = 'inflammatoryGenes')
cd8_bulk_noCC_CAR <-  AddModuleScore(cd8_bulk_noCC_CAR, inflammatoryGenes, name = 'inflammatoryGenes')

VlnPlot(cd8_bulk_noCC_CAR, features = 'activatedGenes1', group.by = 'hypoxia', pt.size = 0)+ theme(legend.position = 'none',
                                                                                                       axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle('CD8 Late Cytotoxic Genes')+
  labs(caption = "HOPX, ENTPD1, LAG3, HAVCR3, GNLY")+
  theme(plot.caption = element_text(color = "orange"))


T_cells <-  AddModuleScore(T_cells, restingMemGenes, name = 'restingMemGenes')
cd4 = subset(T_cells, CD_pred == 'CD4' & CAR_pred == 1)
cd8 = subset(T_cells, CD_pred == 'CD8' & CAR_pred == 1)
cd8Day7 = subset(cd8, day== 'D7')
cd8Day7 = AddModuleScore(cd8Day7, MemoryScore, name = 'MemoryScore')
VlnPlot(cd4Day7, features = 'MT-ND5',  group.by = 'CAR', pt.size = 0) +theme(legend.position = 'none')


#Further zoomed in umaps
cd4_bulk_noCC <- subset(cd4_bulk_noCC, CAR_pred == 1)
cd4_bulk_noCC <- prepUMAP(cd4_bulk_noCC, dimNeighbors = 25, regressCC = 'phase')
DimPlot(cd4_bulk_noCC, group.by = 'day')

cluster_counts <- FetchData(cd4_bulk_noCC, layer="counts", vars=rownames(cd4_bulk_noCC))

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(t(cluster_counts),
                                colData = cd4_bulk_noCC[[]],
                                design = ~Phase + donor_id + CAR + hypoxia + day)
dds <- DESeq(dds)
resultsNames(dds)
resLFC <- lfcShrink(dds, coef="day_D13_vs_D7", type="apeglm")
resOrdered <- resLFC[order(resLFC$padj),]
importantGenes <- resOrdered %>% as.data.frame() %>% filter(log2FoldChange < 0) %>% head(n=100) %>% rownames()

gostres <- gost(query = importantGenes, 
                organism = "hsapiens", ordered_query = FALSE, 
                multi_query = FALSE, significant = TRUE,
                user_threshold = 0.05, evcodes = TRUE)
gostres$result
gostplot(gostres, capped = TRUE, interactive = TRUE)

#visualize clusters
cd8_bulk_noCC_CAR <- prepUMAP(cd8_bulk_noCC_CAR)
DimPlot(cd8_bulk_noCC_CAR, group.by = 'day')
cd8_bulk_noCC_CAR <-  AddModuleScore(cd8_bulk_noCC_CAR, cytotoxicScore, name = 'cytotoxicScore')
VlnPlot(cd8_bulk_noCC_CAR, features = 'cytotoxicScore1', group.by = 'day')
DotPlot(object = cd8_bulk_noCC_CAR, features = 'cytotoxicScore1', group.by = 'hypoxia')

#Look at Deng paper list of genes associated with recovery
DengOutcomeDat <- read.csv('/pfs/stor10/users/home/m/mb223/mystore/cartdata/data/DengCRvsPR.csv', sep = ';')
worseOutcome <- subset(DengOutcomeDat, avg_logFC < 0)
worseOutcome <- list(worseOutcome[,1])
betterOutcome <- subset(DengOutcomeDat, avg_logFC > 0)
betterOutcome <- list(betterOutcome[,1])

cd4_bulk_noCC_CAR <-  AddModuleScore(cd4_bulk_noCC_CAR, worseOutcome, name = 'worseOutcome')
cd8_bulk_noCC_CAR <-  AddModuleScore(cd8_bulk_noCC_CAR, worseOutcome, name = 'worseOutcome')
cd4_bulk_noCC_CAR <-  AddModuleScore(cd4_bulk_noCC_CAR, betterOutcome, name = 'betterOutcome')
cd8_bulk_noCC_CAR <-  AddModuleScore(cd8_bulk_noCC_CAR, betterOutcome, name = 'betterOutcome')

FeaturePlot(cd4_bulk_noCC_CAR, features = 'betterOutcome1')
FeaturePlot(cd8_bulk_noCC_CAR, features = 'betterOutcome1')
DimPlot(cd8_bulk_noCC_CAR, group.by = 'day')
T_cells %>% subset(CAR_pred == 1 & CD_pred == 'CD8') %>% DotPlot(object = ., features = 'betterOutcome1', group.by = 'Phase')

T_cells <-  AddModuleScore(T_cells, worseOutcome, name = 'worseOutcome')
T_cells <-  AddModuleScore(T_cells, betterOutcome, name = 'betterOutcome')
FeaturePlot(T_cells, features = 'worseOutcome1')
DimPlot(T_cells, group.by = 'CD_pred')

#DEGs between CAR groups
cd4_bulk_noCC_CAR = subset(cd4_bulk_noCC, CAR_pred == 1 & CAR != 'untransduced')
cd8_bulk_noCC_CAR = subset(cd8_bulk_noCC, CAR_pred == 1 & CAR != 'untransduced')

#MBBz has very different ratios of cd4/cd8 cells so we want to make it the reference level
cd4_bulk_7 <- subset(cd4_bulk_noCC_CAR, day == 'D7')
cd4_bulk_7$CAR <- relevel(cd4_bulk_7$CAR, 'MBBz')
cd4_bulk_13 <- subset(cd4_bulk_noCC_CAR, day == 'D13')
cd4_bulk_13$CAR <- relevel(cd4_bulk_13$CAR, 'MBBz')
cd8_bulk_7 <- subset(cd8_bulk_noCC_CAR, day == 'D7')
cd8_bulk_7$CAR <- relevel(cd8_bulk_7$CAR, 'MBBz')
cd8_bulk_13 <- subset(cd8_bulk_noCC_CAR, day == 'D13')
cd8_bulk_13$CAR <- relevel(cd8_bulk_13$CAR, 'MBBz')

counts_cd4_d7 <- FetchData(cd4_bulk_7, layer="counts", vars=rownames(cd4_bulk_7))
counts_cd4_d13 <- FetchData(cd4_bulk_13, layer="counts", vars=rownames(cd4_bulk_13))
counts_cd8_d7 <- FetchData(cd8_bulk_7, layer="counts", vars=rownames(cd8_bulk_7))
counts_cd8_d13 <- FetchData(cd8_bulk_13, layer="counts", vars=rownames(cd8_bulk_13))

dds_cd4_7 <- DESeqDataSetFromMatrix(t(counts_cd4_d7),
                              colData = cd4_bulk_7[[]],
                              design = ~hypoxia + donor_id + CAR)
dds_cd4_13 <- DESeqDataSetFromMatrix(t(counts_cd4_d13),
                                    colData = cd4_bulk_13[[]],
                                    design = ~hypoxia + donor_id + CAR)
dds_cd8_7 <- DESeqDataSetFromMatrix(t(counts_cd8_d7),
                                     colData = cd8_bulk_7[[]],
                                     design = ~hypoxia + donor_id + CAR)
dds_cd8_13 <- DESeqDataSetFromMatrix(t(counts_cd8_d13),
                                    colData = cd8_bulk_13[[]],
                                    design = ~hypoxia + donor_id + CAR)

dds_cd4_7 <- DESeq(dds_cd4_7)
dds_cd4_13 <- DESeq(dds_cd4_13)
dds_cd8_7 <- DESeq(dds_cd8_7)
dds_cd8_13 <- DESeq(dds_cd8_13)

vsd_cd4_7 <- vst(dds_cd4_7, blind=FALSE)
vsd_cd4_13 <- vst(dds_cd4_13, blind=FALSE)
vsd_cd8_7 <- vst(dds_cd8_7, blind=FALSE)
vsd_cd8_13 <- vst(dds_cd8_13, blind=FALSE)

plotPCA(vsd_cd4_7, intgroup=c('CAR')) + ggtitle('CD4 day 7, CAR+')
plotPCA(vsd_cd4_13, intgroup=c('CAR')) + ggtitle('CD4 day 13, CAR+')
plotPCA(vsd_cd8_7, intgroup=c('CAR')) + ggtitle('CD8 day 7, CAR+')
plotPCA(vsd_cd8_13, intgroup=c('CAR')) + ggtitle('CD8 day 13, CAR+')

resultsNames(dds_cd8_7)

res_cd4_7 <- lfcShrink(dds_cd4_7, coef="CAR_M28z_vs_MBBz", type="apeglm")
res_cd4_13 <- lfcShrink(dds_cd4_13, coef="CAR_M28z_vs_MBBz", type="apeglm")
res_cd8_7 <- lfcShrink(dds_cd8_7, coef="CAR_M28z_vs_MBBz", type="apeglm")
res_cd8_13 <- lfcShrink(dds_cd8_13, coef="CAR_M28z_vs_MBBz", type="apeglm")

res_cd4_7[order(res_cd4_7$padj),]  %>% as.data.frame() %>% head(n = 500) %>% 
  write.csv(., './mystore/cartdata/data/MBBzDEGLists/CD4_day7_M28zvsMBBz.csv')
res_cd4_13[order(res_cd4_13$padj),]  %>% as.data.frame() %>% head(n = 500) %>% 
  write.csv(., './mystore/cartdata/data/MBBzDEGLists/CD4_day13_M28zvsMBBz.csv')
res_cd8_7[order(res_cd8_7$padj),]  %>% as.data.frame() %>% head(n = 500) %>% 
  write.csv(., './mystore/cartdata/data/MBBzDEGLists/CD8_day7_M28zvsMBBz.csv')
res_cd8_13[order(res_cd8_13$padj),]  %>% as.data.frame() %>% head(n = 500) %>% 
  write.csv(., './mystore/cartdata/data/MBBzDEGLists/CD8_day13_M28zvsMBBz.csv')



#Compare deg lists
cd4Day7Genes <- res_cd4_7[order(res_cd4_7$padj),]  %>% as.data.frame() %>% head(n = 1000) %>% rownames()
cd8Day7Genes <- res_cd8_7[order(res_cd8_7$padj),]  %>% as.data.frame() %>% head(n = 1000) %>% rownames()
geneComps <- c(cd4Day7Genes, cd8Day7Genes)
geneCompCD4 <- res_cd4_7[rownames(res_cd4_7) %in% geneComps,] %>% as.data.frame() %>% select(log2FoldChange) %>% 
  rename(log2FoldChange_CD4 = log2FoldChange) %>% rownames_to_column()
geneCompCD8 <- res_cd8_7[rownames(res_cd8_7) %in% geneComps,] %>% as.data.frame() %>% select(log2FoldChange) %>% 
  rename(log2FoldChange_CD8 = log2FoldChange) %>% rownames_to_column()
cd8Vscd4 <- left_join(geneCompCD4, geneCompCD8, by = 'rowname')
differences <- cd8Vscd4 %>% mutate(diff = log2FoldChange_CD8 - log2FoldChange_CD4) %>% arrange(desc(abs(diff)))
ggplot(cd8Vscd4, aes(x = log2FoldChange_CD4, y = log2FoldChange_CD8))+
  geom_point() +
  geom_abline(slope=1) +
  xlab('CD4 log2 fold change') +
  ylab('CD8 log2 fold change') +
  ggtitle('Day 7 fold change comparison') +
  labs(caption = 'Black line has slope of 1') +
  geom_text(data=subset(cd8Vscd4, rowname %in% head(differences, n = 5)$rowname),
            aes(log2FoldChange_CD4,log2FoldChange_CD8,label=rowname), colour = "red", vjust=2)

plotCounts(dds_cd4_7, gene='ALK', intgroup="CAR")

cd4Day13Genes <- res_cd4_13[order(res_cd4_13$padj),]  %>% as.data.frame() %>% head(n = 1000) %>% rownames()
cd8Day13Genes <- res_cd8_13[order(res_cd8_13$padj),]  %>% as.data.frame() %>% head(n = 1000) %>% rownames()
geneComps <- c(cd4Day13Genes, cd8Day13Genes)
geneCompCD4 <- res_cd4_13[rownames(res_cd4_13) %in% geneComps,] %>% as.data.frame() %>% select(log2FoldChange) %>% 
  rename(log2FoldChange_CD4 = log2FoldChange) %>% rownames_to_column()
geneCompCD8 <- res_cd8_13[rownames(res_cd8_13) %in% geneComps,] %>% as.data.frame() %>% select(log2FoldChange) %>% 
  rename(log2FoldChange_CD8 = log2FoldChange) %>% rownames_to_column()
cd8Vscd4_13 <- left_join(geneCompCD4, geneCompCD8, by = 'rowname')

differences13 <- cd8Vscd4_13 %>% mutate(diff = log2FoldChange_CD8 - log2FoldChange_CD4) %>% arrange(desc(abs(diff)))

ggplot(cd8Vscd4_13, aes(x = log2FoldChange_CD4, y = log2FoldChange_CD8))+
  geom_point() +
  geom_abline(slope=1) +
  xlab('CD4 log2 fold change') +
  ylab('CD8 log2 fold change') +
  ggtitle('Day 13 fold change comparison') +
  labs(caption = 'Black line has slope of 1')+
  geom_text(data=subset(cd8Vscd4_13, rowname %in% head(differences13, n = 5)$rowname),
            aes(log2FoldChange_CD4,log2FoldChange_CD8,label=rowname), colour = "red", 
            position=position_jitter(height=1), size = 2)

#HLA differences
hla4Genes <- rownames(cd4_bulk_noCC_CAR[['RNA']]$scale.data)[grep(pattern = 'HLA', x = rownames(cd4_bulk_noCC_CAR[['RNA']]$scale.data))]
DoHeatmap(object = cd4_bulk_noCC_CAR, features = hla4Genes, group.by = 'CAR')+ ggtitle('CD4 HLA gene expression') + 
  guides(colour=FALSE)
hla8Genes <- rownames(cd8_bulk_noCC_CAR[['RNA']]$scale.data)[grep(pattern = 'HLA', x = rownames(cd8_bulk_noCC_CAR[['RNA']]$scale.data))]
DoHeatmap(object = cd8_bulk_noCC_CAR, features = hla8Genes, group.by = 'CAR')+ ggtitle('CD8 HLA gene expression')+
  guides(colour=FALSE)

hlaGenes <- rownames(T_cells[['RNA']]$data)[grep(pattern = 'HLA', x = rownames(T_cells[['RNA']]$data))]

T_cells_plot <- subset(T_cells, CAR_pred == 1 & CAR != 'untransduced')
cd4 <- subset(T_cells_plot, CD_pred == 'CD4')
hla4Genes <- rownames(cd4[['RNA']]$data)[grep(pattern = 'HLA', x = rownames(cd4[['RNA']]$data))]
DoHeatmap(object = cd4, features = hla4Genes, group.by = 'CAR')+ ggtitle('CD4 HLA gene expression')+
  guides(colour=FALSE)

#Look at hypoxia separately in CAR groups
cd4_bulk_noCC_MBBz_7 <- getHypoxiaPCA(cd4_bulk_noCC, whichDay = 'D7', 'MBBz')
cd4_bulk_noCC_M28z_7 = getHypoxiaPCA(cd4_bulk_noCC, whichDay = 'D7', 'M28z')
cd4_bulk_noCC_M1XX_7 = getHypoxiaPCA(cd4_bulk_noCC, whichDay = 'D7', 'M1XX')


plotPCA(cd4_bulk_noCC_MBBz_7, intgroup=c('hypoxia')) + ggtitle('CD4 day 7, MBBz CAR+')
plotPCA(cd4_bulk_noCC_M28z_7, intgroup=c('hypoxia')) + ggtitle('CD4 day 7, M28z CAR+')
plotPCA(cd4_bulk_noCC_M1XX_7, intgroup=c('hypoxia')) + ggtitle('CD4 day 7, M1XX CAR+')

cd4_bulk_noCC_MBBz_13 <- getHypoxiaPCA(cd4_bulk_noCC, whichDay = 'D13', 'MBBz')
cd4_bulk_noCC_M28z_13 = getHypoxiaPCA(cd4_bulk_noCC, whichDay = 'D13', 'M28z')
cd4_bulk_noCC_M1XX_13 = getHypoxiaPCA(cd4_bulk_noCC, whichDay = 'D13', 'M1XX')


plotPCA(cd4_bulk_noCC_MBBz_13, intgroup=c('hypoxia')) + ggtitle('CD4 day 13, MBBz CAR+')
plotPCA(cd4_bulk_noCC_M28z_13, intgroup=c('hypoxia')) + ggtitle('CD4 day 13, M28z CAR+')
plotPCA(cd4_bulk_noCC_M1XX_13, intgroup=c('hypoxia')) + ggtitle('CD4 day 13, M1XX CAR+')

cd8_bulk_noCC_MBBz_7 <- getHypoxiaPCA(cd8_bulk_noCC, whichDay = 'D7', 'MBBz')
cd8_bulk_noCC_M28z_7 = getHypoxiaPCA(cd8_bulk_noCC, whichDay = 'D7', 'M28z')
cd8_bulk_noCC_M1XX_7 = getHypoxiaPCA(cd8_bulk_noCC, whichDay = 'D7', 'M1XX')

plotPCA(cd8_bulk_noCC_MBBz_7, intgroup=c('hypoxia')) + ggtitle('CD8 day 7, MBBz CAR+')
plotPCA(cd8_bulk_noCC_M28z_7, intgroup=c('hypoxia')) + ggtitle('CD8 day 7, M28z CAR+')
plotPCA(cd8_bulk_noCC_M1XX_7, intgroup=c('hypoxia')) + ggtitle('CD8 day 7, M1XX CAR+')

cd8_bulk_noCC_MBBz_13 <- getHypoxiaPCA(cd8_bulk_noCC, whichDay = 'D13', 'MBBz')
cd8_bulk_noCC_M28z_13 = getHypoxiaPCA(cd8_bulk_noCC, whichDay = 'D13', 'M28z')
cd8_bulk_noCC_M1XX_13 = getHypoxiaPCA(cd8_bulk_noCC, whichDay = 'D13', 'M1XX')

plotPCA(cd8_bulk_noCC_MBBz_13, intgroup=c('hypoxia')) + ggtitle('CD8 day 13, MBBz CAR+')
plotPCA(cd8_bulk_noCC_M28z_13, intgroup=c('hypoxia')) + ggtitle('CD8 day 13, M28z CAR+')
plotPCA(cd8_bulk_noCC_M1XX_13, intgroup=c('hypoxia')) + ggtitle('CD8 day 13, M1XX CAR+')

getHypoxiaPCA <- function(dat, whichDay, car){
  subsetDat = subset(dat, CAR_pred == 1 & CAR == car & day == whichDay)
  counts <- FetchData(subsetDat, layer="counts", vars=rownames(subsetDat))
  
  dds <- DESeqDataSetFromMatrix(t(counts),
                                           colData = subsetDat[[]],
                                           design = ~hypoxia + donor_id)
  dds <- DESeq(dds)
  vsd_dds <- vst(dds, blind=FALSE)
}


#Look at response markers from https://aacrjournals.org/clincancerres/article/29/20/4139/729415/Single-Cell-RNA-Analysis-Reveals-Cell-Intrinsic
DotPlot(T_cells_noCC, features = c('CCL3', 'LAG3', 'IFNG', 'KLRG1', 'CCL4', 'HIF1A'), group.by = 'Phase')
DotPlot(T_cells_noCC, features = c('CD52', 'CXCR4', 'TNFSF8', 'CXCR6', 'JUNB'), group.by = 'hypoxia')
DotPlot(T_cells_noCC, features = c('CXCL2', 'CCL4', 'CCL5'), group.by = 'CAR')

#Looking at genes from https://www.nature.com/articles/s41467-025-59904-x/figures/2
cytokineActivity <- c('IL1A', 'TGFB1', 'TGFB3', 'TNFA', 'IL6', 'IL15', 'CXCL12', 'IL13', 'IFN1',
                      'IL1B', 'IL2', 'IL3', 'IL22', 'IL4', 'IL10', 'IFNL', 'IL36', 'IL17A', 'IL12', 'IL21', 
                      'IFNG', 'IL27')
cytoxicGenes <- c('GZMA', 'GZMB', 'GZMH', 'GZMM', 'GZMK', 'PRF1', 'GNLY', 'NKG7')
DotPlot(T_cells_noCC, features = cytokineActivity, group.by = 'hypoxia')
hypoAvg <- AverageExpression(T_cells_noCC, group.by = 'hypoxia', return.seurat = T)
DoHeatmap(carAvg, features = cytoxicGenes, slot = 'data', draw.lines = FALSE)

#Prep pseudobulk for shiny
cd4_bulk = subset(T_cells_noCC, CD_pred == "CD4" & CAR_pred == 1)
cd8_bulk = subset(T_cells_noCC, CD_pred == "CD8" & CAR_pred == 1)

cd4_bulk_noCC <- getPseudoBulkObject(cd4_bulk, c("CAR","hypoxia", "donor_id", "day"))
cd4_bulk_noCC <- DESeq(cd4_bulk_noCC)
comparisons <- resultsNames(cd4_bulk_noCC)[grepl('hypoxia|CAR',resultsNames(cd4_bulk_noCC))] 
comparisons <- unlist(lapply(comparisons, FUN = function(x){
  gsub('_', ' ', x)
}))

lfcShrink(cd4_bulk_noCC, coef="CAR_M28z_vs_M1XX", type="apeglm")


