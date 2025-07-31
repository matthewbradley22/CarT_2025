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

#### Regress out cell cycle ####
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


#Create pseudobulk UMAPs for fig 2
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

#Plotting for figure 2
varsOfInterest <- c("CAR","hypoxia", "Phase", "donor_id", "carExpression", "day")
titles = c('CD4 CAR UMAP', 'CD4 Hypoxia UMAP', 'CD4 Phase UMAP', 'CD4 Donor UMAP', 'CD4 CAR Expression UMAP', 'CD4 Day UMAP')
titles_cd8 = c('CD8 CAR UMAP', 'CD8 Hypoxia UMAP', 'CD8 Phase UMAP', 'CD8 Donor UMAP', 'CD8 CAR Expression UMAP', 'CD8 Day UMAP')

for(i in 1:length(varsOfInterest)){
  var = varsOfInterest[i]
  pdf(file = paste0('/pfs/stor10/users/home/m/mb223/mystore/cartdata/Plots/pseudobulk/umaps/CD4', var, '.pdf'), width = 6, height = 6)
  print(DimPlot(cd4_bulk_noCC, group.by = var) + ggtitle(titles[i]) +
          theme(plot.title = element_text(size = 28, face = "bold"),
                legend.text=element_text(size=18)))
  dev.off()
}

for(i in 1:length(varsOfInterest)){
  var = varsOfInterest[i]
  pdf(file = paste0('/pfs/stor10/users/home/m/mb223/mystore/cartdata/Plots/pseudobulk/umaps/CD8', var, '.pdf'), width = 6, height = 6)
  print(DimPlot(cd8_bulk_noCC, group.by = var) + ggtitle(titles_cd8[i])+
          theme(plot.title = element_text(size = 28, face = "bold"),
                legend.text=element_text(size=18))) 
  dev.off()
}


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


#### Look at various gene lists of relevant T-Cell genes across groups ####
#Can start here if loading in cd bulk data directly
cd4_bulk_noCC_CAR = subset(cd4_bulk_noCC, carExpression == 'carPos')
cd8_bulk_noCC_CAR = subset(cd8_bulk_noCC,carExpression == 'carPos')
cd4_bulk_noCC_CAR[[]]<- cd4_bulk_noCC_CAR[[]] %>% mutate(specificGroups = paste(day, CAR, hypoxia, sep='_'))
cd8_bulk_noCC_CAR[[]]<- cd8_bulk_noCC_CAR[[]] %>% mutate(specificGroups = paste(day, CAR, hypoxia, sep='_'))
T_cells[[]]<- T_cells[[]] %>% mutate(specificGroups = paste(day, CAR, hypoxia, sep='_'))

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

VlnPlot(cd4_bulk_noCC_CAR, features = 'inflammatoryGenes1', group.by = 'hypoxia', pt.size = 0)+ theme(legend.position = 'none',
                                                                                                       axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle('CD8 Inflammatory Genes')+
  labs(caption = paste0(unlist(inflammatoryGenes), collapse = " "))+
  theme(plot.caption = element_text(color = "orange"))


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

