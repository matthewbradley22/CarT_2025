library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(tidyr)
library(ggrepel)

#Source fucntions
source("mystore/cartdata/scripts/CarT_project_functions.R")

#Make plotting work properly
options(bitmapType="cairo")

#Load in integrated data (can start around line 290 with this data to skip preprocessing)
dengDatInt <- LoadSeuratRds("mystore/cartdata/data/publicDengDatIntegrated.rds")

#Read in cap id dat
dengDatCapInt <- LoadSeuratRds('mystore/cartdata/publicData/dengDatCapInt.rds')

#Read in data and metadata from Characteristics of anti-CD19 CAR T cell infusion products associated with
#efficacy and toxicity in patients with large B cell lymphomas https://www.nature.com/articles/s41591-020-1061-7#data-availability

dengMeta <- read_csv("mystore/cartdata/publicData/dengMetaData.csv")
colnames(dengMeta) = dengMeta[2,]
dengMeta <- dengMeta[-c(1,2),]



DengFiles <- paste0('mystore/cartdata/publicData/singleCellData/', list.files('mystore/cartdata/publicData/singleCellData/'))

dengDat <- lapply(DengFiles, FUN = function(x){
  Read10X(data.dir = x)
})

dengDat <- lapply(dengDat, FUN = function(x){
  CreateSeuratObject(counts = x, project = "deng", min.cells = 3, min.features = 200)
})

dengDat <- lapply(dengDat, FUN=function(x){
  x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^MT-")
  x
})

metaCols <- colnames(dengMeta)[4:17]
for(i in 1:length(dengDat)){
  for(j in 1:length(metaCols)){
    name = metaCols[j]
    dengDat[[i]][[name]] = c(dengMeta[i,name])
  }
}

VlnPlot(dengDat[[3]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Can require 15 cores on hpc2n depending on filtering (integration step later in pipeline is where issues arise with memory).
dengDat <- lapply(dengDat, FUN = function(x){
  x <- subset(x, subset = nFeature_RNA > 50 & nFeature_RNA < 4000 & percent.mt < 10)
})

dengDat <- lapply(dengDat, FUN = function(x){
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x)
  x <- ScaleData(x)
  x <- RunPCA(x)
  x
})

# dengMetaDf <- dengDat[[1]][[]]
# 
# for(i in 2:length(dengDat)){
#   newMeta <- dengDat[[i]][[]]
#   dengMetaDf <- rbind(dengMetaDf, newMeta)
# }

dengDatMerge <- merge(x = dengDat[[1]], y = dengDat[-1], add.cell.ids = as.character(seq(1,24)),
                      merge.data = TRUE)

dengDatMerge <- NormalizeData(dengDatMerge) 
dengDatMerge <- FindVariableFeatures(dengDatMerge)
dengDatMerge <- ScaleData(dengDatMerge)
dengDatMerge <- RunPCA(dengDatMerge)


dengDatInt <- IntegrateLayers(object = dengDatMerge, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
                        verbose = FALSE)

#re-join layers after integration
dengDatInt[["RNA"]] <- JoinLayers(dengDatInt[["RNA"]])

ElbowPlot(dengDatInt, ndims = 30)
dengDatInt <- FindNeighbors(dengDatInt, reduction = "integrated.cca", dims = 1:20)
dengDatInt <- FindClusters(dengDatInt, resolution = 1)
dengDatInt <- RunUMAP(dengDatInt, dims = 1:20, reduction = "integrated.cca")
DimPlot(dengDatInt)
#Saving data runnning UMAP
#SaveSeuratRds(dengDatInt, "mystore/cartdata/data/publicDengDatIntegrated.rds")

#Subset to 24 samples (other 16 are just used for IACs validation)

#Follow same steps with capid data
DengCapFiles <- paste0('mystore/cartdata/publicData/capIdData/', list.files('mystore/cartdata/publicData/capIdData/'))
DengCapFiles <- DengCapFiles[1:24]

dengCapDat <- lapply(DengCapFiles, FUN = function(x){
  Read10X(data.dir = x)
})

#Cannot have min features be too high
dengCapDat <- lapply(dengCapDat, FUN = function(x){
  CreateSeuratObject(counts = x, project = "dengCap", min.cells = 1, min.features = 1)
})

dengCapDat <- lapply(dengCapDat, FUN=function(x){
  x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^MT-")
  x
})

VlnPlot(dengCapDat[[3]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

dengCapDat <- lapply(dengCapDat, FUN = function(x){
  x <- subset(x, subset = nFeature_RNA > 0 & percent.mt < 50)
})

dengCapDat <- lapply(dengCapDat, FUN = function(x){
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x)
  x <- ScaleData(x)
  x <- RunPCA(x)
  x
})


dengCapDatMerge <- merge(x = dengCapDat[[1]], y = dengCapDat[-1], add.cell.ids = as.character(seq(1,24)),
                      merge.data = TRUE)

dengCapDatMerge <- NormalizeData(dengCapDatMerge) 
dengCapDatMerge <- FindVariableFeatures(dengCapDatMerge)
dengCapDatMerge <- ScaleData(dengCapDatMerge)
dengCapDatMerge <- RunPCA(dengCapDatMerge)

dengDatCapInt <- IntegrateLayers(object = dengCapDatMerge, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
                              verbose = FALSE)

#re-join layers after integration
dengDatCapInt[["RNA"]] <- JoinLayers(dengDatCapInt[["RNA"]])

ElbowPlot(dengDatCapInt, ndims = 30)
dengDatCapInt <- FindNeighbors(dengDatCapInt, reduction = "integrated.cca", dims = 1:20)
dengDatCapInt <- FindClusters(dengDatCapInt, resolution = 1)
dengDatCapInt <- RunUMAP(dengDatCapInt, dims = 1:20, reduction = "integrated.cca")

DimPlot(dengDatCapInt, reduction = "umap")

####Hierarchical cell labeling ####

#list of vectors to iterate through for assigning cell types
# gene markers from supplement of Feldman Defining T Cell States Associated with Response to Checkpoint Immunotherapy in Melanoma
cellTypeMarkers <- list('T_cells' = list(required = c('CD3E'), optional = c('CD2'), cannotHave = c('NCR1', 'NCAM1', 'FOXP3')),
                        'CD4_T' = list(required = c('CD4', 'CD3E'), cannotHave = c('NCR1', 'NCAM1', 'FOXP3')),
                        'CD8_T' = list(optional = c('CD8A', 'CD8B'), cannotHave = c('NCR1', 'NCAM1', 'FOXP3')),
                        'CD4_reg' = list(required = c('FOXP3', 'CD3E', 'CD4'), optional = c('CTLA4','IL2RA'), cannotHave = c('CD8A', 'CD8B')),
                        'CD8_reg' = list(required = c('FOXP3', 'CD3E'), optional = c('CD8A','CD8B'), 
                                         optional2 = c('CTLA4','IL2RA'), cannotHave = c('CD4')),
                        'Activated_T_cells' = list(required = c('CD3E', 'CD2', 'CD28'), 
                                                   optional = c('IL2RA', 'CD69', 'ICOS', 'TNFRSF4', 'TNFRSF9','CD27'), 
                                         optional2 = c('IL2', 'TNF', 'IFNG'), 
                                         optional3 = c('CD8A','CD8B', 'CD4'), cannotHave = c('NCR1', 'NCAM1')),
                        'Exhausted_T_cells' =  list(required = c('CD3E', 'CD2'), 
                                                    optional = c('PDCD1','CTLA4', 'BTLA',
                                                                 'KIR3DL1','LAG3','HAVCR2','ADORA2A','HAVCR1'), 
                                                    optional2 = c('CD8A','CD8B','CD4'), cannotHave = c('NCR1', 'NCAM1')),
                        'Memory_T_cells' = list(required = c('CD3E', 'SELL', 'CCR7', 'CD28'), cannotHave = c('FOXP3', 'CD4', 'CD8A', 'CD8B')),
                        'Memory_CD4_cells' = list(required = c('CD3E', 'SELL', 'CCR7', 'CD28', 'CD4'), cannotHave = c('FOXP3', 'CD8A', 'CD8B')),
                        'Memory_CD8_cells' = list(required = c('CD3E', 'SELL', 'CCR7', 'CD28'), optional = c('CD8A','CD8B'), cannotHave = c('FOXP3', 'CD4')))

T_cell_markers <- cellTypeMarkers[[1]]
CD4_markers <- cellTypeMarkers[[2]]
CD8_markers <- cellTypeMarkers[[3]]
CD4_reg_markers <- cellTypeMarkers[[4]]
CD8_reg_markers <- cellTypeMarkers[[5]]
Activated_markers <- cellTypeMarkers[[6]]
Exhausted_markers <- cellTypeMarkers[[7]]
Memory_markers <- cellTypeMarkers[[8]]
Memory_markers_cd4 <- cellTypeMarkers[[9]]
Memory_markers_cd8 <- cellTypeMarkers[[10]]

#Function that takes required and optional genes (and genes that cannot be present) 
#as defined in Feldman supplemental excel sheet 'predefined cell markers fig1E-G'
#and labels cells that fit the description
labelCells <- function(cellMarkers,
                       cellData,
                       cellName) {
  cellData[[cellName]] <- 'no'
  #Drop = false so even if returned data is one row we can subset in the same way in for loop
  required = GetAssayData(object = cellData,
                          assay = "RNA",
                          slot = "data")[cellMarkers$required, , drop = FALSE]
  cannotHave = GetAssayData(object = cellData,
                            assay = "RNA",
                            slot = "data")[cellMarkers$cannotHave, , drop = FALSE]
  optional = GetAssayData(object = cellData,
                          assay = "RNA",
                          slot = "data")[cellMarkers$optional, , drop = FALSE]
  
  optional2 = GetAssayData(object = cellData,
                          assay = "RNA",
                          slot = "data")[cellMarkers$optional2, , drop = FALSE]
  
  optional3 = GetAssayData(object = cellData,
                           assay = "RNA",
                           slot = "data")[cellMarkers$optional3, , drop = FALSE]
  
  for (i in 1:dim(cellData)[2]) {
    cell <- colnames(cellData)[i]
    if (length(optional) == 0 &
        length(required) > 0 & length(cannotHave) > 0) {
      if (all(required[, cell] > 0) & sum(cannotHave[, cell]) == 0) {
        cellData[[cellName]][i,] = 'yes'
      }
    } else if (length(required) == 0) {
      if (sum(cannotHave[, cell]) == 0 & sum(optional[, cell]) > 0) {
        cellData[[cellName]][i,] ='yes'
      }
    } else if(length(optional3) > 0){
      if(all(required[, cell] > 0) & sum(cannotHave[, cell]) == 0 & sum(optional[, cell]) > 0
         & sum(optional2[, cell]) > 0 & sum(optional3[, cell]) > 0){
        cellData[[cellName]][i,] ='yes'
      }
    } else if(length(optional2) > 0){
      if(all(required[, cell] > 0) & sum(cannotHave[, cell]) == 0 & sum(optional[, cell]) > 0
         & sum(optional2[, cell]) > 0){
        cellData[[cellName]][i,] ='yes'
      }
    }
    else if (all(required[, cell] > 0) &
             sum(cannotHave[, cell]) == 0 &
             sum(optional[, cell]) > 0) {
      cellData[[cellName]][i,] = 'yes'
    }
    percentDone <- (i / dim(cellData)[2]) * 100
    print(paste(as.character(percentDone), 'percent done'))
  }
  cellData
}


#Use labelCells function to label cells, if cells are ambiguous they are named indistinguishable
#Could loop through these but they take a while to run
dengDatCapInt <- labelCells(T_cell_markers, cellData = dengDatCapInt, cellName = 'T_cell')
subset(dengDatCapInt, T_cell == '0')[['RNA']]$data[c('CD3E', 'CD2', 'NCR1', 'NCAM1', 'FOXP3'),]

dengDatCapInt <- labelCells(CD4_markers,  cellData = dengDatCapInt, cellName = 'CD4_T_cells')
dengDatCapInt <- labelCells(CD8_markers, cellData = dengDatCapInt, cellName = 'CD8_T_cells')
dengDatCapInt <- labelCells(CD4_reg_markers, cellData = dengDatCapInt, cellName = 'CD4_reg')
dengDatCapInt <- labelCells(CD8_reg_markers, cellData = dengDatCapInt, cellName = 'CD8_reg')
dengDatCapInt <- labelCells(Activated_markers, cellData = dengDatCapInt, cellName = 'Activated_t_cells')
dengDatCapInt <- labelCells(Exhausted_markers, cellData = dengDatCapInt, cellName = 'Exhausted_t_cells')
dengDatCapInt <- labelCells(Memory_markers, cellData = dengDatCapInt, cellName = 'Memory_t_cells')
dengDatCapInt <- labelCells(Memory_markers_cd4, cellData = dengDatCapInt, cellName = 'Memory_CD4_t_cells')
dengDatCapInt <- labelCells(Memory_markers_cd8, cellData = dengDatCapInt, cellName = 'Memory_CD8_t_cells')

#This is long because some columns exclude eachother (cannot be both activated and exhausted), but
#a cell could be both a cd4 cell and an activated cell, so I first rule out cells that have 
#multiple categories excluding eachother by summing number of mutually exclusive 'yes's 
cellAssignments <- dengDatCapInt[[]] %>% select(CD4_T_cells, CD8_T_cells, CD4_reg, CD8_reg, Activated_t_cells, Exhausted_t_cells, 
                                                Memory_t_cells, Memory_CD4_t_cells, Memory_CD8_t_cells) %>% 
  mutate(numAssignments = rowSums(select(.,CD4_reg:Memory_CD8_t_cells) == 'yes'),
         cellType = case_when(numAssignments == 0 & CD4_T_cells == 'yes' & CD8_T_cells == 'no' ~ 'CD4_t_cell',
                              numAssignments == 0 & CD8_T_cells == 'yes' & CD4_T_cells == 'no' ~ 'CD8_t_cell',
                              numAssignments == 0 & CD8_T_cells == 'yes' & CD4_T_cells == 'yes' ~ 'indistinguishable',
                              numAssignments == 0 & CD8_T_cells == 'no' & CD4_T_cells == 'no' ~ 'unknown',
                              numAssignments >= 2 ~ 'indistinguishable',
                              numAssignments == 1 & CD4_reg == 'yes' ~ 'CD4_reg',
                              numAssignments == 1 & CD8_reg == 'yes' ~ 'CD8_reg',
                              numAssignments == 1 & Activated_t_cells == 'yes' ~ 'Activated_t_cell',
                              numAssignments == 1 & Exhausted_t_cells == 'yes' ~ 'Exhausted_t_cell',
                              numAssignments == 1 & Memory_t_cells == 'yes' ~ 'Memory_t_cell',
                              numAssignments == 1 & Memory_CD4_t_cells == 'yes' ~ 'Memory_CD4_t_cell',
                              numAssignments == 1 & Memory_CD8_t_cells == 'yes' ~ 'Memory_CD8_t_cell'))
dengDatCapInt$cellType <- cellAssignments$cellType

#Add labels back to original deng data
dengDatCapInt$cell <- rownames(dengDatCapInt[[]])
dengDatCapInt_dat <- dengDatCapInt[[]][c('cellType', 'cell')]

dengDatInt$cellType <- NULL
dengDatInt$cell <- rownames(dengDatInt[[]])
updatedMeta <- left_join(dengDatInt[[]], dengDatCapInt_dat, by = 'cell')
dengDatInt[[]] <- updatedMeta
DimPlot(dengDatInt, group.by = 'cellType')

#### Can start here if loading in integrated data ####
DimPlot(dengDatInt, reduction = "umap", group.by = c("seurat_clusters"), label = T)

table(dengDatInt$X3mo.PET.CT)
DimPlot(dengDatInt, reduction = "umap", group.by = c("X3mo.PET.CT"), label = T)

VlnPlot(dengDatInt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)

FeaturePlot(dengDatInt, features = 'CD4') 
FeaturePlot(dengDatInt, features = 'CD8A')
VlnPlot(dengDatInt, features = 'CD4', pt.size = 0)
VlnPlot(dengDatInt, features = 'CD8A', pt.size = 0)



#### Begin integration of our data w Deng ####
T_cells_noCC <- LoadSeuratRds('mystore/cartdata/data/T_cells_noCC.rds')

#Transfer labels
cellType.anchors <- FindTransferAnchors(reference = dengDatInt, query = T_cells_noCC, dims = 1:30)

#May need 20+ cores for this 
predictions <- TransferData(anchorset = cellType.anchors, refdata = dengDatInt$cellType, dims = 1:30)
T_cells_noCC <- AddMetaData(T_cells_noCC, metadata = predictions)

DimPlot(T_cells_noCC, group.by = 'predicted.id')
table(T_cells_noCC$predicted.id)
table(subset(T_cells_noCC, CD_pred == 'CD4')$predicted.id)

table(subset(T_cells_noCC, hypoxia == 'HH')$predicted.id)

T_cells_noCC[[]] %>% group_by(day, hypoxia, predicted.id) %>% summarise(count = n()) %>% 
  mutate(dayHyp = paste(day, hypoxia, sep = '_')) %>% ggplot(aes(x = dayHyp, y = count, fill = predicted.id)) +
  geom_bar(stat = 'identity')
#Reduce object sizes for actual integration
dengDatInt <- subset(dengDatInt, cellType != 'unknown' & cellType != 'indistinguishable')
T_cells_noCC <- DietSeurat(T_cells_noCC, 'counts')
dengDatInt <- DietSeurat(dengDatInt, 'counts')

#Integrate our data
allDatMerge <- merge(x = dengDatInt, y = T_cells_noCC, add.cell.ids = as.character(seq(1,2)),
                      merge.data = TRUE)

allDatMerge <- NormalizeData(allDatMerge)
allDatMerge <- FindVariableFeatures(allDatMerge)
allDatMerge <- ScaleData(allDatMerge)
allDatMerge <- RunPCA(allDatMerge)

#setting idents to all stops VlnPlot from grouping violins by seurt cluster
Idents(allDatMerge) <- 'all'
VlnPlot(allDatMerge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)

allDatInt <- IntegrateLayers(object = allDatMerge, method = CCAIntegration, orig.reduction = "pca",
                             new.reduction = "integrated.cca", verbose = FALSE)

# re-join layers after integration
allDatInt[["RNA"]] <- JoinLayers(allDatInt[["RNA"]])

ElbowPlot(allDatInt, ndims = 30)
allDatInt <- FindNeighbors(allDatInt, reduction = "integrated.cca", dims = 1:20)
allDatInt <- FindClusters(allDatInt, resolution = 1)
allDatInt <- RunUMAP(allDatInt, dims = 1:20, reduction = "integrated.cca")

DimPlot(allDatInt, reduction = "umap", group.by = c("orig.ident")) +
  scale_color_manual(labels = c("Our Data", "Deng Data"), values = c("#4ce6d3", "#d16251"))

DimPlot(allDatInt, reduction = "umap", group.by = c("X3mo.PET.CT"))

####DEG comparison between groups####

#Look at proportion of celltype in each group
table(dengDatInt$cellType, dengDatInt$X3mo.PET.CT) %>% as.data.frame() %>% ggplot(aes(x = Var2, y = Freq, fill = Var1))+
  geom_bar(stat = 'identity')

#Look at p or q values and compare to Deng
table(dengDatInt$cellType, dengDatInt$X3mo.PET.CT)%>% as.data.frame() %>% mutate(CD8_Mem = case_when(Var1 == 'Memory_CD8_t_cell' ~ 'CD8Mem',
                                                                                                     Var1 != 'Memory_CD8_t_cell' ~ 'other')) %>% 
  group_by(CD8_Mem, Var2) %>% summarise(count = sum(Freq)) %>% filter(Var2 == 'CR' | Var2 == 'PD') %>% 
  pivot_wider(names_from = Var2, values_from = count) %>% ungroup() %>% select(c(CR, PD)) %>% as.matrix() %>% 
  fisher.test()

#Look at DEGs from responders compared to our data
#Doesn't seem like anything is significant in pseudobulk
dengCD4 <- subset(dengDatInt, cellType == 'CD4_t_cell' | cellType == 'CD4_reg' | cellType == 'Memory_CD4_t_cell')
dengCD4 <- subset(dengCD4, X3mo.PET.CT != 'NE')
dengCD4[[]] <- dengCD4[[]] %>% mutate(CR_vs = ifelse(X3mo.PET.CT == 'CR', 'CR', 'PRPD'))

#Can skip this to seurat DEGs 
dengCD4_bulk <- getPseudoBulkObject(dengCD4, c('Gender', 'Stage', 'CR_vs'))
dengCD4_bulk <- DESeq(dengCD4_bulk)
cr_vs_prpd_cd4 <- lfcShrink(dengCD4_bulk, coef ='CR_vs_PRPD_vs_CR')
cr_vs_prpd_cd4[order(cr_vs_prpd_cd4$padj),]  %>% as.data.frame() %>% head(n = 50) 
cr_vs_prpd_cd4$gene <- rownames(cr_vs_prpd_cd4)

#Seurat comparison
Idents(dengCD4) <- 'CR_vs'
markersCD4 <- FindMarkers(dengCD4, ident.1 = 'CR')
topCD4Markers <- head(markersCD4, n = 60)
topCD4Markers$gene <- rownames(topCD4Markers)
#only keep genes present in t cell data
topCD4Markers <- topCD4Markers[topCD4Markers$gene %in% rownames(T_cells_noCC[['RNA']]$counts[]),]

#CD8
dengCD8 <- subset(dengDatInt, cellType == 'CD8_t_cell' | cellType == 'CD8_reg' | cellType == 'Memory_CD8_t_cell')
dengCD8 <- subset(dengCD8, X3mo.PET.CT != 'NE')
dengCD8[[]] <- dengCD8[[]] %>% mutate(CR_vs = ifelse(X3mo.PET.CT == 'CR', 'CR', 'PRPD'))
Idents(dengCD8) <- 'CR_vs'

markersCD8 <- FindMarkers(dengCD8, ident.1 = 'CR')
topCD8Markers <- head(markersCD8, n = 60)
topCD8Markers$gene <- rownames(topCD8Markers)

#only keep genes present in t cell data
topCD8Markers <- topCD8Markers[topCD8Markers$gene %in% rownames(T_cells_noCC[['RNA']]$counts[]),]

# dengCD4_extra <- subset(dengCD4, cellType == 'CD4_t_cell')
# Idents(dengCD4_extra) <- 'CR_vs'
# markerTest <- FindMarkers(dengCD4_extra, ident.1 = 'CR')
# topDengMarkers <- head(markerTest, n = 50)
# topDengMarkers$gene <- rownames(topDengMarkers)

T_cell_CD4 <- subset(T_cells_noCC, CD_pred == 'CD4')
Idents(T_cell_CD4) <- 'hypoxia'
hypoxiaMarkers <- FindMarkers(T_cell_CD4, ident.1 = 'HH', ident.2 = 'NN', features = topCD4Markers$gene)
hypoxiaMarkers$gene <- rownames(hypoxiaMarkers)

combinedMarkers <- left_join(topCD4Markers, hypoxiaMarkers, by = 'gene')
combinedMarkers <- filter(combinedMarkers, !is.na(avg_log2FC.y))
ggplot(combinedMarkers, aes(x = avg_log2FC.x, y = avg_log2FC.y)) +
  geom_point() +
  xlab('Deng Data CR vs PRPD LFC')+
  ylab('Our Data HH vs NN Log Fold Change') +
  geom_text_repel(aes(label=gene))

T_cell_CD8 <- subset(T_cells_noCC, CD_pred == 'CD8')
Idents(T_cell_CD8) <- 'hypoxia'
hypoxiaMarkers <- FindMarkers(T_cell_CD8, ident.1 = 'HH', ident.2 = 'NN', features = topCD8Markers$gene)
hypoxiaMarkers$gene <- rownames(hypoxiaMarkers)

combinedMarkers <- left_join(topCD8Markers, hypoxiaMarkers, by = 'gene')
combinedMarkers <- filter(combinedMarkers, !is.na(avg_log2FC.y))
ggplot(combinedMarkers, aes(x = avg_log2FC.x, y = avg_log2FC.y)) +
  geom_point() +
  xlab('Deng Data CR vs PRPD LFC')+
  ylab('Our Data HH vs NN Log Fold Change') +
  geom_text_repel(aes(label=gene)) + 
  geom_abline(intercept = 0, slope = 1)


CD4_T_cell_bulk <- getPseudoBulkObject(T_cell_CD4, c("hypoxia", "donor_id", 'CAR', 'day', 'carExpression'))
CD4_T_cell_bulk <- DESeq(CD4_T_cell_bulk)
NN_vs_HH_cd4 <- lfcShrink(CD4_T_cell_bulk, coef ='hypoxia_NN_vs_HH')
NN_vs_HH_cd4$gene <- rownames(NN_vs_HH_cd4)

combinedMarkers <- left_join(as.data.frame(cr_vs_prpd_cd4), as.data.frame(NN_vs_HH_cd4), by = 'gene')


#Cross reference DEG's from our data with Deng's genes comparing response and reactions
CD4_CR_vs_PRPD_deng <- read_csv("mystore/cartdata/publicData/CD4_CR_vs_PRPD_deng.csv")
colnames(CD4_CR_vs_PRPD_deng) <- as.character(CD4_CR_vs_PRPD_deng[3,])
CD4_CR_vs_PRPD_deng <- CD4_CR_vs_PRPD_deng[4:nrow(CD4_CR_vs_PRPD_deng),]

CD8_CR_vs_PRPD_deng <- read_csv("mystore/cartdata/publicData/CD8_CR_vs_PRPD_deng.csv")
colnames(CD8_CR_vs_PRPD_deng) <- as.character(CD8_CR_vs_PRPD_deng[3,])
CD8_CR_vs_PRPD_deng <- CD8_CR_vs_PRPD_deng[4:nrow(CD8_CR_vs_PRPD_deng),]

FeaturePlot(allDatInt, features = 'PGK1')
#Run findCommonGenes function from carVscarnegVsuntrasnduced.R, line to get mostCommonHHNNGenes is right beneath
intersect(names(mostCommonHHNNGenes), CD4_CR_vs_PRPD_deng$Gene)
intersect(names(mostCommonHHNNGenes), CD8_CR_vs_PRPD_deng$Gene)
cd4Deng <- subset(dengDatInt, CD4 > 0 & CD8A == 0)
cd8Deng <- subset(dengDatInt, CD4 == 0 & CD8A > 0)

#Look at intersection genes vs response
VlnPlot(dengDatInt, group.by = 'X3mo.PET.CT', features = 'PGK1', pt.size = 0, assay = 'RNA') +
  ggtitle('PGK1 Deng Data')
AverageExpression(cd4Deng, group.by = 'X3mo.PET.CT', features = 'PGK1')
AverageExpression(cd8Deng, group.by = 'X3mo.PET.CT', features = 'PGK1')

VlnPlot(dengDatInt, group.by = 'X3mo.PET.CT', features = 'GPI', pt.size = 0, assay = 'RNA') +
  ggtitle('GPI Deng Data')

#Look at ICANS and CRS associated genes
ICANS_deng <- read_csv("mystore/cartdata/publicData/IcansAssociatedGenesDeng.csv")
colnames(ICANS_deng) <- as.character(ICANS_deng[3,])
ICANS_deng <- ICANS_deng[4:nrow(ICANS_deng),]

intersect(names(mostCommonHHNNGenes), ICANS_deng$Gene)

CRS_deng <- read_csv("mystore/cartdata/publicData/CRSHighVsLowGenesDeng.csv")
colnames(CRS_deng) <- as.character(CRS_deng[3,])
CRS_deng <- CRS_deng[4:nrow(CRS_deng),]

intersect(names(mostCommonHHNNGenes), CRS_deng$Gene)

#Use t_cell_nocc from cellCycleRegressed.R
cd4Cells <- subset(T_cells_noCC, CD_pred=='CD4')
bulk_CD4 <- getPseudoBulkObject(cd4Cells, c("CAR","hypoxia", "Phase", "donor_id", "carExpression", "day"))
bulk_CD4 <- DESeq(bulk_CD4)

bulk_CD4_hypoxia <- lfcShrink(bulk_CD4, coef = "hypoxia_NN_vs_HH")
CD4_CR_vs_PRPD_genes <- CD4_CR_vs_PRPD_deng$Gene[!is.na(CD4_CR_vs_PRPD_deng$Gene)]
CD4_CR_vs_PRPD_gene_lfcs <- as.data.frame(bulk_CD4_hypoxia)[CD4_CR_vs_PRPD_genes,]
ggplot(CD4_CR_vs_PRPD_gene_lfcs, aes(x = log2FoldChange, y = -log10(padj)))+
  geom_point()

#Need to subset to cd4 cells
# dengDatInt$Age <- as.numeric(dengDatInt$Age)
# 
# dengDatInt[[]] <- dengDatInt[[]] %>% mutate(responseSummarized = case_when(X3mo.PET.CT == 'CR'~ 'CR',
#                                                          X3mo.PET.CT == 'PD'~ 'PRPD',
#                                                          X3mo.PET.CT == 'PR'~ 'PRPD',
#                                                          X3mo.PET.CT == 'NE'~ 'NE'))
# dengBulk <- getPseudoBulkObject(dengDatInt, c("Histology","Stage", "responseSummarized"))
# dengBulk <- DESeq(dengBulk)
# dengBulk_CR_PRPD <- lfcShrink(dengBulk, coef = "responseSummarized_PRPD_vs_CR")
# Deng_CR_vs_PRPD_gene_lfcs <- as.data.frame(dengBulk_CR_PRPD)[CD4_CR_vs_PRPD_genes,]
# ggplot(Deng_CR_vs_PRPD_gene_lfcs, aes(x = log2FoldChange, y = -log10(padj)))+
#   geom_point()

#Check that activated cells cross over completely w other cell types
rna <- dengDatCapInt[['RNA']]$data
activated <- rna[c('CD3E', 'CD2', 'CD28', 'IL2RA','CD69','ICOS','TNFRSF4',
                   'TNFRSF9','CD27','IL2','TNF','IFNG', 'CD8A','CD8B','CD4',
                   'NCR1', 'NCAM1'),]
activated_df <- as.data.frame(t(as.data.frame(activated)))
activated_df <- activated_df %>% mutate(optional = (IL2RA + CD69 + ICOS + TNFRSF4 + TNFRSF9 + CD27),
                                        optional2 = (IL2 + TNF + IFNG),
                                        optional3 = (CD8A + CD8B + CD4))
activeCells <- rownames(activated_df[activated_df$CD3E >0 & activated_df$CD2 > 0 & activated_df$CD28 > 0 & activated_df$optional > 0 &
               activated_df$optional2 > 0 & activated_df$optional3 > 0 & activated_df$NCR1 == 0 & activated_df$NCAM1 == 0,])

activeCellsAll <- dengDatCapInt[['RNA']]$data[,activeCells]

#CD4 and 8 reg genes
regGenes <- activeCellsAll[c('FOXP3', 'CD3E', 'CD4', 'CTLA4','IL2RA', 'CD8A', 'CD8B'),]

