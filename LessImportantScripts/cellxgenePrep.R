#Prepare pseudobulk object for cellxgene app on beagle
library(dplyr)
library(SeuratDisk)
library(reticulate)
library(sceasy)
library(celldex)

#Source fucntions
source("mystore/cartdata/scripts/CarT_project_functions.R")

#Make UMAPs plot properly
options(bitmapType="cairo") 

T_cells_noCC <- LoadSeuratRds('mystore/cartdata/data/T_cells_noCC.rds')
T_cells <-  LoadSeuratRds("./mystore/cartdata/data/tcell_deseqdat_singlets.rds")

bulk <- AggregateExpression(
  T_cells,
  return.seurat = T,
  assays = "RNA",
  group.by = c("CAR","hypoxia", "Phase", "donor_id", "day", "CD_pred", "CAR_pred"))

# Add in number of cells by sample and celltype
n_cells <- T_cells[[]] %>% dplyr::select(all_of(c("CAR","hypoxia", "Phase", "donor_id", "day", "CD_pred", "CAR_pred"))) %>% 
  dplyr::group_by(across(all_of(c("CAR","hypoxia", "Phase", "donor_id", "day", "CD_pred", "CAR_pred")))) %>% 
  dplyr::count() 

bulk[[]] <- n_cells

bulk <- prepUMAP(bulk)
DimPlot(bulk)

#There are some issue with seurat versions and SeuratDisk, the following code fixes  HDF5-API Errors:
bulk[["RNA3"]] <- as(object = bulk[["RNA"]], Class = "Assay")
DefaultAssay(bulk) <- "RNA3"
bulk[["RNA"]] <- NULL
bulk[["RNA"]] <- bulk[["RNA3"]]
DefaultAssay(bulk) <- "RNA"
bulk[["RNA3"]] <- NULL

bulk@assays$RNA$scale.data <- NULL

SaveH5Seurat(bulk, filename = "mystore/cartdata/data/bulkAll.h5Seurat")
Convert("mystore/cartdata/data/bulkAll.h5Seurat", dest = "h5ad")



#Prep sc data, rather than pseudobulk
#There are some issue with seurat versions and SeuratDisk, the following code fixes  HDF5-API Errors:
T_cells[["RNA3"]] <- as(object = T_cells[["RNA"]], Class = "Assay")
DefaultAssay(T_cells) <- "RNA3"
T_cells[["RNA"]] <- NULL
T_cells[["RNA"]] <- T_cells[["RNA3"]]
DefaultAssay(T_cells) <- "RNA"
T_cells[["RNA3"]] <- NULL

T_cells@assays$RNA$scale.data <- NULL

SaveH5Seurat(T_cells, filename = "mystore/cartdata/data/TCellAll.h5Seurat")
Convert("mystore/cartdata/data/TCellAll.h5Seurat", dest = "h5ad")

#Prepare t cells (with cell cycle regressed out) for anndata and scArches

#Start with just a subset of cells to test scArches
T_cells_noCC_sub <- subset(T_cells_noCC, CAR_pred == 1 & CD_pred == 'CD4' & day == 'D7' & hypoxia == 'HH')

T_cells_noCC_sub[["RNA3"]] <- as(object = T_cells_noCC_sub[["RNA"]], Class = "Assay")
DefaultAssay(T_cells_noCC_sub) <- "RNA3"
T_cells_noCC_sub[["RNA"]] <- NULL
T_cells_noCC_sub[["RNA"]] <- T_cells_noCC_sub[["RNA3"]]
DefaultAssay(T_cells_noCC_sub) <- "RNA"
T_cells_noCC_sub[["RNA3"]] <- NULL

T_cells_noCC_sub@assays$RNA$scale.data <- NULL

SaveH5Seurat(T_cells_noCC_sub, filename = "mystore/cartdata/data/TCell_subset.h5Seurat")
Convert("mystore/cartdata/data/TCell_subset.h5Seurat", dest = "h5ad")

#All CD4 D7 cells for celltypist
T_cells_noCC_cd4d7 <- subset(T_cells_noCC, CAR_pred == 1 & CD_pred == 'CD4' & day == 'D7')

T_cells_noCC_cd4d7[["RNA3"]] <- as(object = T_cells_noCC_cd4d7[["RNA"]], Class = "Assay")
DefaultAssay(T_cells_noCC_cd4d7) <- "RNA3"
T_cells_noCC_cd4d7[["RNA"]] <- NULL
T_cells_noCC_cd4d7[["RNA"]] <- T_cells_noCC_cd4d7[["RNA3"]]
DefaultAssay(T_cells_noCC_cd4d7) <- "RNA"
T_cells_noCC_cd4d7[["RNA3"]] <- NULL

T_cells_noCC_cd4d7@assays$RNA$scale.data <- NULL

SaveH5Seurat(T_cells_noCC_cd4d7, filename = "mystore/cartdata/data/TCell_cd4d7.h5Seurat")
Convert("mystore/cartdata/data/TCell_cd4d7.h5Seurat", dest = "h5ad")



#Prep reference immune data for scarches. Doesn't have raw counts so not good for expimap?
azimuthData <- readRDS('mystore/cartdata/data/azimuthHumanPBMC.Rds')
azimuthData[[]]
table(azimuthData$celltype.l2)

azimuthData[["RNA3"]] <- as(object = azimuthData[["RNA"]], Class = "Assay")
DefaultAssay(azimuthData) <- "RNA3"
azimuthData[["RNA"]] <- NULL
azimuthData[["RNA"]] <- azimuthData[["RNA3"]]
DefaultAssay(azimuthData) <- "RNA"
azimuthData[["RNA3"]] <- NULL

azimuthData@assays$RNA$scale.data <- NULL

SaveH5Seurat(azimuthData, filename = "mystore/cartdata/data/TCell_subset.h5Seurat")
Convert("mystore/cartdata/data/TCell_subset.h5Seurat", dest = "h5ad")

#Get CAR Positives for scarches
CarPositives <- subset(T_cells_noCC, CAR != 'untransduced' & CAR_pred == 1)

CarPositives[["RNA3"]] <- as(object = CarPositives[["RNA"]], Class = "Assay")
DefaultAssay(CarPositives) <- "RNA3"
CarPositives[["RNA"]] <- NULL
CarPositives[["RNA"]] <- CarPositives[["RNA3"]]
DefaultAssay(CarPositives) <- "RNA"
CarPositives[["RNA3"]] <- NULL

CarPositives@assays$RNA$scale.data <- NULL

SaveH5Seurat(CarPositives, filename = "mystore/cartdata/data/TCell_CARPos.h5Seurat")
Convert("mystore/cartdata/data/TCell_CARPos.h5Seurat", dest = "h5ad")


# T_cell_deseq <- getPseudoBulkObject(T_cells, c("CAR","hypoxia", "Phase", "donor_id", "day", "CD_pred", "CAR_pred"))
# 
# T_cell_deseq <- vst(T_cell_deseq)
# getPCAGrid(T_cell_deseq, c("hypoxia", "CAR", "Phase", "donor_id"), rowNum = 2)
# 
# pseudoBulk_seurat <- CreateSeuratObject(counts = counts(T_cell_deseq, normalized = FALSE), project = "T_cell_deseq")
# 
# pseudoBulk_seurat <- AddMetaData(pseudoBulk_seurat, colData(T_cell_deseq)$CAR, col.name = "CAR")
# pseudoBulk_seurat <- AddMetaData(pseudoBulk_seurat, colData(T_cell_deseq)$hypoxia, col.name = "hypoxia")
# pseudoBulk_seurat <- AddMetaData(pseudoBulk_seurat, colData(T_cell_deseq)$day, col.name = "day")
# pseudoBulk_seurat <- AddMetaData(pseudoBulk_seurat, colData(T_cell_deseq)$Phase, col.name = "CD_pred")
# pseudoBulk_seurat <- AddMetaData(pseudoBulk_seurat, colData(T_cell_deseq)$Phase, col.name = "CAR_pred")
# 
# pseudoBulk_seurat <- NormalizeData(pseudoBulk_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
# pseudoBulk_seurat <- FindVariableFeatures(pseudoBulk_seurat, selection.method = "vst", nfeatures = 2000)
# all.genes <- rownames(pseudoBulk_seurat)
# pseudoBulk_seurat <- ScaleData(pseudoBulk_seurat, features = all.genes)
# pseudoBulk_seurat <- RunPCA(pseudoBulk_seurat, features = VariableFeatures(object = pseudoBulk_seurat))
# ElbowPlot(pseudoBulk_seurat, ndims = 30)
# 
# pseudoBulk_seurat <- FindNeighbors(pseudoBulk_seurat, dims = 1:25)
# pseudoBulk_seurat <- FindClusters(pseudoBulk_seurat, resolution = 0.5)
# pseudoBulk_seurat <- RunUMAP(pseudoBulk_seurat, dims = 1:25)
# DimPlot(pseudoBulk_seurat, group.by = "Phase")
