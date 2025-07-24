#Load packages
library(Seurat)
library(dplyr)
library(SingleCellExperiment)
library(scDblFinder)
library(ggplot2)
library(DESeq2)
library(pheatmap)
library(tibble)
library(apeglm)
library(tidyr)
library(gridExtra)

#Make plots work properly
options(bitmapType="cairo") 

#Load in some functions for project
source("mystore/cartdata/scripts/CarT_project_functions.R")

#Load here to run everything from the beginning
adata <- LoadSeuratRds("/proj/nobackup/hpc2nstor2024-027/cart/allcount_andvirus.Rdata")
adata$cell = rownames(adata[[]])

#Load here to read in t cells with doublets still present
tCell_adata <- LoadSeuratRds("./mystore/cartdata/data/tCell_adata_with_doublets.rds")

#Load here for t cells with only singlets  (vireo and scdblfinder used), and  no d0 hypoxia groups. 
deseqDat <- LoadSeuratRds("./mystore/cartdata/data/tcell_deseqdat_singlets.rds")

#### Label Donors ####
#Load in donor id for four donors
DonorIDs_four <- read.delim("/pfs/stor10/users/home/m/mb223/mystore/dataset/241211_cart_parse/cellSNPFiles/vireoOutFour/donor_ids.tsv", 
                       row.names=NULL)

#Match donor ids with R row names (shouldn't be necessary for actual data
DonorIDs_four = DonorIDs_four[,c("cell", "donor_id")]
DonorIDs_four$donor_id = factor(DonorIDs_four$donor_id) 

#Merge donor ids to seurat object
adata[[]] = left_join(adata[[]], DonorIDs_four, by = "cell")

adata <- subset(adata, subset = nFeature_RNA < 10000 & nCount_RNA < 50000 & percent.mt < 20)

# Data preprocessing
adata <- NormalizeData(adata, normalization.method = "LogNormalize")
adata <- FindVariableFeatures(adata, selection.method = "vst", nfeatures = 2000)

adata <- ScaleData(adata)
adata <- RunPCA(adata)

DimHeatmap(adata, dims = 1:3, cells = 500, balanced = TRUE)
ElbowPlot(adata, ndims = 30)

adata <- FindNeighbors(adata, dims = 1:20)
adata <- FindClusters(adata, resolution = 0.5)
adata <- RunUMAP(adata, dims = 1:20)

#Load cell type from TCR_Analysis script
#adata[[]] = left_join(adata[[]], cell_types, by = "cell")
DimPlot(adata, reduction = "umap", label = T)
#table(adata$cell_type)

#Subset to just t cells
cd3GeneList = list(c("CD3G", "CD3D", "CD3E"))
adata <- AddModuleScore(object = adata, features = cd3GeneList, name = "cd3Genes") 
VlnPlot(adata, features = "cd3Genes1", pt.size = 0)+
  stat_summary(fun.y = median, geom='point', size = 2)+
  geom_hline(yintercept = 0)

tCell_adata <- subset(adata, seurat_clusters %in% c(0:7, 9:10))

# T cell preprocessing
tCell_adata <- NormalizeData(tCell_adata, normalization.method = "LogNormalize")
tCell_adata <- FindVariableFeatures(tCell_adata, selection.method = "vst", nfeatures = 2000)

#Look at highly variable features for covariate analysis later
top300 <- head(VariableFeatures(tCell_adata), 300)

tCell_adata <- ScaleData(tCell_adata)
tCell_adata <- RunPCA(tCell_adata)

DimHeatmap(tCell_adata, dims = 1:3, cells = 500, balanced = TRUE)
ElbowPlot(tCell_adata, ndims = 30)

tCell_adata <- FindNeighbors(tCell_adata, dims = 1:25)
tCell_adata <- FindClusters(tCell_adata, resolution = 0.5)
tCell_adata <- RunUMAP(tCell_adata, dims = 1:20)

DimPlot(tCell_adata, reduction = "umap", label = T)

#Look for doublets.
#Sample is used to indicate specific captures in multiplexed setups

#Label cells based on capture
tCell_adata$capture <- factor(substr(tCell_adata$cell, nchar(tCell_adata$cell)-1,
                                     nchar(tCell_adata$cell)))

sce = as.SingleCellExperiment(tCell_adata)
sce <- scDblFinder(sce,
                   nfeatures = 750,
                   samples = "capture",
                   dbr.sd=1)

table(sce$scDblFinder.class)
tCell_adata$DblFinderDoublets = sce$scDblFinder.class
table(tCell_adata$donor_id, tCell_adata$DblFinderDoublets)

tCell_adata[[]] <- tCell_adata[[]] %>% mutate(doubletAgreement = case_when(
  donor_id == "doublet" & DblFinderDoublets == "doublet" ~ "agreed_doublet",
  donor_id != "doublet" & DblFinderDoublets == "doublet" ~ "scdblDoublet",
  donor_id == "doublet" & DblFinderDoublets != "doublet" ~ "vireoDoublet",
  donor_id != "doublet" & DblFinderDoublets != "doublet" ~ "singlet"))
ggplot(tCell_adata[[]], aes(x = factor(doubletAgreement, levels = c("agreed_doublet",
                                                                    "scdblDoublet",
                                                                    "vireoDoublet",
                                                                    "singlet")),fill = donor_id))+
  geom_bar()+
  xlab("Group")
VlnPlot(tCell_adata, features = "nCount_RNA", group.by = "doubletAgreement", pt.size = 0) 
tCell_adata$doubletAgreement = factor(tCell_adata$doubletAgreement, levels = c("agreed_doublet",
                                                                               "scdblDoublet",
                                                                               "vireoDoublet",
                                                                               "singlet"))
tCell_adata[[]] %>% group_by(seurat_clusters, doubletAgreement) %>% summarise(counts = n()) %>% 
  ggplot(aes(x = seurat_clusters, y = counts, fill = doubletAgreement))+
  geom_bar(stat = "identity")

VlnPlot(tCell_adata, features = "nCount_RNA")
tCell_adata[[]] %>% group_by(seurat_clusters, doubletAgreement) %>% summarise(counts = n())
#Remove anything labelled a doublet by vireo or scdblfinder 
tCell_adata = subset(tCell_adata, doubletAgreement == "singlet" & donor_id != "unassigned")

#Deseq2 
deseqDat = subset(tCell_adata, hypoxia != "d0" & hypoxia != "d0+il2")

#### Start here with deseqDat loaded above ####
#Relook at umap structure
deseqDat <- FindVariableFeatures(deseqDat)
deseqDat <- ScaleData(deseqDat)
deseqDat <- RunPCA(deseqDat)
DimHeatmap(deseqDat, dims = 1:3, cells = 500, balanced = TRUE)
ElbowPlot(deseqDat, ndims = 30)
deseqDat <- FindNeighbors(deseqDat, dims = 1:25)
deseqDat <- FindClusters(deseqDat, resolution = 0.5)
deseqDat <- RunUMAP(deseqDat, dims = 1:20)
DimPlot(deseqDat, group.by = "hypoxia")

deseqDat$hypoxia = factor(deseqDat$hypoxia)
deseqDat$CAR = factor(deseqDat$CAR)
deseqDat$donor_id = factor(deseqDat$donor_id)
deseqDat$day = factor(deseqDat$day)



####This was an early attempt to see overall effects of covariates on data using MANOVA. ####
#Attempt to see effects of covariates with linear model
#Pull out expression data. Subset for figuring out workflow

#Loop through manovas and plot distribution of pillai's
#Try parametric bootstrap with normal distributions
# datForModel <- t(tCell_adata[["RNA"]]$data[])
# bootstrapManova <- function(){
#   keepCols <- sample.int(ncol(datForModel), 300)
#   datForModel1 <- datForModel[,keepCols] #Start with 300 genes
#   datForModel1 <- datForModel1[,(!colSums(datForModel1) == 0)]
#   linearDependentCols <- caret::findLinearCombos(as.matrix(datForModel1))
#   if(!is.null(linearDependentCols$remove)){
#     datForModel1 <- datForModel1[,-linearDependentCols$remove]
#   }
#   manModel <- manova(as.matrix(datForModel1) ~ tcellPool1$Phase + tcellPool1$hypoxia)
#   testStats = summary(manModel)$stats[,"Pillai"] #Differs from car::Anova sometimes
#   data.frame(groups = c("Phase", "hypoxia"), testStat = testStats[-length(testStats)])
# }
# modelRuns <- lapply(seq_len(100), function(x) bootstrapManova())
# modelRuns = do.call(rbind, modelRuns)
# ggplot(modelRuns, aes(x = groups, y = testStat))+
#   geom_boxplot()
# 
# ggplot(modelRuns, aes(testStat, fill = groups)) + geom_density(alpha = 0.2)

#Paper about annotating t cells https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2023.1306169/full

