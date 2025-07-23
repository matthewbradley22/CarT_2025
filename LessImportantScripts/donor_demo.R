#example script for mapping donor data on to UMAP.
#Only using 1/8 of cells

library(Seurat)
library(dplyr)
library(SingleCellExperiment)
#Have to install this with remotes::install_github("cvarrichio/Matrix.utils")
library(Matrix.utils)
library(data.table)
library(DESeq2)
library(pheatmap)
library(edgeR)
library(car)
library(tibble)
library(ggplot2)
library(EnhancedVolcano)
library(cowplot)
library(tidyr)
library(caret)

#Make UMAPs plot properly
options(bitmapType="cairo") 

#Load single cell data
adata <- LoadSeuratRds("/pfs/stor10/users/home/a/awallin/ondemand/wallin/adata.Rdata")
adata$cell = rownames(adata[[]])
#Load sample donor id data
sampleDonorIDs <- read.delim("/pfs/stor10/users/home/m/mb223/mystore/cartdata/sampleDonorIDs.tsv", 
                             row.names=NULL)

#Match donor ids with R row names (shouldn't be necessary for actual data)
sampleDonorIDs$cell = paste0(sampleDonorIDs$cell, "__s1")
sampleDonorIDs = sampleDonorIDs[,c("cell", "donor_id")]
sampleDonorIDs$donor_id = factor(sampleDonorIDs$donor_id) 

#Merge donor ids to seurat object
adata[[]] = left_join(adata[[]], sampleDonorIDs, by = "cell")

#Subset seurat data to only samples from pool 1
DimPlot(pool1, reduction = "umap", label=F, group.by="donor_id")
DimPlot(pool1, reduction = "umap", label=F)

#Subset to t cells
tcellPool1 <- pool1[,Idents(pool1) %in% 0:2]

#Reprocess data
tcellPool1 <- RunPCA(tcellPool1)
ElbowPlot(tcellPool1)
tcellPool1 <- FindNeighbors(tcellPool1, dims = 1:20)
tcellPool1 <- RunUMAP(tcellPool1, dims = 1:20)
DimPlot(tcellPool1, reduction = "umap", label=F)

#Quick look at variables vs clustering
DimPlot(tcellPool1, reduction = "umap", label=F, split.by="donor_id", group.by = "donor_id")
DimPlot(tcellPool1, reduction = "umap", label=F, split.by="hypoxia", group.by = "hypoxia")
DimPlot(tcellPool1, reduction = "umap", label=F, split.by="Phase", group.by = "Phase")
DimPlot(tcellPool1, reduction = "umap", label=F, group.by="CAR", split.by = "CAR")

SaveSeuratRds(tcellPool1, "/pfs/stor10/users/home/m/mb223/mystore/cartdata/tCellPool1.Rdata")

DimHeatmap(tcellPool1, dims = 1:6, cells = 500, balanced = TRUE)

####Pseudobulk analysis from ####
#https://hbctraining.github.io/scRNA-seq_online/lessons/pseudobulk_DESeq2_scrnaseq.html

#https://hbctraining.github.io/Pseudobulk-for-scRNAseq/schedule/self-learning.html might be 
#more updated? Will try both

#One idea: Look at genes most correlated with PCA dimensions and see which
#variables have largest effect on these genes

#Prepare data for analysis
tCellSCE =   as.SingleCellExperiment(tcellPool1)
colData(tCellSCE)$CAR = factor(colData(tCellSCE)$CAR)
cluster_names <- levels(colData(tCellSCE)$donor_id)

#Columns of interest
groups <- colData(tCellSCE)[, c("donor_id", "CAR")]

#Get counts of CAR and donor comparisons
aggr_counts <- aggregate.Matrix(t(counts(tCellSCE)), 
                                groupings = groups, fun = "sum") 

aggr_counts <- t(aggr_counts)

#Create list of matrices of car/donor columns and gene rows
counts_ls <- list()

for (i in 1:length(cluster_names)) {
  
  ## Extract indexes of columns in the global matrix that match a given cluster
  column_idx <- which(tstrsplit(colnames(aggr_counts), "_")[[1]] == cluster_names[i])
  
  ## Store corresponding sub-matrix as one element of a list
  counts_ls[[i]] <- aggr_counts[, column_idx]
  names(counts_ls)[i] <- cluster_names[i]
  
}

str(counts_ls)

#Adding metadata to new data
metadata <- colData(tCellSCE) %>% 
  as.data.frame() %>% 
  dplyr::select(donor_id,CAR) 
names(metadata) = c("sample_id", "cluster_id")
metadata <- metadata[!duplicated(metadata), ]

t <- table(colData(tCellSCE)$CAR,
           colData(tCellSCE)$donor_id)

#Append cellcounts to metadata table
metadata_ls <- list()

for (i in 1:length(counts_ls)) {
  
  ## Initiate a data frame for cluster i with one row per sample (matching column names in the counts matrix)
  df <- data.frame(cluster_sample_id = colnames(counts_ls[[i]]))
  
  ## Use tstrsplit() to separate cluster (cell type) and sample IDs
  df$cluster_id <- tstrsplit(df$cluster_sample_id, "_")[[1]]
  df$sample_id  <- tstrsplit(df$cluster_sample_id, "_")[[2]]
  
  
  ## Retrieve cell count information for this cluster from global cell count table
  idx <- which(colnames(t) == unique(df$cluster_id))
  cell_counts <- t[, idx]
  
  ## Remove samples with zero cell contributing to the cluster
  cell_counts <- cell_counts[cell_counts > 0]
  
  ## Match order of cell_counts and sample_ids
  sample_order <- match(df$sample_id, names(cell_counts))
  cell_counts <- cell_counts[sample_order]
  
  ## Append cell_counts to data frame
  df$cell_count <- cell_counts
  
  
  ## Join data frame (capturing metadata specific to cluster) to generic metadata
  df <- plyr::join(df, metadata, 
                   by = intersect(names(df), names(metadata)))
  
  ## Update rownames of metadata to match colnames of count matrix, as needed later for DE
  rownames(df) <- df$cluster_sample_id
  
  ## Store complete metadata for cluster i in list
  metadata_ls[[i]] <- df
  names(metadata_ls)[i] <- unique(df$cluster_id)
  
}

# Explore the different components of the list
str(metadata_ls)

#DEG analysis
idx <- which(names(counts_ls) == "donor0")
cluster_counts <- counts_ls[[idx]]
cluster_metadata <- metadata_ls[[idx]]

all(colnames(cluster_counts) == rownames(cluster_metadata))

dds <- DESeqDataSetFromMatrix(cluster_counts, 
                              colData = cluster_metadata, 
                              design = ~ sample_id)

#PCA plot
rld <- rlog(dds, blind=TRUE)
DESeq2::plotPCA(rld, ntop = 500, intgroup = "sample_id")


#heatmap
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)
pheatmap(rld_cor, annotation = cluster_metadata[, c("sample_id"), drop=F])

#Version 2 using deseq2
#Create metadata
meta_columns = c("hypoxia", "CAR", "Phase")
meta <- tcellPool1[[]] %>%
  select(meta_columns) %>%
  unique() %>%
  remove_rownames()

####Starting point for 2nd version of pseudobulk####
#Aggregate counts based on variables of interest
#Switching to factors makes it easier to add cell count below
#https://hbctraining.github.io/Pseudobulk-for-scRNAseq/schedule/self-learning.html


tcellPool1$hypoxia = factor(tcellPool1$hypoxia)
tcellPool1$CAR = factor(tcellPool1$CAR)
tcellPool1$Phase = factor(tcellPool1$Phase)

bulk <- AggregateExpression(
  tcellPool1,
  return.seurat = T,
  assays = "RNA",
  group.by = c("hypoxia", "Phase")
)

# Add in number of cells by sample and celltype
n_cells <- tcellPool1[[]] %>% 
  dplyr::count(hypoxia, Phase) %>% 
  rename("n_cells"="n")

bulk[[]] <- n_cells

#Explore cell counts
ggplot(bulk[[]], aes(x=hypoxia, y=n_cells, fill=Phase)) +
  geom_bar(stat="identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  labs(x="Sample name", y="Number of cells")


# Get count matrix
cluster_counts <- FetchData(bulk, layer="counts", vars=rownames(bulk))

# Create DESeq2 object

dds <- DESeqDataSetFromMatrix(t(cluster_counts),
                              colData = bulk[[]],
                              design = ~Phase+hypoxia)

# Transform counts for data visualization
vstRes <- vst(dds, blind=TRUE)

#For saving all pca plots to compare
pcaVars = colnames(vstRes@colData)[-c(1,5)]
# for(i in 1:length(pcaVars)){
#   outputDir = paste0("./mystore/cartdata/examplePCAOutputs/")
#   png(filename=paste0(outputDir,pcaVars[i], ".png"))
#   plotPCA(vstRes, intgroup=c(pcaVars[i])) + theme_classic()+ ggtitle(pcaVars[i])
#   dev.off()
# }

plotPCA(vstRes, intgroup=c(pcaVars[2])) + theme_classic()+ ggtitle(pcaVars[2])

#Correlation heatmap
# Calculate sample correlation
rld_mat <- assay(vstRes)
rld_cor <- cor(rld_mat)

# Plot heatmap

#Correlations are  high. Expected because most genes are not DEGs.
pheatmap(rld_cor)

# Run DESeq2 differential expression analysis
dds <- DESeq(dds)

#Check dispersion
plotDispEsts(dds)

resultsNames(dds)

res <- results(dds,
               contrast=c("Phase", "G2M", "G1"),
               alpha = 0.05)
head(res)
mcols(res)

# Set thresholds
padj.cutoff <- 0.05

# Turn the results object into a tibble for use with tidyverse functions
res_tbl <- res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

# Subset the significant results
sig_res <- dplyr::filter(res_tbl, 
                         padj < padj.cutoff)

sig_res %>% head()

#Visualize DEGs
#Volcano plot

EnhancedVolcano(sig_res,
                sig_res$gene,
                x="log2FoldChange",
                y="padj"
)

#Heat map

# Extract normalized expression for significant genes from the samples
normalized_counts <- counts(dds, normalized=TRUE) %>% as.data.frame()
norm_sig <- normalized_counts %>% 
  dplyr::filter(row.names(normalized_counts) %in% sig_res$gene)


# Run pheatmap using the metadata data frame for the annotation
anno <- colData(dds) %>% 
  as.data.frame() %>% 
  select(Phase, hypoxia)
pheatmap(norm_sig,
         cluster_rows = TRUE,
         show_rownames = FALSE,
         annotation = anno,
         border_color = NA,
         fontsize = 10,
         scale = "row", 
         fontsize_row = 10, 
         height = 20)

#Look at top DEGs
genes <- sig_res %>% 
  arrange(padj) %>% 
  subset(abs(log2FoldChange) > 0.6) %>% 
  head(6)
genes <- genes$gene
genes

plot_list <- list()
#Plot genes
for (gene in genes) {
  # Save plotcounts to a data frame object
  d <- plotCounts(dds, gene=gene, intgroup="Phase", returnData=TRUE)
  
  # Plot the normalized counts for each sample
  p <- ggplot(d, aes(x = Phase, y = count, color = Phase)) + 
    geom_point(position=position_jitter(w = 0.1,h = 0)) +
    theme_bw() +
    ggtitle(gene) +
    theme(plot.title = element_text(hjust = 0.5)) +
    NoLegend()
  plot_list[[gene]] <- p
}

plot_grid(plotlist=plot_list)

#Now look at single-cell level
DefaultAssay(tcellPool1) <- "RNA"
Idents(tcellPool1) <- "Phase"
VlnPlot(tcellPool1, genes)

#UMAPs based on analysis
# Grab the umap coordinates and condition information for each cell
FeaturePlot(tcellPool1, genes, ncol=3)

#Loop through function to get DEGs from all pairwise comparisons. Then
#pivot data to have each row be a DEG with effects of various variables
deseqLoopFun = function(comparison){
  res <- results(dds,
                 name=comparison, #Using name instead of contrast 
                 alpha = 0.05)
  
  # Set thresholds
  padj.cutoff <- 0.05
  
  # Turn the results object into a tibble for use with tidyverse functions
  res_tbl <- res %>%
    data.frame() %>%
    rownames_to_column(var="gene") %>% 
    as_tibble()
  
  # Subset the significant results
  sig_res <- dplyr::filter(res_tbl, 
                           padj < padj.cutoff) %>% 
    mutate(comp = comparison)
  sig_res
}

#Get degs for all comps
deseqOuts <-  lapply(resultsNames(dds)[-1], deseqLoopFun)
deseqOuts <-  do.call(rbind, deseqOuts)

#See how many degs per comparison
totalDEGsPer <- table(deseqOuts$comp) %>% sort(decreasing = T) %>% as.data.frame()
colnames(totalDEGsPer)[1]="comp"

deseqOuts %>% group_by(comp) %>% summarise(meanLog = mean(log2FoldChange),
                                           meanPadj = mean(padj)) %>% 
  arrange((meanPadj)) %>% left_join(totalDEGsPer)

deseqOuts <-  deseqOuts %>% select(c(gene, log2FoldChange, padj, comp)) %>% 
  pivot_wider(names_from = comp,
              values_from = c(log2FoldChange,padj ))


#Attempt to see effects of covariates with linear model
#Pull out expression data. Subset for figuring out workflow

#Loop through manovas and plot distribution of pillai's
#Try parametric bootstrap with normal distributions
datForModel <- t(tcellPool1[["RNA"]]$data[])
bootstrapManova <- function(){
  keepCols <- sample.int(ncol(datForModel), 300)
  datForModel1 <- datForModel[,keepCols] #Start with 300 genes
  datForModel1 <- datForModel1[,(!colSums(datForModel1) == 0)]
  linearDependentCols <- caret::findLinearCombos(as.matrix(datForModel1))
  if(!is.null(linearDependentCols$remove)){
    datForModel1 <- datForModel1[,-linearDependentCols$remove]
  }
  manModel <- manova(as.matrix(datForModel1) ~ tcellPool1$Phase + tcellPool1$hypoxia)
  testStats = summary(manModel)$stats[,"Pillai"] #Differs from car::Anova sometimes
  data.frame(groups = c("Phase", "hypoxia"), testStat = testStats[-length(testStats)])
}
modelRuns <- lapply(seq_len(100), function(x) bootstrapManova())
modelRuns = do.call(rbind, modelRuns)
ggplot(modelRuns, aes(x = groups, y = testStat))+
  geom_boxplot()

ggplot(modelRuns, aes(testStat, fill = groups)) + geom_density(alpha = 0.2)
####Messing around with Limma####
# pool1Counts = tcellPool1[["RNA"]]$counts
# 
# #Turn to DEGlist object
# d0 <- DGEList(pool1Counts)
# d0 <- calcNormFactors(d0)
# 
# #Cut off low expression genes
# cutoff <- 1
# drop <- which(apply(cpm(d0), 1, max) < cutoff)
# d <- d0[-drop,] 
# dim(d)
# d_names = data.frame(cell = rownames(d$samples))
# donors = left_join(d_names, sampleDonorIDs, by ="cell") 
# donors = donors$donor_id
# 
# mm <- model.matrix(~0 + donors)
# 
# #Very weird plot, probably something is wrong
# y <- voom(d, mm, plot = T)
# 
# fit <- lmFit(y, mm)

