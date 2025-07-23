library(Seurat)
library(gridExtra)
library(DESeq2)
library(ggplot2)
library(dplyr)
#Turn single cell object into a pseudobulk rna object based on variables of interest
getPseudoBulkObject <- function(dat, designVars, intercept = "include", return.Seurat = FALSE){
  bulk <- AggregateExpression(
    dat,
    return.seurat = T,
    assays = "RNA",
    group.by = designVars)
  
  # Add in number of cells by sample and celltype
  n_cells <- dat[[]] %>% dplyr::select(all_of(designVars)) %>% 
    dplyr::group_by(across(all_of(designVars))) %>% 
    dplyr::count() 
  
  bulk[[]] <- n_cells
  
  if(return.Seurat){
    return(bulk)
  }
  # Get count matrix
  cluster_counts <- FetchData(bulk, layer="counts", vars=rownames(bulk))
  
  # Create DESeq2 object
  if(intercept == "include"){
    designForm <- reformulate(termlabels = designVars)
    dds <- DESeqDataSetFromMatrix(t(cluster_counts),
                                  colData = bulk[[]],
                                  design = designForm)
  }else{
    designForm <- reformulate(termlabels = c(0,designVars))
    dds <- DESeqDataSetFromMatrix(t(cluster_counts),
                                  colData = bulk[[]],
                                  design = designForm)
  }
  
  dds
}

#Plot multiple PCAs in one grid
#provide deseq data, the variable set you're interested in, and number of rows
#to plot
getPCAGrid <- function(dat, varSet, rowNum){
  vstRes <- vst(dat, blind=TRUE)
  plotList = list()
  for(i in 1:length(varSet)){
    plotList[[i]] = plotPCA(vstRes, intgroup = varSet[i])
  }
  grid.arrange(grobs = plotList, nrow = rowNum)
}

#Prepare seurat data for UMAP
prepUMAP <- function(dat, regressCC = FALSE, return.elbow = FALSE){
  # Data preprocessing
  dat <- NormalizeData(dat, normalization.method = "LogNormalize")
  dat <- FindVariableFeatures(dat, selection.method = "vst", nfeatures = 2000)
  if(regressCC == 'CC'){
    dat <- ScaleData(dat, vars.to.regress = c("S.Score", "G2M.Score", "Phase"), features = VariableFeatures(dat)) 
  }else if(regressCC == 'phase'){
    dat <- ScaleData(dat, vars.to.regress = c("Phase"), features = VariableFeatures(dat)) 
  }else{
    dat <- ScaleData(dat) 
  }
  dat <- RunPCA(dat)
  return(dat)
}

finishUMAP <- function(scaledDat, dimNeighbors = 25){
  dat <- FindNeighbors(scaledDat, dims = 1:dimNeighbors)
  dat <- FindClusters(dat, resolution = 0.5)
  dat <- RunUMAP(dat, dims = 1:dimNeighbors)
  dat
}

#Compare 1 deseq group to the mean of other groups of same variable
deseqOneVsMean <- function(dat, oneVar, meanVars){
   results(dat, 
          contrast =list(c(oneVar), 
                 c(meanVars)),
           listValues=c(1, -1/length(meanVars)))
}

#Create a pseudobulk UMAP from a single cell object
getPseudoBulkUMAP <- function(dat, vars){
  dat <- getPseudoBulkObject(dat, vars)
  dat <- vst(dat)
  normalized_counts <- assay(dat) %>%
    t() 
  
  dat_UMAP <- umap::umap(normalized_counts)
  
  umap_plot_df <- data.frame(dat_UMAP$layout) %>%
    tibble::rownames_to_column("orig.ident") %>%
    # Add the metadata into this data frame; match by sample IDs
    dplyr::inner_join(data.frame(colData(dat)), by = "orig.ident")
 
  umap_plot_df
}

#Take a deseq object and get a list of pairwise comparison gene lists
getTopGenes <- function(dat, geneList = NA, pos = TRUE){
  comparisons = list(rep(0, length(resultsNames(dat))-1))
  for(i in 1: length(resultsNames(dat))-1){
    comparisons[i] = results(dat, name = resultsNames(dat)[i+1])
  }
  names(comparisons) = resultsNames(dat)[-1]
  importantGeneDat <- lapply(comparisons, FUN = function(x){
    dat = as.data.frame(x)
    if(pos == TRUE){
      dat = filter(dat, log2FoldChange > 0)
    }
    dat = filter(dat, padj < 0.01) %>% arrange(padj)
    if(!all(is.na(geneList))){
      dat = dat[rownames(dat) %in% geneList$gene,]
      if(nrow(dat)>0){
        dat$gene = rownames(dat)
        dat <- dat %>% left_join(CD4Genes, by = "gene")
      }
    }
    dat
  })
  for(i in 1:length(importantGeneDat)){
    if(nrow(importantGeneDat[[i]])>0){
      importantGeneDat[[i]]$comparison = names(importantGeneDat)[[i]]
    }
  }
  
  importantGeneDat
}

#Plot output from above function (getTopGenes)
plotTopGenes <- function(dat){
  dat <- do.call(rbind, dat)
  dat %>% group_by(comparison, type) %>% summarize(total = n()) %>% 
    arrange(desc(total)) %>% filter(total > 1 & !is.na(type)) %>% ggplot(aes(x = comparison, y = total,
                                                                             fill = type))+
    geom_bar(stat = "identity")+
    theme(axis.text.x = element_text(angle = 45, vjust = 0.6))+
    ylab("Number of DEGs")
}

getPCALoadings <- function(dat, numGenes = 500, numCom = 2){
  dat_vst <- vst(dat)
  rv <- rowVars(assay(dat_vst)) #Variance per row
  select <- order(rv, decreasing=TRUE)[seq_len(min(numGenes, length(rv)))]
  pca <- prcomp(t(assay(dat_vst)[select,]))
  loadings <- pca$rotation[, seq_len(numCom)]
  loadings
}


prepTopGO <- function(dat, PC = 1){
  entrez <- mapIds(org.Hs.eg.db, keys =  names(dat[,PC]), 
                   keytype = "SYMBOL", column="ENTREZID")
  genes <- dat[,PC]
  names(genes) <- unname(entrez)
  genes <- genes[!is.na(names(genes))]
  genes
  
}

getFisherEnrichment <- function(dat, goID){
  goID <- goID
  gene.universe <- genes(dat)
  go.genes <- genesInTerm(dat, goID)[[1]]
  sig.genes <- sigGenes(dat)
  
  my.group <- new("classicCount", testStatistic = GOFisherTest, name = "fisher",
                  allMembers = gene.universe, groupMembers = go.genes,
                  sigMembers = sig.genes)
  contTable(my.group)
}
