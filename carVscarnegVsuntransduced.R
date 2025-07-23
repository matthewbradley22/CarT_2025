#Source functions and packages
source("mystore/cartdata/scripts/CarT_project_functions.R")
library(EnhancedVolcano)
library(readr)
library(gprofiler2)
library(stringr)

#Make UMAPs plot properly
options(bitmapType="cairo")

#Run cellCycleRegressed.R script to get T_cells_noCC 
#or load in from
T_cells_noCC <- LoadSeuratRds('./mystore/cartdata/data/T_cells_noCC.rds')
DimPlot(T_cells_noCC, group.by = "Phase")

#Add column to help with analysis of car neg vs untransduced
T_cells_noCC[[]] <- T_cells_noCC[[]] %>% mutate('carGroup' = case_when(CAR == 'untransduced' ~ 'untransduced',
                                                     CAR!='untransduced' & CAR_pred == 1 ~ 'CarPos',
                                                     CAR!='untransduced' & CAR_pred == 0 ~ 'CarNeg' ))
#Function to make a bunch of gene lists and volcano plots comparing car + and car -
compareCarNegAndPos <- function(subsetParams){
  dat <- subset(T_cells_noCC, day == subsetParams[1] & CD_pred == subsetParams[2] & CAR == subsetParams[3]
                & hypoxia == subsetParams[4])
  dat_bulk <- getPseudoBulkObject(dat, c("donor_id", 'CAR_pred'))
  dat_bulk <- DESeq(dat_bulk)
  
  res_dat_bulk <- lfcShrink(dat_bulk, coef="CAR_pred_1_vs_0", type="apeglm")
  pathToSave <- './mystore/cartdata/data/carNegVsPosComps/'
  addOn <- paste0(subsetParams, collapse = '_')
  res_dat_bulk[order(res_dat_bulk$padj),]  %>% as.data.frame() %>% head(n = 1000) %>% 
    write.csv(., paste0(pathToSave, addOn, '_PosNeg.csv'))
  
  pdf(file = paste0('./mystore/cartdata/Plots/CarNegPosPlots/', addOn, '.pdf'), width = 8, height = 8)
  print(EnhancedVolcano(res_dat_bulk,
                  lab = rownames(res_dat_bulk),
                  x = 'log2FoldChange',
                  y = 'pvalue',
                  title = paste0(addOn, 'Car+ vs Car-')))
  dev.off()
}

#get all combinations of car, hypoxia, day, and cd
days = as.character(unique(T_cells_noCC$day))
CD_pred = as.character(unique(T_cells_noCC$CD_pred))
CAR = as.character(unique(T_cells_noCC$CAR))
hypoxia = as.character(unique(T_cells_noCC$hypoxia))

allCombos <- expand.grid(days, CD_pred, CAR, hypoxia)
allCombos <- subset(allCombos, Var1 != 'D7' | Var4 != 'NH')
allCombos <- subset(allCombos, Var3 != 'untransduced')

for(i in 1:nrow(allCombos)){
  params = as.character(unname(unlist(allCombos[i,])))
  compareCarNegAndPos(params)
}

#Read in gene lists
path = 'mystore/cartdata/data/carNegVsPosComps/'
files = list.files(path)

#Create GO tables from gene lists
for(i in 1:length(files)){
  name <-  gsub('.csv', '', files[i]) 
  dat <- read_csv(paste0("mystore/cartdata/data/carNegVsPosComps/", files[i]))
  genes <- unname(unlist(dat[1:100,1]))
  #Get gene ontology table
  gostres <- gost(query = genes, 
                  organism = "hsapiens", ordered_query = FALSE, 
                  multi_query = FALSE, significant = TRUE,
                  user_threshold = 1E-3, evcodes = TRUE,
                  sources = c("GO"))
  
  if(class(gostres) == "list"){
    table <- publish_gosttable(gostres,
                               use_colors = TRUE, 
                               show_columns = c("source", "term_name", "term_size", "intersection_size"),
                               filename = NULL) + 
      ggtitle(name)
    ggsave(plot = table, filename = paste0('mystore/cartdata/Plots/goTablesCarNegPos/', name, '.pdf'),
           width = 20, height = 35)
  }

}
  
datList <- list()
for(i in 1:length(files)){
  dat <- read_csv(paste0("mystore/cartdata/data/carNegVsPosComps/", files[i]))
  genes <- unname(unlist(dat[1:50,1]))
  datList[[i]] = data.frame(place = 1:50, geneName = genes, listNum = i)
}

datList = do.call('rbind', datList)
table(datList$geneName) %>% sort(decreasing = T) %>% subset(.>10) %>% names()
  


#CAR neg vs untransduced
#Function to make a bunch of gene lists and volcano plots comparing car + and car -
compareCarNegAndUntransduced <- function(subsetParams){
  dat <- subset(T_cells_noCC, day == subsetParams[1] & CD_pred == subsetParams[2] & (CAR == subsetParams[3]|CAR == 'untransduced')
                & hypoxia == subsetParams[4] & carGroup != 'CarPos')
  print(table(dat$CAR, dat$CAR_pred))
  dat_bulk <- getPseudoBulkObject(dat, c("donor_id", 'CAR'))
  dat_bulk <- DESeq(dat_bulk)
  
  res_dat_bulk <- lfcShrink(dat_bulk, coef=paste0("CAR_untransduced_vs_", subsetParams[3]), type="apeglm")
  #pathToSave <- './mystore/cartdata/data/carNegVsUntransducedComps/'
  addOn <- paste0(subsetParams, collapse = '_')
  #res_dat_bulk[order(res_dat_bulk$padj),]  %>% as.data.frame() %>% head(n = 1000) %>% 
  # write.csv(., paste0(pathToSave, addOn, '_PosNeg.csv'))
  
  pdf(file = paste0('./mystore/cartdata/Plots/CarNegUntransducedPlots/', addOn, '.pdf'), width = 10, height = 10)
  print(EnhancedVolcano(res_dat_bulk,
                        lab = rownames(res_dat_bulk),
                        x = 'log2FoldChange',
                        y = 'pvalue',
                        title = paste0(addOn, ': Untransduced vs Car-'),
                        max.overlaps =  20,
                        drawConnectors = TRUE))
  dev.off()
}

for(i in 1:nrow(allCombos)){
  params = as.character(unname(unlist(allCombos[i,])))
  compareCarNegAndUntransduced(params)
}

#Check genes from volcano plots
dat <- subset(T_cells_noCC, day == 'D13' & CD_pred == 'CD4' & (CAR == 'M28z'|CAR == 'untransduced')
              & hypoxia == 'NH' & carGroup != 'CarPos')
dat_bulk <- getPseudoBulkObject(dat, c("donor_id", 'CAR'))
dat_bulk <- DESeq(dat_bulk)

res_dat_bulk <- lfcShrink(dat_bulk, coef=paste0("CAR_untransduced_vs_", 'M28z'), type="apeglm")
res_dat_bulk[order(res_dat_bulk$padj),]  %>% as.data.frame() %>% head(n = 20)
plotCounts(dat_bulk, gene='CCR6', intgroup="CAR")

EnhancedVolcano(res_dat_bulk,
                      lab = rownames(res_dat_bulk),
                      x = 'log2FoldChange',
                      y = 'pvalue',
                      title = paste0( ': Untransduced vs Car-'),
                      max.overlaps =  15,
                      drawConnectors = TRUE)

DotPlot(T_cells_noCC, features = 'CCR6', group.by = 'CAR')


#Volcanoes comparing hypoxia groups

compareNNHH <- function(subsetParams){
  dat <- subset(T_cells_noCC, day == subsetParams[1] & CD_pred == subsetParams[2] & CAR == subsetParams[3]
                & CAR_pred == subsetParams[4])
  print(table(dat$CAR, dat$hypoxia))
  dat_bulk <- getPseudoBulkObject(dat, c("donor_id", 'hypoxia'))
  dat_bulk <- DESeq(dat_bulk)
  
  res_dat_bulk <- lfcShrink(dat_bulk, coef="hypoxia_NN_vs_HH", type="apeglm")
  pathToSave <- './mystore/cartdata/data/HHvsNNComps/'
  addOn <- paste0(subsetParams, collapse = '_')
  res_dat_bulk[order(res_dat_bulk$padj),]  %>% as.data.frame() %>% head(n = 1000) %>% 
   write.csv(., paste0(pathToSave, addOn, '_HHNN.csv'))
  
  #pdf(file = paste0('./mystore/cartdata/Plots/CarNNvsHH/', addOn, '.pdf'), width = 10, height = 10)
  # print(EnhancedVolcano(res_dat_bulk,
  #                       lab = rownames(res_dat_bulk),
  #                       x = 'log2FoldChange',
  #                       y = 'pvalue',
  #                       title = paste0(addOn, ': NN vs HH'),
  #                       max.overlaps =  20,
  #                       drawConnectors = TRUE))
  # dev.off()
  
}

#get all combinations of car, hypoxia, day, and cd
days = as.character(unique(T_cells_noCC$day))
CD_pred = as.character(unique(T_cells_noCC$CD_pred))
CAR = as.character(unique(T_cells_noCC$CAR))
CAR_pred= as.character(unique(T_cells_noCC$CAR_pred))

allCombos <- expand.grid(days, CD_pred, CAR, CAR_pred)
allCombos <- subset(allCombos, Var3 != 'untransduced')

for(i in 1:nrow(allCombos)){
  params = as.character(unname(unlist(allCombos[i,])))
  compareNNHH(params)
}


#Compare MBBz to M28z
compareMBBzM28z <- function(subsetParams){
  dat <- subset(T_cells_noCC, day == subsetParams[1] & CD_pred == subsetParams[2] & (CAR == 'MBBz' |CAR == 'M28z')
                & CAR_pred == subsetParams[3] & hypoxia == subsetParams[4])
  print(table(dat$CAR, dat$CD_pred))
  dat_bulk <- getPseudoBulkObject(dat, c("donor_id", 'CAR'))
  dat_bulk <- DESeq(dat_bulk)
  
  res_dat_bulk <- lfcShrink(dat_bulk, coef="CAR_MBBz_vs_M28z", type="apeglm")
  #pathToSave <- './mystore/cartdata/data/carNegVsUntransducedComps/'
  addOn <- paste0(subsetParams, collapse = '_')
  #res_dat_bulk[order(res_dat_bulk$padj),]  %>% as.data.frame() %>% head(n = 1000) %>% 
  # write.csv(., paste0(pathToSave, addOn, '_PosNeg.csv'))
  
  pdf(file = paste0('./mystore/cartdata/Plots/MBBzVsM28zVolcanoes/', addOn, '.pdf'), width = 10, height = 10)
  print(EnhancedVolcano(res_dat_bulk,
                        lab = rownames(res_dat_bulk),
                        x = 'log2FoldChange',
                        y = 'pvalue',
                        title = paste0( addOn, ': MBBz vs M28z'),
                        max.overlaps =  20,
                        drawConnectors = TRUE))
  dev.off()
}

#get all combinations of car, hypoxia, day, and cd
days <-  as.character(unique(T_cells_noCC$day))
CD_pred <-  as.character(unique(T_cells_noCC$CD_pred))
CAR_pred <- as.character(unique(T_cells_noCC$CAR_pred))
hypoxia <- as.character(unique(T_cells_noCC$hypoxia))

allCombos <- expand.grid(days, CD_pred, CAR_pred, hypoxia)
allCombos <- subset(allCombos, Var1 != 'D7' | Var4 != 'NH')

for(i in 1:nrow(allCombos)){
  params = as.character(unname(unlist(allCombos[i,])))
  compareMBBzM28z(params)
}

#Find most often genes in car pos neg tables
carNegPosFiles = list.files('mystore/cartdata/data/carNegVsPosComps/')
carNegPosFiles = paste0('mystore/cartdata/data/carNegVsPosComps/', carNegPosFiles)
carNegPosTables = lapply(carNegPosFiles, read.csv)

carNegPosNames = list.files('mystore/cartdata/data/carNegVsPosComps/')
names(carNegPosTables) <- carNegPosNames
carNegPosTablesSub <- lapply(carNegPosTables, FUN = function(x){
  head(x, n = 30)
})

for(i in 1:length(carNegPosTablesSub)){
  carNegPosTablesSub[[i]]$comp = gsub('.csv','',carNegPosNames[[i]])
}

carNegPosTablesSub <- do.call(rbind, carNegPosTablesSub)


sort(table(carNegPosTablesSub$X))
carNegPosTablesSub[carNegPosTablesSub$X == 'GZMB',]$comp

#Take a list of csvs with DEGs and find most common genes in top n genes of each list
findCommonGenes <- function(filePath, numGenes, topGenes = TRUE){
  geneFiles = list.files(filePath)
  geneFiles = paste0(filePath, geneFiles)
  geneTables = lapply(geneFiles, read.csv)
  
  geneListNames = list.files(filePath)
  names(geneTables) <- geneListNames
  
  if(topGenes){
    geneTables <- lapply(geneTables, FUN = function(x){
      head(x, n = numGenes)
    })
  }

  for(i in 1:length(geneTables)){
    geneTables[[i]]$comp = gsub('.csv','',geneListNames[[i]])
  }
  
  geneTables <- do.call(rbind, geneTables)
  sort(table(geneTables$X))
}

mostCommonHHNNGenes <- findCommonGenes('mystore/cartdata/data/HHvsNNComps/', n = 30, top = TRUE) %>%
  tail(n = 40)

findGeneInTables <- function(filePath, numGenes, geneName){
  geneFiles = list.files(filePath)
  geneFiles = paste0(filePath, geneFiles)
  geneTables = lapply(geneFiles, read.csv)
  
  geneListNames = list.files(filePath)
  names(geneTables) <- geneListNames
  
  
  geneTables <- lapply(geneTables, FUN = function(x){
    head(x, n = numGenes)
  })
  
  present <- lapply(geneTables, FUN = function(x){
    geneName %in% x$X
  })
  
  unname(unlist(present))
  presenceDF <- data.frame(names = names(present), presence =  unname(unlist(present)))
  presenceDF <- presenceDF[presenceDF$presence == TRUE,]
  str_remove(presenceDF$names, '.csv')
}

as.data.frame(findGeneInTables('mystore/cartdata/data/HHvsNNComps/', 30, '7SK'))
