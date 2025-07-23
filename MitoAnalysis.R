#Load libraries
library(scran)
library(Hmisc)

#Source fucntions 
source("mystore/cartdata/scripts/CarT_project_functions.R")

#Make plotting work 
options(bitmapType="cairo") 

#Load data
T_cells_noCC <- LoadSeuratRds('mystore/cartdata/data/T_cells_noCC.rds')

MT_genes <- rownames(T_cells_noCC[['RNA']]$counts)[grep('^MT-', rownames(T_cells_noCC[['RNA']]$counts))]
genesOfInterest <- c(VariableFeatures(T_cells_noCC))

#Look at variably expressed genes vs mt genes
normalCounts <- T_cells_noCC[['RNA']]$data[c(MT_genes, genesOfInterest),]
comparisons_gene1 = rep(MT_genes, each = 2000)
comparions_gene2 = rep(genesOfInterest, 37)
compMatrix <- matrix(c(comparisons_gene1, comparions_gene2), ncol = 2)

geneCors <- correlatePairs(normalCounts, pairings = compMatrix)

geneCors <- data.frame(geneCors@listData$gene1, geneCors@listData$gene2, geneCors@listData$rho, geneCors@listData$FDR)
geneCors <- geneCors[geneCors$geneCors.listData.FDR < 0.01,]
topGeneCors <- geneCors[order(geneCors$geneCors.listData.FDR),] %>% head(n = 500)
topByGroup <- topGeneCors %>% group_by(geneCors.listData.gene1) %>% top_n(n = 20)
table(topByGroup$geneCors.listData.gene2) %>% sort()

#Look at genes vs percent.mt
pecentMt <- as.matrix(T_cells_noCC$percent.mt)
colnames(pecentMt) <- 'mitoPct'
normalCounts <- T_cells_noCC[['RNA']]$data[genesOfInterest,]
normalCountsMat <- t(as.matrix(normalCounts))

#This is slow because it calculates all pairwise correlations and we only
#care about correlations with mitochondrial score. But a better function than base R cor
mitoCor <- Hmisc::rcorr(pecentMt, normalCountsMat)
mitoCor_df <- data.frame(cor = mitoCor$r[,'mitoPct'], pVal = mitoCor$P[,'mitoPct'])
mitoCor_df <- mitoCor_df[order(mitoCor_df$cor, decreasing = TRUE),]

#Bonferonni correction
mitoCor_df$pVal <- mitoCor_df$pVal * nrow(mitoCor_df)
mitoCor_df <- mitoCor_df[-1,]
mitoCor_df$gene <- rownames(mitoCor_df)
mitoCor_df$absoluteCor <- abs(mitoCor_df$cor)
mitoCor_df[order(mitoCor_df$absoluteCor, decreasing = TRUE),] %>% write.csv(file = 'mystore/cartdata/data/mitoCors.csv')

#Look at most prevalent mt genes
MT_genes <- rownames(T_cells_noCC[['RNA']]$counts)[grep('^MT-', rownames(T_cells_noCC[['RNA']]$counts))]
VlnPlot(T_cells_noCC, features = MT_genes)
DotPlot(T_cells_noCC, features = MT_genes, group.by = 'hypoxia', assay = 'RNA') +
  theme(axis.text.x = element_text(angle = 90))

#Look at percent.mt vs other variables (hypoxia etc...)
ggplot(T_cells_noCC[[]], aes(x = hypoxia, y = percent.mt)) + 
  geom_boxplot()
ggplot(T_cells_noCC[[]], aes(x = CAR, y = percent.mt)) + 
  geom_boxplot()
ggplot(T_cells_noCC[[]], aes(x = CD_pred, y = percent.mt)) + 
  geom_boxplot()
ggplot(T_cells_noCC[[]], aes(x = carExpression, y = percent.mt)) + 
  geom_boxplot()
ggplot(T_cells_noCC[[]], aes(x = hypoxia, y = percent.mt, fill = CAR)) + 
   geom_boxplot()
ggplot(T_cells_noCC[[]], aes(x = day, y = percent.mt, fill = CAR)) + 
  geom_boxplot()

