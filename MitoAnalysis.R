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
genesOfInterest <- VariableFeatures(T_cells_noCC)


#Look at genes vs percent.mt
pecentMt <- as.matrix(T_cells_noCC$percent.mt)
colnames(pecentMt) <- 'mitoPct'

normalCounts <- T_cells_noCC[['RNA']]$data[genesOfInterest,]
normalCountsMat <- t(as.matrix(normalCounts))

#This is slow because it calculates all pairwise correlations and we only
#care about correlations with mitochondrial score. But a better function than base R cor maybe
mitoCor <- Hmisc::rcorr(pecentMt, normalCountsMat)
mitoCor_df <- data.frame(cor = mitoCor$r[,'mitoPct'], pVal = mitoCor$P[,'mitoPct'])
mitoCor_df <- mitoCor_df[order(mitoCor_df$cor, decreasing = TRUE),]

#Bonferonni correction
mitoCor_df$pVal <- mitoCor_df$pVal * nrow(mitoCor_df)
mitoCor_df <- mitoCor_df[-1,]
mitoCor_df$gene <- rownames(mitoCor_df)
mitoCor_df$absoluteCor <- abs(mitoCor_df$cor)
mitoCor_df[order(mitoCor_df$absoluteCor, decreasing = TRUE),] %>% write.csv(file = 'mystore/cartdata/data/mitoCors.csv')

#Set gene to factor to determine plotting order
#Saving plot 7x3
topCorGenes <- mitoCor_df[order(mitoCor_df$absoluteCor, decreasing = TRUE),] %>% head(n = 20) %>% arrange(cor)
genes <- factor(topCorGenes$gene, levels = topCorGenes$gene)
mitoCor_df[order(mitoCor_df$absoluteCor, decreasing = TRUE),] %>% head(n = 20) %>% arrange(cor) %>% ggplot(aes(x = factor(gene, levels = genes), y = 1, fill = cor))+
  geom_tile() + theme(axis.text.x = element_text(angle = 90),
                      axis.title.y=element_blank(),
                      axis.title.x=element_blank(),
                      axis.text.y=element_blank(),
                      axis.ticks.y=element_blank(),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank(),
                      axis.text=element_text(size=16)) + 
  scale_y_continuous(breaks = seq(0, 1, by = 1)) +
  scale_fill_gradientn(colours = c("orange", "white", "blue")) +
  ggtitle('Percent.mt Top 20 Correlated Genes')

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

#Look at pairwise comparisons (useful if you want to explore individual mt genes vs other genes. But slow to run)
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
