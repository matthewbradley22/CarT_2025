#Source functions and packages
source("mystore/cartdata/scripts/CarT_project_functions.R") #If this fails can check the IfLoadingSeuratCrashesR script

#Make UMAPs plot properly
options(bitmapType="cairo") 

#Load t cells with cell cycle regressed out'
T_cells_noCC <- LoadSeuratRds('mystore/cartdata/data/T_cells_noCC.rds')
T_cells_noCC[[]] <- T_cells_noCC[[]] %>% mutate(car_hypox = paste(CAR, hypoxia, sep = '_'))
#Create subsets
cd4_day7 <- subset(T_cells_noCC, CD_pred == 'CD4' & day == 'D7')
cd8_day7 <- subset(T_cells_noCC, CD_pred == 'CD8' & day == 'D7')
cd4_day13 <- subset(T_cells_noCC, CD_pred == 'CD4' & day == 'D13')
cd8_day13 <- subset(T_cells_noCC, CD_pred == 'CD8' & day == 'D13')

#Gene lists from https://www.nature.com/articles/s41467-025-59904-x#Sec1
cytokineActivity <- c('IL1A', 'TGFB1', 'TGFB3', 'TNFA', 'IL6', 'IL15', 'CXCL12', 'IL13', 'IFN1',
                      'IL1B', 'IL2', 'IL3', 'IL22', 'IL4', 'IL10', 'IFNL', 'IL36', 'IL17A', 'IL12', 'IL21', 
                      'IFNG', 'IL27')
cytoxicGenes <- c('GZMA', 'GZMB', 'GZMH', 'GZMM', 'GZMK', 'PRF1', 'GNLY', 'NKG7')
negAntiTumor <- c('TGFB1', 'IL6', 'TGFB2', 'TGFB3')
JunGenes <- c('JUN', 'JUNB', 'JUND')
#Average expression heatmap
createAvgExpHeatmap <- function(data, grouping, geneFeats){
  avgExp <- AverageExpression(data, group.by = grouping, return.seurat = T)
  DoHeatmap(avgExp, features = geneFeats, slot = 'data', draw.lines = FALSE)
}

createAvgExpHeatmap(cd8_day7, 'car_hypox', cytokineActivity)


