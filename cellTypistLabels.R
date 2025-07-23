#Source functions and packages
source("mystore/cartdata/scripts/CarT_project_functions.R")
library(readr)


#Make UMAPs plot properly
options(bitmapType="cairo")


#Run cellCycleRegressed.R script to get T_cells_noCC or read in:
T_cells_noCC <- LoadSeuratRds('mystore/cartdata/data/T_cells_noCC.rds')
T_cells_noCC_cd4d7 <- subset(T_cells_noCC, CAR_pred == 1 & CD_pred == 'CD4' & day == 'D7')

cd4d7Labels <- read_csv("mystore/cartdata/data/predicted_labels._cd4d7csv.csv")

T_cells_noCC_cd4d7$cellTypistLabels <- cd4d7Labels$majority_voting
T_cells_noCC_cd4d7 <- prepUMAP(T_cells_noCC_cd4d7)
ElbowPlot(T_cells_noCC_cd4d7)
T_cells_noCC_cd4d7 <- finishUMAP(T_cells_noCC_cd4d7, dimNeighbors = 15)

DimPlot(T_cells_noCC_cd4d7, group.by = 'cellTypistLabels')
