#Source functions and packages
source("mystore/cartdata/scripts/CarT_project_functions.R")
library(BioNERO)
#Make UMAPs plot properly
options(bitmapType="cairo")

#Run cellCycleRegressed.R script to get T_cells_noCC or read in:
T_cells_noCC <- LoadSeuratRds('mystore/cartdata/data/T_cells_noCC.rds')
T_cells_noCC_Car <- subset(T_cells_noCC, CAR_pred == 1)
T_cells_noCC_Car.sce <- as.SingleCellExperiment(T_cells_noCC_Car)

SummarizedExperiment::colData(T_cells_noCC_Car.sce) <- SummarizedExperiment::colData(T_cells_noCC_Car.sce)['CAR']
T_cells_noCC_Car.sce <- filter_by_variance(T_cells_noCC_Car.sce, n = 2000)
T_cells_noCC_Car.sce <- PC_correction(T_cells_noCC_Car.sce)

plot_heatmap(T_cells_noCC_Car.sce, type = "samplecor", show_rownames = FALSE)

