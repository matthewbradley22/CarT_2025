library(readr)
library(dplyr)
library(Seurat)

#Source fucntions
source("mystore/cartdata/scripts/CarT_project_functions.R")

#Load data
T_cells <- LoadSeuratRds("./mystore/cartdata/data/tcell_deseqdat_singlets.rds")

#reportFiles <- list.files("./mystore/cartdata/Trust4Output/", pattern = "*sorted_report*")
barcodeFiles <-  list.files("./mystore/cartdata/Trust4Output/", pattern = "*sorted_barcode*")
#reportFiles <- lapply(paste0("./mystore/cartdata/Trust4Output/",reportFiles), read_delim, delim = "\t",
 #      trim_ws = TRUE)
barcodeFiles <- lapply(paste0("./mystore/cartdata/Trust4Output/",barcodeFiles), read_delim, delim = "\t",
                      trim_ws = TRUE)

#reportFilesComb <- do.call(rbind, reportFiles)
barcodeFilesComb <- do.call(rbind, barcodeFiles)
cell_types = barcodeFilesComb[,1:2]
colnames(cell_types) = c("cell", "cell_type")
T_cells[[]] <- left_join(T_cells[[]], cell_types, by = "cell")

DimPlot(T_cells, split.by = "cell_type", group.by = "cell_type")
table(T_cells$cell_type, useNA="always")

rownames(adata[[]]) %in% cell_types$cell
cell_types$cell %in% rownames(adata[[]])

rownames(adat[[]])[grep("^17_33.*", rownames(adat[[]]))]
