#Sometimes HPC2N RStudio crashes when Seurat is loaded - can do this to fix
db <- available.packages() 
deps <- tools::package_dependencies("Seurat", db)$Seurat 
install.packages(deps[21:45]) 
library(Seurat) 

#Can do this if seuratobject cannot load
install.packages("spatstat.random")
install.packages("reticulate")
install.packages("igraph")
install.packages("Rcpp")
install.packages("Matrix")
