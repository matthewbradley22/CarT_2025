# 2025 Car-T Project

## Script descriptions 

A note before running anything: CarT_project_functions is a script of functions used throughout the analyses. It is sourced in
most other scripts.

### Loose order of how one could run through scripts

1. donor_analysis.R: This script shows how data was originally read in and processed, including donor labelling
from vireo output and doublet labeling.

2. cellPhaseAnalysis.R: Uses data processed in donor_analysis.R to create cell phase plots in Fig. 2, as well
as some other cell phase plots

3. cellCycleRegressed.R: This script was used to regress out cell cycle, and the resulting data
called T_cells_noCC is used in all downstream analyses. Also has Figure 2 UMAP code

4. DEGPathAnalyses.R: Has code to create heatmaps of 50 top important genes in Fig 2

Other:

dengPublicData: Create scatter plots of our data lfc versus deng public data response groups. This script is where
public data from Deng paper https://www.nature.com/articles/s41591-020-1061-7#Sec8 is read in and processed

TODO: Add how we did machine learning labeling
