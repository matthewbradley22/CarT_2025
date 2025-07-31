# 2025 Car-T Project
Code can be found on hpc2n at: /pfs/stor10/users/home/m/mb223/mystore/cartdata/scripts
## Script descriptions 
### (in a loose order of how one could run through the code)

A note before running anything: <em>CarT_project_functions.R</em> is a script of functions used throughout the analyses. It is sourced in
most other scripts.

1. <em>donor_analysis.R</em>: This script shows how data was originally read in and processed, including donor labelling
from vireo output and doublet labeling. Can also make Fig. 1 UMAPs by plotting deseqDat object loaded in here.

2. <em>Split_donor_analysis.R</em>: The beginning of this script loads in the xgboost classifier used for ML labelling
and labels CD4 and CD8 cells. To see the xgboost model creation look at <em>CD4CD8_classifier-copy.R</em>. May combine
the classification from Split_donor_analysis and add it to donor_analysis to reduce number of scripts.
   
3. <em>cellPhaseAnalysis.R</em>: Uses data processed in donor_analysis.R to create cell phase plots in Fig. 2, as well
as some other cell phase plots

4. <em>cellCycleRegressed.R</em>: This script was used to regress out cell cycle, and the resulting data
called T_cells_noCC is used in all downstream analyses. Also has Figure 2 UMAP code

5. <em>DEGPathAnalyses.R</em>: Has code to create heatmaps of 50 top important genes in Fig 2

6. <em>MitoAnalysis.R</em> includes mitochondrial analyses including creating correlation heatmap with other genes, dotplot of individual mitochondrial genes, and violin plots of variables vs percent.mt
   
7. <em>dengPublicData.R</em>: Create scatter plots of our data lfc versus deng public data response groups. This script is where
public data from Deng paper https://www.nature.com/articles/s41591-020-1061-7#Sec8 is read in and processed, and further analyses are done.
### Other:



<em>hypoxiaAnalysis.R</em>: This script is mostly exploratory analysis, UMAPS, volcano plots etc...

<em> CAR_classifier-copy.R</em> and <em>CD4CD8_classifier-copy.R</em> are scripts where the machine learning was done for labelling cells.
I also saved the xgboost models for reloading (my understanding is xgboost is stochastic, so rerunning these scripts will label cells slightly
differently each time. Which is why I saved the models themselves)



