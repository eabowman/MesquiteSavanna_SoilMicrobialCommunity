The folders here contain R scripts and data files for the article "Soil microbial assembly mirrors savanna phase separation amid plant invasion and soil disturbance.". Written by Elizabeth A. Bowman, A. Peyton Smith, Aaron C. Rhodes, Robert M. Plowes, and Lawrence E. Gilbert.


Elizabeth A. Bowman is the author for correspondence on this article. 
Address: 2907 Lake Austin Blvd.
Brackenridge Field Laboratory
University of Texas at Austin
Austin, TX 78703 USA. 
Email: eabowman@utexas.edu

Elizabeth A. Bowman wrote the scripts for data analysis.

------------------------------------------------------------------------------------------

Explanation of folders:
1. The data folder contains all data files used in the R script for analyses with
explanations of columns.
2. The figures folder is an output folder where figures generated in the R script will be
output.
3. The results folder is an output folder for results tables generated in the R script.
4. The script folder contains all scripts organized by type of analysis.

The data folder contains a 'Read me' file with column and file descriptions.

Explanation of each script is below.
------------------------------------------------------------------------------------------

Loadlibraries.R: Loads all libraries used in analyses. The installation of each library is commented out.

0.analysisPlantCommunityClimate.R: Analyses associated with the plant community and climate variables (mean annual precipitation, MAP, and mean annual temperature, MAT) downloaded from WorldClim.

1.analysisWoodlandGrassland.R: Analyses across the savanna including both grasslands and woody patches. (Fig. 1)

2.analysisGrasslands.R: Analyses conducted on the grassland subset of data. (Fig. 2 and fig. 4a)

3.analysisWoodyPatches.R: Analyses conducted on the woody patch subset of data. (Fig. 3 and fig. 4b)

4.analysisTaxonomy.R: Taxonomic analyses of bacterial and fungal soil communities across the savanna. (Fig.5

CoverageBasedRarefaction.R: Code for doing coverage based rarefaction of CF Illumina data. Originally written by Komei Kadowaki, modified by E.A. Bowman.

varpart2.MEM.R: Function for distance based Moran's Eigenvector (dbMEM) for geographical distance. Written by Pierre Legendre. 
