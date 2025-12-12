## Script to rarefy based on sample coverage.
## Created July 8, 2021
## Script is originally from Komei Kadowaki at 
## Ecological Society of Japan meeting 2016, and modified by Yusuke Okazaki.
## Translated by Shuzo Oita
## Modified by Elizabeth Bowman, eabowman@utexas.edu

# This code was probably used in the following article for example.
# https://www.nature.com/articles/ismej201789#supplementary-information

# Read in libraries, (maybe first two is not necessary)
library(tidyverse); library(vegan) 

# 95% sequence similarity read data (Abundance table. Row: samples, Column: OTUs)
# Fungal mock community
# otu.data <- read_csv('G022.data/SequenceData/Fungal/MockCommunity/OTUtables/maxEE0.5_len250/SitexSpecies_95sim_clean.csv')

# Fungal real community 95% sequence similarity
# otu.data <- read.csv('G022.data/Fungal_SitexSpecies_95sim.csv')

# Fungal real community 90% sequence similarity
# otu.data <- read.csv('G022.data/Fungal_SitexSpecies_90sim.csv')

# Bacterial community 99% sequence similarity
otu.data <- read.csv('data/Bacterial.data/Bacterial_SitexSpecies_99sim_50occ.csv')

## double check before running to isolate community data and exclude meta data
comm.data <- dplyr::select(otu.data, starts_with('Otu'))
comm.data <- comm.data[colSums(comm.data) >= 100] # bacterial 100, fungal
rownames(comm.data) <- otu.data$sample

##############################################
#### Rarefy based on the slope (coverage) ####
##############################################

# Export as a list the slopes when sampled read by read for each samples 
# (Need long time!! It took ~6-8 hrs for me)
rareslopelist<-list()
for(i in 1:nrow(comm.data)){
  rareslopelist[[i]]<-rareslope(comm.data[i,],1:(sum(comm.data[i,])-1))
}
# We reduced 1 from the number of total reads because we cannot calculate 
# the slope for the last read ( = 0 ).

# Search the sample with the highest slope (smallest coverage) and pick up 
# the smallest slope value of the sample
getmincov<-c()
for(i in 1:nrow(comm.data)){
  getmincov[i]<-rareslopelist[[i]][length(rareslopelist[[i]])]
}

# Check the coverage of the sample 
(1-max(getmincov))*100

# Pick up the read number of each sample where its slope reached the coverage
cvrfun <- function(x){min(which(x<=max(getmincov)))+1} # Set the function
cvrrare <- unlist(lapply(rareslopelist,cvrfun))　# Get the value by lapply+unlist

set.seed(123) # To get the replicability, choose any number to fix the random number
OTU_covrared <- rrarefy(comm.data,cvrrare) # Subsampling based on the number of rarefying
row.names(OTU_covrared) <- otu.data$Samples

# 99% sequence similarity: Bacterial community


# 95% sequence similarity: fungal  community
# write.csv(OTU_covrared,
#           'G022.data/Fungal_SitexSpecies_95sim_rarefied.csv',
#           row.names = T)  # Export as CSV

# 90% sequence similarity: fungal  community
# write.csv(OTU_covrared,
#           'G022.data/Fungal_SitexSpecies_90sim_rarefied.csv',
#           row.names = T)  # Export as CSV
