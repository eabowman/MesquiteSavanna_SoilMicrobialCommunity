### Data associated with script 0.analysisPlantCommunityClimate.R. ###
field.data_climate.csv
----------------------
plot: Plot sample taken from
block: Site
disturbance: disturbance of site where sample was collected
disturbance.date: Estimated data of last disturbance based on aerial photos
vegetation.type: vegetation type of site where sample was collected
invasion: invasion of site where sample was collected
elevation.m: elevation above sea level in meters
pasture: Pasture name
evidence.of.cattle: Evidence of cattle activity at sample site
Notes: Site notes
MAT: Mean annual temperature, WorldClim
MAP: Mean annual precipitaiton, WorldClim

plant.comm.sxs.csv
------------------
sample: Sample name
treatment: Combination of disturbance, vegetation, and invasion; D = Disturbed, U = undisturbed, G = grassland, M = woody patch, N = Invaded, L = Native
disturbance: disturbance of site where sample was collected
vegetation.type: vegetation type of site where sample was collected
invasion: invasion of site where sample was collected
replicate: replicate of soil community quadrat sampling
block: Site
Species names of plants within quadrat that was sampled. Abundance corresponds to % cover.

soil.data.clean.2023.csv
------------------------
sample: Sample name
block: Site
Rep: replicate
Plot_ind: Internal name of soil analyses
Treatments: Combination of disturbance, vegetation, and invasion; D = Disturbed, U = undisturbed, G = grassland, M = woody patch, N = Invaded, L = Native
disturbance: disturbance of site where sample was collected
vegetation.type: vegetation type of site where sample was collected
invasion: invasion of site where sample was collected
Beta_umolgSOC.1h.1: beta-galactosidase
Nag_umolgSOC.1h.1: NAGase
Phos_umolgSOC.1h.1: Phosphatase
Beta.Nag: Ratio of enzymes
Beta.Phos: Ratio of enzymes
Nag.Phos: Ratio of enzymes
SOC.percent: Percent soil organic carbon
Beta_Prop
Nag_Prop
Phos_Prop
SOC.percent: Percent soil organic carbon
N.percent: Percent nitrogen
H.percent: Percent hydrogen
S.percentIR: Percent sulfur
C.N: Carbon:nitrogen ratio
C.P: Carbon:phosphorus ratio
N.P: Nitrogen: phosphorus ratio
pH
ActiveC.to.TotalC: Ratio of active carbon to total carbon
Active_C_mgkg.1
C.mgperkg: Carbon
N.mgperkg: Nitrogen
P.mgperkg: Phosphorus
Mg.mgperkg: Magnesium
Mn.mgperkg: Manganese
Fe.mgperkg: Iron
K.mgperkg:Potassium
Ca.mgperkg: Calcium
Ca_NoHIGH.mgperkg
Clay.percent
Sand.percent
Silt.percent
Texture.class


### Data associated with script 1.analysisWoodlandGrassland.R. ###
data.diversity.Overall.csv
--------------------------
block: site
invasion: invasion of site where sample was collected
disturbance: disturbance of site where sample was collected
vegetation.type: vegetation type of site where sample was collected
sample: Sample name
Sequencing.name: Internal sequencing name
plot: Plot sample taken from
fungal.spec.richness: OTU richness
fungal.fisher.alpha: Fisher's alpha
bacterial.spec.richness: OTU richness
bacterial.fisher.alpha: Fisher's alpha
plant.shannon: Plant shannon's diversity
plant.invsimpson: Plant Simpson's diversity, inverse of
plant.spec.richness: Plant species richness
Beta_umolgSOC.1h.1: beta-galactosidase
Nag_umolgSOC.1h.1: NAGase
Phos_umolgSOC.1h.1: Phosphatase

data.fungal.Overall.csv and data.bacterial.Overall.csv
------------------------------------------------------
sample: Sample name
block: Site
Rep: replicate
Treatments: Combination of disturbance, vegetation, and invasion; D = Disturbed, U = undisturbed, G = grassland, M = woody patch, N = Invaded, L = Native
disturbance: disturbance of site where sample was collected
vegetation.type: vegetation type of site where sample was collected
invasion: invasion of site where sample was collected
Beta_umolgSOC.1h.1: beta-galactosidase
Nag_umolgSOC.1h.1: NAGase
Phos_umolgSOC.1h.1: Phosphatase
Beta.Nag: Ratio of enzymes
Beta.Phos: Ratio of enzymes
Nag.Phos: Ratio of enzymes
SOC.percent: Percent soil organic carbon
N.percent: Percent nitrogen
H.percent: Percent hydrogen
S.percentIR: Percent sulfur
C.N: Carbon:nitrogen ratio
C.P: Carbon:phosphorus ratio
N.P: Nitrogen: phosphorus ratio
pH
ActiveC.to.TotalC: Ratio of active carbon to total carbon
Active_C_mgkg.1
C.mgperkg: Carbon
N.mgperkg: Nitrogen
P.mgperkg: Phosphorus
Mg.mgperkg: Magnesium
Mn.mgperkg: Manganese
Fe.mgperkg: Iron
K.mgperkg:Potassium
Ca.mgperkg: Calcium
Clay.percent
Sand.percent
Silt.percent
Texture.class
log.beta: Log tranformed Beta_umolgSOC.1h.1
cube.root.nag: Cube root transformed Nag_umolgSOC.1h.1
cube.root.phos: Cube root transformed Phos_umolgSOC.1h.1:
log.C
log.N
log.P
log.Mg
log.K
log.Ca
log.Fe
Sequencing.name: Internal sequencing name
plot: Plot sample taken from
plant.shannon: Plant shannon's diversity
plant.spec.richness: Plant species richness
PCNM1: Geographic distance calculated with Moran's eigenvector; First axis
PCNM2: Geographic distance calculated with Moran's eigenvector; Second axis
Otu###: 95% sequence similarity clustered OTUs with read abundance. These have been rarified with coverage based rarefaction.

### Data associated with script 2.analysisGrasslands.R. They have a similar structure. ###
data.diversity.Grasslands.csv
data.fungal.Grasslands.csv
data.bacterial.Grasslands.csv
----------------------------
Same structure as for Overall files

data.SEM.Grassland.csv
------------------------------
sample: Sample name
invasion: invasion of site where sample was collected
disturbance: disturbance of site where sample was collected
PC1.soil: PC1 of soil properties
PC2.soil: PC2 of soil properties
PCNM1: Geographic distance calculated with p
fungal.spec.richness: OTU richness
bacterial.fisher.alpha: Fisher's alpha
MAP: Mean annual precipitation, WorldClim
plant.shannon: Plant shannon's diversity
jacc.f.nmds1: Fungal community NMDS ordination axis 1 calculated with Jaccard index
horn.b.nmds1: Bacterial community NMDS ordination axis 1 calculated with Morisita-horn index


### Data associated with script 3.analysisWoodyPatches.R. They have a similar structure. ###
data.diversity.WoodyPatches.csv
data.fungal.WoodyPatches.csv
data.bacterial.WoodyPatches.csv
------------------------------
Same structure as for Overall files


data.SEM.Woodypatch.csv
-----------------------
sample: Sample name
PC1.soil: PC1 of soil properties
PC2.soil: PC2 of soil properties
PCNM1: Geographic distance calculated with Moran's eigenvector; First axis
log.beta: Log transformed beta-galactosidase enzymes
log.nag: Log transformed NAGase enzyme
log.phos: Log transformed phosphatase enzyme
fungal.spec.richness: OTU richness
bacterial.fisher.alpha: Fisher's alpha
bacterial.spec.richness: OTU richness
MAP: Mean annual precipitation, WorldClim
plant.shannon: Plant shannon's diversity
jacc.f.nmds1: Fungal community NMDS ordination axis 1 calculated with Jaccard index
horn.b.nmds1: Bacterial community NMDS ordination axis 1 calculated with Morisita-horn index

### Data associated with script 4.analysisTaxonomy.R. ###
TaxonomyFungi.csv
------------------
Otu: OTU clustered at 95% sequence similarity
phylum: Fungal phylum classification
subphylum: Fungal subphylum classification
class: Fungal class classification
order: Fungal order classification
family: Fungal family classification
genus: Fungal genus classification
species: Fungal species classification
trophic.mode: Trophic mode of fungal group from FunGuild
guild: Guild from FunGuild

BacterialTaxonomy_Qiime.csv
---------------------------
Otu: OTU clustered at 95% sequence similarity
Taxon: Overall taxonomic classification
Confidence: Confidence score assigned by Qiime
domain: Bacterial domain classification
phylum: Bacterial phylum classification
class: Bacterial class classification
order: Bacterial order classification
family: Bacterial family classification
genus: Bacterial genus classification
species: Bacterial species classification