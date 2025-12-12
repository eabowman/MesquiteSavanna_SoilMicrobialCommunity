## Script written to assess variation in soil traits (nutrients, fungal 
## communities, and bacterial communities) in mottes.
## Written by Dr. Elizabeth Bowman, May 19, 2025
## eabowman@utexas.edu
## 
## Here, I am assessing how the distinct soil characteristics and plant
## communities found in mottes structures soil microbial communities
## (fungi and bacteria) and to what extent does invasion by Guinea grass
## and disturbance alter these patterns.

## 0. Collinearity
## 1. Soil characteristics
## 2. Fungi
## 3. Bacteria
## 4. Structural equation modeling


# Import metadata ----
site.data <- read_csv('data/field.data_climate.csv')

# Import soil data ----
soil.data <- read.csv('G022.data/Soil.data/soil.data.clean.2023.csv')

filter(soil.data, vegetation.type == 'Motte',
       sample != '3UMN.4') -> soil.m

# Import diversity data ----
div.data <- read_csv('G022.data/Diversity.data_Complete.csv')

filter(div.data, vegetation.type == 'motte',
       sample != '3UMN.4') -> div.m

## Combine soil data and div data, add PC1 of soil data ----
# combine soil data with diversity data
soil.m %>%
  dplyr::select(sample, Active_C_mgkg.1, C.mgperkg, N.mgperkg, P.mgperkg,
                Mg.mgperkg, Mn.mgperkg, Fe.mgperkg, K.mgperkg, Ca.mgperkg, pH,
                Beta_umolgSOC.1h.1, Nag_umolgSOC.1h.1,Phos_umolgSOC.1h.1) %>%
  mutate(log.beta = log(Beta_umolgSOC.1h.1),
         log.nag = log(Nag_umolgSOC.1h.1)) %>%
  right_join(., div.m) -> div.soil.m

# Add MAP to div.soil.g
site.data %>%
  dplyr::select(plot, MAP) %>%
  right_join(div.soil.m) -> div.soil.m

### Create PCA of soil traits (same traits as used in the Woodland-grassland analysis)
div.soil.m %>%
  mutate(log.C.mgperkg = log(C.mgperkg),
         log.N.mgperkg = log(N.mgperkg),
         log.P.mgperkg = log(P.mgperkg), # not normal even with transf.
         log.Mg.mgperkg = log(Mg.mgperkg),
         log.Mn.mgperkg = log(Mn.mgperkg),
         log.Fe.mgperkg = log(Fe.mgperkg),
         log.K.mgperkg = log(K.mgperkg),
         log.Ca.mgperkg = log(Ca.mgperkg),
         log.pH = log(pH)) %>%
  dplyr::select(Active_C_mgkg.1, log.C.mgperkg, log.N.mgperkg, log.P.mgperkg,
                log.Mg.mgperkg, log.Mn.mgperkg, log.Fe.mgperkg, log.K.mgperkg,
                log.Ca.mgperkg, log.pH) -> PCA.soil

# Standardize data using Z-scores (default method is "standardize")
standardized_data <- decostand(PCA.soil, method = "standardize")

# PCA
PCA.m.soil <- prcomp(standardized_data,
                     center = F, 
                     scale. = F)
summary(PCA.m.soil)
biplot(PCA.m.soil)

# add PC1 to div.soil.g
div.soil.m$PC1.soil <- PCA.m.soil$x[,1]
div.soil.m$PC2.soil <- PCA.m.soil$x[,2]

# factor invasion, disturbance, and plant community
div.soil.m$invasion <- factor(div.soil.m$invasion,
                              levels = c('invaded','uninvaded'),
                              labels = c('Invaded', 'Uninvaded'))
div.soil.m$disturbance <- factor(div.soil.m$disturbance,
                                 levels = c('disturbed', 'undisturbed'),
                                 labels = c('Disturbed', 'Undisturbed'))
div.soil.m$vegetation.type <- factor(div.soil.m$vegetation.type,
                                     levels = c('grassland', 'motte'),
                                     labels = c('Grassland', 'Motte'))

# Import fungal data ----
fung.data <- read.csv('G022.data/Fungal.data/Fungal_SitexSpecies_95sim_rarefied.csv')

# Separate out the Motte data
filter(fung.data, vegetation.type == 'motte',
       sample != '3UMN.4') -> fung.m

## Standardize soil data and distance data (dbMEM)----
div.soil.m %>%
  dplyr::select(sample, PC1.soil, PC2.soil, log.beta, log.nag,
                Phos_umolgSOC.1h.1,
                fungal.fisher.alpha, fungal.spec.richness,
                bacterial.fisher.alpha, bacterial.spec.richness,
                MAP) %>%
  right_join(fung.m) -> fung.m

### distance based Moran's eigenvector (dbMEM) -----
# Load varpart2.MEM.R function file
# varpart2.MEM.R was written by Legendre et al. 2012 'Variation partitioning
# involving orthogonal spatial eigenfunction submodels'
source('G022.scripts/varpart2.MEM.R')

# isolate location data
fung.xy <- fung.m[c('lat','long')]

### PCNM
# Create distance matrix of spatial data
fung.dist <- dist(fung.xy, method = 'euclidean')

# transform spatial data 
fung.pcnm <- pcnm(fung.dist)

# access eigenvectors with scores() function

# map spatial data
op <- par(mfrow = c(1,2))
ordisurf(fung.xy, scores(fung.pcnm, choi=1), bubble = 4, main = "PCNM 1")
ordisurf(fung.xy, scores(fung.pcnm, choi=2), bubble = 4, main = "PCNM 1")
par(op)

fung.m$PCNM1 <- scores(fung.pcnm, choices = 1)
fung.m$PCNM2 <- scores(fung.pcnm, choices = 2)

fung.m[c(1:21, 1373:1374, 22:1372)] -> fung.m

# make block a categorical variable
fung.m$block <- factor(fung.m$block)

# factor invasion, disturbance, and plant community
fung.m$invasion <- factor(fung.m$invasion,
                          levels = c('invaded','uninvaded'),
                          labels = c('Invaded', 'Uninvaded'))
fung.m$disturbance <- factor(fung.m$disturbance,
                             levels = c('disturbed', 'undisturbed'),
                             labels = c('Disturbed', 'Undisturbed'))
fung.m$vegetation.type <- factor(fung.m$vegetation.type,
                                 levels = c('grassland', 'motte'),
                                 labels = c('Grassland', 'Motte'))

# Import bacterial data ----
bac.data <- read.csv('G022.data/Bacterial.data/Bacterial_SitexSpecies_99sim_50occ.csv')

# Separate out the Motte data
filter(bac.data, vegetation.type == 'motte',
       sample != '3UMN.4') -> bac.m

## Standardize soil data and distance data (dbMEM)----
div.soil.m %>%
  dplyr::select(sample, PC1.soil, PC2.soil, log.beta, log.nag,
                Phos_umolgSOC.1h.1,
                fungal.fisher.alpha, fungal.spec.richness,
                bacterial.fisher.alpha, bacterial.spec.richness,
                MAP) %>%
  right_join(bac.m) -> bac.m

### distance based Moran's eigenvector (dbMEM) -----
# Load varpart2.MEM.R function file
# varpart2.MEM.R was written by Legendre et al. 2012 'Variation partitioning
# involving orthogonal spatial eigenfunction submodels'
source('G022.scripts/varpart2.MEM.R')

# isolate location data
bac.xy <- bac.m[c('lat','long')]

### PCNM
# Create distance matrix of spatial data
bac.dist <- dist(bac.xy, method = 'euclidean')

# transform spatial data 
bac.pcnm <- pcnm(bac.dist)

# access eigenvectors with scores() function

# map spatial data
op <- par(mfrow = c(1,2))
ordisurf(bac.xy, scores(bac.pcnm, choi=1), bubble = 4, main = "PCNM 1")
ordisurf(bac.xy, scores(bac.pcnm, choi=2), bubble = 4, main = "PCNM 1")
par(op)

bac.m$PCNM1 <- scores(bac.pcnm, choices = 1)
bac.m$PCNM2 <- scores(bac.pcnm, choices = 2)

bac.m[c(1:21, 5684:5685, 22:5683)]  -> bac.m

# make block a categorical variable
bac.m$block <- factor(bac.m$block)

# factor invasion, disturbance, and plant community
bac.m$invasion <- factor(bac.m$invasion,
                         levels = c('invaded','uninvaded'),
                         labels = c('Invaded', 'Uninvaded'))
bac.m$disturbance <- factor(bac.m$disturbance,
                            levels = c('disturbed', 'undisturbed'),
                            labels = c('Disturbed', 'Undisturbed'))
bac.m$vegetation.type <- factor(bac.m$vegetation.type,
                                levels = c('grassland', 'motte'),
                                labels = c('Grassland', 'Motte'))

#--------------------------------------------------------------#
# Pre-made data ----
#--------------------------------------------------------------#
fung.m <- read.csv('data/data.fungal.WoodyPatches.csv')
bac.m <- read.csv('data/data.bacterial.WoodyPatches.csv')
div.soil.m <- read.csv('data/data.diversity.WoodyPatches.csv')

#--------------------------------------------------------------#
# 0. Assess collinearity of data ----
#--------------------------------------------------------------#

## rcorr ----
cor.data <- fung.m[, c("PC1.soil", "PC2.soil", "MAP", "plant.shannon",
                       "fungal.fisher.alpha", "bacterial.fisher.alpha")]
Hmisc::rcorr(as.matrix(cor.data),
             type = "spearman")

#                           PC1.soil PC2.soil MAP   plant.shannon fungal.fisher.alpha bacterial.fisher.alpha
# PC1.soil                   1.00     0.00 -0.17          0.06               -0.07                  -0.02
# PC2.soil                   0.00     1.00  0.23          0.37                0.22                   0.51
# MAP                       -0.17     0.23  1.00         -0.06                0.08                   0.06
# plant.shannon              0.06     0.37 -0.06          1.00                0.11                   0.14
# fungal.fisher.alpha       -0.07     0.22  0.08          0.11                1.00                   0.48
# bacterial.fisher.alpha    -0.02     0.51  0.06          0.14                0.48                   1.00
# 
# n= 59 
# 
# 
# P
#                         PC1.soil PC2.soil MAP    plant.shannon fungal.fisher.alpha bacterial.fisher.alpha
# PC1.soil                        0.9737   0.2031 0.6254        0.6053              0.8550                
# PC2.soil               0.9737            0.0816 0.0036        0.0895              0.0000                
# MAP                    0.2031   0.0816          0.6592        0.5351              0.6379                
# plant.shannon          0.6254   0.0036   0.6592               0.4194              0.3038                
# fungal.fisher.alpha    0.6053   0.0895   0.5351 0.4194                            0.0001                
# bacterial.fisher.alpha 0.8550   0.0000   0.6379 0.3038        0.0001   

#--------------------------------------------------------------#
# 1. Soil characteristics----
#--------------------------------------------------------------#

## Create PCA of soil traits (same traits as used in the Woodland-grassland analysis) ----
div.soil.m %>%
  mutate(log.C.mgperkg = log(C.mgperkg),
         log.N.mgperkg = log(N.mgperkg),
         log.P.mgperkg = log(P.mgperkg), # not normal even with transf.
         log.Mg.mgperkg = log(Mg.mgperkg),
         log.Mn.mgperkg = log(Mn.mgperkg),
         log.Fe.mgperkg = log(Fe.mgperkg),
         log.K.mgperkg = log(K.mgperkg),
         log.Ca.mgperkg = log(Ca.mgperkg),
         log.pH = log(pH)) %>%
  dplyr::select(Active_C_mgkg.1, log.C.mgperkg, log.N.mgperkg, log.P.mgperkg,
                log.Mg.mgperkg, log.Mn.mgperkg, log.Fe.mgperkg, log.K.mgperkg,
                log.Ca.mgperkg, log.pH) -> PCA.soil

# Standardize data using Z-scores (default method is "standardize")
standardized_data <- decostand(PCA.soil, method = "standardize")

# PCA
PCA.m.soil <- prcomp(standardized_data,
                     center = F, 
                     scale. = F)
summary(PCA.m.soil)
biplot(PCA.m.soil)

# Extract loadings for PC1 and PC2
loadings <- as.data.frame(PCA.m.soil$rotation[, 1:2])
loadings$trait <- rownames(loadings)  # Add variable names for labeling
loadings$trait <- c("ActiveC", "C", "N", "P", "Mg", "Mn", "Fe", "K", "Ca", "pH") # Shorten names

# Scale loadings (adjust the multiplier if arrows are too short or too long)
arrow_multiplier <- 5
loadings <- loadings %>%
  mutate(PC1 = PC1 * arrow_multiplier, PC2 = PC2 * arrow_multiplier)

# Function to calculate confidence ellipse
calc_ellipse <- function(x, y, conf = 0.95) {
  ellipse::ellipse(cbind(x, y), level = conf)
}

### Plot: Vegetation type ----
ggplot(fung.m,
       aes(x = PC1.soil, y = PC2.soil,
           fill = plant.shannon,
           shape = invasion)) +
  geom_point(size = 4, color = "black") +  # black border
  xlab("PC1 (60.3%)") +
  ylab("PC2 (13.9%)") +
  # Add arrows for soil traits
  geom_segment(data = loadings,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "black", inherit.aes = F) +
  # Add labels for soil traits
  ggrepel::geom_text_repel(data = loadings,
                           aes(x = PC1, y = PC2, label = trait),
                           color = "black", size = 5,
                           max.overlaps = 20,
                           box.padding = 0.5,
                           point.padding = 0.2,
                           segment.color = "darkred",
                           inherit.aes = FALSE) +
  scale_fill_gradient(low = "#F1BB83", high = "grey50") +
  scale_shape_manual(values = c(21, 24)) +  # both support fill + border
  labs(fill = "Plant Shannon's\ndiversity",
       shape = "") +
  theme_classic() +
  theme(axis.text = element_text(size = 12, color = 'black'),
        axis.title = element_text(size = 14),
        strip.text.x = element_text(size = 12),
        legend.position = 'none',
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave('figures/Fig3a.tiff', plot = last_plot(),
       device = 'tiff', width = 5, height = 4)

## PC1 as a function of MAP and Plant diversity ----
PC1.lm <- lm(PC1.soil ~ MAP + plant.shannon,
             # random = ~ 1 | block,
             data = fung.m)
summary(PC1.lm)
plot(PC1.lm$residuals)
shapiro.test(PC1.lm$residuals)

## PC2 as a function of MAP and Plant diversity ----
PC2.lm <- lm(PC2.soil ~ MAP + plant.shannon,
              # random = ~ 1 | block,
              data = fung.m)
summary(PC2.lm)
plot(PC2.lm$residuals)
shapiro.test(PC2.lm$residuals)

# Quadratic
lm_poly <- lm(PC2.soil ~ poly(plant.shannon, 2) + MAP, data = fung.m)
summary(lm_poly)
plot(lm_poly$residuals)
shapiro.test(lm_poly$residuals)

ggplot(fung.m, aes(x = plant.shannon,
                   y = PC2.soil)) +
  geom_point(size = 3, shape = 21, fill = '#F1BB83', color = 'black') +
  geom_smooth(method = "loess", se = FALSE, color = 'darkgrey') +
  labs(y = 'Soil PC2 (13.9%)',
       x = "Plant Shannon's diversity") +
  theme_classic() +
  theme(axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 12,
                                   color = 'black'),
        axis.text.x = element_text(size = 12,
                                   color = 'black'),
        legend.position = 'right',
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave('figures/Fig3b.tiff', plot = last_plot(),
       device = 'tiff', width = 6, height = 4)

## Enzymatic activity ----

### Beta galactosidase ----
beta.lm <- lm(log.beta ~ fungal.fisher.alpha + bacterial.fisher.alpha +
                 MAP + plant.shannon + PC1.soil,
               # random = ~ 1 | block,
               data = fung.m)
summary(beta.lm)
plot(beta.lm$residuals)
shapiro.test(beta.lm$residuals)

ggplot(fung.m, aes(x = bacterial.fisher.alpha,
                   y = log.beta)) +
  geom_point(aes(fill = MAP), size = 3, shape = 21, color = 'black') +
  scale_fill_gradient(low = "white", high = "#f0a354") +  # change colors here
  geom_smooth(method = "lm", se = FALSE, color = 'darkgrey') +
  ylab(expression(paste("Log ",beta, "-glucosidase (", mu, "mol g",
                        SOC^-1, hr^-1, ")"))) +
  labs(x = "Bacterial Fisher's Alpha",
       fill = "MAP", ) +
  theme_classic() +
  theme(axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 12,
                                   color = 'black'),
        axis.text.x = element_text(size = 12,
                                   color = 'black'),
        legend.position = 'right',
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave('figures/WoodyBetaGlucoside.tiff', plot = last_plot(),
       device = 'tiff', width = 6, height = 4)

### Nag----
# remove sample with -Inf
fung.m %>%
  filter(log.nag != '-Inf') -> fung.m.nag

nag.lm <- lm(log.nag ~ fungal.fisher.alpha + bacterial.fisher.alpha +
                MAP + plant.shannon + PC1.soil,
              #random = ~ 1 | block,
              data = fung.m.nag)
summary(nag.lm)
plot(nag.lm$residuals)
shapiro.test(nag.lm$residuals)

ggplot(fung.m.nag, aes(x = MAP,
                   y = log.nag)) +
  geom_point(size = 3, shape = 21, fill = '#f0a354', color = 'black') +
  geom_smooth(method = "lm", se = FALSE, color = 'darkgrey') +
  ylab(expression(paste("Log NAGase (",mu, "mol g",
                        SOC^-1, hr^-1, ")"))) +
  labs(x = "Mean annual precipitaiton (MAP)") +
  theme_classic() +
  theme(axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 12,
                                   color = 'black'),
        axis.text.x = element_text(size = 12,
                                   color = 'black'),
        legend.position = 'right',
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave('figures/WoodyNAGase.tiff', plot = last_plot(),
       device = 'tiff', width = 6, height = 4)

### Phos----
fung.m %>%
  mutate(log.phos = log(Phos_umolgSOC.1h.1)) %>%
  filter(log.phos > 6) -> fung.m.phos

phos.lm <- lm(log.phos ~ fungal.fisher.alpha + bacterial.fisher.alpha +
                 MAP + plant.shannon + PC1.soil,
                # random = ~ 1 | block,
               data = fung.m.phos)
summary(phos.lm)
plot(phos.lm$residuals)
shapiro.test(phos.lm$residuals)

ggplot(fung.m.phos, aes(x = bacterial.fisher.alpha,
                        y = log.phos)) +
  geom_point(aes(fill = PC1.soil), size = 3, shape = 21, color = 'black') +
  geom_smooth(method = "lm", se = FALSE, color = 'darkgrey') +
  scale_fill_gradient(low = "white", high = "#f0a354") +  # change colors here
  ylab(expression(paste("Log Phosphatase (",mu, "mol g",
                        SOC^-1, hr^-1, ")"))) +
  labs(x = "Bacterial Fisher's alpha",
       fill = "Soil PC1") +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 12,
                                   color = 'black'),
        axis.text.x = element_text(size = 12,
                                   color = 'black'),
        legend.position = 'right',
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave('figures/Fig3c.tiff', plot = last_plot(),
       device = 'tiff', width = 6, height = 4)

#--------------------------------------------------------------#
# 2. Fungi ----
#--------------------------------------------------------------#

## Species richness and diversity ----

# Check normality
hist(div.soil.m$fungal.spec.richness)
shapiro.test(div.soil.m$fungal.spec.richness)  
hist(div.soil.m$fungal.fisher.alpha)
shapiro.test(div.soil.m$fungal.fisher.alpha)  

### Species richness: normal ----
### Fine
fung.sr.lm <- lm(fungal.spec.richness ~ plant.shannon + PC2.soil + MAP,
                 #random = ~ 1 | block,
                 data = fung.m)
summary(fung.sr.lm)
hist(fung.sr.lm$residuals)
shapiro.test(fung.sr.lm$residuals)

### Fisher's alpha: normal ----
### Fine
fung.div.lm <- lm(fungal.fisher.alpha ~ plant.shannon + PC2.soil + MAP,
                  #random = ~ 1 | block,
                  data = fung.m)
summary(fung.div.lm)
hist(fung.div.lm$residuals)
shapiro.test(fung.div.lm$residuals)

## NMDS and PERMANOVA: Fungi ----
# isolate fungal community data
fung.comm <- dplyr::select(fung.m, starts_with('Otu'))
# Remove singletons
fung.comm <- fung.comm[colSums(fung.comm) > 1]

### PERMANOVA: Jaccard----
# Constrains permutations to blocks (aka pastures)
perm <- how(nperm = 199)
setBlocks(perm) <- with (fung.m, block)

# fine
fung.jac.adonis <- adonis2(fung.comm ~ PC1.soil * plant.shannon * MAP + PCNM1,
                           method = 'jaccard',
                           by = 'terms',
                           data = fung.m,
                           permutations = perm)

fung.jac.adonis.df <- as.data.frame(fung.jac.adonis)
write.csv(fung.jac.adonis.df,
          'results/Table3_Fungal_JaccardPERMANOVA_Grasslands.csv',
          row.names = T)

### NMDS: Jaccard ----
# create distance matrix
jacc.dist <- vegdist(fung.comm, method = 'jaccard', binary = T)

# GGplot: Jaccard
jacc.mds <- metaMDS(jacc.dist, dist = 'bray',
                    try = 1000, trymax = 1000)
jacc.stress <- jacc.mds$stress

# format data for plot
data.scores <- data.frame(NMDS1 = jacc.mds$points[,1],
                          NMDS2 = jacc.mds$points[,2],
                          inv = fung.m$invasion,
                          dist = fung.m$disturbance,
                          MAP = fung.m$MAP,
                          soil = fung.m$PC1.soil,
                          plant.div = fung.m$plant.shannon,
                          geo.dist = fung.m$PCNM1)
# install.packages("viridis")  # Install
library("viridis")  

jacc.plot <- ggplot() + 
  geom_point(data = data.scores, aes(x = NMDS1,
                                     y = NMDS2,
                                     color = soil,
                                     shape = inv),
             size = 4) +
  #scale_fill_manual(values = c('darkblue', 'lightblue')) + # for invasion
  scale_color_viridis() +
  coord_equal() +
  theme_classic() +
  labs(shape = "",
       color = "Soil PC1") +  
  theme(axis.text = element_text(size = 12,
                                 colour = 'black'),
        axis.title = element_text(size = 14),
        legend.position = 'none',
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave('figures/Fig3e.tiff',
       device = 'tiff',
       plot = jacc.plot,
       width = 6, height = 5, units = 'in')

# Isolate outlier invaded communitys with a high soil.1 value
data.scores %>%
  filter(inv == 'Uninvaded' & NMDS1 > 0.5) -> odd.uninvaded

div.soil.m[div.soil.m$PC1.soil %in% odd.uninvaded$soil,] -> odd.uninvaded.id

### PERMANOVA: Morisita-Horn ----
# Constrains permutations to blocks (aka pastures)
perm <- how(nperm = 199)
setBlocks(perm) <- with (fung.m, block)

# fine
fung.horn.adonis <- adonis2(fung.comm ~ PC1.soil * plant.shannon * MAP + PCNM1,
                            method = 'horn',
                            by = 'terms',
                            data = fung.m,
                            permutations = perm)

fung.horn.adonis.df <- as.data.frame(fung.horn.adonis)
write.csv(fung.horn.adonis.df, 'results/Table3_Fungal_MorisitaHornPERMANOVA_Mottes.csv',
          row.names = T)

### NMDS: Morisita-Horn ----
# create distance matrix
horn.dist <- vegdist(fung.comm, method = 'horn', binary = F)

# GGplot: Jaccard
horn.mds <- metaMDS(horn.dist, dist = 'bray',
                    try = 1000, trymax = 1000)
horn.stress <- horn.mds$stress

# format data for plot
data.scores <-  data.frame(NMDS1 = jacc.mds$points[,1],
                           NMDS2 = jacc.mds$points[,2],
                           inv = fung.m$invasion,
                           dist = fung.m$disturbance,
                           MAP = fung.m$MAP,
                           soil = fung.m$PC1.soil,
                           plant.div = fung.m$plant.shannon,
                           geo.dist = fung.m$PCNM1)

horn.plot <- ggplot() + 
  geom_point(data = data.scores, aes(x = NMDS1,
                                     y = NMDS2,
                                     color = soil,
                                     shape = inv),
             size = 4) +
  #scale_fill_manual(values = c('darkblue', 'lightblue')) + # for invasion
  scale_color_viridis() +
  coord_equal() +
  labs(color = "Soil PC1",
       shape = "") +
  theme_classic() +
  theme(axis.text = element_text(size = 14,
                                 colour = 'black'),
        axis.title = element_text(size = 16),
        legend.position = 'bottom',
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave('figures/SupplementaryFigS6b.jpeg',
       device = 'jpeg',
       plot = horn.plot,
       width = 6, height = 5, units = 'in')

#--------------------------------------------------------------#
# 3. Bacteria ----
#--------------------------------------------------------------#

## Be sure to run the section "Combine soil data and div data, add PC1
## of soil data" above under 2. Fungi section

## Species richness and diversity ----

# Check normality
hist((div.soil.m$bacterial.spec.richness)^2)
shapiro.test((div.soil.m$bacterial.spec.richness)^2)
hist(div.soil.m$bacterial.fisher.alpha)
shapiro.test(div.soil.m$bacterial.fisher.alpha)  

div.soil.m %>% mutate(sq.root.bacterial.spec.richness = (div.soil.m$bacterial.spec.richness)^2) -> div.soil.m

### Species richness ----

### Fine
bac.sr.lm <- lm(sq.root.bacterial.spec.richness ~ plant.shannon + PC2.soil + MAP,
                #random = ~ 1 | block,
                data = div.soil.m)
summary(bac.sr.lm)
hist(bac.sr.lm$residuals)
shapiro.test(bac.sr.lm$residuals)

ggplot(div.soil.m, aes(x = PC2.soil,
                        y = bacterial.spec.richness)) +
  geom_point(size = 3, shape = 21, color = 'black', fill = '#f0a354') +
  geom_smooth(method = "lm", se = FALSE, color = 'darkgrey') +
  # scale_fill_gradient(low = "white", high = "#f0a354") +  # change colors here
  ylab(expression(paste("Log Phosphatase (",mu, "mol g",
                        SOC^-1, hr^-1, ")"))) +
  labs(x = "Soil PC2",
       y = "Bacterial OTU richness") +
  theme_classic() +
  theme(axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 12,
                                   color = 'black'),
        axis.text.x = element_text(size = 12,
                                   color = 'black'),
        legend.position = 'right',
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave('figures/Woody_BacOTUrich.tiff', plot = last_plot(),
       device = 'tiff', width = 6, height = 4)

### Fisher's alpha ----
### Fine
bac.div.lm <- lm(bacterial.fisher.alpha ~ plant.spec.richness + PC2.soil + MAP,
                 #random = ~ 1 | block,
                 data = div.soil.m)
summary(bac.div.lm)
hist(bac.div.lm$residuals)
shapiro.test(bac.div.lm$residuals)

ggplot(div.soil.m, aes(x = PC2.soil,
                       y = bacterial.fisher.alpha)) +
  geom_point(size = 3, shape = 21, color = 'black', fill = '#f0a354') +
  geom_smooth(method = "lm", se = FALSE, color = 'darkgrey') +
  # scale_fill_gradient(low = "white", high = "#f0a354") +  # change colors here
  ylab(expression(paste("Log Phosphatase (",mu, "mol g",
                        SOC^-1, hr^-1, ")"))) +
  labs(x = "Soil PC2",
       y = "Bacterial Fisher's alpha") +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 12,
                                   color = 'black'),
        axis.text.x = element_text(size = 12,
                                   color = 'black'),
        legend.position = 'right',
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave('figures/Fig3d.tiff', plot = last_plot(),
       device = 'tiff', width = 6, height = 4)

## NMDS and PERMANOVA: Bacteria ----
# isolate fungal community data
bac.comm <- dplyr::select(bac.m, starts_with('Otu'))
# Remove singletons
bac.comm <- bac.comm[colSums(bac.comm) > 1]

### PERMANOVA: Jaccard----
# Constrains permutations to blocks (aka pastures)
perm <- how(nperm = 199)
setBlocks(perm) <- with(bac.m, block)

# fine
bac.jac.adonis <- adonis2(bac.comm ~ PC1.soil * plant.shannon * MAP + PCNM1,
                          method = 'jaccard',
                          by = 'terms',
                          data = bac.m,
                          permutations = perm)

bac.jac.adonis.df <- as.data.frame(bac.jac.adonis)
write.csv(bac.jac.adonis.df,
          'results/Table3_Bacterial_JaccardPERMANOVA_Grasslands.csv',
          row.names = T)

### NMDS: Jaccard ----
# create distance matrix
jacc.dist <- vegdist(bac.comm, method = 'jaccard', binary = T)

# GGplot: Jaccard
jacc.mds <- metaMDS(jacc.dist, dist = 'bray',
                    try = 1000, trymax = 1000)
jacc.stress <- jacc.mds$stress

# format data for plot
data.scores <- data.frame(NMDS1 = jacc.mds$points[,1],
                          NMDS2 = jacc.mds$points[,2],
                          inv = fung.m$invasion,
                          dist = fung.m$disturbance,
                          MAP = fung.m$MAP,
                          soil = fung.m$PC1.soil,
                          plant.div = fung.m$plant.shannon,
                          geo.dist = fung.m$PCNM1)

jacc.plot <- ggplot() + 
  geom_point(data = data.scores, aes(x = NMDS1,
                                     y = NMDS2,
                                     color = soil,
                                     shape = inv),
             size = 5) +
  #scale_fill_manual(values = c('darkblue', 'lightblue')) + # for invasion
  scale_color_viridis() +
  coord_equal() +
  theme_classic() +
  labs(color = 'Soil PC1',
       shape = "") +
  theme(axis.text = element_text(size = 14,
                                 colour = 'black'),
        axis.title = element_text(size = 16),
        legend.position = 'right',
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave('figures/SupplementaryFig6d.tiff',
       device = 'tiff',
       plot = jacc.plot,
       width = 6, height = 5, units = 'in')

### PERMANOVA: Morisita-Horn ----
# Constrains permutations to blocks (aka pastures)
perm <- how(nperm = 199)
setBlocks(perm) <- with(bac.g, block)

# fine
bac.horn.adonis <- adonis2(bac.comm ~  PC1.soil * plant.shannon * MAP + PCNM1,
                           method = 'horn',
                           by = 'terms',
                           data = bac.m,
                           permutations = perm)

bac.horn.adonis.df <- as.data.frame(bac.horn.adonis)
write.csv(bac.horn.adonis.df,
          'results/Table3_Bacterial_MorisitaHornPERMANOVA_Mottes.csv',
          row.names = T)

### NMDS: Morisita-Horn ----
# create distance matrix
horn.dist <- vegdist(bac.comm, method = 'horn', binary = F)

# GGplot: Morisita-Horn
horn.mds <- metaMDS(horn.dist, dist = 'bray',
                    try = 1000, trymax = 1000)
horn.stress <- horn.mds$stress

# format data for plot
data.scores <- data.frame(NMDS1 = horn.mds$points[,1],
                          NMDS2 = horn.mds$points[,2],
                          inv = fung.m$invasion,
                          dist = fung.m$disturbance,
                          MAP = fung.m$MAP,
                          soil = fung.m$PC1.soil,
                          plant.div = fung.m$plant.shannon,
                          geo.dist = fung.m$PCNM1)

horn.plot <- ggplot() + 
  geom_point(data = data.scores, aes(x = NMDS1,
                                     y = NMDS2,
                                     color = soil,
                                     shape = inv),
             size = 4) +
  #scale_fill_manual(values = c('darkblue', 'lightblue')) + # for invasion
  scale_color_viridis() +
  coord_equal() +
  theme_classic() +
  labs(color = 'Soil PC1',
       shape = "") +
  theme(axis.text = element_text(size = 12,
                                 colour = 'black'),
        axis.title = element_text(size = 14),
        legend.position = 'none',
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave('figures/Fig3f.tiff',
       device = 'tiff',
       plot = horn.plot,
       width = 6, height = 5, units = 'in')


#--------------------------------------------------------------#
# 4. Structural equation modeling ----
#--------------------------------------------------------------#

## Read in data from above ----
sem.data <- read.csv('data/data.SEM.Woodypatch.csv')

sem.data[-1] -> sem.data
sem.data.std <- as.data.frame(scale(sem.data))

## Define the model ----
sem_model_woody <- '
  # Soil axes as function of predictors (optional if not latent)
  PC1.soil ~ a * MAP + b * plant.shannon + c * PCNM1
  PC2.soil ~ d * MAP + e * plant.shannon + f * PCNM1

  # Enzyme activity driven by both PCs and MAP
  log.phos ~ g * MAP + h * PC1.soil + i * PC2.soil

  # Microbial diversity
  bacterial.fisher.alpha ~ j * PC1.soil + k * PC2.soil + l * log.phos
  bacterial.spec.richness ~ m * PC1.soil + n * PC2.soil

  # Community dissimilarity
  jacc.f.nmds1 ~ o * PC1.soil + p * PC2.soil + q * plant.shannon + r * PCNM1 + s * MAP
  horn.b.nmds1 ~ t * PC1.soil + u * PC2.soil + v * plant.shannon + w * PCNM1 + x * MAP

  # Indirect effects from MAP
  indirect_MAP_PC2_bact_alpha := d * k    # MAP → PC2 → Bacterial Fisher
  indirect_MAP_phos_bact_alpha := g * l   # MAP → Phosphatase → Bacterial Fisher
  indirect_MAP_total_bact_alpha := indirect_MAP_PC2_bact_alpha + indirect_MAP_phos_bact_alpha

  # Indirect effects from PC2
  indirect_PC2_phos_bact_alpha := i * l   # PC2 → Phosphatase → Bacterial Fisher

  # Indirect effects from PC2 → diversity
  indirect_PC2_to_richness := n   # direct path is significant, no mediator needed

  # Indirect effects from PC2 → community structure
  indirect_PC2_to_jacc := p             # direct
  indirect_PC2_to_horn := u             # direct
'

## Fit the model ----
fit <- sem(sem_model_woody, data = sem.data.std, missing = "ML")
summary(fit, fit.measures = TRUE, standardized = TRUE, rsquare = TRUE)

## Visualize SEM ----

# Get standardized parameter estimates
estimates <- parameterEstimates(fit, standardized = TRUE)
# Filter for significant pathways
sig_paths <- subset(estimates, op == "~" & pvalue < 0.05)
# Create a named vector of edge colors
edge_colors <- ifelse(sig_paths$std.all < 0, "red", "blue")
names(edge_colors) <- paste0(sig_paths$lhs, "~", sig_paths$rhs)

labels = list(#invasion = "Invasion",
  plant.shannon = "Plant diversity",
  PC1.soil = "Soil PC1",
  PC2.soil = "Soil PC2",
  log.phos = 'Phosphatase',
  MAP = "MAP",
  PCNM1 = "Geographical\ndistance",
  # fungal.spec.richness = "Fungal OTU richenss",
  bacterial.fisher.alpha = "Bacterial\ndiversity",
  bacterial.spec.richness = "Bacterial\nOTU richness",
  jacc.f.nmds1 = "Fungal\ncomm. composition",
  horn.b.nmds1 = "Bacterial\ncomm. composition")

lavaanPlot(model = fit,
           labels = labels,
           stand = T,
           coefs = T,
           covs = F,
           sig = 0.05,
           node_options = list(shape = "box",
                               fontname = "Helvetica"),
           edge_options = list(color = "darkgrey", penwidth = 1),
           stars = c('regress'))
