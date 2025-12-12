## Script written to assess variation in soil traits (nutrients, fungal 
## communities, and bacterial communities) in mottes.
## Written by Dr. Elizabeth Bowman, May 19, 2025
## eabowman@utexas.edu
## 
## Here, I am assessing how the distinct soil characteristics and plant
## communities found in grasslands structures soil microbial communities
## (fungi and bacteria) and to what extent does invasion by Guinea grass
## and disturbance alter these patterns.

## 0. Collinearity of explanatory variables
## 1. Soil characteristics
## 2. Fungi
## 3. Bacteria
## 4. Structural equation model

fung.g <- read.csv('data/data.fungal.Grasslands.csv')
bac.g <- read.csv('data/data.bacterial.Grasslands.csv')
div.soil.g <- read.csv('data/data.diversity.Grasslands.csv')

#--------------------------------------------------------------#
# 0. Assess collinearity of data ----
#--------------------------------------------------------------#

## rcorr ----
cor.data <- fung.g[, c("PC1.soil", "PC2.soil", "MAP", "plant.shannon",
                       "fungal.fisher.alpha", "bacterial.fisher.alpha")]
Hmisc::rcorr(as.matrix(cor.data),
             type = "spearman")

#                        PC1.soil PC2.soil MAP    plant.shannon fungal.fisher.alpha bacterial.fisher.alpha
# PC1.soil                        0.2393   0.5713 0.0000        0.0237              0.0020                
# PC2.soil               0.2393            0.1812 0.7480        0.1352              0.0904                
# MAP                    0.5713   0.1812          0.8060        0.6390              0.0002                
# plant.shannon          0.0000   0.7480   0.8060               0.0189              0.0704                
# fungal.fisher.alpha    0.0237   0.1352   0.6390 0.0189                            0.0560                
# bacterial.fisher.alpha 0.0020   0.0904   0.0002 0.0704        0.0560                       

## Calculate residuals ----

### Fungi
pc1.div.lm <- glm(plant.shannon ~ PC1.soil, data = fung.g)
summary(pc1.div.lm)

shapiro.test(pc1.div.lm$residuals)
hist(pc1.div.lm$residuals)

fung.g$plant.shannon.residuals <- pc1.div.lm$residuals

fung.g[c(1:23, 1375, 24:1374)] -> fung.g

### Bacteria
pc1.div.b.lm <- glm(plant.shannon ~ PC1.soil, data = bac.g)
summary(pc1.div.b.lm)

shapiro.test(pc1.div.b.lm$residuals)
hist(pc1.div.b.lm$residuals)

bac.g$plant.shannon.residuals <- pc1.div.b.lm$residuals

bac.g[c(1:23, 5686, 24:5685)] -> bac.g

#--------------------------------------------------------------#
# 1. Soil characteristics----
#--------------------------------------------------------------#

## PCA of soil traits with PC1 and PC2 ----
### Create PCA of soil traits (same traits as used in the Woodland-grassland analysis)
div.soil.g %>%
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
PCA.g.soil <- prcomp(standardized_data,
                     center = F, 
                     scale. = F)
summary(PCA.g.soil)
biplot(PCA.g.soil)

### Plot: Vegetation type ----
# Extract loadings for PC1 and PC2
loadings <- as.data.frame(PCA.g.soil$rotation[, 1:2])
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

ggplot(fung.g,
       aes(x = PC1.soil, y = PC2.soil,
           fill = plant.shannon,
           shape = invasion)) +
  geom_point(size = 4, color = "black") +  # black border
  xlab("PC1 (43.4%)") +
  ylab("PC2 (15.8%)") +
  # Add arrows for soil traits
  geom_segment(data = loadings,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "black", inherit.aes = FALSE) +
  # Add trait labels
  ggrepel::geom_text_repel(data = loadings,
                           aes(x = PC1, y = PC2, label = trait),
                           color = "black", size = 5,
                           max.overlaps = 20,
                           box.padding = 0.5,
                           point.padding = 0.2,
                           segment.color = "darkred",
                           inherit.aes = FALSE) +
  scale_fill_gradient(low = "white", high = "#905971") +
  scale_shape_manual(values = c(21, 24)) +  # both support fill + border
  labs(fill = "Plant Shannon's\ndiversity",
       shape = "") +
  theme_classic() +
  theme(axis.text = element_text(size = 12, color = 'black'),
        axis.title = element_text(size = 14),
        strip.text.x = element_text(size = 12),
        #legend.position = 'none',
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave('figures/Fig2a.tiff', plot = last_plot(),
       device = 'tiff', width = 5, height = 4)

## PC1 of soil traits ----
PC1.lm <- lm(PC1.soil ~ MAP + plant.shannon,
              #random = ~ 1 | block,
              data = fung.g)
summary(PC1.lm)
plot(PC1.lm$residuals)
shapiro.test(PC1.lm$residuals)

ggplot(fung.g, aes(x = plant.shannon,
                   y = PC1.soil)) +
  geom_point(size = 3, shape = 21, color = 'black', fill = '#905971') +
  geom_smooth(method = "lm", se = FALSE, color = 'darkgrey') +
  labs(y = 'Soil PC1 (43.4%)',
       x = "Plant Shannon's diversity" ) +
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

ggsave('figures/Fig2b.tiff', plot = last_plot(),
       device = 'tiff', width = 6, height = 4)

## PC2 of soil traits ----
PC2.lm <- lm(PC2.soil ~ MAP + plant.shannon,
             #random = ~ 1 | block,
             data = fung.g)
summary(PC2.lm)
plot(PC2.lm)
shapiro.test(PC2.lm$residuals)

## Enzymatic activity ----

### Beta galactosidase ----
beta.lm <- lm(log.beta ~ fungal.fisher.alpha + bacterial.fisher.alpha +
                 MAP + plant.shannon.residuals + PC1.soil,
               # random = ~ 1 | block,
               data = fung.g)
summary(beta.lm)
plot(beta.lm$residuals)
shapiro.test(beta.lm$residuals)

### Nag----
nag.lm <- lm(log.nag ~ fungal.fisher.alpha + bacterial.fisher.alpha +
                MAP + invasion + PC1.soil,
              # random = ~ 1 | block,
              data = fung.g)
summary(nag.lm)
plot(nag.lm$residuals)
shapiro.test(nag.lm$residuals)

### Phos----
phos.lm <- lme(Phos_umolgSOC.1h.1 ~ fungal.fisher.alpha + bacterial.fisher.alpha +
                 MAP + plant.shannon.residuals + PC1.soil,
    random = ~ 1 | block,
    data = fung.g)
summary(phos.lm)
plot(phos.lm$residuals)
shapiro.test(phos.lm$residuals)

#--------------------------------------------------------------#
# 2. Fungi ----
#--------------------------------------------------------------#

## Species richness and diversity ----

# Check normality
hist(div.soil.g$fungal.spec.richness)
shapiro.test(div.soil.g$fungal.spec.richness)  
hist(div.soil.g$fungal.fisher.alpha)
shapiro.test(div.soil.g$fungal.fisher.alpha)  

### OTU richness: normal ----
### Fine
fung.sr.lm <- lm(fungal.spec.richness ~ plant.shannon.residuals + PC1.soil + MAP,
                 #random = ~ 1 | block,
                 data = fung.g)
summary(fung.sr.lm)
hist(fung.sr.lm$residuals)
shapiro.test(fung.sr.lm$residuals)

ggplot(fung.g, aes(x = PC1.soil,
                   y = fungal.spec.richness)) +
  geom_point(size = 3, shape = 21, fill = '#905971', color = 'black') +
  geom_smooth(method = "lm", se = FALSE, color = 'darkgrey') +
  ylab("Fungal OTU richness") +
  xlab("PC1 soil (43.4%)") +
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

ggsave('figures/Fig2c.tiff', plot = last_plot(),
       device = 'tiff', width = 6, height = 4)

### Fisher's alpha: normal ----
### Fine
fung.div.lm <- lm(fungal.fisher.alpha ~ plant.shannon.residuals + PC1.soil + MAP,
                   #random = ~ 1 | block,
                   data = fung.g)
summary(fung.div.lm)
hist(fung.div.lm$residuals)
shapiro.test(fung.div.lm$residuals)

## NMDS and PERMANOVA: Fungi ----
# isolate fungal community data
fung.comm <- dplyr::select(fung.g, starts_with('Otu'))
# Remove singletons
fung.comm <- fung.comm[colSums(fung.comm) > 1]

### PERMANOVA: Jaccard----
# Constrains permutations to blocks (aka pastures)
perm <- how(nperm = 199)
setBlocks(perm) <- with (fung.g, block)

# fine
fung.jac.adonis <- adonis2(fung.comm ~ PC1.soil * plant.shannon.residuals * MAP + PCNM1,
                          method = 'jaccard',
                          by = 'terms',
                          data = fung.g,
                          permutations = perm)

fung.jac.adonis.df <- as.data.frame(fung.jac.adonis)
write.csv(fung.jac.adonis.df, 'results/Table3_Fungal_JaccardPERMANOVA_Grasslands.csv',
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
                          inv = fung.g$invasion,
                          dist = fung.g$disturbance,
                          MAP = fung.g$MAP,
                          soil = fung.g$PC1.soil,
                          plant.div = fung.g$plant.shannon.residuals,
                          geo.dist = fung.g$PCNM1)
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

ggsave('figures/Fig2e.tiff',
       device = 'tiff',
       plot = jacc.plot,
       width = 5, height = 4, units = 'in')

 ### PERMANOVA: Morisita-Horn ----
# Constrains permutations to blocks (aka pastures)
perm <- how(nperm = 199)
setBlocks(perm) <- with (fung.g, block)

# fine
fung.horn.adonis <- adonis2(fung.comm ~ PC1.soil * plant.shannon.residuals * MAP + PCNM1,
                          method = 'horn',
                          by = 'terms',
                          data = fung.g,
                          permutations = perm)

fung.horn.adonis.df <- as.data.frame(fung.horn.adonis)
write.csv(fung.horn.adonis.df,
          'results/Table3_Fungal_MorisitaHornPERMANOVA_Grasslands.csv',
          row.names = T)

### NMDS: Morisita-Horn ----
# create distance matrix
horn.dist <- vegdist(fung.comm, method = 'horn', binary = F)

# GGplot: Jaccard
horn.mds <- metaMDS(horn.dist, dist = 'bray',
                    try = 1000, trymax = 1000)
horn.stress <- horn.mds$stress

# format data for plot
data.scores <- data.frame(NMDS1 = horn.mds$points[,1],
                          NMDS2 = horn.mds$points[,2],
                          inv = fung.g$invasion,
                          dist = fung.g$disturbance,
                          MAP = fung.g$MAP,
                          soil = fung.g$PC1.soil,
                          plant.div = fung.g$plant.shannon.residuals,
                          geo.dist = fung.g$PCNM1)

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
        legend.position = 'right')

ggsave('figures/SupplementaryFigS6a.jpeg',
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
hist(div.soil.g$bacterial.spec.richness)
shapiro.test(div.soil.g$bacterial.spec.richness)  
hist(div.soil.g$bacterial.fisher.alpha)
shapiro.test(div.soil.g$bacterial.fisher.alpha)  

### Species richness ----
### Fine
bac.sr.lm <- lm(bacterial.spec.richness ~ plant.shannon.residuals + PC1.soil + MAP,
                #random = ~ 1 | block,
                data = bac.g)
summary(bac.sr.lm)
hist(bac.sr.lm$residuals)
shapiro.test(bac.sr.lm$residuals)

### Fisher's alpha ----
### Fine
bac.div.lm <- lm(bacterial.fisher.alpha ~ plant.shannon.residuals + PC1.soil + MAP,
                 #random = ~ 1 | block,
                 data = bac.g)
summary(bac.div.lm)
hist(bac.div.lm$residuals)
shapiro.test(bac.div.lm$residuals)

ggplot(bac.g, aes(x = PC1.soil,
                  y = bacterial.fisher.alpha)) +
  geom_point(aes(fill = MAP), size = 3, shape = 21, color = 'black') +
  geom_smooth(method = "lm", se = FALSE, color = 'darkgrey') +
  scale_fill_gradient(low = "white", high = "#905971") +  # change colors here
  ylab("Bacterial Fisher's alpha") +
  xlab("PC1 soil (43.4%)") +
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
        panel.grid.minor = element_blank()
  )

ggsave('figures/Fig2d.tiff', plot = last_plot(),
       device = 'tiff', width = 6, height = 4)

## NMDS and PERMANOVA: Bacteria ----
# isolate fungal community data
bac.comm <- dplyr::select(bac.g, starts_with('Otu'))
# Remove singletons
bac.comm <- bac.comm[colSums(bac.comm) > 1]

### PERMANOVA: Jaccard----
# Constrains permutations to blocks (aka pastures)
perm <- how(nperm = 199)
setBlocks(perm) <- with (bac.g, block)

# fine
bac.jac.adonis <- adonis2(bac.comm ~ PC1.soil * plant.shannon.residuals * MAP + PCNM1,
                          method = 'jaccard',
                          by = 'terms',
                          data = bac.g,
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
                          inv = bac.g$invasion,
                          dist = bac.g$disturbance,
                          soil = bac.g$PC1.soil,
                          MAP = bac.g$MAP,
                          plant.div = bac.g$plant.shannon,
                          geo.dist = bac.g$PCNM1)

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
  theme(axis.text = element_text(size = 12,
                                 colour = 'black'),
        axis.title = element_text(size = 14),
        legend.position = 'right',
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  )

ggsave('figures/SupplementaryFigS6c.jpeg',
       device = 'jpeg',
       plot = jacc.plot,
       width = 6, height = 5, units = 'in')

### PERMANOVA: Morisita-Horn ----
# Constrains permutations to blocks (aka pastures)
perm <- how(nperm = 199)
setBlocks(perm) <- with(bac.g, block)

# fine
bac.horn.adonis <- adonis2(bac.comm ~ PC1.soil * plant.shannon.residuals * MAP + PCNM1,
                          method = 'horn',
                          by = 'terms',
                          data = bac.g,
                          permutations = perm)

bac.horn.adonis.df <- as.data.frame(bac.horn.adonis)
write.csv(bac.horn.adonis.df,
          'results/Table3_Bacterial_MorisitaHornPERMANOVA_Grasslands.csv',
          row.names = T)

### NMDS: Morisita-Horn ----
# create distance matrix
horn.dist <- vegdist(bac.comm, method = 'horn', binary = F)

# GGplot: Jaccard
horn.mds <- metaMDS(horn.dist, dist = 'bray',
                    try = 1000, trymax = 1000)
horn.stress <- horn.mds$stress

# format data for plot
data.scores <- data.frame(NMDS1 = horn.mds$points[,1],
                          NMDS2 = horn.mds$points[,2],
                          sample = bac.g$sample,
                          inv = bac.g$invasion,
                          dist = bac.g$disturbance,
                          MAP = bac.g$MAP,
                          soil = bac.g$PC1.soil,
                          plant.div = bac.g$plant.shannon,
                          geo.dist = bac.g$PCNM1)

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
  theme(axis.text = element_text(size = 12,
                                 colour = 'black'),
        axis.title = element_text(size = 14),
        legend.position = 'none',
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  )

ggsave('figures/Fig2f.tiff',
       device = 'tiff',
       plot = horn.plot,
       width = 5, height = 4, units = 'in')

#--------------------------------------------------------------#
# 4. Structural equation model ----
#--------------------------------------------------------------#

## Read in data from above ----
sem.data <- read.csv('data/data.SEM.Grassland.csv')

## Define the model ----
sem_model_grassland <- '
  # Soil affected by plants
  PC1.soil ~ a * plant.shannon

  # Microbial diversity
  fungal.spec.richness ~ b * PC1.soil
  bacterial.fisher.alpha ~ c * PC1.soil + d * MAP

  # Community dissimilarity (NMDS1)
  jacc.f.nmds1 ~ e * PC1.soil + f * MAP + g * plant.shannon + h * PCNM1
  horn.b.nmds1 ~ i * PC1.soil + j * MAP + k * plant.shannon + l * PCNM1

  # -----------------------
  # Nuanced Indirect Effects
  # -----------------------

  # Effects from plant.shannon
  indirect_plant_on_fungal_div := a * b     # plant → soil → fungal richness
  indirect_plant_on_bact_div := a * c       # plant → soil → bacterial diversity
  indirect_plant_on_jacc := a * e + g       # plant → soil + direct → Jaccard
  indirect_plant_on_horn := a * i + k       # plant → soil + direct → Horn

  # Effects from MAP
  indirect_MAP_on_bact_div := d               # direct only
  indirect_MAP_on_jacc := f                   # direct only
  indirect_MAP_on_horn := j                   # direct only

  # Spatial effects (if meaningful to track)
  spatial_on_jacc := h
  spatial_on_horn := l
'

## Fit the model ----
fit <- sem(sem_model_grassland, data = sem.data)
summary(fit, fit.measures = TRUE, standardized = TRUE, rsquare = TRUE)

# inspect residuals
residuals(fit, type = "standardized")

# Check modification indices
modindices(fit, sort. = TRUE, minimum.value = 3.84)

## Visualize SEM ----
graph_sem(fit,
          shape = T,
          stand = T)

labels = list(#invasion = "Invasion",
           plant.shannon = "Plant diversity",
           PC1.soil = "Soil PC1",
           MAP = "MAP",
           PCNM1 = "Geographical\ndistance",
           fungal.spec.richness = "Fungal OTU richenss",
           bacterial.fisher.alpha = "Bacterial diversity",
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
           edge_options = list(color = "darkgrey"),
           stars = c('regress'))

