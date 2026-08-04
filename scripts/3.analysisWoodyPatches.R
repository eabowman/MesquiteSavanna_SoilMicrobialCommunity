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

fung.m <- read.csv('data/data.fungal.WoodyPatches.csv')
bac.m <- read.csv('data/data.bacterial.WoodyPatches.csv')
div.soil.m <- read.csv('data/data.diversity.WoodyPatches.csv')

# Change plot from block + treatment to just a plot number
## Diversity data
# Remove the first character (block number)
plot_code <- substring(div.soil.m$plot, 2)

# Convert treatments to numbers
div.soil.m$plot_num <- as.integer(factor(plot_code,
                                       levels = c("DMN", "UMN",
                                                  "DML", "UML")))

## Fungal data
# Remove the first character (block number)
plot_code <- substring(fung.m$plot, 2)

# Convert treatments to numbers
fung.m$plot_num <- as.integer(factor(plot_code,
                                       levels = c("DMN", "UMN",
                                                  "DML", "UML")))

## Bacterial data
# Remove the first character (block number)
plot_code <- substring(bac.m$plot, 2)

# Convert treatments to numbers
bac.m$plot_num <- as.integer(factor(plot_code,
                                       levels = c("DMN", "UMN",
                                                  "DML", "UML")))

# 0. Assess collinearity of data ----

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

# 1. Soil characteristics----

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
  geom_point(size = 3, color = "black") +  # black border
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
  theme(axis.text = element_text(size = 10, color = 'black'),
        axis.title = element_text(size = 12),
        strip.text.x = element_text(size = 10),
        legend.position = 'right',
        legend.title = element_text(size = 10),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) -> Fig.4a

## PC1 as a function of MAP and Plant diversity ----
PC1.lmer <- lmer(PC1.soil ~ plant.shannon + (1 | block/plot),
             data = fung.m)
summary(PC1.lmer)

plot(PC1.lmer)
shapiro.test(residuals(PC1.lmer))

## PC2 as a function of MAP and Plant diversity ----
PC2.lmer <- lmer(PC2.soil ~ plant.shannon + (1 | block/plot),
              data = fung.m)
summary(PC2.lmer)

plot(PC2.lmer)
shapiro.test(residuals(PC2.lmer))

## Enzymatic activity ----

### Beta galactosidase ----
beta.lmer <- lmer(log.beta ~ fungal.fisher.alpha + bacterial.fisher.alpha +
              plant.shannon + PC1.soil + (1 | block/plot),
              data = fung.m)
summary(beta.lmer)

plot(beta.lm$residuals)
shapiro.test(beta.lm$residuals)

ggplot(fung.m, aes(x = bacterial.fisher.alpha,
                   y = log.beta)) +
  geom_point(size = 3, shape = 21, color = 'black', fill = '#f0a354') +
  geom_smooth(method = "lm", se = FALSE, color = 'darkgrey') +
  ylab(expression(paste("Log ",beta, "-glucosidase (", mu, "mol g",
                        SOC^-1, hr^-1, ")"))) +
  labs(x = "Bacterial Fisher's Alpha",
       fill = "MAP", ) +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 10,
                                   color = 'black'),
        axis.text.x = element_text(size = 10,
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

nag.lmer <- lmer(log.nag ~ fungal.fisher.alpha + bacterial.fisher.alpha +
               plant.shannon + PC1.soil + (1 | block/plot),
             data = fung.m.nag)
summary(nag.lmer)

plot(nag.lmer)
shapiro.test(residuals(nag.lmer))

### Phos----
fung.m %>%
  mutate(log.phos = log(Phos_umolgSOC.1h.1)) %>%
  filter(log.phos > 6) -> fung.m.phos

phos.lmer <- lmer(log.phos ~ fungal.fisher.alpha + bacterial.fisher.alpha +
                plant.shannon + PC1.soil + (1 | block/plot),
                data = fung.m.phos)
summary(phos.lmer)

plot(phos.lmer)
shapiro.test(residuals(phos.lmer))

ggplot(fung.m.phos, aes(x = bacterial.fisher.alpha,
                        y = log.phos)) +
  geom_point(aes(fill = PC1.soil), size = 3, shape = 21, color = 'black') +
  geom_smooth(method = "lm", se = FALSE, color = 'darkgrey') +
  scale_fill_gradient(low = "white", high = "#c96a04") +  # change colors here
  ylab(expression(paste("Log Phosphatase (",mu, "mol g",
                        SOC^-1, hr^-1, ")"))) +
  labs(x = "Bacterial Fisher's alpha",
       fill = "Soil PC1") +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 10,
                                   color = 'black'),
        axis.text.x = element_text(size = 10,
                                   color = 'black'),
        legend.position = 'right',
        legend.title = element_text(size = 10),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) -> Fig.4b

# 2. Fungi ----

## Species richness and diversity ----

# Check normality
hist(div.soil.m$fungal.spec.richness)
shapiro.test(div.soil.m$fungal.spec.richness)  
hist(div.soil.m$fungal.fisher.alpha)
shapiro.test(div.soil.m$fungal.fisher.alpha)  

fung.sr.lmer <- lmer(fungal.spec.richness ~ plant.shannon + PC2.soil + (1 | block/plot),
                 data = fung.m)
summary(fung.sr.lmer)

plot(fung.sr.lmer)
shapiro.test(residuals(fung.sr.lmer))

### Fisher's alpha: normal ----
### Fine
fung.div.lmer <- lmer(fungal.fisher.alpha ~ plant.shannon + PC2.soil + (1 | block/plot),
                  data = fung.m)
summary(fung.div.lmer)

plot(fung.div.lmer)
shapiro.test(residuals(fung.div.lmer))

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
fung.jac.adonis <- adonis2(fung.comm ~ PC1.soil * plant.shannon + PCNM1,
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

Fig.4d <- ggplot() + 
  geom_point(data = data.scores, aes(x = NMDS1,
                                     y = NMDS2,
                                     fill = soil,
                                     shape = inv),
             size = 3, color = 'black') +
  scale_fill_gradient(low = "white", high = "#c96a04") +  # change colors here
  # scale_color_viridis() +
  coord_equal() +
  scale_shape_manual(values = c(21, 24)) +
  theme_classic() +
  labs(shape = "",
       fill = "Soil PC1") +  
  guides(shape = 'none') +
  theme(axis.text = element_text(size = 10,
                                 colour = 'black'),
        axis.title = element_text(size = 12),
        legend.position = 'right',
        legend.title = element_text(size = 10),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# Isolate outlier invaded communitys with a high soil.1 value
data.scores %>%
  filter(inv == 'Uninvaded' & NMDS1 > 0.5) -> odd.uninvaded

div.soil.m[div.soil.m$PC1.soil %in% odd.uninvaded$soil,] -> odd.uninvaded.id

### PERMANOVA: Morisita-Horn ----
# Constrains permutations to blocks (aka pastures)
perm <- how(nperm = 199)
setBlocks(perm) <- with (fung.m, block)

# fine
fung.horn.adonis <- adonis2(fung.comm ~ PC1.soil * plant.shannon + PCNM1,
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
data.scores <-  data.frame(NMDS1 = horn.mds$points[,1],
                           NMDS2 = horn.mds$points[,2],
                           inv = fung.m$invasion,
                           dist = fung.m$disturbance,
                           MAP = fung.m$MAP,
                           soil = fung.m$PC1.soil,
                           plant.div = fung.m$plant.shannon,
                           geo.dist = fung.m$PCNM1)

Fig.S6b <- ggplot() + 
  geom_point(data = data.scores, aes(x = NMDS1,
                                     y = NMDS2,
                                     color = soil,
                                     shape = inv),
             size = 3) +
  #scale_fill_manual(values = c('darkblue', 'lightblue')) + # for invasion
  scale_color_viridis() +
  coord_equal() +
  labs(color = "Soil PC1",
       shape = "") +
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = 'black'),
        axis.title = element_text(size = 12),
        legend.position = 'bottom',
        legend.title = element_text(size = 10),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# 3. Bacteria ----

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

bac.sr.lmer <- lmer(sq.root.bacterial.spec.richness ~ plant.shannon + PC2.soil + (1 | block/plot),
                data = div.soil.m)
summary(bac.sr.lmer)

plot(bac.sr.lmer)
shapiro.test(residuals(bac.sr.lmer))

### Fisher's alpha ----
bac.div.lmer <- lmer(bacterial.fisher.alpha ~ plant.spec.richness + PC2.soil + (1 | block/plot),
                 data = div.soil.m)
summary(bac.div.lmer)

plot(bac.div.lmer)
shapiro.test(residuals(bac.div.lmer))

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
  theme(axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 9,
                                   color = 'black'),
        axis.text.x = element_text(size = 10,
                                   color = 'black'),
        legend.position = 'right',
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) -> Fig.4c

## NMDS and PERMANOVA: Bacteria ----
# isolate fungal community data
bac.comm <- dplyr::select(bac.m, starts_with('Otu'))
# Remove singletons
bac.comm <- bac.comm[colSums(bac.comm) > 10]

### PERMANOVA: Jaccard----
# Constrains permutations to blocks (aka pastures)
perm <- how(nperm = 199)
setBlocks(perm) <- with(bac.m, block)

# fine
bac.jac.adonis <- adonis2(bac.comm ~ PC1.soil * plant.shannon + PCNM1,
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

Fig.S6d <- ggplot() + 
  geom_point(data = data.scores, aes(x = NMDS1,
                                     y = NMDS2,
                                     color = soil,
                                     shape = inv),
             size = 3) +
  #scale_fill_manual(values = c('darkblue', 'lightblue')) + # for invasion
  scale_color_viridis() +
  coord_equal() +
  theme_classic() +
  labs(color = 'Soil PC1',
       shape = "") +
  theme(axis.text = element_text(size = 10,
                                 colour = 'black'),
        axis.title = element_text(size = 12),
        legend.position = 'none',
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# Generate Appendix S2: Figure S6 ----
((Fig.S6a | Fig.S6b) /
  (Fig.S6c | Fig.S6d)) +
    plot_layout(
      heights = c(1.1,1),
      # width = c(1,1.2,1,1.2),
      guides = 'collect') + 
  plot_annotation(tag_levels = "a") -> Fig.S6

ggsave("figures/FigureS6.jpg",
       plot = Fig.S6,
       width = 10, height = 8,
       device = 'jpg',
       dpi = 600)


### PERMANOVA: Morisita-Horn ----
# Constrains permutations to blocks (aka pastures)
perm <- how(nperm = 199)
setBlocks(perm) <- with(bac.g, block)

# fine
bac.horn.adonis <- adonis2(bac.comm ~  PC1.soil * plant.shannon + PCNM1,
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

Fig.4e <- ggplot() + 
  geom_point(data = data.scores, aes(x = NMDS1,
                                     y = NMDS2,
                                     fill = soil,
                                     shape = inv),
             size = 3, color = 'black') +
  scale_fill_gradient(low = "white", high = "#c96a04") +  
  coord_equal() +
  scale_shape_manual(values = c(21,24)) +
  theme_classic() +
  labs(color = 'Soil PC1',
       shape = "") +
  theme(axis.text = element_text(size = 10,
                                 colour = 'black'),
        axis.title = element_text(size = 12),
        legend.position = 'none',
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# Generate Figure 4 ----
((Fig.4a) / 
  (Fig.4b | Fig.4c) /
  (Fig.4d | Fig.4e)) +
    plot_layout(heights = c(1.3, 1, 1),
                guides = 'collect') + 
  plot_annotation(tag_levels = "a") -> Fig.4

ggsave("figures/Figure4.jpg",
       plot = Fig.4,
       width = 10, height = 10,
       device = 'jpg',
       dpi = 600)

# 4. Structural equation modeling ----

## Read in data from above ----
sem.data <- read.csv('data/data.SEM.Woodypatch.csv')

# Scale continuous variables prior to analysis
sem.data.scaled <- sem.data

vars <- c("plant.shannon",  
          "PC1.soil",  "PC2.soil",  "log.phos", 
          "bacterial.fisher.alpha",
          "jacc.f.nmds1","horn.b.nmds1",
          "PCNM1")

sem.data.scaled[vars] <-
    lapply(sem.data[vars], scale)

## Define the model ----
sem_model_woody.v2 <- '
# Soil
PC2.soil ~ a * plant.shannon

# Enzyme
log.phos ~ b * PC2.soil

# Diversity
bacterial.fisher.alpha ~ c * PC2.soil +
                         d * log.phos

# Community composition
jacc.f.nmds1 ~ e * PC1.soil +
               f * PC2.soil

horn.b.nmds1 ~ g * PC1.soil +
               h * PCNM1
               

# residual covariance
jacc.f.nmds1 ~~ horn.b.nmds1

# mediation
indirect_phos := b * d
'

## Fit the model ----
fit <- sem(sem_model_woody.v2, data = sem.data.scaled, missing = "ML")
summary(fit, fit.measures = TRUE, standardized = TRUE, rsquare = TRUE)

## Visualize SEM ----
graph_sem(fit,
          shape = T,
          stand = T)

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
  bacterial.fisher.alpha = "Bacterial\ndiversity",
  jacc.f.nmds1 = "Fungal\ncomm. composition",
  horn.b.nmds1 = "Bacterial\ncomm. composition",
  PCNM1 = 'Geo. distance')

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

# remove all objects from the R environment
rm(list = ls())
