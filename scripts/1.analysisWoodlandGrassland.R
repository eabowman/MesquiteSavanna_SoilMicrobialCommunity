## Script written to assess variation in soil traits (nutrients, fungal 
## communities, and bacterial communities) as a function of vegetation type 
## in South Texas savannas. 
## The two vegetation types explored are woodlands dominated by Prosopis
## glandulosa, Honey mesquite, as well as other woody plant species and 
## grasslands dominated by herbaceous plants.
## Written by Dr. Elizabeth Bowman, Dec. 9, 2024
## eabowman@utexas.edu
## 
## Here I am testing the hypothesis that soil traits and microbial communities 
## vary most strongly between vegetation types. 
## The explanatory variables here are the three broad factors: vegetation type,
## invasion, and disturbance.

## 1. Soil characteristics
## 2. Fungi
## 3. Bacteria

# Import data ----
## Import metadata ----
site.data <- read.csv('data/field.data_climate.csv')

## Import soil data ----
soil.data <- read.csv('data/soil.data.clean.2023.csv')

# change 0 values to 0.1 for transformation
soil.data[soil.data$Nag_umolgSOC.1h.1 == 0, 'Nag_umolgSOC.1h.1'] <- 0.000000001
soil.data[soil.data$Phos_umolgSOC.1h.1 == 0, 'Phos_umolgSOC.1h.1'] <- 0.000000001

# Transform soil data
soil.data %>%
  mutate(log.beta = log(Beta_umolgSOC.1h.1),
         cube.root.nag = (Nag_umolgSOC.1h.1)^(1/3),
         cube.root.phos = (Phos_umolgSOC.1h.1)^(1/3),
         log.C = log(C.mgperkg),
         log.N = log(N.mgperkg),
         log.P = log(P.mgperkg),
         log.Mg = log(Mg.mgperkg),
         log.K = log(K.mgperkg),
         log.Ca = log(Ca.mgperkg),
         log.Fe = log(Fe.mgperkg)) -> soil.data

# make disturbance, invasion, and vegetation type columns all lower case 
tolower(soil.data$invasion) -> soil.data$invasion
tolower(soil.data$disturbance) -> soil.data$disturbance
tolower(soil.data$vegetation.type) -> soil.data$vegetation.type

## Import plant community data ----
plant.data <- read.csv('data/plant.comm.sxs.csv')

# Isolate community data
plant.comm <- plant.data[8:length(plant.data)]
rownames(plant.comm) <- plant.data$sample

# Remove Thatch and Bare from community data
plant.comm[!colnames(plant.comm) %in% c('Thatch', 'Bare')] -> plant.comm

# Calculate species richness and diversity
plant.data$spec.richness <- specnumber(plant.comm)
plant.data$shannon.div <- diversity(plant.comm, index = 'shannon')
plant.data$invsimpson.div <- diversity(plant.comm, index = 'invsimpson')

## Import diversity data ----
div.data <- read.csv('data/data.diversity.Overall.csv')

## Import fungal data ----
fung.data <- read.csv('data/data.fungal.Overall.csv')

## Import bacterial data ----
bac.data <- read.csv('data/data.bacterial.Overall.csv')

#--------------------------------------------------------------#
# 1. Soil traits ----
#--------------------------------------------------------------#

# make block a categorical variable
soil.data$block <- factor(soil.data$block)

# factor invasion, disturbance, and plant community
soil.data$invasion <- factor(soil.data$invasion,
                             levels = c('invaded','uninvaded'),
                             labels = c('Invaded', 'Uninvaded'))
soil.data$disturbance <- factor(soil.data$disturbance,
                                levels = c('disturbed', 'undisturbed'),
                                labels = c('Disturbed', 'Undisturbed'))
soil.data$vegetation.type <- factor(soil.data$vegetation.type,
                                    levels = c('grassland', 'motte'),
                                    labels = c('Grassland', 'Woody patch'))

## Enzymatic activity ----

### Beta glucosidase ----
soil.data %>%
  filter(Beta_umolgSOC.1h.1 < 39900) -> soil.beta

# assess normality
shapiro.test(soil.beta$log.beta)

# ANOVA
beta.lm <- lm(log.beta ~ vegetation.type * invasion * disturbance,
                data = soil.beta)
beta.anova <- anova(beta.lm)

# plot
beta.plot <- ggplot(soil.beta, aes(x = vegetation.type,
                                   y = Beta_umolgSOC.1h.1,
                                   fill = vegetation.type)) +
  geom_boxplot(width = 0.25) +
  geom_point(size = 1, color = 'grey') +
  # facet_grid(. ~ Invasion) +
  xlab('') +
  ylab(expression(paste(beta, "-glucosidase (", mu, "mol g",
                        SOC^-1, hr^-1, ")"))) +
  scale_fill_manual(values = c('#905971', '#F1BB83')) +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 8,
                                 color = 'black'),
        axis.text.x = element_text(size = 12,
                                 color = 'black',
                                 angle = 45,
                                 hjust = 1),
        legend.position = 'none',
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  )

### Nag ----
mutate(soil.data, cube.root.nag = (Nag_umolgSOC.1h.1)^(1/3)) -> nag.soil

# assess normality
shapiro.test(nag.soil$cube.root.nag)
hist(nag.soil$cube.root.nag)

# ANOVA
nag.lm <- lm(cube.root.nag ~ vegetation.type * invasion * disturbance,
                data = nag.soil)
nag.anova <- anova(nag.lm)

# plot
nag.plot <- ggplot(nag.soil, aes(x = vegetation.type,
                                   y = Nag_umolgSOC.1h.1,
                                 fill = vegetation.type)) +
  geom_boxplot(width = 0.25) +
  geom_point(size = 1, color = 'grey') +
  # ylim(0, 1500) +
  xlab('') +
  ylab(expression(paste("NAGase (",mu, "mol g",
                        SOC^-1, hr^-1, ")"))) +
  scale_fill_manual(values = c('#905971', '#F1BB83')) +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 8,
                                 color = 'black'),
        axis.text.x = element_text(size = 12,
                                 color = 'black',
                                 angle = 45,
                                 hjust = 1),
        legend.position = 'none',
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  )

### Phosphatase ----
mutate(soil.data,
       cube.root.phos = (Phos_umolgSOC.1h.1)^(1/3)) -> phos.soil

# assess normality
shapiro.test(phos.soil$cube.root.phos)
hist(phos.soil$cube.root.phos)

phos.lm <- lm(cube.root.phos ~ vegetation.type * invasion * disturbance,
                data = phos.soil)
phos.anova <- anova(phos.lm)

plot(phos.lm)

# plot
phos.plot <- ggplot(phos.soil, aes(x = vegetation.type,
                                   y = Phos_umolgSOC.1h.1,
                                   fill = vegetation.type)) +
  geom_boxplot(width = 0.25) +
  geom_point(size = 1, color = 'grey') +
  xlab('') +
  ylab(expression(paste("Phosphatase (",mu, "mol g",
                        SOC^-1, hr^-1, ")"))) +
  scale_fill_manual(values = c('#905971', '#F1BB83')) +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 8,
                                 color = 'black'),
        axis.text.x = element_text(size = 12,
                                 color = 'black',
                                 angle = 45,
                                 hjust = 1),
        legend.position = 'none',
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  )

veg.plot <- gridExtra::arrangeGrob(beta.plot, nag.plot, phos.plot,
                                   nrow = 1, ncol = 3)
ggsave('figures/Fig1btod.tiff', plot = veg.plot,
       device = 'tiff', width = 6, height = 5, units = 'in')

## Nutrients and pH -----

### PCA of soil traits ----
# isolate soil data traits of interest
soil.data %>%
  dplyr::select(Active_C_mgkg.1, C.mgperkg, N.mgperkg, P.mgperkg, Mg.mgperkg,
         Mn.mgperkg, Fe.mgperkg, K.mgperkg, Ca.mgperkg, pH) -> soil.pca.data

# Assess data
soil.pca.data %>%
  pivot_longer(cols = everything(), names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value)) +
  geom_histogram(binwidth = 1, fill = "blue", color = "black", alpha = 0.7) +
  facet_wrap(~ variable, scales = "free") +
  theme_minimal() +
  labs(x = "Value", y = "Frequency", title = "Histograms for Variables")

# Standardize data using Z-scores (default method is "standardize")
standardized_data <- decostand(soil.pca.data, method = "standardize")

# PCA
pca.soil <- prcomp(standardized_data,
                   center = F, 
                   scale. = F)

summary(pca.soil)

# Extract principal components 
pca.soil.plot <- as.data.frame(pca.soil$x)

# Add treatments
cbind(sample = soil.data$sample,
      pca.soil.plot, disturbance = soil.data$disturbance,
      invasion = soil.data$invasion,
      vegetation.type = soil.data$vegetation.type) -> pca.soil.plot

# assess how soil traits are influenced by veg. type, inv., and dist.
# PC1
lm.pc1 <- lm(PC1 ~ vegetation.type * invasion * disturbance, data = pca.soil.plot)
anova(lm.pc1)

# PC2
lm.pc2 <- lm(PC2 ~ vegetation.type * invasion * disturbance, data = pca.soil.plot)
anova(lm.pc2)

# Extract loadings for PC1 and PC2
loadings <- as.data.frame(pca.soil$rotation[, 1:2])
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

#### Plot: Vegetation type ----
ggplot(pca.soil.plot,
       aes(x = PC1, y = PC2, fill = vegetation.type)) +
  geom_point(shape = 21, size = 3, color = 'black') +
  stat_ellipse(data = pca.soil.plot,
               aes(x = PC1, y = PC2, fill = vegetation.type),
               geom = "polygon",
               type = "norm",
               level = 0.95,
               alpha = 0.3) +
  xlab("PC1 (55.9%)") +
  ylab("PC2 (13.6%)") +
  scale_fill_manual(values = c("#905971", "#F1BB83")) +
  # Add arrows for soil traits
  geom_segment(data = loadings,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "black", inherit.aes = F) +
  # Add labels for soil traits
  # geom_text(data = loadings,
  #           aes(x = PC1, y = PC2, label = trait),
  #           color = "black", size = 5,
  #           vjust = -0.5, hjust = 0.5,
  #           inherit.aes = F) +
  ggrepel::geom_text_repel(data = loadings,
            aes(x = PC1, y = PC2, label = trait),
            color = "black", size = 5,
            max.overlaps = 20,
            box.padding = 0.5, 
            point.padding = 0.2,
            segment.color = "darkred",
            inherit.aes = F) +
  theme_classic() +
  theme(axis.text = element_text(size = 18, color = 'black'),
        axis.title = element_text(size = 20),
        strip.text.x = element_text(size = 18),
        legend.position = 'none',
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  )

ggsave('figures/Fig1a.tiff', plot = last_plot(),
       device = 'tiff', width = 6, height = 5)

#### Plot: Invasion ----
ggplot(pca.soil.plot,
       aes(x = PC1, y = PC2,
           fill = invasion)) +
  stat_ellipse(data = pca.soil.plot,
               aes(x = PC1, y = PC2,
                   fill = invasion),
               geom = "polygon",
               type = "norm",
               level = 0.95,
               alpha = 0.3) +
  geom_point(shape = 21, size = 3, color = 'black') +
  # facet_grid(. ~ vegetation.type) +
  scale_fill_manual(values = c("#4E84C4", "#293352")) +
  # Add arrows for soil traits
  geom_segment(data = loadings,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "black", inherit.aes = F) +
  # Add labels for soil traits
  # geom_text(data = loadings,
  #           aes(x = PC1, y = PC2, label = trait),
  #           color = "black", size = 5,
  #           vjust = -0.5, inherit.aes = F) +
  ggrepel::geom_text_repel(data = loadings,
                           aes(x = PC1, y = PC2, label = trait),
                           color = "black", size = 5,
                           max.overlaps = 20,
                           box.padding = 0.5, 
                           point.padding = 0.2,
                           segment.color = "darkred",
                           inherit.aes = F) +
  theme_classic() +
  xlab("PC1 (55.9%)") +
  ylab("PC2 (13.6%)") +
  theme(axis.text = element_text(size = 18, color = 'black'),
        axis.title = element_text(size = 20),
        strip.text.x = element_text(size = 18),
        strip.background = element_rect(color = 'white',
                                        fill = 'white'),
        legend.position = 'none')

ggsave('figures/SupplementaryFigS4a.tiff', plot = last_plot(),
       device = 'tiff', width = 6, height = 5)

#### Plot: Disturbance ----
ggplot(pca.soil.plot,
       aes(x = PC1, y = PC2,
           fill = disturbance)) +
  stat_ellipse(data = pca.soil.plot,
               aes(x = PC1, y = PC2,
                   fill = disturbance),
               geom = "polygon",
               type = "norm",
               level = 0.95,
               alpha = 0.2) +
  geom_point(shape = 21, size = 3, color = 'black') +
  # facet_grid(. ~ vegetation.type) +
  scale_fill_manual(values = c("#D1D0DE","#636D97")) +
  # Add arrows for soil traits
  geom_segment(data = loadings,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "black", inherit.aes = F) +
  # Add labels for soil traits
  # geom_text(data = loadings,
  #           aes(x = PC1, y = PC2, label = trait),
  #           color = "black", size = 5,
  #           vjust = -0.5, inherit.aes = F) +
  ggrepel::geom_text_repel(data = loadings,
                           aes(x = PC1, y = PC2, label = trait),
                           color = "black", size = 5,
                           max.overlaps = 20,
                           box.padding = 0.5, 
                           point.padding = 0.2,
                           segment.color = "darkred",
                           inherit.aes = F) +
  theme_classic() +
  xlab("PC1 (55.9%)") +
  ylab("PC2 (13.6%)") +
  theme(axis.text = element_text(size = 18, color = 'black'),
        axis.title = element_text(size = 20),
        strip.text.x = element_text(size = 18),
        strip.background = element_rect(color = 'white',
                                        fill = 'white'),
        legend.position = 'none')

ggsave('figures/SupplementaryFigS4b.tiff', plot = last_plot(),
       device = 'tiff', width = 6, height = 5)

#### Plot: boxplot ----
pca.soil.plot %>%
  dplyr::select(PC1, PC2, disturbance, invasion, vegetation.type) -> pca.soil.plot

# Create the box plot
veg.plot <- ggplot(pca.soil.plot, aes(x = vegetation.type,
                          y = PC1,
                          fill = vegetation.type)) +
  geom_boxplot(width = 0.5) +
  theme_classic() +
  scale_fill_manual(values = c("#905971", "#F1BB83")) +
  # facet_grid(invasion ~ disturbance) +
  xlab("") +
  ylab("PC1 scores") +
  theme(axis.text.y = element_text(size = 18, color = 'black'),
        axis.text.x = element_text(size = 18, color = 'black',
                                   angle = 45, hjust = 1),
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 18),
        strip.background = element_rect(color = 'white',
                                        fill = 'white'),
        legend.position = 'none')

inv.plot <- ggplot(pca.soil.plot, aes(x = invasion,
                                      y = PC1,
                                      fill = invasion)) +
  geom_boxplot(width = 0.5) +
  theme_classic() +
  scale_fill_manual(values = c("#4E84C4", "#293352")) +
  # facet_grid(invasion ~ disturbance) +
  xlab("") +
  ylab("PC1 scores") +
  theme(axis.text.y = element_text(size = 18, color = 'black'),
        axis.text.x = element_text(size = 18, color = 'black',
                                   angle = 45, hjust = 1),
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 18),
        strip.background = element_rect(color = 'white',
                                        fill = 'white'),
        legend.position = 'none')

dist.plot <- ggplot(pca.soil.plot, aes(x = disturbance,
                                      y = PC1,
                                      fill = disturbance)) +
  geom_boxplot(width = 0.5) +
  theme_classic() +
  scale_fill_manual(values = c("#D1D0DE","#636D97")) +
  # facet_grid(invasion ~ disturbance) +
  xlab("") +
  ylab("PC1 scores") +
  theme(axis.text.y = element_text(size = 18, color = 'black'),
        axis.text.x = element_text(size = 18, color = 'black',
                                   angle = 45, hjust = 1),
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 18),
        strip.background = element_rect(color = 'white',
                                        fill = 'white'),
        legend.position = 'none')

box.plot <- gridExtra::arrangeGrob(veg.plot, inv.plot, dist.plot,
                                   nrow = 1, ncol = 3)
ggsave('figures/SupplementaryFigS4c.tiff', plot = box.plot,
       device = 'tiff', width = 6, height = 5, units = 'in')

#--------------------------------------------------------------#
# 2. Fungi ----
#--------------------------------------------------------------#

## Species richness and diversity ----

### Summary data ----
div.data %>%
  group_by(.) %>%
  summarise(fung.sr.mean = mean(fungal.spec.richness),
            fung.sr.sd = sd(fungal.spec.richness),
            fung.div.mean = mean(fungal.fisher.alpha),
            fung.div.sd = sd(fungal.fisher.alpha)) -> overall.fung

div.data %>%
  group_by(vegetation.type) %>%
  summarise(fung.sr.mean = mean(fungal.spec.richness),
            fung.sr.sd = sd(fungal.spec.richness),
            fung.div.mean = mean(fungal.fisher.alpha),
            fung.div.sd = sd(fungal.fisher.alpha)) -> veg.fung

div.data %>%
  group_by(invasion) %>%
  summarise(fung.sr.mean = mean(fungal.spec.richness),
            fung.sr.sd = sd(fungal.spec.richness),
            fung.div.mean = mean(fungal.fisher.alpha),
            fung.div.sd = sd(fungal.fisher.alpha)) -> inv.fung

div.data %>%
  group_by(disturbance) %>%
  summarise(fung.sr.mean = mean(fungal.spec.richness),
            fung.sr.sd = sd(fungal.spec.richness),
            fung.div.mean = mean(fungal.fisher.alpha),
            fung.div.sd = sd(fungal.fisher.alpha)) -> dist.fung

### Species richness: normal ----
# test normality
hist(div.data$fungal.spec.richness)
shapiro.test(div.data$fungal.spec.richness)

fung.sr.lme <- lm(fungal.spec.richness ~ vegetation.type * invasion * disturbance,
                  # random = ~ 1 | block,
                 data = div.data)
anova(fung.sr.lme)
hist(fung.sr.lme$residuals)

ggplot(div.data,
       aes(x = vegetation.type,
           y = fungal.spec.richness,
           fill = vegetation.type)) +
  geom_boxplot(width = 0.25) +
  facet_grid(. ~ invasion) +
  scale_fill_manual(values = c("#905971", "#F1BB83")) +
  xlab("") +
  #facet_grid(. ~ disturbance) +
  ylab("Fungal OTU richness") +
  theme_classic() +
  theme(axis.text.y = element_text(size = 12, color = 'black'),
        axis.text.x = element_text(size = 12, color = 'black',
                                   angle = 45, hjust = 1),
        axis.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        strip.background = element_blank(),
        legend.position = 'none')

ggsave('figures/SupplementaryFigS5a_SR.tiff',
       device = 'tiff',
       plot = last_plot(),
       width = 4, height = 5, units = 'in')

### Fisher's alpha: normal ----
# test normality
hist(div.data$fungal.fisher.alpha)
shapiro.test(div.data$fungal.fisher.alpha)

fung.div.lme <- lm(fungal.fisher.alpha ~ vegetation.type * invasion * disturbance,
                  # random = ~ 1 | block,
                 data = div.data)
anova(fung.div.lme)
hist(fung.div.lme$residuals)

ggplot(div.data,
       aes(x = vegetation.type,
           y = fungal.fisher.alpha,
           fill = vegetation.type)) +
  geom_boxplot(width = 0.25) +
  facet_grid(. ~ invasion) +
  scale_fill_manual(values = c("#905971", "#F1BB83")) +
  xlab("") +
  #facet_grid(. ~ disturbance) +
  ylab("Fungal Fisher's alpha") +
  theme_classic() +
  theme(axis.text.y = element_text(size = 12, color = 'black'),
        axis.text.x = element_text(size = 12, color = 'black',
                                   angle = 45, hjust = 1),
        axis.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        strip.background = element_blank(),
        legend.position = 'none')

ggsave('figures/SupplementaryFigS5c_div.tiff',
       device = 'tiff',
       plot = last_plot(),
       width = 3, height = 5, units = 'in')

## Species richness and diversity as a function of soil traits ----
# How to do this as all soil traits are correlated?
# Should I use a PCA of a few or all traits

div.data %>%
  # dplyr::select(sample, Beta_umolgSOC.1h.1, Phos_umolgSOC.1h.1, Nag_umolgSOC.1h.1) %>%
  # full_join(div.data) %>%
  mutate(log.beta = log(Beta_umolgSOC.1h.1),
         cube.root.phos = (Phos_umolgSOC.1h.1)^(1/3),
         cube.root.nag = (Nag_umolgSOC.1h.1)^(1/3)) -> div.enzyme

# Species richness
enzyme.lm <- lm(fungal.spec.richness ~ log.beta * vegetation.type, 
                 # random = ~ 1 | block,
                 data = div.enzyme)

summary(enzyme.lm)
plot(enzyme.lm)
hist(enzyme.lm$residuals)

shapiro.test(enzyme.lm$residuals)

# Fisher's alpha
# Check this out; originally lme and perumtations randomixed by block as abov.
enzyme.lm <- lm(fungal.fisher.alpha ~ cube.root.nag,
                data = div.enzyme)

summary(enzyme.lm)
plot(enzyme.lm)

## NMDS and PERMANOVA: Fungi ----
# isolate fungal community data
fung.comm <- dplyr::select(fung.data, starts_with('Otu'))
# Remove singletons
fung.comm <- fung.comm[colSums(fung.comm) > 1]

### ANOSIM: Jaccard ----
# jaccard distance matrix
jac.dist <- vegdist(fung.comm,
                     method = "jaccard",
                     binary = T)

# check homoegeniety of variance
jac.disper <- betadisper(jac.dist, group = fung.data$vegetation.type)
anova(jac.disper) # normal

# ANOSIM
anosim(fung.comm, grouping = fung.data$vegetation.type,
       distance = 'jaccard', strata = fung.data$block)

### NMDS: Jaccard ----
# create distance matrix
jacc.dist <- vegdist(fung.comm, method = 'jaccard', binary = T)

# GGplot: Jaccard
jacc.mds <- metaMDS(jacc.dist, dist = 'bray',
                    try = 1000, trymax = 1000)
jacc.stress <- jacc.mds$stress
# rownames(jacc.mds$points) <- fung.data$sample

# format data for plot

data.scores <- data.frame(NMDS1 = jacc.mds$points[,1],
                          NMDS2 = jacc.mds$points[,2],
                          invasion = fung.data$invasion,
                          disturbance = fung.data$disturbance,
                          vegetation.type = fung.data$vegetation.type)

jacc.plot <- ggplot() + 
  geom_point(data = data.scores, aes(x = NMDS1,
                                     y = NMDS2,
                                     fill = vegetation.type),
             size = 3,
             shape = 21) +
  scale_fill_manual(values = c("#905971", "#F1BB83")) +
  coord_equal() +
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = 'black'),
        axis.title = element_text(size = 12),
        legend.position = 'none',
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  )

ggsave('figures/Fig1e.tiff',
       device = 'tiff',
       plot = jacc.plot,
       width = 4, height = 3, units = 'in')

### ANOSIM: Morisita-Horn, plant community structure ---
# Morisita-Horn distance matrix
horn.dist <- vegdist(fung.comm,
                     method = "horn",
                     binary = F)

# check homoegeniety of variance
horn.disper <- betadisper(horn.dist, group = fung.data$vegetation.type)
anova(horn.disper) # normal

### ANOSIM: Morisita-Horn ----
anosim(fung.comm, grouping = fung.data$vegetation.type,
       distance = 'horn', strata = fung.data$block)

### NMDS: Morisita-Horn ----
# create distance matrix
horn.dist <- vegdist(fung.comm, method = 'horn', binary = F)

# GGplot: Jaccard
horn.mds <- metaMDS(horn.dist, dist = 'bray',
                    try = 1000, trymax = 1000)
horn.stress <- horn.mds$stress
# rownames(horn.mds$points) <- fung.data$sample

# format data for plot
data.scores <- data.frame(NMDS1 = horn.mds$points[,1],
                          NMDS2 = horn.mds$points[,2],
                          invasion = fung.data$invasion,
                          disturbance = fung.data$disturbance,
                          vegetation.type = fung.data$vegetation.type)

horn.plot <- ggplot() + 
  geom_point(data = data.scores, aes(x = NMDS1,
                                     y = NMDS2,
                                     fill = vegetation.type),
             size = 4,
             shape = 21) +
  scale_fill_manual(values = c("#905971", "#F1BB83")) +
  coord_equal() +
  theme_classic() +
  theme(axis.text = element_text(size = 14,
                                 colour = 'Black'),
        axis.title = element_text(size = 16),
        legend.position = 'none')

ggsave('figures/SupplementaryFigS5d.tiff',
       device = 'tiff',
       plot = horn.plot,
       width = 5, height = 3, units = 'in')


#--------------------------------------------------------------#
# 3. Bacteria ----
#--------------------------------------------------------------#

## Species richness and diversity ----
### Summary data ----
div.data %>%
  group_by(.) %>%
  summarise(bac.sr.mean = mean(bacterial.spec.richness),
            bac.sr.sd = sd(bacterial.spec.richness),
            bac.div.mean = mean(bacterial.fisher.alpha),
            bac.div.sd = sd(bacterial.fisher.alpha)) -> overall.bac

div.data %>%
  group_by(vegetation.type) %>%
  summarise(bac.sr.mean = mean(bacterial.spec.richness),
            bac.sr.sd = sd(bacterial.spec.richness),
            bac.div.mean = mean(bacterial.fisher.alpha),
            bac.div.sd = sd(bacterial.fisher.alpha)) -> veg.bac

div.data %>%
  group_by(invasion) %>%
  summarise(bac.sr.mean = mean(bacterial.spec.richness),
            bac.sr.sd = sd(bacterial.spec.richness),
            bac.div.mean = mean(bacterial.fisher.alpha),
            bac.div.sd = sd(bacterial.fisher.alpha)) -> inv.bac

div.data %>%
  group_by(disturbance) %>%
  summarise(bac.sr.mean = mean(bacterial.spec.richness),
            bac.sr.sd = sd(bacterial.spec.richness),
            bac.div.mean = mean(bacterial.fisher.alpha),
            bac.div.sd = sd(bacterial.fisher.alpha)) -> dist.bac

filter(div.data, bacterial.spec.richness > 500) -> bac.sr.data

hist(bac.sr.data$bacterial.spec.richness)
shapiro.test(bac.sr.data$bacterial.spec.richness)
  
### Species richness: normal ----
bac.sr.lme <- lm(bacterial.spec.richness ~ vegetation.type * invasion * disturbance,
                  # random = ~ 1 | block,
                 data = bac.sr.data)
anova(bac.sr.lme)
 
### Fisher's alpha ----
hist(bac.sr.data$bacterial.fisher.alpha)
shapiro.test(bac.sr.data$bacterial.fisher.alpha)

bac.div.lme <- lm(bacterial.fisher.alpha ~ vegetation.type * invasion * disturbance,
                  # random = ~ 1 | block,
                 data = bac.sr.data)
anova(bac.div.lme)

# check residuals
histogram(residuals(bac.div.lme))

# plot
bac.div.plot <- ggplot(div.data,
                       aes(x = vegetation.type,
                           y = bacterial.fisher.alpha,
                           fill = vegetation.type)) +
  geom_boxplot(width = 0.25) +
  theme_classic() +
  scale_fill_manual(values = c("#905971", "#F1BB83")) +
  xlab("") +
  #facet_grid(. ~ disturbance) +
  ylab("Bacterial Fisher's alpha") +
  theme(axis.text.y = element_text(size = 12, color = 'black'),
        axis.text.x = element_text(size = 12, color = 'black',
                                   angle = 45, hjust = 1),
        axis.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        strip.background = element_blank(),
        legend.position = 'none')

ggsave('figures/SupplementaryFigS5b.tiff',
       device = 'tiff',
       plot = bac.div.plot,
       width = 2, height = 5, units = 'in')

bac.div.plot <- ggplot(bac.sr.data,
                       aes(x = invasion,
                           y = bacterial.fisher.alpha,
                           fill = invasion)) +
  geom_boxplot(width = 0.25) +
  theme_classic() +
  scale_fill_manual(values = c("#8DB6AB", "#EDE6DE")) +
  xlab("") +
  facet_grid(. ~ disturbance) +
  ylab("Bacterial Fisher's alpha") +
  theme(axis.text.y = element_text(size = 12, color = 'black'),
        axis.text.x = element_text(size = 12, color = 'black',
                                   angle = 45, hjust = 1),
        axis.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        strip.background = element_blank(),
        legend.position = 'none')

ggsave('figures/SupplementaryFigS5c.tiff',
       device = 'tiff',
       plot = bac.div.plot,
       width = 3, height = 5, units = 'in')

## Species richness and diversity as a function of soil traits ----
# Species richness
enzyme.lm <- lme(bacterial.spec.richness ~ log.beta * vegetation.type, 
                 random = ~ 1 | block,
                 data = div.enzyme)

summary(enzyme.lm)
plot(enzyme.lm)

# Fisher's alpha
enzyme.lm <- lme(bacterial.fisher.alpha ~ log.beta * vegetation.type,
                 random = ~ 1 | block,
                 data = div.enzyme)

summary(enzyme.lm)
plot(enzyme.lm)

## NMDS and PERMANOVA: Bacteria ----
# isolate bacterial community data
bac.comm <- dplyr::select(bac.data, starts_with('Otu'))
# Remove singletons
bac.comm <- bac.comm[colSums(bac.comm) > 10]

# check homoegeniety of variance
# jaccard distance matrix
jac.dist <- vegdist(bac.comm,
                    method = "jaccard",
                    binary = T)

jac.disper <- betadisper(jac.dist, group = bac.data$vegetation.type)
anova(jac.disper) # not normal

### PERMANOVA: Jaccard ----
# Constrains permutations to blocks (aka pastures)
perm <- how(nperm = 199)
setBlocks(perm) <- with (bac.data, block)

bac.jac.adonis <- adonis2(bac.comm ~ vegetation.type,
                          method = 'jaccard',
                          data = bac.data,
                          permutations = perm)
bac.jac.adonis

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
                          invasion = bac.data$invasion,
                          disturbance = bac.data$disturbance,
                          vegetation.type = bac.data$vegetation.type)

jacc.plot <- ggplot() + 
  geom_point(data = data.scores, aes(x = NMDS1,
                                     y = NMDS2,
                                     fill = vegetation.type),
             size = 3,
             shape = 21) +
  scale_fill_manual(values = c("#905971", "#F1BB83")) +
  coord_equal() +
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = 'black'),
        axis.title = element_text(size = 12),
        legend.position = 'none',
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  )

ggsave('figures/SupplementaryFigS5e.tiff',
       device = 'tiff',
       plot = jacc.plot,
       width = 4, height = 3, units = 'in')

### PERMANOVA: Morisita-Horn  ----

# jaccard distance matrix
horn.dist <- vegdist(bac.comm,
                     method = "horn",
                     binary = T)

# check homoegeniety of variance
horn.disper <- betadisper(horn.dist,
                          group = bac.data$vegetation.type)
anova(horn.disper) # not normal

# Constrains permutations to blocks (aka pastures)
perm <- how(nperm = 199)
setBlocks(perm) <- with (bac.data, block)

bac.horn.adonis <- adonis2(bac.comm ~ vegetation.type,
                           method = 'horn',
                           data = bac.data,
                           permutations = perm)

bac.horn.adonis

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
                          invasion = bac.data$invasion,
                          disturbance = bac.data$disturbance,
                          vegetation.type = bac.data$vegetation.type)

horn.plot <- ggplot() + 
  geom_point(data = data.scores, aes(x = NMDS1,
                                     y = NMDS2,
                                     fill = vegetation.type),
             size = 3,
             shape = 21) +
  scale_fill_manual(values = c("#905971", "#F1BB83")) +
  coord_equal() +
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = 'Black'),
        axis.title = element_text(size = 12),
        legend.position = 'none',
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  )

ggsave('figures/Fig1f.tiff',
       device = 'tiff',
       plot = horn.plot,
       width = 4, height = 3, units = 'in')

#--------------------------------------------------------------#
# 4. Comparison of two communities ----
#--------------------------------------------------------------#

# Distance matrices
jac.f.dist <- vegdist(fung.comm,
                    method = "jaccard",
                    binary = T)

horn.f.dist <- vegdist(fung.comm,
                      method = "horn",
                      binary = F)

jac.b.dist <- vegdist(bac.comm,
                      method = "jaccard",
                      binary = T)

horn.b.dist <- vegdist(bac.comm,
                       method = "horn",
                       binary = F)

## Combine all distance matrices into a single dataframe ----

as_long_df <- function(d, sample_names, dist_name = "Distance") {
  # Assign sample names to the dist object
  attr(d, "Labels") <- sample_names
  
  # Convert to full matrix
  mat <- as.matrix(d)
  
  # Get lower triangle only
  inds <- lower.tri(mat)
  
  # Build data frame
  data.frame(
    Sample1 = rownames(mat)[row(mat)[inds]],
    Sample2 = colnames(mat)[col(mat)[inds]],
    Distance = mat[inds]
  ) |> 
    dplyr::rename(!!dist_name := Distance)
}

jac.f <- as_long_df(jac.f.dist,
                    fung.data$sample,
                    "Fungi.Jaccard")
horn.f <- as_long_df(horn.f.dist,
                    fung.data$sample,
                    "Fungi.MorisitaHorn")
jac.b <- as_long_df(jac.b.dist,
                    bac.data$sample,
                    "Bacteria.Jaccard")
horn.b <- as_long_df(horn.b.dist,
                     bac.data$sample,
                     "Bacteria.MorisitaHorn")

# Merge all dataframes
merged.dist <- jac.f %>%
  inner_join(horn.f, by = c("Sample1", "Sample2")) %>%
  inner_join(jac.b, by = c("Sample1", "Sample2")) %>%
  inner_join(horn.b, by = c("Sample1", "Sample2"))

# Add metadata
veg.type.df <- fung.data %>%
  dplyr::select(sample, vegetation.type)
merged.dist %>%
  left_join(veg.type.df, by = c("Sample1" = "sample")) %>%
  rename(Veg1 = vegetation.type) %>%
  left_join(veg.type.df, by = c("Sample2" = "sample")) %>%
  rename(Veg2 = vegetation.type) %>%
  mutate(Group = if_else(Veg1 == Veg2, "Within", "Between")) -> merged.dist

## Plot ----
merged.dist$Group <- factor(merged.dist$Group,
                            levels = c('Between', 'Within'))

merged.dist %>%
  gather(key = "Distance.type",
         value = "Distance",
         -Sample1, -Sample2, -Veg1, -Veg2, -Group) %>%
  separate(Distance.type, into = c("Microbial.group", "Index"),
           sep = "\\.") %>%
  mutate(Microbial.group = factor(Microbial.group,
                                  levels = c("Fungi", "Bacteria"))) %>%
  ggplot(aes(x = Group, y = Distance)) +
  geom_boxplot(width = 0.25) +
  ggpubr::stat_compare_means(method = "wilcox.test",
                             label = "p.signif",
                             label.y = 1.0) +
  xlab('') +
  ylab('Dissimilarity') +
  facet_grid(. ~ Microbial.group + Index) +
  theme_classic() +
  theme(axis.text = element_text(size = 12, colour = 'black'),
        axis.title = element_text(size = 13),
        legend.position = 'none',
        strip.background = element_blank(),
        strip.text = element_text(size = 12, color = 'black'),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_cartesian(clip = "off")

ggsave('figures/Fig1g.tiff',
       device = 'tiff',
       plot = last_plot(),
       width = 7, height = 3, units = 'in')

## Wilcox-test ----
# Make dataframe
merged.dist %>%
  gather(key = "Distance.type",
         value = "Distance",
         -Sample1, -Sample2, -Veg1, -Veg2, -Group) %>%
  separate(Distance.type, into = c("Microbial.group", "Index"),
           sep = "\\.") -> merged.dist.long
# Separate out each group
# Jaccard Bacteria
merged.dist.long %>%
  filter(Index == 'Jaccard',
         Microbial.group == 'Bacteria') -> merged.dist.long.bac.jac
wilcox.test(Distance ~ Group, data = merged.dist.long.bac.jac)

# Morisita-Horn Bacteria
merged.dist.long %>%
  filter(Index == 'MorisitaHorn',
         Microbial.group == 'Bacteria') -> merged.dist.long.bac.horn
wilcox.test(Distance ~ Group, data = merged.dist.long.bac.horn)

# Jaccard Fungi
merged.dist.long %>%
  filter(Index == 'Jaccard',
         Microbial.group == 'Fungi') -> merged.dist.long.fun.jac
wilcox.test(Distance ~ Group, data = merged.dist.long.fun.jac)

# Morisita-Horn Fungi
merged.dist.long %>%
  filter(Index == 'MorisitaHorn',
         Microbial.group == 'Fungi') -> merged.dist.long.fun.horn
wilcox.test(Distance ~ Group, data = merged.dist.long.fun.horn)
