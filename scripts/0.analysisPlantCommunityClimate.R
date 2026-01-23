## Script written to assess in plant community and climate across our
## sampled sites.
## Written by Dr. Elizabeth Bowman, Dec. 9, 2024
## eabowman@utexas.edu
## 
## For the plant community, I focused mainly on invaded areas as uninvaded areas
## had zero Guinea grass.
# 1. Climate
# 2. Guinea grass cover
# 3. Plant community

pl.data <- read.csv('data/plant.comm.sxs.csv')
field.data <- read.csv('data/field.data_climate.csv')

# Restructure plant community data to isolate Megathyrsis maximus cover ----
pl.data %>%
  dplyr::select(disturbance, vegetation, invasion,
                replicate, block, Megathyrsus.maximus) %>%
  group_by(disturbance, vegetation, invasion, replicate,  block) %>%
  summarise(Meg.max = Megathyrsus.maximus) -> mm.data

mm.data$invasion <- factor(mm.data$invasion,
                           levels = c('invaded', 'uninvaded'),
                           labels = c('Invaded', 'Uninvaded'))

mm.data$disturbance <- factor(mm.data$disturbance,
                           levels = c('disturbed', 'undisturbed'),
                           labels = c('Disturbed', 'Undisturbed'))

mm.data$vegetation <- factor(mm.data$vegetation,
                           levels = c('grassland', 'motte'),
                           labels = c('Grassland', 'Woody patch'))

# Factor plant community data ----
pl.data$invasion <- factor(pl.data$invasion,
                           levels = c('invaded', 'uninvaded'),
                           labels = c('Invaded', 'Uninvaded'))

pl.data$disturbance <- factor(pl.data$disturbance,
                              levels = c('disturbed', 'undisturbed'),
                              labels = c('Disturbed', 'Undisturbed'))

pl.data$vegetation <- factor(pl.data$vegetation,
                             levels = c('grassland', 'motte'),
                             labels = c('Grassland', 'Woody patch'))

#--------------------------------------------------------------#
# 1. Climate ----
#--------------------------------------------------------------#

## MAP ----
map.lm <- lm(MAP ~ block, data = field.data)
anova(map.lm)

max(field.data$MAP) - min(field.data$MAP)

# Read in raster files (if using multiple files, read them into a stack or brick)
# This file was not included in the repository as it is publicly available and very large. 
# To run this, please download the file at https://www.worldclim.org/data/worldclim21.html
raster_files <- raster("../G022.data/WorldClim/wc2.1_30s_bio/wc2.1_30s_bio_12.tif")

# Define the extent manually using extent()
subset_extent <- extent(c(-98.0, -97.85, 27.15, 27.25))

# Subset the raster using crop()
subset_raster <- crop(raster_files, subset_extent)

# Convert subset raster to a data frame
subset_df <- as.data.frame(subset_raster, xy = TRUE)

# isolate plot data
field.data %>%
  dplyr::select(block, lat, long) %>%
  rename(x = long, y = lat) -> plot.data

# rename MAP data in subset_df
subset_df %>%
  rename(MAP = wc2.1_30s_bio_12) -> subset_df

# Plotting subsetted data frame
ggplot() +
  geom_raster(data = subset_df, aes(x = x, y = y, fill = MAP)) +
  scale_fill_viridis_c() +
  labs(x = "Longitude", y = "Latitude", fill = "Mean Annual\nPrec. (mm)") +
  theme_minimal() +
  geom_point(data = plot.data, aes(x = x, y = y), size = 4) +
  theme(axis.title = element_text(size = 18, color = 'black'),
        axis.text = element_text(size = 16, color = 'black'),
        legend.title = element_text(size = 14,
                                    margin = margin(b = 15)),
        legend.text = element_text(size = 12))

ggsave('figures/SupplementaryFigS1.tiff', plot = last_plot(),
       width = 10, height = 10, units = 'in', device = 'tiff')

## MAT ----
mat.lm <- lm(MAT ~ block, data = field.data)
anova(mat.lm)

max(field.data$MAT) - min(field.data$MAT)

#--------------------------------------------------------------#
# 2. Guinea grass cover ----
#--------------------------------------------------------------#

## A. Distribution of Guinea grass cover ----
ggplot(mm.data, aes(x = Meg.max,
                    fill = invasion)) +
  geom_histogram(bins = 15, position = 'identity', alpha = 1.0, color = 'black') +
  ylab('Count') +
  xlab('Percent cover of M. maximus') +
  theme_classic() +
  scale_fill_manual(values = c('#8DB6AB', '#EDE6DE')) +
  labs(fill = 'Invasion') +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12,
                                 color = 'black'),
        legend.position = 'right') -> Fig.S1a

## B. Invaded sites only ----
# Isolate invaded site data
filter(mm.data, invasion == 'Invaded') -> mm.inv.data

### Mean Guinea grass cover ------
mm.inv.data %>%
  group_by(invasion) %>%
  summarise(mean.GGcover = mean(mm.inv.data$Meg.max),
            sd.GGcover = sd(mm.inv.data$Meg.max))

### assess cover across blocks to see if there is any bias ----
mm.inv.data$block <- factor(mm.inv.data$block)

kruskal.test(Meg.max ~ block, data = mm.inv.data) 
# no sig. difference in percent cover of Guinea grass between blocks

ggplot(mm.inv.data, aes(x = block,
                    y = Meg.max)) +
  geom_point(size = 0.1, color = 'darkgrey') +
  geom_boxplot(width = 0.25) +
  ylab('Percent cover') +
  xlab('') +
  theme_classic() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14,
                                 color = 'black'))

### assess cover as a function of disturbance and plant community to assess variation  ----
GG.lm <- lme(Meg.max ~ disturbance * vegetation, 
             random = ~ 1| block,
             data = mm.inv.data)
GG.anova <- anova(GG.lm)
# Overall higher percent cover of Guinea grass in grasslands
mm.inv.data %>%
  group_by(vegetation) %>%
  summarise(mean.Meg.max = mean(Meg.max),
            sd.Meg.max = sd(Meg.max))

ggplot(mm.inv.data, aes(x = vegetation,
                    y = Meg.max,
                    fill = vegetation)) +
  geom_point(size = 0.1, color = 'darkgrey') +
  geom_boxplot(width = 0.25) +
  scale_fill_manual(values = c("#905971", "#F1BB83")) +
  ylab(expression("Percent cover of " * italic ("M. maximus"))) +
  labs(fill = 'Plant community') +
  xlab('') +
  theme_classic() +
  theme(axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12,
                                 color = 'black'),
        axis.text.x = element_text(size = 12,
                                   color = 'black'),
        legend.position = 'right') -> Fig.S1b

# Generate Appendix S2: Fig. S1
Fig.S1a +
  Fig.S1b +
    plot_layout(
    nrow = 1,
    ncol = 2,
    guides = 'collect',
    heights = c(1,1)
  ) + 
  plot_annotation(tag_levels = "a") -> Fig.S1

ggsave("figures/FigureS1.jpg",
       plot = Fig.S1,
       width = 8, height = 5,
       device = 'jpg',
       dpi = 600)

#--------------------------------------------------------------#
# 3. Plant community ----
#--------------------------------------------------------------#
## Plant species richness ------
# isolate plant community data
pl.comm <- pl.data[8:length(pl.data)]

# remove bare and thatch columns
dplyr::select(pl.comm, -Bare, -Thatch) -> pl.comm

# calculate Species richness 
pl.data$plant.sr <- specnumber(pl.comm)

sr.lme <- lme(plant.sr ~ invasion * disturbance * vegetation, 
               random = ~ 1| block,
               data = pl.data)
anova(sr.lme)

# assess residuals of model
plot(sr.lme)

# Mean species richness based on invasion
pl.data %>%
  group_by(invasion) %>%
  summarise(mean.sr = mean(plant.sr),
            sd.sr = sd(plant.sr))

# Mean species richness based on vegetation
pl.data %>%
  group_by(vegetation) %>%
  summarise(mean.sr = mean(plant.sr),
            sd.sr = sd(plant.sr))

# Mean species richness based on invasion * disturbance
pl.data %>%
  group_by(invasion, disturbance) %>%
  summarise(mean.sr = mean(plant.sr),
            sd.sr = sd(plant.sr))

# Plot
ggplot(pl.data, aes(x = invasion,
                    y = plant.sr,
                    fill = invasion)) +
  geom_boxplot(width= 0.5) +
  facet_grid(. ~ vegetation) +
  scale_fill_manual(values = c("#8DB6AB", "#EDE6DE")) +
  ylab("Plant species richness") +
  xlab('') +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 10,
                                   color = 'black'),
        axis.text.x = element_text(size = 10,
                                   color = 'black',
                                   angle = 45,
                                   hjust = 1),
        legend.position = 'none',
        strip.text = element_text(size = 10,
                                  color = 'black'),
        strip.background = element_blank()) -> Fig.S3a


# Plot
ggplot(pl.data, aes(x = invasion,
                    y = plant.sr,
                    fill = invasion)) +
  geom_boxplot(width= 0.5) +
  facet_grid(. ~ disturbance) +
  scale_fill_manual(values = c("#8DB6AB", "#EDE6DE")) +
  ylab("Plant species richness") +
  xlab('') +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 10,
                                   color = 'black'),
        axis.text.x = element_text(size = 10,
                                   color = 'black',
                                   angle = 45,
                                   hjust = 1),
        legend.position = 'none',
        strip.text = element_text(size = 10,
                                  color = 'black'),
        strip.background = element_blank()) -> Fig.S3c


## Plant diversity ----

# calculate Shannon diversity
pl.data$plant.shannon <- diversity(pl.comm, index = 'shannon')
log1p(pl.data$plant.shannon) -> pl.data$log1p.plant.shannon

div.lme <- lme(log1p.plant.shannon ~ invasion * disturbance * vegetation, 
               random = ~ 1| block,
               data = pl.data)
anova(div.lme)

# Assess model residuals
plot(div.lme)

# Mean shannon based on invasion
pl.data %>%
  group_by(invasion) %>%
  summarise(mean.sh = mean(plant.shannon),
            sd.sh = sd(plant.shannon))

# Mean shannon based on vegetation
pl.data %>%
  group_by(vegetation) %>%
  summarise(mean.sh = mean(plant.shannon),
            sd.sh = sd(plant.shannon))

# Plot
ggplot(pl.data, aes(x = invasion,
                    y = plant.shannon,
                    fill = invasion)) +
  geom_boxplot(width= 0.5) +
  facet_grid(. ~ vegetation) +
  scale_fill_manual(values = c("#8DB6AB", "#EDE6DE")) +
  ylab("Plant Shannon's diversity") +
  xlab('') +
  labs(fill = 'Invasion') +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 10,
                                 color = 'black'),
        axis.text.x = element_text(size = 10,
                                   color = 'black',
                                   angle = 45,
                                   hjust = 1),
        legend.position = 'right',
        strip.text = element_text(size = 10,
                                  color = 'black'),
        strip.background = element_blank()) -> Fig.S3b

# Generate Appendix S2: Fig. S1
Fig.S3a +
  Fig.S3b +
  Fig.S3c +
    plot_layout(
    nrow = 1,
    ncol = 3,
    guides = 'collect',
    heights = c(1,1,1)
  ) + 
  plot_annotation(tag_levels = "a") -> Fig.S3

ggsave("figures/FigureS3.jpg",
       plot = Fig.S3,
       width = 10, height = 5,
       device = 'jpg',
       dpi = 600)

## Multivariate dispersion ----
# isolate plant community data
pl.comm <- pl.data[8:length(pl.data)]

# remove bare and thatch columns
dplyr::select(pl.comm, -Bare, -Thatch) -> pl.comm

# distance matrix 
pl.dist <- vegdist(pl.comm, method = 'bray')

# add a plot column
pl.data$plot_id <- with(pl.data,
  interaction(vegetation, invasion, disturbance, block, drop = TRUE)
)

### Quantify within-plot dispersion ----
disp_within <- betadisper(pl.dist, group = pl.data$plot_id)

within_vals <- disp_within$distances

### Quantify between-plot dispersion ----
disp_between <- betadisper(pl.dist, group = pl.data$block)


tapply(disp_within$distances, pl.data$plot_id, mean)
tapply(disp_between$distances, pl.data$block, mean)

permutest(disp_within)
