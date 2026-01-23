## Script written to assess variation in community taxonomy.
## The two vegetation types explored are woodlands dominated by Prosopis
## glandulosa, Honey mesquite, as well as other woody plant species and 
## grasslands dominated by herbaceous plants.
## Written by Dr. Elizabeth Bowman, July 3, 2025
## eabowman@utexas.edu

library(tidyverse)

fung.soil.data <- read.csv('data/data.fungal.Overall.csv')
bac.soil.data <- read.csv('data/data.bacterial.Overall.csv')

## Taxonomy: Fungi ----
# read in taxonomic data
tax.data <- read.csv('data/TaxonomyFungi.csv')

## Taxonomy: Bacteria ----
# read in taxonomic data
qiime.data <- read.csv('data/BacterialTaxonomy_Qiime.csv',
                       as.is = T)

# ncbi.data <- read.csv(paste('G022.results/BacterialTaxonomy/',
#                             'BacterialNCBI_zotu99_TaxonomySummary_genus.csv'))


#--------------------------------------------------------------#
# Fungi ----
#--------------------------------------------------------------#

### Create taxonomic dataframe
fung.soil.data %>%
  dplyr::select(block, disturbance, vegetation.type, invasion, pH, P.mgperkg,
                C.N, Mn.mgperkg, plant.shannon, plant.spec.richness,
                contains('Otu')) %>%
  gather(key = Otu, value = read.abund, -block, -disturbance, -vegetation.type,
         -invasion, -pH, -P.mgperkg, -C.N, -Mn.mgperkg, -plant.shannon,
         -plant.spec.richness) %>%
  filter(read.abund > 1) %>%
  left_join(., tax.data, by = 'Otu') -> tax.full

# Change vegetation type labels from Grassland and Woody patch to G and W
tax.full$vegetation.type <- factor(tax.full$vegetation.type,
                                   levels = c('Grassland', 'Woody patch'),
                                   labels = c('G', 'W'))


# summary
tax.full %>%
  group_by(class, family, vegetation.type) %>%
  summarise(otu.count = n(),
            read.abund = sum(read.abund)) %>%
  filter(read.abund > 1000) -> summary

# Species richness and abundance dataframe
fung.soil.data %>%
  dplyr::select(sample, block, disturbance, vegetation.type, invasion, pH, P.mgperkg,
                C.N, Mn.mgperkg, plant.shannon, plant.spec.richness,
                contains('Otu')) %>%
  gather(key = Otu, value = read.abund, -sample, -block, -disturbance, -vegetation.type,
         -invasion, -pH, -P.mgperkg, -C.N, -Mn.mgperkg, -plant.shannon,
         -plant.spec.richness) %>%
  filter(read.abund > 1) %>%
  left_join(., tax.data, by = 'Otu') %>%
  group_by(sample, vegetation.type, Otu, phylum, family, order, genus) %>%
  summarise(read.abund = sum(read.abund)) %>%
  mutate(count = 1) -> abund.full

### Phylum overall ----
### Vegetation type
# Chi-square
tax.full %>%
  filter(!phylum %in% c('unknown','Zygomycota')) %>%
  dplyr::select(vegetation.type, phylum) %>%
  mutate(count = 1) %>%
  group_by(vegetation.type, phylum) %>%
  summarise(Otu.count = sum(count)) %>%
  spread(key = phylum, value = Otu.count) -> chi.phylum.veg
# change to data frame and add vegetation type as row names
chi.phylum.veg <- as.data.frame(chi.phylum.veg)
rownames(chi.phylum.veg) <- chi.phylum.veg$vegetation.type
chi.phylum.veg <- chi.phylum.veg[-1]

chisq.test(chi.phylum.veg)

# plot
tax.full %>%
  filter(!phylum %in% c('unknown','Zygomycota')) %>%
  ggplot(aes(x = vegetation.type,
             fill = factor(phylum))) +
  geom_bar(position = 'fill') +
  # facet_grid(. ~ vegetation.type * disturbance) +
  xlab('') +
  ylab('Relative OTU abundance') +
  labs(fill = 'Phylum') +
  scale_fill_brewer(palette = 'Spectral', direction = -1) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 10, color = 'black'),
        axis.text.x = element_text(size = 10, color = 'black',
                                   # angle = 45, 
                                   # hjust = 1
        ),
        axis.title = element_text(size = 12),
        #strip.text.x = element_blank(),
        #legend.position = 'none',
        aspect.ratio = 5/1,
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.key.height = unit(0.35, "cm"),
        legend.key.width  = unit(0.35, "cm"),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  ) -> Fig.5a

### Vegetation type: Ascomycota class -----
# Chi-square
tax.full %>%
  filter(phylum == 'Ascomycota',
         !class %in% c('Candelariomycetes',
                       'Lecanoromycetes',
                       'Orbiliomycetes',
                       'Xylobotryomycetes')) %>%
  dplyr::select(vegetation.type, class) %>%
  mutate(count = 1) %>%
  group_by(vegetation.type, class) %>%
  summarise(Otu.count = sum(count)) %>%
  spread(key = class, value = Otu.count) -> chi.asco.veg
# change to data frame and add vegetation type as row names
chi.asco.veg <- as.data.frame(chi.asco.veg)
rownames(chi.asco.veg) <- chi.asco.veg$vegetation.type
chi.asco.veg <- chi.asco.veg[-1]

chisq.test(chi.asco.veg)

# plot
tax.full %>%
  filter(phylum == 'Ascomycota',
         !class %in% c('Candelariomycetes',
                       'Lecanoromycetes',
                       'Orbiliomycetes',
                       'Xylobotryomycetes')) %>%
  ggplot(aes(x = vegetation.type, 
             fill = factor(class))) +
  geom_bar(position = 'fill') +
  # facet_grid(. ~ vegetation.type * disturbance) +
  xlab('') +
  ylab('Relative OTU abundance') +
  labs(fill = "Ascomycota\nClass") +
  scale_fill_brewer(palette = 'Blues') +
  theme_classic() +
  theme(axis.text.y = element_text(size = 10, color = 'black'),
        axis.text.x = element_text(size = 10, color = 'black',
                                   # angle = 45, hjust = 1
        ),
        axis.title = element_text(size = 12),
        #strip.text.x = element_blank(),
        #legend.position = 'none',
        aspect.ratio = 5/1,
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.key.height = unit(0.35, "cm"),
        legend.key.width  = unit(0.35, "cm"),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) -> Fig.5b

### Vegetation type: Basidiomycota order ----
# Chi-square
tax.full %>%
  filter(class == 'Agaricomycetes') %>%
  dplyr::select(vegetation.type, order) %>%
  mutate(count = 1) %>%
  group_by(vegetation.type, order) %>%
  summarise(Otu.count = sum(count)) %>%
  spread(key = order, value = Otu.count) -> chi.order.veg
# change to data frame and add vegetation type as row names
chi.order.veg <- as.data.frame(chi.order.veg)
rownames(chi.order.veg) <- chi.order.veg$vegetation.type
chi.order.veg <- chi.order.veg[-1]

# change NAs to 0s
chi.order.veg[is.na(chi.order.veg)] <- 0

# filter out low abundance
chi.order.veg[colSums(chi.order.veg) > 10] -> chi.order.veg

chisq.test(chi.order.veg)

# plot
tax.full %>%
  filter(class == 'Agaricomycetes',
         !order %in% c('Atheliales', 'Hymenochaetales', 'Corticiales')) %>%
  mutate(order = replace(order, order == 'Agaricomycetes_unidentified', 'Unidentified')) %>%
  ggplot(aes(x = vegetation.type, 
             fill = factor(order))) +
  geom_bar(position = 'fill') +
  # facet_grid(. ~ vegetation.type * disturbance) +
  xlab('') +
  ylab('Relative OTU abundance') +
  labs(fill = "Basidiomycota\norder") +
  scale_fill_brewer(palette = 'PRGn') +
  theme_classic() +
  theme(axis.text.y = element_text(size = 10, color = 'black'),
        axis.text.x = element_text(size = 10, color = 'black',
                                   # angle = 45, hjust = 1
        ),
        axis.title = element_text(size = 12),
        #strip.text.x = element_blank(),
        #legend.position = 'none',
        aspect.ratio = 5/1,
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.key.height = unit(0.35, "cm"),
        legend.key.width  = unit(0.35, "cm"),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) -> Fig.5c

# summary of Agaricales by genus
tax.full %>%
  filter(order == 'Agaricales') %>%
  group_by(class, genus, vegetation.type) %>%
  summarise(otu.count = n(),
            read.abund = sum(read.abund)) %>%
  filter(read.abund > 1) -> summary


### Glomeromycotina ----
### Vegetation type: Glomeromycotina genus
# Chi-square
tax.full %>%
  filter(subphylum == 'Glomeromycotina') %>%
  dplyr::select(vegetation.type, genus) %>%
  mutate(count = 1) %>%
  group_by(vegetation.type, genus) %>%
  summarise(Otu.count = sum(count)) %>%
  spread(key = genus, value = Otu.count, fill = 0) -> chi.AM.veg
# change to data frame and add vegetation type as row names
chi.AM.veg <- as.data.frame(chi.AM.veg)
rownames(chi.AM.veg) <- chi.AM.veg$vegetation.type
chi.AM.veg <- chi.AM.veg[-c(1)]

chisq.test(chi.AM.veg)

tax.full %>%
  filter(subphylum == 'Glomeromycotina',
         !genus %in% c('unknown', 'Orientoglomus')) %>%
  mutate(genus = replace(genus, genus == 'Archaeosporaceae_unidentified',
                         'Archaeosporaceae'),
         genus = replace(genus, genus == 'Diversisporaceae_unidentified',
                         'Diversisporaceae')) %>%
  ggplot(aes(x = vegetation.type, 
             fill = factor(genus))) +
  geom_bar(position = 'fill') +
  xlab('') +
  ylab('Relative OTU abundance') +
  labs(fill = "Glomeromycotina\ngenus") +
  scale_fill_brewer(palette = 'Oranges') +
  theme_classic() +
  theme(axis.text.y = element_text(size = 10, color = 'black'),
        axis.text.x = element_text(size = 10, color = 'black',
                                   # angle = 45, hjust = 1
        ),
        axis.title = element_text(size = 12),
        #strip.text.x = element_blank(),
        #legend.position = 'none',
        aspect.ratio = 5/1,
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.key.height = unit(0.35, "cm"),
        legend.key.width  = unit(0.35, "cm"),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) -> Fig.5d

### Species richness and read abundance of Glomeromycotina ----
# isolate OTUs identified to Glomeromycotina
tax.full %>%
  filter(subphylum == 'Glomeromycotina',
         !genus %in% c('unknown', 'Orientoglomus')) %>%
  dplyr::select(Otu) %>%
  distinct(.) -> Glomero.Otu

# isolate fungal OTU table
fung.comm <- dplyr::select(fung.soil.data, starts_with('Otu'))
# pull out OTU identified as Glomeromycotina
fung.comm[colnames(fung.comm) %in% Glomero.Otu$Otu] -> glomero.comm

# assess read abundance per OTU
colSums(glomero.comm)

# calculate species richness
glomero.sr <- specnumber(glomero.comm)

# calculate abundance
glomero.abund <- rowSums(glomero.comm)

# bind metadata, species richness, and community data
fung.soil.data %>%
  dplyr::select(sample, vegetation.type, invasion, disturbance) %>%
  mutate(spec.rich = glomero.sr, 
         abundance = glomero.abund) -> glomero.data

#### Species richness ----
# assess normality
shapiro.test(glomero.data$spec.rich)
histogram(glomero.data$spec.rich)

# Kruskal wallis
kruskal.test(glomero.sr ~ vegetation.type,
             data = glomero.data)

# Change labels for vegetation type to G and W as above
glomero.data$vegetation.type <- factor(glomero.data$vegetation.type,
                                   levels = c('Grassland', 'Woody patch'),
                                   labels = c('G', 'W'))

ggplot(glomero.data, aes(x = vegetation.type,
                         y = spec.rich,
                         fill = vegetation.type)) +
  geom_boxplot(width = 0.25) +
  geom_violin(alpha = 0.5) +
  xlab("") +
  ylab("OTU richness") +
  theme_classic() +
  scale_fill_manual(values = c("#905971", "#F1BB83")) +
  theme(axis.text.y = element_text(size = 10, color = 'black'),
        axis.text.x = element_text(size = 10, color = 'black',
                                   # angle = 45, hjust = 1
        ),
        axis.title = element_text(size = 12),
        #strip.text.x = element_blank(),
        legend.position = 'none',
        aspect.ratio = 5/1,
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) -> Fig.5e

#### Abundance ----
# assess normality
shapiro.test(log1p(glomero.data$abundance))
histogram(log1p(glomero.data$abundance))

mutate(glomero.data, log1p.glomero.abund = log1p(abundance)) -> glomero.data

# t-test
t.test(log1p.glomero.abund ~ vegetation.type,
       data = glomero.data)

ggplot(glomero.data, aes(x = vegetation.type,
                         y = abundance,
                         fill = vegetation.type)) +
  geom_boxplot(width = 0.25) +
  geom_violin(alpha = 0.5) +
  xlab("") +
  ylab("Read abundance") +
  ylim(0, 100) +
  scale_fill_manual(values = c("#905971", "#F1BB83")) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 10, color = 'black'),
        axis.text.x = element_text(size = 10, color = 'black',
                                   # angle = 45, hjust = 1
        ),
        axis.title = element_text(size = 12),
        #strip.text.x = element_blank(),
        legend.position = 'none',
        aspect.ratio = 5/1,
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) -> Fig.5f

#--------------------------------------------------------------#
# Bacteria ----
#--------------------------------------------------------------#

### Create taxonomic dataframe
bac.soil.data %>%
  dplyr::select(block, disturbance, vegetation.type, invasion,
                contains('Otu')) %>%
  gather(key = Otu, value = read.abund, -block, -disturbance, -vegetation.type, -invasion) %>%
  filter(read.abund > 10) %>%
  left_join(., qiime.data, by = 'Otu') %>%
  dplyr::select(-Taxon) %>%
  filter(Confidence > 0.7) -> tax.full

# make read.abund column numeric
tax.full$read.abund <- as.numeric(tax.full$read.abund)

# summary
tax.full %>%
  group_by(vegetation.type, class) %>%
  summarise(otu.count = n(),
            read.abund = sum(read.abund)) %>%
  filter(read.abund > 1000) -> summary

# Species richness and abundance dataframe
bac.soil.data %>%
  dplyr::select(sample, block, disturbance, vegetation.type, invasion,
                contains('Otu')) %>%
  gather(key = Otu, value = read.abund, -sample, -block, -disturbance, -vegetation.type,
         -invasion) %>%
  filter(read.abund > 1) %>%
  left_join(., qiime.data, by = 'Otu') %>%
  group_by(sample, vegetation.type, Otu, phylum, family, order, genus) %>%
  summarise(read.abund = sum(read.abund)) %>%
  mutate(count = 1) -> abund.full

### Phylum overall ----
# Isolate phyla with greater than 250 reads, filter out low abundance phyla
tax.full %>% 
  group_by(phylum) %>%
  summarise(read.abund = sum(read.abund)) %>%
  filter(read.abund > 750,
         !phylum %in% c("","Ascomycota")) -> phylum.groups

### Vegetation type
# Chi-square
tax.full %>%
  filter(phylum %in% phylum.groups$phylum,
         phylum != 'WPS-2') %>%
  dplyr::select(vegetation.type, phylum) %>%
  mutate(count = 1) %>%
  group_by(vegetation.type, phylum) %>%
  summarise(Otu.count = sum(count)) %>%
  spread(key = phylum, value = Otu.count) -> chi.phylum.veg
# change to data frame and add vegetation type as row names
chi.phylum.veg <- as.data.frame(chi.phylum.veg)
rownames(chi.phylum.veg) <- chi.phylum.veg$vegetation.type
chi.phylum.veg <- chi.phylum.veg[-1]

chisq.test(chi.phylum.veg)
# create color pallete
mycolors = c(RColorBrewer::brewer.pal(name="BrBG", n = 11),
             RColorBrewer::brewer.pal(name="Greens", n = 5))

# Change labels for vegetation type
tax.full$vegetation.type <- factor(tax.full$vegetation.type,
                                   levels = c('Grassland', 'Woody patch'),
                                   labels = c('G', 'W'))

### Vegetation type
tax.full %>%
  # filter(phylum %in% c('Acidobacteriota','Proteobacteria','Actinobacteriota',
  #                      'Methylomirabilota','Firmicutes','Planctomycetota',
  #                      'Myxococcota','Cyanobacteria','Bacteroidota',
  #                      'Gemmatimonadota','Chloroflexi','Desulfobacterota',
  #                      'Nitrospirota','Verrucomicrobiota','Elusimicrobiota',
  #                      'Abditibacteriota','WPS-2','Thermoplasmatota',
  #                      'Armatimonadota','Entotheonellaeota','Dadabacteria',
  #                      'NB1-j')) %>%
  filter(phylum %in% phylum.groups$phylum,
         phylum != 'WPS-2') %>%
  ggplot(aes(x = vegetation.type,
             fill = factor(phylum))) +
  geom_bar(position = 'fill') +
  # facet_grid(. ~ vegetation.type * disturbance) +
  xlab('') +
  ylab('Relative OTU abundance') +
  labs(fill = "Phylum") +
  scale_fill_manual(values = mycolors) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 10, color = 'black'),
        axis.text.x = element_text(size = 10, color = 'black',
                                   # angle = 45, hjust = 1
        ),
        axis.title = element_text(size = 12),
        strip.text.x = element_blank(),
        legend.position = 'right',
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.key.height = unit(0.35, "cm"),
        legend.key.width  = unit(0.35, "cm"),
        aspect.ratio = 5/1,
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) -> Fig.5g

### Actinobacteriota ----
# Isolate phyla with greater than 250 reads, filter out low abundance phyla
tax.full %>% 
  filter(phylum == "Actinobacteriota") %>%
  group_by(class) %>%
  summarise(read.abund = sum(read.abund)) -> actino.groups

### Vegetation type
# Chi-square
tax.full %>%
  filter(phylum == "Actinobacteriota") %>%
  dplyr::select(vegetation.type, class) %>%
  mutate(count = 1) %>%
  group_by(vegetation.type, class) %>%
  summarise(Otu.count = sum(count)) %>%
  spread(key = class, value = Otu.count) -> chi.class.veg
# change to data frame and add vegetation type as row names
chi.class.veg <- as.data.frame(chi.class.veg)
rownames(chi.class.veg) <- chi.class.veg$vegetation.type
chi.class.veg <- chi.class.veg[-1]

chisq.test(chi.class.veg)

### Vegetation type
tax.full %>%
  filter(phylum == "Actinobacteriota") %>%
  ggplot(aes(x = vegetation.type,
             fill = factor(class))) +
  geom_bar(position = 'fill') +
  # facet_grid(. ~ vegetation.type * disturbance) +
  xlab('') +
  ylab('Relative OTU abundance') +
  labs(fill = "Actinobacteriota\norder") +
  scale_fill_manual(values = mycolors) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 10, color = 'black'),
        axis.text.x = element_text(size = 10, color = 'black',
                                   # angle = 45, hjust = 1
        ),
        axis.title = element_text(size = 12),
        strip.text.x = element_blank(),
        legend.position = 'right',
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.key.height = unit(0.35, "cm"),
        legend.key.width  = unit(0.35, "cm"),
        aspect.ratio = 5/1,
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) -> Fig.5h

### Cyanobacteria ----
### Isolate OTU assigned to phylum Cyanobacteria
tax.full %>%
  filter(phylum == 'Cyanobacteria') %>%
  dplyr::select(Otu) %>%
  distinct() -> Cyanobacteria.otu
# isolate OTU belonging to Cyanbacteria
cyano.otu <-colnames(bac.soil.data[colnames(bac.soil.data) %in% Cyanobacteria.otu$Otu])

bac.soil.data %>%
  dplyr::select(sample, plot, block,
                disturbance, vegetation.type,invasion,
                N.mgperkg, N.percent, Nag_umolgSOC.1h.1,
                all_of(cyano.otu)) -> cyano.data

# isolate community data to calculate species richness and diversity
cyano.comm <- dplyr::select(cyano.data, starts_with('Otu'))

# Chi-square
tax.full %>%
  filter(phylum == 'Cyanobacteria',
         family != 'Chloroplast'
         ) %>%
  dplyr::select(vegetation.type, family) %>%
  mutate(count = 1) %>%
  group_by(vegetation.type, family) %>%
  summarise(Otu.count = sum(count)) %>%
  spread(key = family, value = Otu.count) -> chi.cyano.veg
# change to data frame and add vegetation type as row names
chi.cyano.veg <- as.data.frame(chi.cyano.veg)
rownames(chi.cyano.veg) <- chi.cyano.veg$vegetation.type
chi.cyano.veg <- chi.cyano.veg[-1]

# change NAs to 0s
chi.cyano.veg[is.na(chi.cyano.veg)] <- 0

chisq.test(chi.cyano.veg)

tax.full %>%
  filter(phylum == 'Cyanobacteria',
         family != 'Chloroplast'
         ) %>%
  mutate(family = ifelse(family == "", "Unknown", family),
         family = ifelse(family == 'uncultured', 'Uncultured', family)) %>%
  ggplot(aes(x = vegetation.type,
             fill = factor(family))) +
  geom_bar(position = 'fill') +
  # facet_grid(. ~ vegetation.type * disturbance) +
  xlab('') +
  ylab('Relative OTU abundance') +
  labs(fill = "Cyanobacteria\nfamily") +
  scale_fill_manual(values = mycolors) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 10, color = 'black'),
        axis.text.x = element_text(size = 10, color = 'black',
                                   # angle = 45, hjust = 1
        ),
        axis.title = element_text(size = 12),
        strip.text.x = element_blank(),
        legend.position = 'right',
        aspect.ratio = 5/1,
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.key.height = unit(0.35, "cm"),
        legend.key.width  = unit(0.35, "cm"),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) -> Fig.5i

### Species richness and read abundance of Cyanobacteria ----
# Diversity calculation is odd, instead will focus on species richness and abundance
cyano.data$diversity <- fisher.alpha(cyano.comm) 

cyano.data$spec.richness <- specnumber(cyano.comm)
hist(cyano.data$spec.richness)
shapiro.test(cyano.data$spec.richness) # not normal

cyano.data$abundance <- rowSums(cyano.comm)
hist(cyano.data$abundance)
shapiro.test(cyano.data$abundance) # normal

#### Species richness -----
kruskal.test(spec.richness ~ vegetation.type, data = cyano.data)

# Change labels for vegetation type
cyano.data$vegetation.type <- factor(cyano.data$vegetation.type,
                                     levels = c('Grassland','Woody patch'),
                                     labels = c('G', 'W'))

# plot for Species richness
ggplot(cyano.data, aes(x = vegetation.type,
                       y = spec.richness,
                       fill = vegetation.type)) +
  geom_boxplot() +
  geom_violin(alpha = 0.5) +
  #facet_grid(. ~ invasion * disturbance) +
  xlab('') +
  ylab('OTU richness') +
  scale_fill_manual(values = c("#905971", "#F1BB83")) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 10, color = 'black'),
        axis.text.x = element_text(size = 10, color = 'black',
                                   # angle = 45, hjust = 1
        ),
        axis.title = element_text(size = 12),
        #strip.text.x = element_blank(),
        legend.position = 'none',
        aspect.ratio = 5/1,
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) -> Fig.5j

#### Read abundance -----
kruskal.test(abundance ~ vegetation.type, data = cyano.data)

# plot for read abundance
ggplot(cyano.data, aes(x = vegetation.type, 
                                      y = abundance,
                                      fill = vegetation.type)) +
  geom_boxplot() +
  geom_violin(alpha = 0.5) +
  #facet_grid(. ~ invasion * disturbance) +
  xlab('') +
  ylab('Read abundance') +
  ylim(0, 200) +
  scale_fill_manual(values = c("#905971", "#F1BB83")) +
  theme_classic() +
 theme(axis.text.y = element_text(size = 10, color = 'black'),
        axis.text.x = element_text(size = 10, color = 'black',
                                   # angle = 45, hjust = 1
        ),
        axis.title = element_text(size = 12),
        #strip.text.x = element_blank(),
        legend.position = 'none',
        aspect.ratio = 5/1,
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) -> Fig.5k

# Generate Figure 5
(# ---- Top row (6 plots) ----
  Fig.5a | Fig.5b | Fig.5c | Fig.5d | Fig.5e | Fig.5f) /
(# ---- Bottom row (5 plots + spacer) ----
  Fig.5g | Fig.5h | Fig.5i | Fig.5j | Fig.5k) +
  plot_layout(
    widths = c(1, 1, 1, 1, 0.6, 0.6),
    heights = 0.9,
    guides = "keep"
  ) + 
  plot_annotation(tag_levels = "a") -> Fig.5
  
ggsave("figures/Figure5.jpg",
       plot = Fig.5,
       width = 15, height = 7,
       device = 'jpg',
       dpi = 600)
