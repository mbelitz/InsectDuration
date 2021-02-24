library(tidyverse)
library(sf)
library(rnaturalearth)

# study location map

# read in data
mdf <- read.csv("data/LMM_data/phenoTraits_data_withNumObsData.csv", stringsAsFactors = FALSE) %>% 
  filter(scientificName != "Diabrotica undecimpunctata")
library(phyr)
library(INLA)
insect_tree3 = ape::read.tree("data/phylogeny/insect_test2.tre")
plot(insect_tree3, type = "fan")

insect_tree3$tip.label <- word(insect_tree3$tip.label, start = 1, end = 2, sep = "_") %>% 
  sub(pattern = "_", replacement = " ")


eusocial <- c("Bombus bimaculatus",
              "Bombus fervidus",
              "Bombus griseocollis",
              "Bombus huntii",
              "Bombus impatiens",
              "Bombus melanopygus",
              "Bombus pensylvanicus",
              "Bombus sonorus",
              "Bombus ternarius",
              "Polistes fuscatus")

tree_sp <- insect_tree3$tip.label

mdf_phylo_spp <- mdf %>% 
  filter(scientificName %in% tree_sp)

spp_year_cell <- mdf_phylo_spp %>% 
  group_by(lat, lon) %>% 
  summarise(Combinations = n(), years = length(unique(year)), 
            spp = length(unique(scientificName)))

na2 <- maps::map(regions=c("usa", "mexico", "canada"), plot = FALSE, fill = TRUE) %>%
  st_as_sf() %>% 
  st_simplify(dTolerance = 0.1)
na_equalArea <- st_transform(na2, crs="+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs ")
lake <- rnaturalearth::ne_download(scale = 110, type = 'lakes', category = 'physical')
lake_sf <- st_as_sf(lake)
great_lakes <- lake_sf %>% 
  filter(name_alt == "Great Lakes") %>% 
  st_transform(na2, crs="+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs ")

ggplot() + 
  geom_sf(na_equalArea, mapping = aes(), fill = "grey95") + 
  geom_sf(great_lakes, mapping = aes(), fill = "paleturquoise3", alpha = 0.4) +
  geom_tile(spp_year_cell, 
            mapping = aes(x = lon, y = lat, fill = Combinations)) +
  coord_sf(xlim = c(-2303077, 1900923), ylim = c(-2305719,1149281)) +
  scale_fill_gradient(high = "magenta4", low = "orange", trans = "log", breaks = c(1,7,50)) +
  #scale_fill_viridis_c(option = "inferno", trans = "log", breaks = c(1,7,50), direction = -1) + 
  labs(x = "", y = "", fill = "Pheno-estimates") +
  theme_void()

ggsave(filename = "Figures/studyarea.png", dpi = 400, width = 8, height = 6)  
