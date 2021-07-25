library(tidyverse)
library(raster)

# read in data
# removed D. undecimpunctata from analyses b/c there is no true diapause in this species
# removed L. cyanea b/c it wasn't in phylogeny
mdf <- read.csv("data/LMM_data/phenoTraits_data_withNumObsData.csv", stringsAsFactors = FALSE) %>% 
  filter(scientificName != "Diabrotica undecimpunctata") %>% 
  filter(scientificName != "Libellula cyanea") %>% 
  mutate(numObs = scale(numObs)) # Note the other variables were already scaled

ggplot(mdf, mapping = aes(x = lon, y = lat, fill = onset)) + 
  geom_tile() +
  theme_classic()

r <- raster("data/elevation_raster/elevation_NorthAmerica.tif")
rproj <- projectRaster(r, crs =  "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs ")
rdf <- raster::as.data.frame(rproj, xy = T) %>% 
  rename("lon" = "x", "lat" = "y")

ggplot() + 
  geom_tile(rdf, mapping = aes(x = lon, y = lat, fill = elevation_NorthAmerica)) +
  geom_point(mdf, mapping = aes(x = lon, y = lat), color = "yellow") +
  scale_fill_viridis_c(na.value = "transparent") +
  theme_classic()


mdf2 <- left_join(mdf, rdf)

coordinates(mdf) <- ~ lon + lat
elev <- extract(rproj, mdf)

mdf2 <- as.data.frame(mdf) %>% 
  dplyr::mutate(elevation = elev)

ggplot() + 
  geom_tile(mdf2, mapping = aes(x = lon, y = lat, fill = elevation)) +
  scale_fill_viridis_c(trans = "log", na.value = "transparent") +
  theme_classic()

lowest_cell <- mdf2 %>% 
  filter(lat == min(mdf2$lat))

lowest_cell_lon <- lowest_cell$lon
lowest_cell_lat <- lowest_cell$lat
lowest_cell_elev <- lowest_cell$elevation

mdf3 <- mdf2 %>% 
  mutate(lon_cor = lon - lowest_cell_lon,
         lat_cor = lat - lowest_cell_lat,
         elev_cor = round((elev * 3.28) - (lowest_cell_elev * 3.28)))

# for every m North, phenology delays (4/100000) days
# for every m east, phenology delays (4/500000) days
# for every foot higher in elevation, phenology delays (1/100) days

mdf3 <- mdf3 %>% 
  mutate(hopkins_cor =  (lat_cor * -(4/100000)) + 
                        (lon_cor * -(4/500000)) +
                        (elev_cor * -(1/100)))

ggplot() + 
  geom_tile(mdf3, mapping = aes(x = lon, y = lat, fill = hopkins_cor)) +
  scale_fill_gradient2() +
  theme_classic()


mdf4 <- mutate(mdf3, hopkins_onset = onset + hopkins_cor)

spp_onset <- mdf4 %>% 
  group_by(scientificName) %>% 
  summarise(mean_onset = mean(onset), sd_onset = sd(onset))

ggplot() + 
  geom_histogram(spp_onset, mapping = aes(x = mean_onset), fill = "turquoise") +
  geom_vline(xintercept = c(79, 172, 265))
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() 

## Can't use equinox to group species into traits. let's try quantiles
quantile(spp_onset$mean_onset, c(0.25,0.5, 0.75))

ggplot() + 
  geom_histogram(spp_onset, mapping = aes(x = mean_onset), fill = "turquoise") +
  geom_vline(xintercept = c(123, 167)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() 

spp_onset <- spp_onset %>% 
  mutate(seas2 = case_when(mean_onset <= 123 ~ "Spring",
                          mean_onset > 123 & mean_onset < 167 ~ "Summer",
                          mean_onset >167 ~ "Fall"))

mdf5 <- left_join(mdf4, spp_onset)
nrow(filter(mdf5, seas == seas2))
nrow(filter(mdf5, seas != seas2))
184/2459
t <- filter(mdf5, seas != seas2)
unique(t$scientificName)
head(mdf5)

write.csv(x = mdf5, file = "data/LMM_data/phenoTraits_data_withNumObsData_updatedSeasTrat.csv")
