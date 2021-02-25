library(data.table)
library(dplyr)
library(sf)

setwd("/srv/gspeedq32/mbelitz/insect-pheno-duration/")

library(tidyverse)

# combine all predictors into a single dataframe

files <- list.files("phenesse_outputs_25km_V2/", recursive = TRUE, full.names = TRUE)

total_duration <- do.call(rbind, lapply(files, readr::read_csv))

# make grid cells
na <-  rnaturalearth::ne_countries(country = c("United States of America", "Mexico", "Canada"),
                                   returnclass = "sp")
na_map <-  rnaturalearth::ne_countries(country = c("United States of America", "Mexico", "Canada"),
                                       returnclass = "sf")

# make grid over na
na_map <- st_transform(na_map, crs = "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs ")
grids <-  st_make_grid(na_map, cellsize = c(25000, 25000))
grids_sf <-  mutate(st_sf(geometry = grids), id_cells = 1:n())

coords <- st_centroid(grids) %>% 
  st_coordinates() %>% 
  data.frame()

grid_coords <- bind_cols(grids_sf, coords) %>% 
  st_set_geometry(NULL)

total_duration2 <- left_join(total_duration, grid_coords) %>% 
  mutate(lon = X, lat = Y) %>% 
  mutate(OS = paste(Order, scientificName, sep = ":"))

ggplot(total_duration2) +
  geom_tile(aes(x = lon, y = lat, fill = duration))

write.csv(total_duration2, "total_phenesse_outputs/total_duration_25km_V2.csv")