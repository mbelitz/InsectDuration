library(phyr)
library(ape)
library(tidyverse)
library(stringr)
library(grid)
library(ggplotify)
library(geoR)
library(sf)

# read in data
mdf <- read.csv("data/LMM_data/phenoTraits_data_withNumObsData_updatedSeasTrat.csv", stringsAsFactors = FALSE) %>% 
  filter(scientificName != "Diabrotica undecimpunctata") %>% 
  filter(scientificName != "Libellula cyanea") %>% 
  mutate(numObs = scale(numObs)) # Note the other variables we

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

insect_tree4 <- ape::drop.tip(phy = insect_tree3, tip = eusocial) %>% 
  ape::drop.tip(tip = "Aeshna cyanea")  %>% 
  ape::drop.tip(tip = "Diabrotica undecimpunctata")


## ## Add spatial component to pglmm

# set up spatial correlation covariance matrix
## grab coordinates for each unique cell
cell_id <- dplyr::distinct(mdf_phylo_spp, id_cells, .keep_all = T) %>% 
  dplyr::select(id_cells, lon, lat)
x.coord <- cell_id$lon
y.coord <- cell_id$lat

nsite <- length(x.coord)

# generate matrix
Dist <- matrix(0, nrow = nsite, ncol = nsite)
for (i in 1:nsite)
  for (j in 1:nsite)
    Dist[i, j] <-
  ((x.coord[i] - x.coord[j]) ^ 2 + (y.coord[i] - y.coord[j]) ^ 2) ^ .5

range <- 14446371
sd.space <- 30000 # sd of b0_site
V.space <- sd.space ^ 2 * exp(-Dist / range)
rownames(V.space) <- 1:nsite
colnames(V.space) <- 1:nsite
V.space <- V.space / max(V.space)
V.space <- V.space / exp(determinant(V.space)$modulus[1] / nsite)
row.names(V.space) <- cell_id$id_cells
colnames(V.space) <- cell_id$id_cells

## Fit models

pglmm_onset <- pglmm(onset ~ temp + pop + prec + temp_seas + numObs + 
                       diapause.stage + flights +
                       (1|id_cells) + (1|scientificName__) +
                       (0 + temp | scientificName__) + 
                       (0 + prec | scientificName__) +
                       (0 + temp_seas | scientificName__) +
                       temp:prec + temp:pop + 
                       temp:diapause.stage + prec:flights,
                     data = mdf_phylo_spp, 
                     cov_ranef = list(scientificName = insect_tree4), 
                     bayes = TRUE)


# new PGLMM offset with spatial correlation matrix
pglmm_offset_sp <- pglmm(offset ~ temp_seas + prec +
                           seas2 + diapause.stage + immature.habitat + larval.diet + 
                           temp_seas:immature.habitat +
                           numObs + 
                           (1 | id_cells__) + (1 | scientificName__) + 
                           (0 + temp_seas | scientificName__) + 
                           (0 + prec | scientificName__), 
                         data = mdf_phylo_spp, 
                         cov_ranef = list(scientificName = insect_tree4,
                                          id_cells = V.space), 
                         bayes = TRUE)

summary(pglmm_duration_sp)

## Duration

# new PGLMM offset with spatial correlation matrix
pglmm_duration_sp <- pglmm(duration ~ temp + prec + temp_seas +
                             diapause.stage + 
                             immature.habitat +
                             larval.diet + 
                             numObs + 
                             (1 | id_cells__) + (1 | scientificName__) + 
                             (0 + temp | scientificName__) + 
                             (0 + prec | scientificName__) + 
                             temp:prec + temp:larval.diet + temp:immature.habitat,
                           data = mdf_phylo_spp, 
                           cov_ranef = list(scientificName = insect_tree4,
                                            id_cells = V.space), 
                           bayes = TRUE)

summary(pglmm_duration_sp)

save(pglmm_onset, file = "ModelOutputs/resubmission/pglmm_onset.Rdata")
save(pglmm_offset_sp, file = "ModelOutputs/resubmission/pglmm_offset.Rdata")
save(pglmm_duration_sp, file = "ModelOutputs/resubmission/pglmm_duration.Rdata")

### Model Assumption Checks

# residuals
resids <- residuals(pglmm_onset)
hist(resids)
qqnorm(resids)
qqline(resids)
shapiro.test(resids)

on_r_hist <- ggplot() + 
  geom_histogram(mapping = aes(x = resids)) +
  theme_minimal() + 
  labs(x = "Residuals") +
  scale_y_continuous(expand = c(0,0)) +
  ggtitle(label = "Onset")

gg_qq <- function(x, distribution = "norm", ..., line.estimate = NULL, conf = 0.95,
                  labels = names(x)){
  q.function <- eval(parse(text = paste0("q", distribution)))
  d.function <- eval(parse(text = paste0("d", distribution)))
  x <- na.omit(x)
  ord <- order(x)
  n <- length(x)
  P <- ppoints(length(x))
  df <- data.frame(ord.x = x[ord], z = q.function(P, ...))
  
  if(is.null(line.estimate)){
    Q.x <- quantile(df$ord.x, c(0.25, 0.75))
    Q.z <- q.function(c(0.25, 0.75), ...)
    b <- diff(Q.x)/diff(Q.z)
    coef <- c(Q.x[1] - b * Q.z[1], b)
  } else {
    coef <- coef(line.estimate(ord.x ~ z))
  }
  
  zz <- qnorm(1 - (1 - conf)/2)
  SE <- (coef[2]/d.function(df$z)) * sqrt(P * (1 - P)/n)
  fit.value <- coef[1] + coef[2] * df$z
  df$upper <- fit.value + zz * SE
  df$lower <- fit.value - zz * SE
  
  if(!is.null(labels)){ 
    df$label <- ifelse(df$ord.x > df$upper | df$ord.x < df$lower, labels[ord],"")
  }
  
  p <- ggplot(df, aes(x=z, y=ord.x)) +
    geom_point() + 
    geom_abline(intercept = coef[1], slope = coef[2]) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha=0.2) +
    theme_minimal() + 
    labs(x = "Theoretical Quantiles", y = "Sample Quantiles") 
  if(!is.null(labels)) p <- p + geom_text( aes(label = label))
  return(p)
  coef
}


on_qq <- gg_qq(x = resids) 


# spatial autocorrelation map
na2 <- maps::map(regions=c("usa", "mexico", "canada"), plot = FALSE, fill = TRUE) %>%
  st_as_sf() %>% 
  st_simplify(dTolerance = 0.1)
na_equalArea <- st_transform(na2, crs="+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs ")
lake <- rnaturalearth::ne_download(scale = 110, type = 'lakes', category = 'physical')
lake_sf <- st_as_sf(lake)
great_lakes <- lake_sf %>% 
  filter(name_alt == "Great Lakes") %>% 
  st_transform(na2, crs="+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs ")

resids_df <- mdf_phylo_spp %>% 
  mutate(resids = resids)
resids_df <-  resids_df %>% 
  group_by(lon, lat) %>% 
  summarise(mean_resid = mean(resids))

on_resid_map <- ggplot() + 
  geom_sf(na_equalArea, mapping = aes(), fill = "grey97") + 
  geom_sf(great_lakes, mapping = aes(), fill = "paleturquoise3", alpha = 0.4) +
  geom_tile(resids_df, mapping = aes(x = lon, y = lat, fill = mean_resid)) +
  coord_sf(xlim = c(-2303077, 1900923), ylim = c(-2305719,1149281)) +
  scale_fill_gradient2() +
  labs(fill = "Mean Residual") +
  theme_void()

cp_onset <- cowplot::plot_grid(on_r_hist,on_qq, ncol = 1,
                               labels = c("A", "B"))

cp_onset

## OFFSET TIME

# residuals
resids_off <- residuals(pglmm_offset_sp)
hist(resids_off)
qqnorm(resids_off)
qqline(resids_off)
shapiro.test(resids_off)

off_r_hist <- ggplot() + 
  geom_histogram(mapping = aes(x = resids_off)) +
  theme_minimal() + 
  labs(x = "Residuals") +
  scale_y_continuous(expand = c(0,0)) +
  ggtitle(label = "Offset")

off_qq <- gg_qq(x = resids_off) 


# spatial autocorrelation
resids_df <- mdf_phylo_spp %>% 
  mutate(resids = resids_off)
resids_df <-  resids_df %>% 
  group_by(lon, lat) %>% 
  summarise(mean_resid = mean(resids))

off_resid_map <- ggplot() + 
  geom_sf(na_equalArea, mapping = aes(), fill = "grey97") + 
  geom_sf(great_lakes, mapping = aes(), fill = "paleturquoise3", alpha = 0.4) +
  geom_tile(resids_df, mapping = aes(x = lon, y = lat, fill = mean_resid)) +
  coord_sf(xlim = c(-2303077, 1900923), ylim = c(-2305719,1149281)) +
  scale_fill_gradient2() +
  labs(fill = "Mean Residual") +
  theme_void()

cp_offset <- cowplot::plot_grid(off_r_hist,off_qq, ncol = 1,
                                labels = c("C", "D"))

cp_offset

## Duration TIME

# residuals
resids_dur <- residuals(pglmm_duration_sp)
hist(resids_dur)
qqnorm(resids_dur)
qqline(resids_dur)
shapiro.test(resids_dur)

dur_r_hist <- ggplot() + 
  geom_histogram(mapping = aes(x = resids_dur)) +
  theme_minimal() + 
  labs(x = "Residuals") +
  scale_y_continuous(expand = c(0,0)) +
  ggtitle(label = "Duration")

dur_qq <- gg_qq(x = resids_dur) 

# spatial autocorrelation
resids_df <- mdf_phylo_spp %>% 
  mutate(resids = resids_dur)
resids_df <-  resids_df %>% 
  group_by(lon, lat) %>% 
  summarise(mean_resid = mean(resids))

dur_resid_map <- ggplot() + 
  geom_sf(na_equalArea, mapping = aes(), fill = "grey97") + 
  geom_sf(great_lakes, mapping = aes(), fill = "paleturquoise3", alpha = 0.4) +
  geom_tile(resids_df, mapping = aes(x = lon, y = lat, fill = mean_resid)) +
  coord_sf(xlim = c(-2303077, 1900923), ylim = c(-2305719,1149281)) +
  scale_fill_gradient2() +
  labs(fill = "Mean Residual") +
  theme_void()

cp_duration <- cowplot::plot_grid(dur_r_hist,dur_qq, ncol = 1,
                                  labels = c("E", "F"))

cp_duration


## TOtal cp

cp_total <- cowplot::plot_grid(cp_onset, cp_offset, cp_duration, ncol = 3, nrow = 1)

cp_total

ggsave(plot = cp_total, filename = "Figures/resubmission/model_assumptions.png", width = 8, height = 5)

## Total maps
map_total <- cowplot::plot_grid(on_resid_map, off_resid_map, dur_resid_map,
                                labels = c("A", "B", "C"),
                                ncol = 1, nrow = 3)
                                
ggsave(map_total, filename = "Figures/resubmission/spatialResiduals.png",
                                width = 4.5, height = 8)
