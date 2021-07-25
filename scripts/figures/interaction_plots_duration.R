library(tidyverse)
library(gridExtra)
library(phyr)
library(INLA)

# read in phenology data
mdf <- read.csv("data/LMM_data/phenoTraits_data_withNumObsData_updatedSeasTrat.csv", stringsAsFactors = FALSE) %>% 
  filter(scientificName != "Diabrotica undecimpunctata") %>% 
  mutate(immature.habitat = case_when(immature.habitat == "above ground" ~ "Above Ground",
                                      immature.habitat == "fresh water" ~ "Freshwater",
                                      immature.habitat == "underground" ~ "Underground")) %>% 
  mutate(larval.diet = case_when(larval.diet == "carnivorous" ~ "Carnivorous",
                                      larval.diet == "detritivorous" ~ "Detritivorous",
                                      larval.diet == "live plants" ~ "Herbivorous"
                                      ))


# read in phylogeny
insect_tree3 = ape::read.tree("data/phylogeny/insect_test2.tre")
insect_tree3$tip.label <- word(insect_tree3$tip.label, start = 1, end = 2, sep = "_") %>% 
  sub(pattern = "_", replacement = " ")

# remove eusocial species from dataset
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

# remove species that aren't in phylogeny
mdf_phylo_spp <- mdf %>% 
  filter(scientificName %in% tree_sp)

# remove species that do not diapause
insect_tree4 <- ape::drop.tip(phy = insect_tree3, tip = eusocial) %>% 
  ape::drop.tip(tip = "Aeshna cyanea")  %>% 
  ape::drop.tip(tip = "Diabrotica undecimpunctata")

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

### Function to plot Figure fig5_c
# parameters are:
# model_df = dataframe consisting of data used to fit model
# pred_vals = the values that you want to predict over on the x - axis of your figure

plot_pglmm_inter_5a <- function(model_df, pred_vals){
  
  mdf <- model_df

  # here I create a new data frame based on the all of the columns of the original
  # model_df. NA's are used for any columns that you want to hold constant in the
  # predicted data frame and mean values are used for any predictor variables
  # (because we standardized variables, mean values are 0). The column that you
  # want to be the x-axis of your plot should be assigned to pred_vals.
  # The only value that will chance between these new dataframes and following 
  # new dataframes is the value you want to change of the interactive effect.
  # In this case, that is prec column in this case. For now, we are going to get the
  # predicted values for when prec = 0 (mean), but following dataframes will change
  # the prec values
  new_inter_term1 <- data.frame(X =  rep(NA, 7), onset = rep(NA, 7), onset_low = rep(NA, 7), onset_high = rep(NA,7),
                                year = rep(NA, 7), offset = rep(NA, 7), offset_low = rep(NA, 7), offset_high = rep(NA, 7),
                                duration =  rep(NA, 7), scientificName =  rep(NA, 7), Order =  rep(NA, 7), id_cells =  rep(NA, 7),
                                lon =  rep(NA, 7), lat = rep(NA, 7), OS= rep(NA, 7), 
                                prec = rep(0,7), temp = pred_vals, bio4 = rep(0, 7), bio15 = rep(0, 7), pop = rep(0, 7),
                                numObs = rep(0, 7), prec_seas = rep(0, 7), temp_seas = rep(0, 7), 
                                higher_taxon= rep(NA, 7), flights = rep(NA, 7), development= rep(NA, 7),
                                immature.habitat =  rep(NA, 7), diapause.stage= rep(NA, 7), 
                                larval.diet= rep(NA, 7), seas= rep(NA, 7), elevation = rep(NA, 7),
                                lon_cor= rep(NA, 7), lat_cor= rep(NA, 7), elev_cor= rep(NA, 7),
                                hopkins_cor= rep(NA, 7), hopkins_onset= rep(NA, 7), mean_onset= rep(NA, 7),
                                sd_onset= rep(NA, 7), seas2= rep(NA, 7))
  
  pred.df <- rbind(mdf, new_inter_term1) # combine new data frame for predicting values 
                                         # to original dataframe
  
  # repeat above method but change for high value of interactive effect.
  # In this case that is prec = 1 (1 s.d.)
  new_inter_term.high <- data.frame(X =  rep(NA, 7), onset = rep(NA, 7), onset_low = rep(NA, 7), onset_high = rep(NA,7),
                                    year = rep(NA, 7), offset = rep(NA, 7), offset_low = rep(NA, 7), offset_high = rep(NA, 7),
                                    duration =  rep(NA, 7), scientificName =  rep(NA, 7), Order =  rep(NA, 7), id_cells =  rep(NA, 7),
                                    lon =  rep(NA, 7), lat = rep(NA, 7), OS= rep(NA, 7), 
                                    prec = rep(1,7), temp = pred_vals, bio4 = rep(0, 7), bio15 = rep(0, 7), pop = rep(0, 7),
                                    numObs = rep(0, 7), prec_seas = rep(0, 7), temp_seas = rep(0, 7), 
                                    higher_taxon= rep(NA, 7), flights = rep(NA, 7), development= rep(NA, 7),
                                    immature.habitat =  rep(NA, 7), diapause.stage= rep(NA, 7), 
                                    larval.diet= rep(NA, 7), seas= rep(NA, 7), elevation = rep(NA, 7),
                                    lon_cor= rep(NA, 7), lat_cor= rep(NA, 7), elev_cor= rep(NA, 7),
                                    hopkins_cor= rep(NA, 7), hopkins_onset= rep(NA, 7), mean_onset= rep(NA, 7),
                                    sd_onset= rep(NA, 7), seas2= rep(NA, 7))
  
  pred.df.high <- rbind(mdf, new_inter_term.high)
  
  # repeat above method but change for low value of interactive effect.
  # In this case that is prec = -1 (1 s.d.)
  new_inter_term.low <- data.frame(X =  rep(NA, 7), onset = rep(NA, 7), onset_low = rep(NA, 7), onset_high = rep(NA,7),
                                   year = rep(NA, 7), offset = rep(NA, 7), offset_low = rep(NA, 7), offset_high = rep(NA, 7),
                                   duration =  rep(NA, 7), scientificName =  rep(NA, 7), Order =  rep(NA, 7), id_cells =  rep(NA, 7),
                                   lon =  rep(NA, 7), lat = rep(NA, 7), OS= rep(NA, 7), 
                                   prec = rep(-1,7), temp = pred_vals, bio4 = rep(0, 7), bio15 = rep(0, 7), pop = rep(0, 7),
                                   numObs = rep(0, 7), prec_seas = rep(0, 7), temp_seas = rep(0, 7), 
                                   higher_taxon= rep(NA, 7), flights = rep(NA, 7), development= rep(NA, 7),
                                   immature.habitat =  rep(NA, 7), diapause.stage= rep(NA, 7), 
                                   larval.diet= rep(NA, 7), seas= rep(NA, 7), elevation = rep(NA, 7),
                                   lon_cor= rep(NA, 7), lat_cor= rep(NA, 7), elev_cor= rep(NA, 7),
                                   hopkins_cor= rep(NA, 7), hopkins_onset= rep(NA, 7), mean_onset= rep(NA, 7),
                                   sd_onset= rep(NA, 7), seas2= rep(NA, 7))
  
  pred.df.low <- rbind(mdf, new_inter_term.low)
  
  # fit model with combined data frame that has mean value for interactive effect that
  # you are changing in your legend (Prec in this example). 
  m <-  pglmm(duration ~ temp + prec + temp_seas +
                diapause.stage + 
                immature.habitat +
                larval.diet + 
                numObs + 
                (1 | id_cells__) + (1 | scientificName__) + 
                (0 + temp | scientificName__) + 
                (0 + prec | scientificName__) + 
                temp:prec + temp:larval.diet + temp:immature.habitat,
              data = pred.df, 
              cov_ranef = list(scientificName = insect_tree4,
                               id_cells = V.space), 
              bayes = TRUE)
  
  # fit model with data frame that has high value for interactive effect that
  # you are changing in your legend (Prec in this example). 
  m.high <- pglmm(duration ~ temp + prec + temp_seas +
                    diapause.stage + 
                    immature.habitat +
                    larval.diet + 
                    numObs + 
                    (1 | id_cells__) + (1 | scientificName__) + 
                    (0 + temp | scientificName__) + 
                    (0 + prec | scientificName__) + 
                    temp:prec + temp:larval.diet + temp:immature.habitat,
                  data = pred.df.high, 
                  cov_ranef = list(scientificName = insect_tree4,
                                   id_cells = V.space), 
                  bayes = TRUE)
  
  # fit model with data frame that has low value for interactive effect that
  # you are changing in your legend (Prec in this example). 
  m.low <- pglmm(duration ~ temp + prec + temp_seas +
                   diapause.stage + 
                   immature.habitat +
                   larval.diet + 
                   numObs + 
                   (1 | id_cells__) + (1 | scientificName__) + 
                   (0 + temp | scientificName__) + 
                   (0 + prec | scientificName__) + 
                   temp:prec + temp:larval.diet + temp:immature.habitat,
                 data = pred.df.low, 
                 cov_ranef = list(scientificName = insect_tree4,
                                  id_cells = V.space), 
                 bayes = TRUE)
  
  # pull out model predictions for the last number of rows in your dataframe
  # that are predicted not fit values.
  # Here I also add the characters that will be used in the Figure legend and 
  # name them inter_eff2
  a <- m$inla.model$summary.linear.predictor[2644:2650,] %>% 
    mutate(inter_eff1 = pred_vals,
           inter_eff2 = "Mid",
           inter_eff2_vals = 0)
  a.low <- m.low$inla.model$summary.linear.predictor[2644:2650,] %>% 
    mutate(inter_eff1 = pred_vals,
           inter_eff2 = "Low",
           inter_eff2_vals = -1)
  a.high <- m.high$inla.model$summary.linear.predictor[2644:2650,] %>% 
    mutate(inter_eff1 = pred_vals,
           inter_eff2 = "High",
           inter_eff2_vals = 1)
  
  a.total <- rbind(a, a.low, a.high) # combine predicted results
  a.total$inter_eff2 <- factor(a.total$inter_eff2, levels = c("Low", "Mid", "High"))
  
  mdf_phylo_spp <- mdf_phylo_spp %>% 
    mutate(inter_eff2 = case_when(prec >=1 ~ "High",
                                  prec <= -1 ~ "Low",
                                  prec > -1 & prec <1 ~ "Mid"))
  mdf_phylo_spp$inter_eff2 <- factor(mdf_phylo_spp$inter_eff2, levels = c("Low", "Mid", "High"))
  
  
  # plot results 
  p <-  ggplot(a.total, mapping = aes(x = inter_eff1, y = mean)) +
    geom_point(mdf_phylo_spp, mapping = aes(x = temp, y = duration, color = inter_eff2), # adds original datapoints
               alpha = 0.08) + 
    geom_ribbon(mapping = aes(ymin = `0.025quant`, ymax = `0.975quant`, # add credible intervals
                              fill = inter_eff2), 
                alpha = 0.15) +
    geom_line(mapping = aes(color = inter_eff2), size = 1.05) +
    labs(x = "Temperature", y = "Duration", fill = "Precipitation", color = "Precipitation") + 
    scale_color_manual(values = c("brown", "cyan", "Blue")) +
    scale_fill_manual(values = c("brown", "cyan", "Blue")) +
    guides(fill = FALSE) +
    theme_bw()
  
  return(p)
  
}

# plot onset interactions
fig5_a <- plot_pglmm_inter_5a(model_df = mdf_phylo_spp, pred_vals = -3:3)
fig5_a

#5b temperatuer:immature.habitat
# function to plot credible intervals of model fit for Figure 5b
# note that I change the 7 to be length of pred_vals -- this ensures 
# function won't error out because of mismatched lengths.
# but you'll still have to manually change many lines of code that tailor the
# figure to the data and specific interaction that you are looking at. 

# This example is plotting an intercation between temperature (x-axis) and
# categorical data (immature habitat).
plot_pglmm_inter_5b <- function(model_df, pred_vals){
  
  mdf <- model_df
  
  new_inter_term_ag <- data.frame(X =  rep(NA, length(pred_vals)), onset = rep(NA, length(pred_vals)), onset_low = rep(NA, length(pred_vals)), onset_high = rep(NA,length(pred_vals)),
                                  year = rep(NA, length(pred_vals)), offset = rep(NA, length(pred_vals)), offset_low = rep(NA, length(pred_vals)), offset_high = rep(NA, length(pred_vals)),
                                  duration =  rep(NA, length(pred_vals)), scientificName =  rep(NA, length(pred_vals)), Order =  rep(NA, length(pred_vals)), id_cells =  rep(NA, length(pred_vals)),
                                  lon =  rep(NA, length(pred_vals)), lat = rep(NA, length(pred_vals)), OS= rep(NA, length(pred_vals)), 
                                  prec = rep(0,length(pred_vals)), temp = pred_vals, bio4 = rep(0, length(pred_vals)), bio15 = rep(0, length(pred_vals)), pop = rep(0, length(pred_vals)),
                                  numObs = rep(0, length(pred_vals)), prec_seas = rep(0, length(pred_vals)), temp_seas = rep(0, length(pred_vals)), 
                                  higher_taxon= rep(NA, length(pred_vals)), flights = rep(NA, length(pred_vals)), development= rep(NA, length(pred_vals)),
                                  immature.habitat =  rep("Above Ground"), diapause.stage= rep(NA, length(pred_vals)), 
                                  larval.diet= rep(NA, length(pred_vals)), seas= rep(NA, length(pred_vals)), elevation = rep(NA, length(pred_vals)),
                                  lon_cor= rep(NA, length(pred_vals)), lat_cor= rep(NA, length(pred_vals)), elev_cor= rep(NA, length(pred_vals)),
                                  hopkins_cor= rep(NA, length(pred_vals)), hopkins_onset= rep(NA, length(pred_vals)), mean_onset= rep(NA, length(pred_vals)),
                                  sd_onset= rep(NA, length(pred_vals)), seas2= rep(NA, length(pred_vals)))
  
  pred.df.ag <- rbind(mdf, new_inter_term_ag)
  
  new_inter_term.fw <- data.frame(X =  rep(NA, length(pred_vals)), onset = rep(NA, length(pred_vals)), onset_low = rep(NA, length(pred_vals)), onset_high = rep(NA,length(pred_vals)),
                                  year = rep(NA, length(pred_vals)), offset = rep(NA, length(pred_vals)), offset_low = rep(NA, length(pred_vals)), offset_high = rep(NA, length(pred_vals)),
                                  duration =  rep(NA, length(pred_vals)), scientificName =  rep(NA, length(pred_vals)), Order =  rep(NA, length(pred_vals)), id_cells =  rep(NA, length(pred_vals)),
                                  lon =  rep(NA, length(pred_vals)), lat = rep(NA, length(pred_vals)), OS= rep(NA, length(pred_vals)), 
                                  prec = rep(0,length(pred_vals)), temp = pred_vals, bio4 = rep(0, length(pred_vals)), bio15 = rep(0, length(pred_vals)), pop = rep(0, length(pred_vals)),
                                  numObs = rep(0, length(pred_vals)), prec_seas = rep(0, length(pred_vals)), temp_seas = rep(0, length(pred_vals)), 
                                  higher_taxon= rep(NA, length(pred_vals)), flights = rep(NA, length(pred_vals)), development= rep(NA, length(pred_vals)),
                                  immature.habitat =  rep("Freshwater", length(pred_vals)), diapause.stage= rep(NA, length(pred_vals)), 
                                  larval.diet= rep(NA, length(pred_vals)), seas= rep(NA, length(pred_vals)), elevation = rep(NA, length(pred_vals)),
                                  lon_cor= rep(NA, length(pred_vals)), lat_cor= rep(NA, length(pred_vals)), elev_cor= rep(NA, length(pred_vals)),
                                  hopkins_cor= rep(NA, length(pred_vals)), hopkins_onset= rep(NA, length(pred_vals)), mean_onset= rep(NA, length(pred_vals)),
                                  sd_onset= rep(NA, length(pred_vals)), seas2= rep(NA, length(pred_vals)))
  
  pred.df.fw <- rbind(mdf, new_inter_term.fw)
  
  new_inter_term.ug <- data.frame(X =  rep(NA, length(pred_vals)), onset = rep(NA, length(pred_vals)), onset_low = rep(NA, length(pred_vals)), onset_high = rep(NA,length(pred_vals)),
                                  year = rep(NA, length(pred_vals)), offset = rep(NA, length(pred_vals)), offset_low = rep(NA, length(pred_vals)), offset_high = rep(NA, length(pred_vals)),
                                  duration =  rep(NA, length(pred_vals)), scientificName =  rep(NA, length(pred_vals)), Order =  rep(NA, length(pred_vals)), id_cells =  rep(NA, length(pred_vals)),
                                  lon =  rep(NA, length(pred_vals)), lat = rep(NA, length(pred_vals)), OS= rep(NA, length(pred_vals)), 
                                  prec = rep(0,length(pred_vals)), temp = pred_vals, bio4 = rep(0, length(pred_vals)), bio15 = rep(0, length(pred_vals)), pop = rep(0, length(pred_vals)),
                                  numObs = rep(0, length(pred_vals)), prec_seas = rep(0, length(pred_vals)), temp_seas = rep(0, length(pred_vals)), 
                                  higher_taxon= rep(NA, length(pred_vals)), flights = rep(NA, length(pred_vals)), development= rep(NA, length(pred_vals)),
                                  immature.habitat =  rep("Underground", length(pred_vals)), diapause.stage= rep(NA, length(pred_vals)), 
                                  larval.diet= rep(NA, length(pred_vals)), seas= rep(NA, length(pred_vals)), elevation = rep(NA, length(pred_vals)),
                                  lon_cor= rep(NA, length(pred_vals)), lat_cor= rep(NA, length(pred_vals)), elev_cor= rep(NA, length(pred_vals)),
                                  hopkins_cor= rep(NA, length(pred_vals)), hopkins_onset= rep(NA, length(pred_vals)), mean_onset= rep(NA, length(pred_vals)),
                                  sd_onset= rep(NA, length(pred_vals)), seas2= rep(NA, length(pred_vals)))
  
  pred.df.ug <- rbind(mdf, new_inter_term.ug)
  
  # fit models
  m.ag <-  pglmm(duration ~ temp + prec + temp_seas +
                   diapause.stage + 
                   immature.habitat +
                   larval.diet + 
                   numObs + 
                   (1 | id_cells__) + (1 | scientificName__) + 
                   (0 + temp | scientificName__) + 
                   (0 + prec | scientificName__) + 
                   temp:prec + temp:larval.diet + temp:immature.habitat,
                 data = pred.df.ag, 
                 cov_ranef = list(scientificName = insect_tree4,
                                  id_cells = V.space), 
                 bayes = TRUE)
  
  m.fw <-  pglmm(duration ~ temp + prec + temp_seas +
                   diapause.stage + 
                   immature.habitat +
                   larval.diet + 
                   numObs + 
                   (1 | id_cells__) + (1 | scientificName__) + 
                   (0 + temp | scientificName__) + 
                   (0 + prec | scientificName__) + 
                   temp:prec + temp:larval.diet + temp:immature.habitat,
                 data = pred.df.fw, 
                 cov_ranef = list(scientificName = insect_tree4,
                                  id_cells = V.space), 
                 bayes = TRUE)
  
  m.ug <-  pglmm(duration ~ temp + prec + temp_seas +
                   diapause.stage + 
                   immature.habitat +
                   larval.diet + 
                   numObs + 
                   (1 | id_cells__) + (1 | scientificName__) + 
                   (0 + temp | scientificName__) + 
                   (0 + prec | scientificName__) + 
                   temp:prec + temp:larval.diet + temp:immature.habitat,
                 data = pred.df.ug, 
                 cov_ranef = list(scientificName = insect_tree4,
                                  id_cells = V.space), 
                 bayes = TRUE)
  
  # make dataframe of predicted results 
  a.ag <- m.ag$inla.model$summary.linear.predictor[2644:2650,] %>% 
    mutate(inter_eff1 = pred_vals,
           inter_eff2 = "Above Ground",
           inter_eff2_vals = 0)
  a.fw <- m.fw$inla.model$summary.linear.predictor[2644:2650,] %>% 
    mutate(inter_eff1 = pred_vals,
           inter_eff2 = "Freshwater",
           inter_eff2_vals = -1)
  a.ug <- m.ug$inla.model$summary.linear.predictor[2644:2650,] %>% 
    mutate(inter_eff1 = pred_vals,
           inter_eff2 = "Underground",
           inter_eff2_vals = 1)
  
  a.total <- rbind(a.ag, a.fw, a.ug)
  
  p <-  ggplot(a.total, mapping = aes(x = inter_eff1, y = mean)) +
   geom_point(mdf_phylo_spp, mapping = aes(x = temp, y = duration, color = immature.habitat),
               alpha = 0.08) + 
    geom_ribbon(mapping = aes(ymin = `0.025quant`, ymax = `0.975quant`,
                              fill = inter_eff2), 
                alpha = 0.15) +
    geom_line(mapping = aes(color = inter_eff2), size = 1.05) +
    labs(x = "Temperature", y = "Duration", fill = "Immature Habitat", color = "Immature Habitat") + 
    scale_color_manual(values = c("olivedrab4", "deepskyblue2", "orangered4")) +
    scale_fill_manual(values = c("olivedrab4", "deepskyblue2", "orangered4")) +
    guides(fill = FALSE) +
    theme_bw()
  
  return(p)
  
}

## run function 
fig5b <- plot_pglmm_inter_5b(model_df = mdf_phylo_spp, pred_vals = -3:3)
fig5b 


#5c temperature:larval.diet
# function to plot credible intervals of model fit for Figure 5c which is 
# examining the interactive effects between temperature and larval diet. 

plot_pglmm_inter_5c <- function(model_df, pred_vals){
  
  mdf <- model_df
  
  new_inter_term_ag <- data.frame(X =  rep(NA, length(pred_vals)), onset = rep(NA, length(pred_vals)), onset_low = rep(NA, length(pred_vals)), onset_high = rep(NA,length(pred_vals)),
                                  year = rep(NA, length(pred_vals)), offset = rep(NA, length(pred_vals)), offset_low = rep(NA, length(pred_vals)), offset_high = rep(NA, length(pred_vals)),
                                  duration =  rep(NA, length(pred_vals)), scientificName =  rep(NA, length(pred_vals)), Order =  rep(NA, length(pred_vals)), id_cells =  rep(NA, length(pred_vals)),
                                  lon =  rep(NA, length(pred_vals)), lat = rep(NA, length(pred_vals)), OS= rep(NA, length(pred_vals)), 
                                  prec = rep(0,length(pred_vals)), temp = pred_vals, bio4 = rep(0, length(pred_vals)), bio15 = rep(0, length(pred_vals)), pop = rep(0, length(pred_vals)),
                                  numObs = rep(0, length(pred_vals)), prec_seas = rep(0, length(pred_vals)), temp_seas = rep(0, length(pred_vals)), 
                                  higher_taxon= rep(NA, length(pred_vals)), flights = rep(NA, length(pred_vals)), development= rep(NA, length(pred_vals)),
                                  immature.habitat =  rep(NA, length(pred_vals)), diapause.stage= rep(NA, length(pred_vals)), 
                                  larval.diet= rep("Carnivorous", length(pred_vals)), seas= rep(NA, length(pred_vals)), elevation = rep(NA, length(pred_vals)),
                                  lon_cor= rep(NA, length(pred_vals)), lat_cor= rep(NA, length(pred_vals)), elev_cor= rep(NA, length(pred_vals)),
                                  hopkins_cor= rep(NA, length(pred_vals)), hopkins_onset= rep(NA, length(pred_vals)), mean_onset= rep(NA, length(pred_vals)),
                                  sd_onset= rep(NA, length(pred_vals)), seas2= rep(NA, length(pred_vals)))
  
  pred.df.ag <- rbind(mdf, new_inter_term_ag)
  
  new_inter_term.fw <- data.frame(X =  rep(NA, length(pred_vals)), onset = rep(NA, length(pred_vals)), onset_low = rep(NA, length(pred_vals)), onset_high = rep(NA,length(pred_vals)),
                                  year = rep(NA, length(pred_vals)), offset = rep(NA, length(pred_vals)), offset_low = rep(NA, length(pred_vals)), offset_high = rep(NA, length(pred_vals)),
                                  duration =  rep(NA, length(pred_vals)), scientificName =  rep(NA, length(pred_vals)), Order =  rep(NA, length(pred_vals)), id_cells =  rep(NA, length(pred_vals)),
                                  lon =  rep(NA, length(pred_vals)), lat = rep(NA, length(pred_vals)), OS= rep(NA, length(pred_vals)), 
                                  prec = rep(0,length(pred_vals)), temp = pred_vals, bio4 = rep(0, length(pred_vals)), bio15 = rep(0, length(pred_vals)), pop = rep(0, length(pred_vals)),
                                  numObs = rep(0, length(pred_vals)), prec_seas = rep(0, length(pred_vals)), temp_seas = rep(0, length(pred_vals)), 
                                  higher_taxon= rep(NA, length(pred_vals)), flights = rep(NA, length(pred_vals)), development= rep(NA, length(pred_vals)),
                                  immature.habitat =  rep(NA, length(pred_vals)), diapause.stage= rep(NA, length(pred_vals)), 
                                  larval.diet= rep("Detritivorous", length(pred_vals)), seas= rep(NA, length(pred_vals)), elevation = rep(NA, length(pred_vals)),
                                  lon_cor= rep(NA, length(pred_vals)), lat_cor= rep(NA, length(pred_vals)), elev_cor= rep(NA, length(pred_vals)),
                                  hopkins_cor= rep(NA, length(pred_vals)), hopkins_onset= rep(NA, length(pred_vals)), mean_onset= rep(NA, length(pred_vals)),
                                  sd_onset= rep(NA, length(pred_vals)), seas2= rep(NA, length(pred_vals)))
  
  pred.df.fw <- rbind(mdf, new_inter_term.fw)
  
  new_inter_term.ug <- data.frame(X =  rep(NA, length(pred_vals)), onset = rep(NA, length(pred_vals)), onset_low = rep(NA, length(pred_vals)), onset_high = rep(NA,length(pred_vals)),
                                  year = rep(NA, length(pred_vals)), offset = rep(NA, length(pred_vals)), offset_low = rep(NA, length(pred_vals)), offset_high = rep(NA, length(pred_vals)),
                                  duration =  rep(NA, length(pred_vals)), scientificName =  rep(NA, length(pred_vals)), Order =  rep(NA, length(pred_vals)), id_cells =  rep(NA, length(pred_vals)),
                                  lon =  rep(NA, length(pred_vals)), lat = rep(NA, length(pred_vals)), OS= rep(NA, length(pred_vals)), 
                                  prec = rep(0,length(pred_vals)), temp = pred_vals, bio4 = rep(0, length(pred_vals)), bio15 = rep(0, length(pred_vals)), pop = rep(0, length(pred_vals)),
                                  numObs = rep(0, length(pred_vals)), prec_seas = rep(0, length(pred_vals)), temp_seas = rep(0, length(pred_vals)), 
                                  higher_taxon= rep(NA, length(pred_vals)), flights = rep(NA, length(pred_vals)), development= rep(NA, length(pred_vals)),
                                  immature.habitat =  rep(NA, length(pred_vals)), diapause.stage= rep(NA, length(pred_vals)), 
                                  larval.diet= rep("Herbivorous", length(pred_vals)), seas= rep(NA, length(pred_vals)), elevation = rep(NA, length(pred_vals)),
                                  lon_cor= rep(NA, length(pred_vals)), lat_cor= rep(NA, length(pred_vals)), elev_cor= rep(NA, length(pred_vals)),
                                  hopkins_cor= rep(NA, length(pred_vals)), hopkins_onset= rep(NA, length(pred_vals)), mean_onset= rep(NA, length(pred_vals)),
                                  sd_onset= rep(NA, length(pred_vals)), seas2= rep(NA, length(pred_vals)))
  
  pred.df.ug <- rbind(mdf, new_inter_term.ug)
  
  # fit models
  m.ag <-  pglmm(duration ~ temp + prec + temp_seas +
                   diapause.stage + 
                   immature.habitat +
                   larval.diet + 
                   numObs + 
                   (1 | id_cells__) + (1 | scientificName__) + 
                   (0 + temp | scientificName__) + 
                   (0 + prec | scientificName__) + 
                   temp:prec + temp:larval.diet + temp:immature.habitat,
                 data = pred.df.ag, 
                 cov_ranef = list(scientificName = insect_tree4,
                                  id_cells = V.space), 
                 bayes = TRUE)
  
  m.fw <-  pglmm(duration ~ temp + prec + temp_seas +
                   diapause.stage + 
                   immature.habitat +
                   larval.diet + 
                   numObs + 
                   (1 | id_cells__) + (1 | scientificName__) + 
                   (0 + temp | scientificName__) + 
                   (0 + prec | scientificName__) + 
                   temp:prec + temp:larval.diet + temp:immature.habitat,
                 data = pred.df.fw, 
                 cov_ranef = list(scientificName = insect_tree4,
                                  id_cells = V.space), 
                 bayes = TRUE)
  
  m.ug <-  pglmm(duration ~ temp + prec + temp_seas +
                   diapause.stage + 
                   immature.habitat +
                   larval.diet + 
                   numObs + 
                   (1 | id_cells__) + (1 | scientificName__) + 
                   (0 + temp | scientificName__) + 
                   (0 + prec | scientificName__) + 
                   temp:prec + temp:larval.diet + temp:immature.habitat,
                 data = pred.df.ug, 
                 cov_ranef = list(scientificName = insect_tree4,
                                  id_cells = V.space), 
                 bayes = TRUE)
  
  # make dataframe of predicted results 
  a.ag <- m.ag$inla.model$summary.linear.predictor[2644:2650,] %>% 
    mutate(inter_eff1 = pred_vals,
           inter_eff2 = "Carnivorous",
           inter_eff2_vals = 0)
  a.fw <- m.fw$inla.model$summary.linear.predictor[2644:2650,] %>% 
    mutate(inter_eff1 = pred_vals,
           inter_eff2 = "Detritivorous",
           inter_eff2_vals = -1)
  a.ug <- m.ug$inla.model$summary.linear.predictor[2644:2650,] %>% 
    mutate(inter_eff1 = pred_vals,
           inter_eff2 = "Herbivorous",
           inter_eff2_vals = 1)
  
  a.total <- rbind(a.ag, a.fw, a.ug)
  
  p <-  ggplot(a.total, mapping = aes(x = inter_eff1, y = mean)) +
    geom_point(mdf_phylo_spp, mapping = aes(x = temp, y = duration, color = larval.diet),
               alpha = 0.08) + 
    geom_ribbon(mapping = aes(ymin = `0.025quant`, ymax = `0.975quant`,
                              fill = inter_eff2), 
                alpha = 0.15) +
    geom_line(mapping = aes(color = inter_eff2), size = 1.05) +
    labs(x = "Temperature", y = "Duration", fill = "Larval Diet", color = "Larval Diet") + 
    scale_color_manual(values = c("#d1495b", "#edae49",  "#66a182")) +
    scale_fill_manual(values = c("#d1495b", "#edae49",  "#66a182")) +
    guides(fill = FALSE) +
    theme_bw()
  
  return(p)
  
}

## Something 
fig5c <- plot_pglmm_inter_5c(model_df = mdf_phylo_spp, pred_vals = -3:3)
fig5c

# Combine and Save Duration plots
duration_interaction_plots <- egg::ggarrange(fig5_a, fig5b, fig5c, 
                                          ncol = 2, labels = c("A", "B", 'C'))

ggsave(duration_interaction_plots, filename = "Figures/resubmission/duration_interactions.png",
       width = 8, height = 4)
