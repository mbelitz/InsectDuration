library(tidyverse)
library(gridExtra)
library(phyr)
library(INLA)

# read in data
mdf <- read.csv("data/LMM_data/phenoTraits_data_withNumObsData_updatedSeasTrat.csv", stringsAsFactors = FALSE) %>% 
  filter(scientificName != "Diabrotica undecimpunctata")  %>% 
  mutate(immature.habitat = case_when(immature.habitat == "above ground" ~ "Above Ground",
                             immature.habitat == "fresh water" ~ "Freshwater",
                             immature.habitat == "underground" ~ "Underground"))

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


### fig3_c
plot_pglmm_inter_4 <- function(model_df, pred_vals){
  
  mdf <- model_df
  n <- names(mdf)
  
  new_inter_term_ag <- data.frame(X =  rep(NA, length(pred_vals)), onset = rep(NA, length(pred_vals)), onset_low = rep(NA, length(pred_vals)), onset_high = rep(NA,length(pred_vals)),
                                     year = rep(NA, length(pred_vals)), offset = rep(NA, length(pred_vals)), offset_low = rep(NA, length(pred_vals)), offset_high = rep(NA, length(pred_vals)),
                                     duration =  rep(NA, length(pred_vals)), scientificName =  rep(NA, length(pred_vals)), Order =  rep(NA, length(pred_vals)), id_cells =  rep(NA, length(pred_vals)),
                                     lon =  rep(NA, length(pred_vals)), lat = rep(NA, length(pred_vals)), OS= rep(NA, length(pred_vals)), 
                                     prec = rep(0,length(pred_vals)), temp = rep(0, length(pred_vals)), bio4 = rep(0, length(pred_vals)), bio15 = rep(0, length(pred_vals)), pop = rep(0, length(pred_vals)),
                                     numObs = rep(0, length(pred_vals)), prec_seas = rep(0, length(pred_vals)), temp_seas = pred_vals, 
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
                                   prec = rep(0,length(pred_vals)), temp = rep(0,length(pred_vals)), bio4 = rep(0, length(pred_vals)), bio15 = rep(0, length(pred_vals)), pop = rep(0, length(pred_vals)),
                                   numObs = rep(0, length(pred_vals)), prec_seas = rep(0, length(pred_vals)), temp_seas = pred_vals, 
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
                                     prec = rep(0,length(pred_vals)), temp = rep(0,length(pred_vals)), bio4 = rep(0, length(pred_vals)), bio15 = rep(0, length(pred_vals)), pop = rep(0, length(pred_vals)),
                                     numObs = rep(0, length(pred_vals)), prec_seas = rep(0, length(pred_vals)), temp_seas = pred_vals, 
                                     higher_taxon= rep(NA, length(pred_vals)), flights = rep(NA, length(pred_vals)), development= rep(NA, length(pred_vals)),
                                     immature.habitat =  rep("Underground", length(pred_vals)), diapause.stage= rep(NA, length(pred_vals)), 
                                     larval.diet= rep(NA, length(pred_vals)), seas= rep(NA, length(pred_vals)), elevation = rep(NA, length(pred_vals)),
                                     lon_cor= rep(NA, length(pred_vals)), lat_cor= rep(NA, length(pred_vals)), elev_cor= rep(NA, length(pred_vals)),
                                     hopkins_cor= rep(NA, length(pred_vals)), hopkins_onset= rep(NA, length(pred_vals)), mean_onset= rep(NA, length(pred_vals)),
                                     sd_onset= rep(NA, length(pred_vals)), seas2= rep(NA, length(pred_vals)))
  
  pred.df.ug <- rbind(mdf, new_inter_term.ug)
  
  # fit models
  m.ag <- pglmm(offset ~ temp_seas + prec +
                  seas2 + diapause.stage + immature.habitat + larval.diet + 
                  temp_seas:immature.habitat +
                  numObs + 
                  (1 | id_cells__) + (1 | scientificName__) + 
                  (0 + temp_seas | scientificName__) + 
                  (0 + prec | scientificName__), 
                data = pred.df.ag, 
                cov_ranef = list(scientificName = insect_tree4,
                                 id_cells = V.space), 
                bayes = TRUE)
  
  m.fw <- pglmm(offset ~ temp_seas + prec +
                  seas2 + diapause.stage + immature.habitat + larval.diet + 
                  temp_seas:immature.habitat +
                  numObs + 
                  (1 | id_cells__) + (1 | scientificName__) + 
                  (0 + temp_seas | scientificName__) + 
                  (0 + prec | scientificName__), 
                data = pred.df.fw, 
                cov_ranef = list(scientificName = insect_tree4,
                                 id_cells = V.space), 
                bayes = TRUE)
  
  m.ug <- pglmm(offset ~ temp_seas + prec +
                  seas2 + diapause.stage + immature.habitat + larval.diet + 
                  temp_seas:immature.habitat +
                  numObs + 
                  (1 | id_cells__) + (1 | scientificName__) + 
                  (0 + temp_seas | scientificName__) + 
                  (0 + prec | scientificName__), 
                data = pred.df.ug, 
                cov_ranef = list(scientificName = insect_tree4,
                                 id_cells = V.space), 
                bayes = TRUE)
  
  # make dataframe of predicted results 
  a.ag <- m.ag$inla.model$summary.linear.predictor[2644:2649,] %>% 
    mutate(inter_eff1 = pred_vals,
           inter_eff2 = "Above Ground",
           inter_eff2_vals = 0)
  a.fw <- m.fw$inla.model$summary.linear.predictor[2644:2649,] %>% 
    mutate(inter_eff1 = pred_vals,
           inter_eff2 = "Freshwater",
           inter_eff2_vals = -1)
  a.ug <- m.ug$inla.model$summary.linear.predictor[2644:2649,] %>% 
    mutate(inter_eff1 = pred_vals,
           inter_eff2 = "Underground",
           inter_eff2_vals = 1)
  
  a.total <- rbind(a.ag, a.fw, a.ug)
  
  p <-  ggplot(a.total, mapping = aes(x = inter_eff1, y = mean)) +
    geom_point(mdf_phylo_spp, mapping = aes(x = temp_seas, y = offset, color = immature.habitat),
               alpha = 0.08) + 
    geom_ribbon(mapping = aes(ymin = `0.025quant`, ymax = `0.975quant`,
                              fill = inter_eff2), 
                alpha = 0.15) +
    geom_line(mapping = aes(color = inter_eff2), size = 1.05) +
    labs(x = "Temperature Seasonality", y = "Termination", fill = "Immature Habitat", color = "Immature Habitat") + 
    scale_color_manual(values = c("olivedrab4", "deepskyblue2", "orangered4")) +
    scale_fill_manual(values = c("olivedrab4", "deepskyblue2", "orangered4")) +
    guides(fill = FALSE) +
    theme_bw()
  
  return(p)
  
}

## Something 
fig4 <- plot_pglmm_inter_4(model_df = mdf_phylo_spp, pred_vals = -3:2)
fig4 

ggsave(plot = fig4, filename = "Figures/resubmission/offset_interaction.png",
       width = 5, height = 3)

