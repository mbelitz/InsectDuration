library(tidyverse)
library(gridExtra)
library(phyr)
library(INLA)

# read in data
mdf <- read.csv("data/LMM_data/phenoTraits_data_withNumObsData_updatedSeasTrat.csv", stringsAsFactors = FALSE) %>% 
  filter(scientificName != "Diabrotica undecimpunctata") %>% 
  mutate(diapause.stage = if_else(diapause.stage == "Adults", true = "Adult",
                                  false = diapause.stage)) %>% 
  mutate(flights = case_when(flights == "not univoltine" ~ "Not Univoltine",
                             flights == "univoltine" ~ "Univoltine"))

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

#'function to plot PGLMM interactions

plot_pglmm_inter <- function(model_df, pred_vals){
  
  mdf <- model_df
  n <- names(mdf)
  
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
  
  pred.df <- rbind(mdf, new_inter_term1)
  
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
  
  m <- pglmm(onset ~ temp + pop + prec + temp_seas + numObs + 
               diapause.stage + flights +
               (1|id_cells) + (1|scientificName__) +
               (0 + temp | scientificName__) + 
               (0 + prec | scientificName__) +
               (0 + temp_seas | scientificName__) +
               temp:prec + temp:pop + 
               temp:diapause.stage + prec:flights,
             data = pred.df, 
             cov_ranef = list(scientificName = insect_tree4), 
             bayes = TRUE)
  
  m.high <- pglmm(onset ~ temp + pop + prec + temp_seas + numObs + 
                    diapause.stage + flights +
                    (1|id_cells) + (1|scientificName__) +
                    (0 + temp | scientificName__) + 
                    (0 + prec | scientificName__) +
                    (0 + temp_seas | scientificName__) +
                    temp:prec + temp:pop + 
                    temp:diapause.stage + prec:flights,
                  data = pred.df.high, 
                  cov_ranef = list(scientificName = insect_tree4), 
                  bayes = TRUE)
  
  m.low <- pglmm(onset ~ temp + pop + prec + temp_seas + numObs + 
                   diapause.stage + flights +
                   (1|id_cells) + (1|scientificName__) +
                   (0 + temp | scientificName__) + 
                   (0 + prec | scientificName__) +
                   (0 + temp_seas | scientificName__) +
                   temp:prec + temp:pop + 
                   temp:diapause.stage + prec:flights,
                 data = pred.df.low, 
                 cov_ranef = list(scientificName = insect_tree4), 
                 bayes = TRUE)
  
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
  
  a.total <- rbind(a, a.low, a.high)
  a.total$inter_eff2 <- factor(a.total$inter_eff2, levels = c("Low", "Mid", "High"))
  
  mdf_phylo_spp <- mdf_phylo_spp %>% 
    mutate(inter_eff2 = case_when(prec >=1 ~ "High",
                                   prec <= -1 ~ "Low",
                                   prec > -1 & prec <1 ~ "Mid"))
  mdf_phylo_spp$inter_eff2 <- factor(mdf_phylo_spp$inter_eff2, levels = c("Low", "Mid", "High"))
  
  
  p <- ggplot(a.total, mapping = aes(x = inter_eff1, y = mean)) +
    geom_point(mdf_phylo_spp, mapping = aes(x = temp, y = onset, color = inter_eff2),
               alpha = 0.08) + 
    geom_ribbon(mapping = aes(ymin = `0.025quant`, ymax = `0.975quant`,
                              fill = inter_eff2), 
                alpha = 0.15) +
    geom_line(mapping = aes(color = inter_eff2), size = 1.05) +
    labs(x = "Temperature", y = "Emergence", fill = "Precipitation", color = "Precipitation") + 
    scale_color_manual(values = c("brown", "cyan", "Blue")) +
    scale_fill_manual(values = c("brown", "cyan", "Blue")) +
    guides(fill = FALSE) +
    theme_bw()
  
  return(p)
  
}

# plot onset interactions
fig2_b <- plot_pglmm_inter(model_df = mdf_phylo_spp, pred_vals = -3:3)
fig2_b

### try for fig2_a
plot_pglmm_inter_2a <- function(model_df, pred_vals){
  
  mdf <- model_df
  n <- names(mdf)
  
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
  
  pred.df <- rbind(mdf, new_inter_term1)
  
  new_inter_term.high <- data.frame(X =  rep(NA, 7), onset = rep(NA, 7), onset_low = rep(NA, 7), onset_high = rep(NA,7),
                                    year = rep(NA, 7), offset = rep(NA, 7), offset_low = rep(NA, 7), offset_high = rep(NA, 7),
                                    duration =  rep(NA, 7), scientificName =  rep(NA, 7), Order =  rep(NA, 7), id_cells =  rep(NA, 7),
                                    lon =  rep(NA, 7), lat = rep(NA, 7), OS= rep(NA, 7), 
                                    prec = rep(0,7), temp = pred_vals, bio4 = rep(0, 7), bio15 = rep(0, 7), pop = rep(1, 7),
                                    numObs = rep(0, 7), prec_seas = rep(0, 7), temp_seas = rep(0, 7), 
                                    higher_taxon= rep(NA, 7), flights = rep(NA, 7), development= rep(NA, 7),
                                    immature.habitat =  rep(NA, 7), diapause.stage= rep(NA, 7), 
                                    larval.diet= rep(NA, 7), seas= rep(NA, 7), elevation = rep(NA, 7),
                                    lon_cor= rep(NA, 7), lat_cor= rep(NA, 7), elev_cor= rep(NA, 7),
                                    hopkins_cor= rep(NA, 7), hopkins_onset= rep(NA, 7), mean_onset= rep(NA, 7),
                                    sd_onset= rep(NA, 7), seas2= rep(NA, 7))
  
  pred.df.high <- rbind(mdf, new_inter_term.high)
  
  new_inter_term.low <- data.frame(X =  rep(NA, 7), onset = rep(NA, 7), onset_low = rep(NA, 7), onset_high = rep(NA,7),
                                   year = rep(NA, 7), offset = rep(NA, 7), offset_low = rep(NA, 7), offset_high = rep(NA, 7),
                                   duration =  rep(NA, 7), scientificName =  rep(NA, 7), Order =  rep(NA, 7), id_cells =  rep(NA, 7),
                                   lon =  rep(NA, 7), lat = rep(NA, 7), OS= rep(NA, 7), 
                                   prec = rep(0,7), temp = pred_vals, bio4 = rep(0, 7), bio15 = rep(0, 7), pop = rep(-1, 7),
                                   numObs = rep(0, 7), prec_seas = rep(0, 7), temp_seas = rep(0, 7), 
                                   higher_taxon= rep(NA, 7), flights = rep(NA, 7), development= rep(NA, 7),
                                   immature.habitat =  rep(NA, 7), diapause.stage= rep(NA, 7), 
                                   larval.diet= rep(NA, 7), seas= rep(NA, 7), elevation = rep(NA, 7),
                                   lon_cor= rep(NA, 7), lat_cor= rep(NA, 7), elev_cor= rep(NA, 7),
                                   hopkins_cor= rep(NA, 7), hopkins_onset= rep(NA, 7), mean_onset= rep(NA, 7),
                                   sd_onset= rep(NA, 7), seas2= rep(NA, 7))
  
  pred.df.low <- rbind(mdf, new_inter_term.low)
  
  m <- pglmm(onset ~ temp + pop + prec + temp_seas + numObs + 
               diapause.stage + flights +
               (1|id_cells) + (1|scientificName__) +
               (0 + temp | scientificName__) + 
               (0 + prec | scientificName__) +
               (0 + temp_seas | scientificName__) +
               temp:prec + temp:pop + 
               temp:diapause.stage + prec:flights,
             data = pred.df, 
             cov_ranef = list(scientificName = insect_tree4), 
             bayes = TRUE)
  
  m.high <- pglmm(onset ~ temp + pop + prec + temp_seas + numObs + 
                    diapause.stage + flights +
                    (1|id_cells) + (1|scientificName__) +
                    (0 + temp | scientificName__) + 
                    (0 + prec | scientificName__) +
                    (0 + temp_seas | scientificName__) +
                    temp:prec + temp:pop + 
                    temp:diapause.stage + prec:flights,
                  data = pred.df.high, 
                  cov_ranef = list(scientificName = insect_tree4), 
                  bayes = TRUE)
  
  m.low <- pglmm(onset ~ temp + pop + prec + temp_seas + numObs + 
                   diapause.stage + flights +
                   (1|id_cells) + (1|scientificName__) +
                   (0 + temp | scientificName__) + 
                   (0 + prec | scientificName__) +
                   (0 + temp_seas | scientificName__) +
                   temp:prec + temp:pop + 
                   temp:diapause.stage + prec:flights,
                 data = pred.df.low, 
                 cov_ranef = list(scientificName = insect_tree4), 
                 bayes = TRUE)
  
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
  
  a.total <- rbind(a, a.low, a.high)
  a.total$inter_eff2 <- factor(a.total$inter_eff2, levels = c("Low", "Mid", "High"))
  
  mdf_phylo_spp <- mdf_phylo_spp %>% 
    mutate(inter_eff2 = case_when(pop >=1 ~ "High",
                                  pop <= -1 ~ "Low",
                                  pop > -1 & pop < 1 ~ "Mid"))
  mdf_phylo_spp$inter_eff2 <- factor(mdf_phylo_spp$inter_eff2, levels = c("Low", "Mid", "High"))
  
  
  
  p <-  ggplot(a.total, mapping = aes(x = inter_eff1, y = mean)) +
    geom_point(mdf_phylo_spp, mapping = aes(x = temp, y = onset, color = inter_eff2),
               alpha = 0.08) + 
    geom_ribbon(mapping = aes(ymin = `0.025quant`, ymax = `0.975quant`,
                              fill = inter_eff2), 
                alpha = 0.15) +
    geom_line(mapping = aes(color = inter_eff2), size = 1.05) +
    labs(x = "Temperature", y = "Emergence", fill = "Population", color = "Population") + 
    scale_color_manual(labels = c("Low", "Mid", "High"), values = c("turquoise4", "grey5", "violetred2")) +
    scale_fill_manual(values = c("turquoise4", "grey5", "violetred2")) +
    guides(fill = FALSE) +
    theme_bw()
  
  return(p)
  
}

fig2_a <- plot_pglmm_inter_2a(model_df = mdf_phylo_spp, pred_vals = -3:3)
fig2_a

### fig3_c
plot_pglmm_inter_3c <- function(model_df, pred_vals){
  
  mdf <- model_df
  n <- names(mdf)
  
  new_inter_term_adult <- data.frame(X =  rep(NA, 7), onset = rep(NA, 7), onset_low = rep(NA, 7), onset_high = rep(NA,7),
                                     year = rep(NA, 7), offset = rep(NA, 7), offset_low = rep(NA, 7), offset_high = rep(NA, 7),
                                     duration =  rep(NA, 7), scientificName =  rep(NA, 7), Order =  rep(NA, 7), id_cells =  rep(NA, 7),
                                     lon =  rep(NA, 7), lat = rep(NA, 7), OS= rep(NA, 7), 
                                     prec = rep(0,7), temp = pred_vals, bio4 = rep(0, 7), bio15 = rep(0, 7), pop = rep(0, 7),
                                     numObs = rep(0, 7), prec_seas = rep(0, 7), temp_seas = rep(0, 7), 
                                     higher_taxon= rep(NA, 7), flights = rep(NA, 7), development= rep(NA, 7),
                                     immature.habitat =  rep(NA, 7), diapause.stage= rep("Adult", 7), 
                                     larval.diet= rep(NA, 7), seas= rep(NA, 7), elevation = rep(NA, 7),
                                     lon_cor= rep(NA, 7), lat_cor= rep(NA, 7), elev_cor= rep(NA, 7),
                                     hopkins_cor= rep(NA, 7), hopkins_onset= rep(NA, 7), mean_onset= rep(NA, 7),
                                     sd_onset= rep(NA, 7), seas2= rep(NA, 7))
  
  pred.df.adult <- rbind(mdf, new_inter_term_adult)
  
  new_inter_term.egg <- data.frame(X =  rep(NA, 7), onset = rep(NA, 7), onset_low = rep(NA, 7), onset_high = rep(NA,7),
                                   year = rep(NA, 7), offset = rep(NA, 7), offset_low = rep(NA, 7), offset_high = rep(NA, 7),
                                   duration =  rep(NA, 7), scientificName =  rep(NA, 7), Order =  rep(NA, 7), id_cells =  rep(NA, 7),
                                   lon =  rep(NA, 7), lat = rep(NA, 7), OS= rep(NA, 7), 
                                   prec = rep(0,7), temp = pred_vals, bio4 = rep(0, 7), bio15 = rep(0, 7), pop = rep(0, 7),
                                   numObs = rep(0, 7), prec_seas = rep(0, 7), temp_seas = rep(0, 7), 
                                   higher_taxon= rep(NA, 7), flights = rep(NA, 7), development= rep(NA, 7),
                                   immature.habitat =  rep(NA, 7), diapause.stage= rep("Egg", 7), 
                                   larval.diet= rep(NA, 7), seas= rep(NA, 7), elevation = rep(NA, 7),
                                   lon_cor= rep(NA, 7), lat_cor= rep(NA, 7), elev_cor= rep(NA, 7),
                                   hopkins_cor= rep(NA, 7), hopkins_onset= rep(NA, 7), mean_onset= rep(NA, 7),
                                   sd_onset= rep(NA, 7), seas2= rep(NA, 7))
  
  pred.df.egg <- rbind(mdf, new_inter_term.egg)
  
  new_inter_term.larva <- data.frame(X =  rep(NA, 7), onset = rep(NA, 7), onset_low = rep(NA, 7), onset_high = rep(NA,7),
                                     year = rep(NA, 7), offset = rep(NA, 7), offset_low = rep(NA, 7), offset_high = rep(NA, 7),
                                     duration =  rep(NA, 7), scientificName =  rep(NA, 7), Order =  rep(NA, 7), id_cells =  rep(NA, 7),
                                     lon =  rep(NA, 7), lat = rep(NA, 7), OS= rep(NA, 7), 
                                     prec = rep(0,7), temp = pred_vals, bio4 = rep(0, 7), bio15 = rep(0, 7), pop = rep(0, 7),
                                     numObs = rep(0, 7), prec_seas = rep(0, 7), temp_seas = rep(0, 7), 
                                     higher_taxon= rep(NA, 7), flights = rep(NA, 7), development= rep(NA, 7),
                                     immature.habitat =  rep(NA, 7), diapause.stage= rep("Larvae", 7), 
                                     larval.diet= rep(NA, 7), seas= rep(NA, 7), elevation = rep(NA, 7),
                                     lon_cor= rep(NA, 7), lat_cor= rep(NA, 7), elev_cor= rep(NA, 7),
                                     hopkins_cor= rep(NA, 7), hopkins_onset= rep(NA, 7), mean_onset= rep(NA, 7),
                                     sd_onset= rep(NA, 7), seas2= rep(NA, 7))
  
  pred.df.larva <- rbind(mdf, new_inter_term.larva)
  
  new_inter_term_pupa <- data.frame(X =  rep(NA, 7), onset = rep(NA, 7), onset_low = rep(NA, 7), onset_high = rep(NA,7),
                                    year = rep(NA, 7), offset = rep(NA, 7), offset_low = rep(NA, 7), offset_high = rep(NA, 7),
                                    duration =  rep(NA, 7), scientificName =  rep(NA, 7), Order =  rep(NA, 7), id_cells =  rep(NA, 7),
                                    lon =  rep(NA, 7), lat = rep(NA, 7), OS= rep(NA, 7), 
                                    prec = rep(0,7), temp = pred_vals, bio4 = rep(0, 7), bio15 = rep(0, 7), pop = rep(0, 7),
                                    numObs = rep(0, 7), prec_seas = rep(0, 7), temp_seas = rep(0, 7), 
                                    higher_taxon= rep(NA, 7), flights = rep(NA, 7), development= rep(NA, 7),
                                    immature.habitat =  rep(NA, 7), diapause.stage= rep("Pupae", 7), 
                                    larval.diet= rep(NA, 7), seas= rep(NA, 7), elevation = rep(NA, 7),
                                    lon_cor= rep(NA, 7), lat_cor= rep(NA, 7), elev_cor= rep(NA, 7),
                                    hopkins_cor= rep(NA, 7), hopkins_onset= rep(NA, 7), mean_onset= rep(NA, 7),
                                    sd_onset= rep(NA, 7), seas2= rep(NA, 7))
  
  pred.df.pupa <- rbind(mdf, new_inter_term_pupa)
  
  
  m.adult <- pglmm(onset ~ temp + pop + prec + temp_seas + numObs + 
                     diapause.stage + flights +
                     (1|id_cells) + (1|scientificName__) +
                     (0 + temp | scientificName__) + 
                     (0 + prec | scientificName__) +
                     (0 + temp_seas | scientificName__) +
                     temp:prec + temp:pop + 
                     temp:diapause.stage + prec:flights,
                   data = pred.df.adult, 
                   cov_ranef = list(scientificName = insect_tree4), 
                   bayes = TRUE)
  
  m.egg <- pglmm(onset ~ temp + pop + prec + temp_seas + numObs + 
                   diapause.stage + flights +
                   (1|id_cells) + (1|scientificName__) +
                   (0 + temp | scientificName__) + 
                   (0 + prec | scientificName__) +
                   (0 + temp_seas | scientificName__) +
                   temp:prec + temp:pop + 
                   temp:diapause.stage + prec:flights,
                 data = pred.df.egg, 
                 cov_ranef = list(scientificName = insect_tree4), 
                 bayes = TRUE)
  
  m.larva <- pglmm(onset ~ temp + pop + prec + temp_seas + numObs + 
                     diapause.stage + flights +
                     (1|id_cells) + (1|scientificName__) +
                     (0 + temp | scientificName__) + 
                     (0 + prec | scientificName__) +
                     (0 + temp_seas | scientificName__) +
                     temp:prec + temp:pop + 
                     temp:diapause.stage + prec:flights,
                   data = pred.df.larva, 
                   cov_ranef = list(scientificName = insect_tree4), 
                   bayes = TRUE)
  
  m.pupa <- pglmm(onset ~ temp + pop + prec + temp_seas + numObs + 
                    diapause.stage + flights +
                    (1|id_cells) + (1|scientificName__) +
                    (0 + temp | scientificName__) + 
                    (0 + prec | scientificName__) +
                    (0 + temp_seas | scientificName__) +
                    temp:prec + temp:pop + 
                    temp:diapause.stage + prec:flights,
                  data = pred.df.pupa, 
                  cov_ranef = list(scientificName = insect_tree4), 
                  bayes = TRUE)
  
  a.adult <- m.adult$inla.model$summary.linear.predictor[2644:2650,] %>% 
    mutate(inter_eff1 = pred_vals,
           inter_eff2 = "Adult",
           inter_eff2_vals = 0)
  a.egg <- m.egg$inla.model$summary.linear.predictor[2644:2650,] %>% 
    mutate(inter_eff1 = pred_vals,
           inter_eff2 = "Egg",
           inter_eff2_vals = -1)
  a.larva <- m.larva$inla.model$summary.linear.predictor[2644:2650,] %>% 
    mutate(inter_eff1 = pred_vals,
           inter_eff2 = "Larvae",
           inter_eff2_vals = 1)
  a.pupa <- m.pupa$inla.model$summary.linear.predictor[2644:2650,] %>% 
    mutate(inter_eff1 = pred_vals,
           inter_eff2 = "Pupae",
           inter_eff2_vals = 1)
  
  a.total <- rbind(a.adult, a.egg, a.larva, a.pupa)
  
  p <-  ggplot(a.total, mapping = aes(x = inter_eff1, y = mean)) +
    geom_point(mdf_phylo_spp, mapping = aes(x = temp, y = onset, color = diapause.stage),
               alpha = 0.08) + 
    geom_ribbon(mapping = aes(ymin = `0.025quant`, ymax = `0.975quant`,
                              fill = inter_eff2), 
                alpha = 0.15) +
    geom_line(mapping = aes(color = inter_eff2), size = 1.05) +
    labs(x = "Temperature", y = "Emergence", fill = "Diapause Stage", color = "Diapause Stage") + 
    scale_color_viridis_d() +
    scale_fill_viridis_d() +
    guides(fill = FALSE) +
    theme_bw()
  
  return(p)
  
}

## Something 
fig3_c <- plot_pglmm_inter_3c(model_df = mdf_phylo_spp, pred_vals = -3:3)
fig3_c

### fig3_d
plot_pglmm_inter_3d <- function(model_df, pred_vals){
  
  mdf <- model_df
  
  new_inter_term_uv <- data.frame(X =  rep(NA, 7), onset = rep(NA, 7), onset_low = rep(NA, 7), onset_high = rep(NA,7),
                                  year = rep(NA, 7), offset = rep(NA, 7), offset_low = rep(NA, 7), offset_high = rep(NA, 7),
                                  duration =  rep(NA, 7), scientificName =  rep(NA, 7), Order =  rep(NA, 7), id_cells =  rep(NA, 7),
                                  lon =  rep(NA, 7), lat = rep(NA, 7), OS= rep(NA, 7), 
                                  prec = pred_vals, temp = rep(0,7), bio4 = rep(0, 7), bio15 = rep(0, 7), pop = rep(0, 7),
                                  numObs = rep(0, 7), prec_seas = rep(0, 7), temp_seas = rep(0, 7), 
                                  higher_taxon= rep(NA, 7), flights = rep("univoltine", 7), development= rep(NA, 7),
                                  immature.habitat =  rep(NA, 7), diapause.stage= rep(NA, 7), 
                                  larval.diet= rep(NA, 7), seas= rep(NA, 7), elevation = rep(NA, 7),
                                  lon_cor= rep(NA, 7), lat_cor= rep(NA, 7), elev_cor= rep(NA, 7),
                                  hopkins_cor= rep(NA, 7), hopkins_onset= rep(NA, 7), mean_onset= rep(NA, 7),
                                  sd_onset= rep(NA, 7), seas2= rep(NA, 7))
  
  pred.df.uv <- rbind(mdf, new_inter_term_uv)
  
  new_inter_term_nuv <- data.frame(X =  rep(NA, 7), onset = rep(NA, 7), onset_low = rep(NA, 7), onset_high = rep(NA,7),
                                   year = rep(NA, 7), offset = rep(NA, 7), offset_low = rep(NA, 7), offset_high = rep(NA, 7),
                                   duration =  rep(NA, 7), scientificName =  rep(NA, 7), Order =  rep(NA, 7), id_cells =  rep(NA, 7),
                                   lon =  rep(NA, 7), lat = rep(NA, 7), OS= rep(NA, 7), 
                                   prec = pred_vals, temp = rep(0,7), bio4 = rep(0, 7), bio15 = rep(0, 7), pop = rep(0, 7),
                                   numObs = rep(0, 7), prec_seas = rep(0, 7), temp_seas = rep(0, 7), 
                                   higher_taxon= rep(NA, 7), flights = rep("not univoltine", 7), development= rep(NA, 7),
                                   immature.habitat =  rep(NA, 7), diapause.stage= rep(NA, 7), 
                                   larval.diet= rep(NA, 7), seas= rep(NA, 7), elevation = rep(NA, 7),
                                   lon_cor= rep(NA, 7), lat_cor= rep(NA, 7), elev_cor= rep(NA, 7),
                                   hopkins_cor= rep(NA, 7), hopkins_onset= rep(NA, 7), mean_onset= rep(NA, 7),
                                   sd_onset= rep(NA, 7), seas2= rep(NA, 7))
  
  pred.df.nuv <- rbind(mdf, new_inter_term_nuv)
  
  
  m.uv <- pglmm(onset ~ temp + pop + prec + temp_seas + numObs + 
                  diapause.stage + flights +
                  (1|id_cells) + (1|scientificName__) +
                  (0 + temp | scientificName__) + 
                  (0 + prec | scientificName__) +
                  (0 + temp_seas | scientificName__) +
                  temp:prec + temp:pop + 
                  temp:diapause.stage + prec:flights,
                data = pred.df.uv, 
                cov_ranef = list(scientificName = insect_tree4), 
                bayes = TRUE)
  
  m.nuv <- pglmm(onset ~ temp + pop + prec + temp_seas + numObs + 
                   diapause.stage + flights +
                   (1|id_cells) + (1|scientificName__) +
                   (0 + temp | scientificName__) + 
                   (0 + prec | scientificName__) +
                   (0 + temp_seas | scientificName__) +
                   temp:prec + temp:pop + 
                   temp:diapause.stage + prec:flights,
                 data = pred.df.nuv, 
                 cov_ranef = list(scientificName = insect_tree4), 
                 bayes = TRUE)
  
  a.uv <- m.uv$inla.model$summary.linear.predictor[2644:2650,] %>% 
    mutate(inter_eff1 = pred_vals,
           inter_eff2 = "Univoltine",
           inter_eff2_vals = 0)
  a.nuv <- m.nuv$inla.model$summary.linear.predictor[2644:2650,] %>% 
    mutate(inter_eff1 = pred_vals,
           inter_eff2 = "Not Univoltine",
           inter_eff2_vals = -1)
  
  a.total <- rbind(a.uv, a.nuv)
  
  p <-  ggplot(a.total, mapping = aes(x = inter_eff1, y = mean)) +
    geom_point(mdf_phylo_spp, mapping = aes(x = prec, y = onset, color = flights),
               alpha = 0.08) + 
    geom_ribbon(mapping = aes(ymin = `0.025quant`, ymax = `0.975quant`,
                              fill = inter_eff2), 
                alpha = 0.15) +
    geom_line(mapping = aes(color = inter_eff2), size = 1.05) +
    labs(x = "Precipitation", y = "Emergence", fill = "Voltinism", color = "Voltinism") + 
    scale_color_manual(values = c("purple", "black")) +
    scale_fill_manual(values = c("purple", "black")) +
    scale_x_continuous(limits = c(-3,3)) +
    guides(fill = FALSE) +
    theme_bw()
  
  return(p)
  
}

## Something 
fig3_d <- plot_pglmm_inter_3d(model_df = mdf_phylo_spp, pred_vals = c(-3:3))
fig3_d 

# Combine and Save Duration plots
onset_interaction_plots <- egg::ggarrange(fig2_a, fig2_b, fig3_c, fig3_d,
                                          ncol = 2, labels = c("A", "B", 'C', "D"))

ggsave(onset_interaction_plots, filename = "Figures/resubmission/onset_interactions.png",
       width = 8, height = 4)
