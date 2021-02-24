library(tidyverse)
library(gridExtra)

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

insect_tree4 <- ape::drop.tip(phy = insect_tree3, tip = eusocial) %>% 
  ape::drop.tip(tip = "Aeshna cyanea")  %>% 
  ape::drop.tip(tip = "Diabrotica undecimpunctata")

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

inla_onset <-  pglmm_onset$inla.model
onset_fe <- inla_onset$summary.fixed %>% 
  dplyr::add_rownames(var = "vars")

# functions to plot interactions

# Onset interaction response curves
tmp.x <- mdf_phylo_spp$temp

## pop interactions
pop_high <- 1
pop <- 0
pop_low <- -1

on_int <- filter(onset_fe, vars == "(Intercept)")$mean
on_tmp <- filter(onset_fe, vars == "temp")$mean
on_pop <- filter(onset_fe, vars == "pop")$mean
on_tmppop <- filter(onset_fe, vars == "temp:pop")$mean

onset.low = on_int + 
  (on_tmp * tmp.x) + 
  (on_pop * pop_low) + 
  (on_tmppop * pop_low * tmp.x)

onset = on_int + 
  (on_tmp * tmp.x) + 
  (on_pop * pop) + 
  (on_tmppop * pop * tmp.x)

onset.high = on_int + 
  (on_tmp * tmp.x) + 
  (on_pop * pop_high) + 
  (on_tmppop * pop_high * tmp.x)

tmpz <- c(tmp.x, tmp.x, tmp.x)
onset_pop <- c(onset.high, onset, onset.low)
popz <- c(rep(.95, 2643), rep(-0.02, 2643), rep(-1, 2643))

pop_temp <- data.frame(temp = tmpz, onset = onset_pop, Pop = as.factor(popz)) 

tpop <- ggplot() + 
  geom_line(pop_temp, mapping = aes(x = temp, y = onset, color = Pop), size = 1.05) +
  labs(x = "Temperature", y = "Emergence") + 
  scale_color_manual(labels = c("Low", "Mid", "High"), values = c("turquoise4", "grey5", "violetred2")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, vjust = 2.5))


tpop


## prec interactions
prec_high <- 1
prec <- 0
prec_low <- -1

on_int <- filter(onset_fe, vars == "(Intercept)")$mean
on_tmp <- filter(onset_fe, vars == "temp")$mean
on_prec <- filter(onset_fe, vars == "prec")$mean
on_tmpprec <- filter(onset_fe, vars == "temp:prec")$mean

onset.low.prec = on_int + 
  (on_tmp * tmp.x) + 
  (on_prec * prec_low) + 
  (on_tmpprec * prec_low * tmp.x)

onset.prec = on_int + 
  (on_tmp * tmp.x) + 
  (on_prec * prec) + 
  (on_tmpprec * prec* tmp.x)

onset.high.prec = on_int + 
  (on_tmp * tmp.x) + 
  (on_prec * prec_high) + 
  (on_tmpprec * prec_high * tmp.x)

tmpz <- c(tmp.x, tmp.x, tmp.x)
onset_prec <- c(onset.high.prec, onset.prec, onset.low.prec)
precz <- c(rep(.99, 2643), rep(0.09, 2643), rep(-0.81, 2643))

prec_temp <- data.frame(temp = tmpz, onset = onset_prec, Prec = as.factor(precz)) 

tprec <- ggplot() + 
  geom_line(prec_temp, mapping = aes(x = temp, y = onset_prec, color = Prec), size = 1.05) +
  labs(x = "Temperature", y = "Emergence") + 
  scale_color_manual(labels = c("Low", "Mid", "High"), values = c("brown", "cyan", "Blue")) +
  theme_bw()

tprec

# temp:diapause stage
on_int <- filter(onset_fe, vars == "(Intercept)")$mean
on_tmp <- filter(onset_fe, vars == "temp")$mean
on_egg <- filter(onset_fe, vars == "diapause.stageEgg")$mean
on_larvae <- filter(onset_fe, vars == "diapause.stageLarvae")$mean
on_pupae <- filter(onset_fe, vars == "diapause.stagePupae")$mean
on_tmpegg <- filter(onset_fe, vars == "temp:diapause.stageEgg")$mean
on_tmplarvae <- filter(onset_fe, vars == "temp:diapause.stageLarvae")$mean
on_tmppupae <- filter(onset_fe, vars == "temp:diapause.stagePupae")$mean

onset.tmp.egg = on_int + 
  (on_tmp * tmp.x) + 
  (on_egg) + 
  (on_tmpegg * tmp.x)

onset.tmp.larvae = on_int + 
  (on_tmp * tmp.x) + 
  (on_larvae) + 
  (on_tmplarvae * tmp.x)

onset.tmp.pupae = on_int + 
  (on_tmp * tmp.x) + 
  (on_pupae) + 
  (on_tmppupae * tmp.x)

onset.tmp.adult= on_int + 
  (on_tmp * tmp.x)

tmpz <- c(tmp.x, tmp.x, tmp.x, tmp.x)
onset_diapausestage_tmp <- c(onset.tmp.egg, onset.tmp.larvae, onset.tmp.pupae, onset.tmp.adult)
Diapause.Stage <- c(rep("Egg", 2643), rep("Larvae", 2643)
                    , rep("Pupae", 2643), rep("Adult", 2643))

diapausestage_temp <- data.frame(temp = tmpz, onset = onset_diapausestage_tmp, 
                                 Diapause.Stage = as.factor(Diapause.Stage)) 

ds_onset <- ggplot() + 
  geom_line(diapausestage_temp, mapping = aes(x = temp, y = onset_diapausestage_tmp, color = Diapause.Stage), size = 1.05) +
  labs(x = "Temperature", y = "Emergence") + 
  scale_color_viridis_d() +
  theme_bw()

ds_onset

# prec:flights
prec.x <- mdf_phylo_spp$prec

on_int <- filter(onset_fe, vars == "(Intercept)")$mean
on_prec <- filter(onset_fe, vars == "prec")$mean
on_univoltine <- filter(onset_fe, vars == "flightsunivoltine")$mean
on_precuni <- filter(onset_fe, vars == "prec:flightsunivoltine")$mean

onset.prec.uni = on_int + 
  (on_prec * prec.x) + 
  (on_univoltine) + 
  (on_precuni * prec.x)

onset.prec.notuni = on_int + 
  (on_prec * prec.x) 

precz <- c(prec.x, prec.x)
onset_prec_voltinism <- c(onset.prec.uni, onset.prec.notuni)
volt <- c(rep("Univoltine", 2643), rep("Not Univoltine", 2643))

voltinism_prec <- data.frame(prec = precz, onset = onset_prec_voltinism, 
                             Voltinism = as.factor(volt)) 

volt_onset <- ggplot() + 
  geom_line(voltinism_prec, mapping = aes(x = prec, y = onset_prec_voltinism, color = volt), size = 1.05) +
  labs(x = "Precipitation", y = "Emergence") + 
  scale_color_manual(values = c("purple", "black")) +
  labs(color = "Voltinism") + 
  theme_bw()

volt_onset

# Combine and Save onset plots
onset_plots <- egg::ggarrange(tpop, tprec,  ds_onset, volt_onset,
                              ncol = 2, labels = c("A", "B", 'C', "D"))

ggsave(onset_plots, filename = "Figures/onset_interactions.png",
       width = 8, height = 4)


## offset interactions
pglmm_offset <- pglmm(offset ~ temp_seas + prec +
                        seas + diapause.stage + immature.habitat + larval.diet + 
                        temp_seas:diapause.stage +
                        numObs + 
                        (1 | id_cells) + (1 | scientificName__) + 
                        (0 + temp_seas | scientificName__) + 
                        (0 + prec | scientificName__), 
                      data = mdf_phylo_spp, 
                      cov_ranef = list(scientificName = insect_tree4), 
                      bayes = TRUE)

inla_offset <-  pglmm_offset$inla.model
offset_fe <- inla_offset$summary.fixed %>% 
  dplyr::add_rownames(var = "vars")

# temp_seas:diapause.stage
off_int <- filter(offset_fe, vars == "(Intercept)")$mean
off_tmp_seas <- filter(offset_fe, vars == "temp_seas")$mean
off_egg <- filter(offset_fe, vars == "diapause.stageEgg")$mean
off_larvae <- filter(offset_fe, vars == "diapause.stageLarvae")$mean
off_pupae <- filter(offset_fe, vars == "diapause.stagePupae")$mean
off_tmpseas_egg <- filter(offset_fe, vars == "temp_seas:diapause.stageEgg")$mean
off_tmpseas_larvae <- filter(offset_fe, vars == "temp_seas:diapause.stageLarvae")$mean
off_tmpseas_pupae <- filter(offset_fe, vars == "temp_seas:diapause.stagePupae")$mean

tmp_seas.x2 <- mdf_phylo_spp$temp_seas

offset.tmp.egg = off_int + 
  (off_tmp_seas * tmp_seas.x2) + 
  (off_egg) + 
  (off_tmpseas_egg * tmp_seas.x2)

offset.tmp.larvae = off_int + 
  (off_tmp_seas * tmp_seas.x2) + 
  (off_larvae) + 
  (off_tmpseas_larvae * tmp_seas.x2)

offset.tmp.pupae = off_int + 
  (off_tmp_seas * tmp_seas.x2) + 
  (off_pupae) + 
  (off_tmpseas_pupae * tmp_seas.x2)

offset.tmp.adult= off_int + 
  (off_tmp_seas * tmp_seas.x2)

tmpzeas2 <- c(tmp_seas.x2, tmp_seas.x2, tmp_seas.x2, tmp_seas.x2)
off_diapause <- c(offset.tmp.egg, offset.tmp.larvae, offset.tmp.pupae, offset.tmp.adult)
diapause_stages <- c(rep("Egg", 2643), rep("Larvae", 2643), 
                 rep("Pupae", 2643),  rep("Adult", 2643))

off_tmp_larval.diet <- data.frame(Temp = tmpzeas2, Offset = off_diapause, Diapause.Stage = as.factor(diapause_stages)) 


off_ds <- ggplot() + 
  geom_line(off_tmp_larval.diet, mapping = aes(x = tmpzeas2 , y = Offset, color = Diapause.Stage), size = 1.05) +
  labs(x = "Temperature Seasonality", y = "Termination") + 
  scale_color_viridis_d() +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, vjust = 2.5))

off_ds

ggsave(plot = off_ds, filename = "Figures/offset_interactions.png",
       width = 5, height = 3)

## Duration
pglmm_duration <- pglmm(duration ~ temp + prec + temp_seas +
                          diapause.stage + 
                          immature.habitat +
                          larval.diet + 
                          numObs + 
                          (1 | id_cells) + (1 | scientificName__) + 
                          (0 + temp | scientificName__) + 
                          (0 + prec | scientificName__) + 
                          temp:prec + temp:larval.diet + temp:immature.habitat,
                        data = mdf_phylo_spp, 
                        cov_ranef = list(scientificName = insect_tree4), 
                        bayes = TRUE)

inla_duration <-  pglmm_duration$inla.model
duration_fe <- inla_duration$summary.fixed %>% 
  dplyr::add_rownames(var = "vars")


## prec interactions
prec_high <- 1
prec <- 0
prec_low <- -1

dur_int <- filter(duration_fe, vars == "(Intercept)")$mean
dur_tmp <- filter(duration_fe, vars == "temp")$mean
dur_prec <- filter(duration_fe, vars == "prec")$mean
dur_tmpprec <- filter(duration_fe, vars == "temp:prec")$mean

duration.low.prec = dur_int + 
  (dur_tmp * tmp.x) + 
  (dur_tmpprec * prec_low * tmp.x)

duration.prec = dur_int + 
  (dur_tmp * tmp.x) + 
  (dur_tmpprec * prec* tmp.x)

duration.high.prec = dur_int + 
  (dur_tmp * tmp.x) + 
  (dur_tmpprec * prec_high * tmp.x)

tmpz <- c(tmp.x, tmp.x, tmp.x)
duration_prec <- c(duration.high.prec, duration.prec, duration.low.prec)
precz <- c(rep(.99, 2643), rep(0.09, 2643), rep(-0.81, 2643))

prec_temp <- data.frame(temp = tmpz, duration = duration_prec, Prec = as.factor(precz)) 

dur_tprec <- ggplot() + 
  geom_line(prec_temp, mapping = aes(x = temp, y = duration_prec, color = Prec), size = 1.05) +
  labs(x = "Temperature", y = "Duration") + 
  scale_color_manual(labels = c("Low", "Mid", "High"), values = c("brown", "cyan", "Blue")) +
  theme_bw()

dur_tprec


#temp:immature.habitat
dur_int<- filter(duration_fe, vars == "(Intercept)")$mean
dur_tmp <- filter(duration_fe, vars == "temp")$mean
dur_freshwater <- filter(duration_fe, vars == "immature.habitatfresh water")$mean
dur_underground <- filter(duration_fe, vars == "immature.habitatunderground")$mean
dur_tmp_fw <- filter(duration_fe, vars == "temp:immature.habitatfresh water")$mean
dur_tmp_under <- filter(duration_fe, vars == "temp:immature.habitatunderground")$mean

dur_fw_tmp <- dur_int + 
  (tmp.x * dur_tmp) + 
  (dur_freshwater) + 
  (dur_tmp_fw* tmp.x)

dur_underground_tmp <- dur_int + 
  (tmp.x * dur_tmp) + 
  (dur_underground) + 
  (dur_tmp_under * tmp.x)

dur_above.ground.vegetation_tmp <- dur_int + 
  (tmp.x * dur_tmp)

tmpz.x <- c(tmp.x, tmp.x, tmp.x)
dur_ih_tmp <- c(dur_above.ground.vegetation_tmp, dur_fw_tmp, dur_underground_tmp)
Immature.Habitat <- c(rep("Above Ground", 2643), rep("Freshwater", 2643)
                      , rep("Underground", 2643))

tmp_immature.habitat <- data.frame(tmpz = tmpz.x, duration = dur_ih_tmp, 
                                   Immature.Habitat = as.factor(Immature.Habitat)) 


dur_tmp_ih <- ggplot() + 
  geom_line(tmp_immature.habitat, mapping = aes(x = tmpz, y = duration, color = Immature.Habitat), size = 1.05) +
  labs(x = "Temperature", y = "Duration") + 
  scale_color_manual(values = c("olivedrab4", "deepskyblue2", "orangered4")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, vjust = 2.5))

dur_tmp_ih


#temp:larval.diet
dur_int<- filter(duration_fe, vars == "(Intercept)")$mean
dur_tmp <- filter(duration_fe, vars == "temp")$mean
dur_detritivorous<- filter(duration_fe, vars == "larval.dietdetritivorous")$mean
dur_herbivorous <- filter(duration_fe, vars == "larval.dietlive plants")$mean
dur_tmp_det <- filter(duration_fe, vars == "temp:larval.dietdetritivorous")$mean
dur_tmp_herb <- filter(duration_fe, vars == "temp:larval.dietlive plants")$mean

dur_herb_tmp <- dur_int + 
  (tmp.x * dur_tmp) + 
  (dur_herbivorous) + 
  (dur_tmp_herb* tmp.x)

dur_det_tmp <- dur_int + 
  (tmp.x * dur_tmp) + 
  (dur_detritivorous) + 
  (dur_tmp_det * tmp.x)

dur_carn_tmp <- dur_int + 
  (tmp.x * dur_tmp)

tmpz.x <- c(tmp.x, tmp.x, tmp.x)
dur_ld_tmp <- c(dur_carn_tmp, dur_det_tmp, dur_herb_tmp)
Larval.Diet <- c(rep("Carnivorous", 2643), rep("Detritivorous", 2643)
                      , rep("Herbivorous", 2643))

tmp_immature.habitat <- data.frame(tmpz = tmpz.x, duration = dur_ld_tmp, 
                                   Larval.Diet = as.factor(Larval.Diet)) 


dur_tmp_ld <- ggplot() + 
  geom_line(tmp_immature.habitat, mapping = aes(x = tmpz, y = duration, color = Larval.Diet), size = 1.05) +
  labs(x = "Temperature", y = "Duration") + 
  scale_color_manual(values = c("#d1495b", "#edae49",  "#66a182"))+
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, vjust = 2.5))

dur_tmp_ld

# Combine and Save Duration plots
duration_plots <- egg::ggarrange(dur_tprec, dur_tmp_ih, dur_tmp_ld,
                              ncol = 2, labels = c("A", "B", 'C'))

ggsave(duration_plots, filename = "Figures/duration_interactions.png",
       width = 8, height = 4)

