library(tidyverse)
library(lmerTest)
library(lme4)
library(car)
library(kableExtra)
library(sjPlot)
library(geoR)

# read in data
# removed D. undecimpunctata from analyses b/c there is no true diapause in this species
# removed L. cyanea b/c it wasn't in phylogeny
mdf <- read.csv("data/LMM_data/phenoTraits_data_withNumObsData_updatedSeasTrat.csv", stringsAsFactors = FALSE) %>% 
  filter(scientificName != "Diabrotica undecimpunctata") %>% 
  filter(scientificName != "Libellula cyanea") %>% 
  mutate(numObs = scale(numObs)) # Note the other variables were already scaled

########## Onset Models ###########
m1on <- lmer(onset ~ temp + pop + prec + prec_seas + temp_seas + temp:pop + temp:prec +
               numObs + 
               (1|id_cells) + (1|scientificName) +
               (0 + temp | scientificName) + 
               (0 + pop | scientificName) + 
               (0 + prec | scientificName) +
               (0 + temp_seas | scientificName) + 
               (0 + prec_seas | scientificName) ,
             data = mdf, REML = FALSE, 
             lmerControl(optimizer = "bobyqa"))

on_step <- step(m1on)
on_step

#Model found:
#  onset ~ temp + pop + prec + temp_seas + numObs + (1 | id_cells) + (1 |scientificName) + 
#  (0 + temp | scientificName) + (0 + prec |  scientificName) + 
#  (0 + temp_seas | scientificName) + (0 + prec_seas | scientificName) + 
#  temp:pop + temp:prec

on_top <- lmer(onset ~ temp + pop + prec + temp_seas + temp:pop + temp:prec +
                 numObs + 
                 (1 | id_cells) + (1 | scientificName) + 
                 (0 + temp | scientificName) + 
                 (0 + prec |scientificName) + 
                 (0 + temp_seas | scientificName),
               data = mdf,
               REML = FALSE,
               lmerControl(optimizer = "bobyqa"))

# make sure model is stable

on_top_s <- step(on_top)
on_top_s 

on_top <- lmer(onset ~ temp + pop + prec + temp_seas + temp:pop + temp:prec +
                 numObs + 
                 (1 | id_cells) + (1 | scientificName) + 
                 (0 + temp | scientificName) + 
                 (0 + prec |scientificName) + 
                 (0 + temp_seas | scientificName),
               data = mdf,
               REML = FALSE,
               lmerControl(optimizer = "bobyqa"))

summary(on_top)
car::vif(on_top)
plot_model(on_top, type = "pred", terms = c("temp", "prec"), ci.lvl = NA)
plot_model(on_top, type = "pred", terms = c("temp", "pop"), ci.lvl = NA)
MuMIn::r.squaredGLMM(on_top) # R2m 0.355  R2c 0.769


## Add traits to model
# onset models with traits

on_traits <- lmer(onset ~ temp + pop + prec + temp_seas + temp:prec + temp:pop +
                    numObs +
                    development + prec:development + pop:development + prec:development + 
                    diapause.stage + temp:diapause.stage + pop:diapause.stage + prec:diapause.stage + temp_seas:diapause.stage +
                    flights + temp:flights + pop:flights + prec:flights + temp_seas:flights +
                    immature.habitat + temp:immature.habitat + pop:immature.habitat + prec:immature.habitat + temp_seas:immature.habitat +
                    larval.diet + temp:larval.diet + pop:larval.diet + prec:larval.diet + 
                    (1 | id_cells) + (1 | scientificName) + 
                    (0 + temp | scientificName) + 
                    (0 + temp_seas | scientificName) +
                    (0 + prec |scientificName),
                  data = mdf,
                  REML = FALSE,
                  lmerControl(optimizer = "bobyqa"))

mon_traits_s <- step(on_traits)
mon_traits_s

#Model found:
#  onset ~ temp + pop + prec + temp_seas + numObs  + diapause.stage + flights + 
#  (1 | id_cells) + (1 | scientificName) + (0 + temp | scientificName) + 
#  (0 + prec | scientificName) + 
# + temp:prec + temp:pop + temp:diapause.stage + prec:flights

on_traits_top <- lmer(onset ~ temp + pop + prec + temp_seas + numObs + 
                        diapause.stage + flights +
                        (1|id_cells) + (1|scientificName) +
                        (0 + temp | scientificName) + 
                        (0 + prec | scientificName) +
                        (0 + temp_seas | scientificName) +
                        temp:prec + temp:pop + 
                        temp:diapause.stage + prec:flights,
                      data = mdf, 
                      REML = FALSE, 
                      lmerControl(optimizer = "bobyqa"))

# double check model stability
on_traits_top_s <- step(on_traits_top)
on_traits_top_s

summary(on_traits_top)
car::vif(on_traits_top)
plot_model(on_traits_top, type = "pred", terms = c("temp", "prec"), ci.lvl = NA)
plot_model(on_traits_top, type = "pred", terms = c("temp", "pop"), ci.lvl = NA)
plot_model(on_traits_top, type = "pred", terms = c("temp", "diapause.stage"), ci.lvl = NA)
plot_model(on_traits_top, type = "pred", terms = c("prec", "flights"), ci.lvl = NA)
MuMIn::r.squaredGLMM(on_traits_top) ## R2m = 0.454 R2c = 0.768

### PLMM #########
# use package phyr and inla to generate phylogenetic linear mixed models

library(phyr)
library(INLA)
# read in phylogeny
insect_tree3 = ape::read.tree("data/phylogeny/insect_test2.tre")

insect_tree3$tip.label <- word(insect_tree3$tip.label, start = 1, end = 2, sep = "_") %>% 
  sub(pattern = "_", replacement = " ")

# Remove eusocial species from phylogeny
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

mdf2 <- read.csv("data/LMM_data/phenoTraits_data_withNumObsData_updatedSeasTrat.csv", stringsAsFactors = FALSE) %>% 
  filter(scientificName != "Diabrotica undecimpunctata") %>% 
  mutate(numObs = scale(numObs))

mdf_phylo_spp <- mdf2 %>% 
  filter(scientificName %in% tree_sp)

insect_tree4 <- ape::drop.tip(phy = insect_tree3, tip = eusocial) %>% 
  ape::drop.tip(tip = "Aeshna cyanea") %>% 
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

summary(pglmm_onset)

re <- phyr::ranef(pglmm_onset) %>% 
  knitr::kable()
fe <- phyr::fixef(pglmm_onset) %>% 
  knitr::kable() 

kableExtra::save_kable(fe, file = "Tables/pglmm_FE_onset.html")
kableExtra::save_kable(re, file = "Tables/pglmm_RE_onset.html")


inla_onset = pglmm_onset$inla.model
summary(inla_onset)
inla_onset$summary.fixed # kld are small, good

## Model Assumption Checks
# residuals
resids <- residuals(pglmm_onset)
hist(resids)
qqnorm(resids)
qqline(resids)
shapiro.test(resids)

# spatial autocorrelation
jitter_dist = 12500 # 12.5 km
resid_onset <- mdf_phylo_spp %>% 
  mutate(resid = resids)

gdf <- as.geodata(resid_onset, coords.col = c("lon", "lat"), data.col = "resid")
gdf <-  geoR::jitterDupCoords(gdf, max = jitter_dist)
v1 <- variog(gdf, trend = "1st")

plot(v1)

## Goodness of fit

rr2::R2(pglmm_onset) #0.747
