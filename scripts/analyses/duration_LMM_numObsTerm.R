library(tidyverse)
library(lmerTest)
library(lme4)
library(car)
library(sjPlot)
library(geoR)

# read in data
# removed D. undecimpunctata from analyses b/c there is no true diapause in this species
# removed L. cyanea b/c it wasn't in phylogeny
mdf <- read.csv("data/LMM_data/phenoTraits_data_withNumObsData.csv", stringsAsFactors = FALSE) %>% 
  filter(scientificName != "Diabrotica undecimpunctata") %>% 
  filter(scientificName != "Libellula cyanea") %>% 
  mutate(numObs = scale(numObs)) # Note the other variables were already scaled

########## Onset Models ###########
m1dur <- lmer(duration ~ temp + pop + prec + prec_seas + temp_seas + temp:pop + temp:prec +
                numObs + 
                (1|id_cells) + (1|scientificName) +
                (0 + temp | scientificName) + 
                (0 + pop | scientificName) + 
                (0 + prec | scientificName) +
                (0 + temp_seas | scientificName) + 
                (0 + prec_seas | scientificName),
              data = mdf, REML = FALSE, 
              lmerControl(optimizer = "bobyqa"))

dur_step <- step(m1dur)
dur_step

#Model found:
#  duration ~ temp + prec + temp_seas + numObs + (1 | id_cells) + (1 | scientificName) + 
#  (0 + temp | scientificName) + (0 + prec | scientificName) + 
#  (0 + temp_seas | scientificName) + (0 + prec_seas | scientificName) + 
#  (0 + temp:prec | scientificName) + temp:prec

dur_top <- lmer(duration ~ temp + prec + temp_seas + temp:prec +
                  numObs + 
                  (1 | id_cells) + (1 | scientificName) + 
                  (0 + temp | scientificName) + 
                  (0 + prec |scientificName) + 
                  (0 + temp_seas | scientificName), 
                data = mdf,
                REML = FALSE, 
                lmerControl(optimizer = "bobyqa"))

# make sure top model is stable
dur_top_s <- step(dur_top)
dur_top_s

dur_top <- lmer(duration ~ temp + prec + temp_seas +  temp:prec +
                  numObs + 
                  (1 | id_cells) + (1 | scientificName) + 
                  (0 + temp | scientificName) + 
                  (0 + prec | scientificName) +
                  (0 + temp_seas | scientificName), 
                data = mdf,
                REML = FALSE, 
                lmerControl(optimizer = "bobyqa"))


car::vif(dur_top)
summary(dur_top)

plot_model(dur_top, type = "pred", terms = c("temp", "prec"), ci.lvl = NA)
MuMIn::r.squaredGLMM(dur_top) #R2m 0.35 #R2c 0.77

# trait LMM duration

m_dur_traits <- lmer(duration ~ temp + prec + temp_seas +  temp:prec +
                       development + temp:development + prec:development + temp_seas:development +
                       seas + temp:seas + prec:seas + temp_seas:seas + 
                       diapause.stage + temp:diapause.stage + prec:diapause.stage + temp_seas:diapause.stage + 
                       flights + temp:flights + prec:flights + temp_seas:flights +
                       immature.habitat + temp:immature.habitat + prec:immature.habitat + temp_seas:immature.habitat +
                       larval.diet + temp:larval.diet + prec:larval.diet + temp_seas:larval.diet + 
                       numObs +
                       (1 | id_cells) + (1 | scientificName) + 
                       (0 + temp | scientificName) + 
                       (0 + prec | scientificName) +
                       (0 + temp_seas | scientificName), 
                     data = mdf,
                     REML = FALSE, 
                     lmerControl(optimizer = "bobyqa"))


# stepwise regression to select model
m_dur_traits_s = step(m_dur_traits, reduce.random = T)
m_dur_traits_s

#Model found:
#  duration ~ temp + prec + temp_seas + diapause.stage + 
#  immature.habitat +  larval.diet + numObs + (1 | id_cells) + (1 | scientificName) + 
#  (0 + temp | scientificName) + (0 + temp_seas | scientificName) + 
#  temp:prec + temp:immature.habitat + temp:larval.diet

dur_traits_top <- lmer(duration ~ temp + prec + temp_seas +
                         diapause.stage + 
                         immature.habitat +
                         larval.diet + 
                         numObs + 
                         (1 | id_cells) + (1 | scientificName) + 
                         (0 + temp | scientificName) + 
                         (0 + prec | scientificName) + 
                         temp:prec + temp:larval.diet + temp:immature.habitat,               
                       data = mdf, 
                       REML = FALSE, 
                       lmerControl(optimizer = "bobyqa"))

# make sure the top model is stable

dur_traits_top_s <- step(dur_traits_top)
dur_traits_top_s

car::vif(dur_traits_top) 
summary(dur_traits_top)
plot_model(dur_traits_top, type = "pred", terms = c("temp", "prec"), ci.lvl = NA)
plot_model(dur_traits_top, type = "pred", terms = c("temp", "immature.habitat"), ci.lvl = NA)
plot_model(dur_traits_top, type = "pred", terms = c("temp", "larval.diet"), ci.lvl = NA)
MuMIn::r.squaredGLMM(dur_traits_top) #0.462 0.762

### PLMM #########

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

mdf2 <- read.csv("data/LMM_data/phenoTraits_data_withNumObsData.csv", stringsAsFactors = FALSE) %>% 
  filter(scientificName != "Diabrotica undecimpunctata") %>% 
  mutate(numObs = scale(numObs))

mdf_phylo_spp <- mdf2 %>% 
  filter(scientificName %in% tree_sp)

insect_tree4 <- ape::drop.tip(phy = insect_tree3, tip = eusocial)%>% 
  ape::drop.tip(tip = "Aeshna cyanea") %>% 
  ape::drop.tip(tip = "Diabrotica undecimpunctata")

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

summary(pglmm_duration)

re <- ranef(pglmm_duration) %>% knitr::kable()
fe <- fixef(pglmm_duration) %>% knitr::kable()
""
kableExtra::save_kable(fe, file = "Tables/pglmm_FE_duration.html")
kableExtra::save_kable(re, file = "Tables/pglmm_RE_duration.html")

inla_duration = pglmm_duration$inla.model
summary(inla_duration)
inla_duration$summary.fixed # kld are small, good

# residuals
resids <- residuals(pglmm_duration)
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
rr2::R2(pglmm_duration) # 0.749
