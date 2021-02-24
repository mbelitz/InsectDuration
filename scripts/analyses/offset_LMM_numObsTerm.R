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

########## offset Models ###########
m1off <- lmer(offset ~ temp + pop + prec + prec_seas + temp_seas + temp:pop + temp:prec +
                numObs + 
                (1|id_cells) + (1|scientificName) +
                (0 + temp | scientificName) + 
                (0 + pop | scientificName) + 
                (0 + prec | scientificName) +
                (0 + temp_seas | scientificName) + 
                (0 + prec_seas | scientificName),
              data = mdf, REML = FALSE, 
              lmerControl(optimizer = "bobyqa"))

off_step <- step(m1off)
off_step

#Model found:
#  offset ~ prec + temp_seas + numObs + (1 | id_cells) + (1 | scientificName) + 
#  (0 + temp | scientificName) + (0 + pop | scientificName) + 
#  (0 + temp_seas | scientificName) + (0 + prec_seas | scientificName) 

off_top <- lmer(offset ~ prec + temp_seas +
                  numObs + 
                  (1|id_cells) + (1|scientificName) +
                  (0 + temp_seas | scientificName) +
                  (0 + prec | scientificName),
                data = mdf, 
                REML = FALSE, 
                lmerControl(optimizer = "bobyqa"))

# make sure top model is stable
off_top_s <- step(off_top)
off_top_s

off_top <- lmer(offset ~ temp_seas + prec + numObs + 
                  (1|id_cells) + (1|scientificName) +
                  (0 + temp_seas | scientificName) +
                  (0 + prec | scientificName),
                data = mdf, 
                REML = FALSE, 
                lmerControl(optimizer = "bobyqa"))

summary(off_top)
car::vif(off_top)
MuMIn::r.squaredGLMM(off_top) #R2m 0.0757 R2c 0.653

## offset LMM with traits

m_off_traits <- lmer(offset ~ temp_seas + prec + 
                       seas + temp_seas:seas + prec:seas +
                       development + temp_seas:development + prec:development +
                       diapause.stage + temp_seas:diapause.stage + prec:diapause.stage +
                       flights +  temp_seas:flights + prec:flights +
                       immature.habitat +  temp_seas:immature.habitat + prec:immature.habitat +
                       larval.diet +  temp_seas:larval.diet + prec:larval.diet + 
                       numObs + 
                       (1|id_cells) + (1|scientificName) +
                       (0 + temp_seas | scientificName) +
                       (0 + prec | scientificName),
                     data = mdf, 
                     REML = FALSE, 
                     lmerControl(optimizer = "bobyqa"))

# stepwise regression to select model
m_off_traits_s = step(m_off_traits, reduce.random = T)
m_off_traits_s

#Model found:
#  offset ~ temp_seas + prec + seas + diapause.stage + immature.habitat + 
#  larval.diet + numObs + (1 | id_cells) + (1 | scientificName) + (0 +temp_seas | scientificName) +
# (0 + prec | scientificName) + 
#  temp_seas: diapause.stage + temp_seas:immature.habitat + temp_seas:larval.diet

off_traits_top <- lmer(offset ~ temp_seas + prec +
                         seas + diapause.stage + immature.habitat + larval.diet + 
                         temp_seas:immature.habitat + temp_seas:larval.diet + temp_seas:diapause.stage +
                         numObs + 
                         (1 | id_cells) + (1 | scientificName) + 
                         (0 + temp_seas | scientificName) + 
                         (0 + prec | scientificName), 
                       data = mdf, 
                       REML = FALSE, 
                       lmerControl(optimizer = "bobyqa"))

car::vif(off_traits_top) 

## Vif >5 w/ temp_seas, so removing temp_seas:immature.habitat

off_traits_top <- lmer(offset ~ temp_seas + prec +
                         seas + diapause.stage + immature.habitat + larval.diet + 
                         temp_seas:larval.diet + temp_seas:diapause.stage +
                         numObs + 
                         (1 | id_cells) + (1 | scientificName) + 
                         (0 + temp_seas | scientificName) + 
                         (0 + prec | scientificName), 
                       data = mdf, 
                       REML = FALSE, 
                       lmerControl(optimizer = "bobyqa"))


# see if top model is stable
off_traits_top_s <- step(off_traits_top)
off_traits_top_s

off_traits_top <- lmer(offset ~ temp_seas + prec +
                         seas + diapause.stage + immature.habitat + larval.diet + 
                         temp_seas:diapause.stage +
                         numObs + 
                         (1 | id_cells) + (1 | scientificName) + 
                         (0 + temp_seas | scientificName) + 
                         (0 + prec | scientificName), 
                       data = mdf, 
                       REML = FALSE, 
                       lmerControl(optimizer = "bobyqa"))


summary(off_traits_top)
car::vif(off_traits_top)
plot_model(off_traits_top, type = "pred", terms = c("temp_seas", "diapause.stage"), ci.lvl = NA) 
MuMIn::r.squaredGLMM(off_traits_top) #R2m 0.189 R2c 0.605

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

summary(pglmm_offset)

re <- ranef(pglmm_offset) %>% knitr::kable()
fe <- fixef(pglmm_offset) %>% knitr::kable()

kableExtra::save_kable(fe, file = "Tables/pglmm_FE_offset.html")
kableExtra::save_kable(re, file = "Tables/pglmm_RE_offset.html")

inla_offset = pglmm_offset$inla.model
summary(inla_offset)
inla_offset$summary.fixed # kld are small, good

# residuals
resids <- residuals(pglmm_offset)
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

rr2::R2(pglmm_offset) # 0.580