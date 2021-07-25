library(tidyverse)
library(lmerTest)
library(lme4)
library(car)
library(sjPlot)
library(geoR)

# read in data
# removed D. undecimpunctata from analyses b/c there is no true diapause in this species
# removed L. cyanea b/c it wasn't in phylogeny
mdf <- read.csv("data/LMM_data/phenoTraits_data_withNumObsData_updatedSeasTrat.csv", stringsAsFactors = FALSE) %>% 
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
                       seas2 + temp_seas:seas2 + prec:seas2 +
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
# temp_seas:immature.habitat + temp_seas:larval.diet

off_traits_top <- lmer(offset ~ temp_seas + prec +
                         seas2 + diapause.stage + immature.habitat + larval.diet + 
                         temp_seas:immature.habitat + temp_seas:larval.diet + 
                         numObs + 
                         (1 | id_cells) + (1 | scientificName) + 
                         (0 + temp_seas | scientificName) + 
                         (0 + prec | scientificName), 
                       data = mdf, 
                       REML = FALSE, 
                       lmerControl(optimizer = "bobyqa"))

car::vif(off_traits_top) 

## Vif >5 w/ temp_seas, so removing temp_seas:immature.habitat

off_traits_top1 <- lmer(offset ~ temp_seas + prec +
                         seas2 + diapause.stage + immature.habitat + larval.diet + 
                         temp_seas:larval.diet +
                         numObs + 
                         (1 | id_cells) + (1 | scientificName) + 
                         (0 + temp_seas | scientificName) + 
                         (0 + prec | scientificName), 
                       data = mdf, 
                       REML = FALSE, 
                       lmerControl(optimizer = "bobyqa"))

off_traits_top2 <- lmer(offset ~ temp_seas + prec +
                          seas2 + diapause.stage + immature.habitat + larval.diet + 
                          temp_seas:immature.habitat +
                          numObs + 
                          (1 | id_cells) + (1 | scientificName) + 
                          (0 + temp_seas | scientificName) + 
                          (0 + prec | scientificName), 
                        data = mdf, 
                        REML = FALSE, 
                        lmerControl(optimizer = "bobyqa"))


car::vif(off_traits_top1)
car::vif(off_traits_top2)

MuMIn::AICc(off_traits_top1, off_traits_top2)

off_traits_top <- off_traits_top2

# see if top model is stable
off_traits_top_s <- step(off_traits_top)
off_traits_top_s

off_traits_top <- lmer(offset ~ temp_seas + prec +
                         seas2 + diapause.stage + immature.habitat + larval.diet + 
                         temp_seas:immature.habitat +
                         numObs + 
                         (1 | id_cells) + (1 | scientificName) + 
                         (0 + temp_seas | scientificName) + 
                         (0 + prec | scientificName), 
                       data = mdf, 
                       REML = FALSE, 
                       lmerControl(optimizer = "bobyqa"))


summary(off_traits_top)
car::vif(off_traits_top)
plot_model(off_traits_top, type = "pred", terms = c("temp_seas", "immature.habitat")) 
MuMIn::r.squaredGLMM(off_traits_top) #R2m 0.185 R2c 0.613

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

mdf_phylo_spp <- mdf %>% 
  filter(scientificName %in% tree_sp)

insect_tree4 <- ape::drop.tip(phy = insect_tree3, tip = eusocial)%>% 
  ape::drop.tip(tip = "Aeshna cyanea") %>% 
  ape::drop.tip(tip = "Diabrotica undecimpunctata")

pglmm_offset <- pglmm(offset ~ temp_seas + prec +
                        seas2 + diapause.stage + immature.habitat + larval.diet + 
                        temp_seas:immature.habitat +
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

kableExtra::save_kable(fe, file = "Tables/resubmission/pglmm_FE_offset.html")
kableExtra::save_kable(re, file = "Tables/resubmission/pglmm_RE_offset.html")

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

## Examine spatial autocorrelation
library(ncf)
rdf <- mutate(mdf_phylo_spp, residuals = resids)

fit3 = correlog(x = rdf$lon, y = rdf$lat, z = rdf$residuals,
                increment = 25000, latlon = F, resamp = 100)

summary(fit3)
plot(fit3)

fit3_df <- data.frame(distance = fit3$mean.of.class, 
                      correlation = fit3$correlation,
                      p.value = if_else(fit3$p < 0.025, 
                                        true = "Sig",
                                        false = "Not Sig")) %>% 
  filter(distance < 500000 & distance > 1)



offset_correlogram_plot <- ggplot() +
  geom_line(fit3_df,
            mapping = aes(x = distance, y = correlation)) +
  geom_point(fit3_df,
             mapping = aes(x = distance, y = correlation, fill = p.value),
             size = 3, shape = 21, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("white", "black")) +
  labs(x = "Distance", y = "Moran's I", fill = "P-Value") +
  theme_classic()

offset_correlogram_plot

save(offset_correlogram_plot, file = "Figures/resubmission/offset_correlogram_plot.Rdata")


## Add spatial component to pglmm

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

summary(pglmm_offset_sp)
rr2::R2(pglmm_offset_sp) # 0.580

resids_sp <- resid(pglmm_offset_sp)
rdf_sp <- mutate(mdf_phylo_spp, residuals = resids)

fit_sp = correlog(x = rdf_sp$lon, y = rdf_sp$lat, z = rdf_sp$residuals,
                increment = 100000, latlon = F, resamp = 100)

summary(fit_sp)
plot(fit_sp)

## get spatial correlation model coef
re <- ranef(pglmm_offset_sp) %>% knitr::kable()
fe <- fixef(pglmm_offset_sp) %>% knitr::kable()

kableExtra::save_kable(fe, file = "Tables/resubmission/pglmm_FE_offset_spatial_correlation.html")
kableExtra::save_kable(re, file = "Tables/resubmission/pglmm_RE_offset_spatial_correlation.html")