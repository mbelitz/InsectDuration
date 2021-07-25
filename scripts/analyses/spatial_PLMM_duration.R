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

########## Onset Models ###########
m1dur <- lmer(duration ~ temp + pop + prec + prec_seas + temp_seas + 
                temp:pop + temp:prec +
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

plot_model(dur_top, type = "pred", terms = c("temp", "prec"))
MuMIn::r.squaredGLMM(dur_top) #R2m 0.35 #R2c 0.77

# trait LMM duration

m_dur_traits <- lmer(duration ~ temp + prec + temp_seas +  temp:prec +
                       development + temp:development + prec:development + temp_seas:development +
                       seas2 + temp:seas2 + prec:seas2 + temp_seas:seas2 + 
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
#  (0 + temp | scientificName) + (0 + prec | scientificName) + 
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
plot_model(dur_traits_top, type = "pred", terms = c("temp", "prec"))
plot_model(dur_traits_top, type = "pred", terms = c("temp", "immature.habitat"))
plot_model(dur_traits_top, type = "pred", terms = c("temp", "larval.diet"))
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

mdf_phylo_spp <- mdf %>% 
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

kableExtra::save_kable(fe, file = "Tables/resubmission/pglmm_FE_duration.html")
kableExtra::save_kable(re, file = "Tables/resubmission/pglmm_RE_duration.html")

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


## Examine spatial autocorrelation
library(ncf)
rdf <- mutate(mdf_phylo_spp, residuals = resids)

fit3 = correlog(x = rdf$lon, y = rdf$lat, z = rdf$residuals,
                increment = 25000, latlon = F, resamp = 100)
plot(fit3)

fit3_df <- data.frame(distance = fit3$mean.of.class, 
                      correlation = fit3$correlation,
                      p.value = if_else(fit3$p < 0.025, 
                                        true = "Sig",
                                        false = "Not Sig")) %>% 
  filter(distance < 500000 & distance > 1)


duration_correlogram_plot <- ggplot() +
  geom_line(fit3_df,
            mapping = aes(x = distance, y = correlation)) +
  geom_point(fit3_df,
             mapping = aes(x = distance, y = correlation, fill = p.value),
             size = 3, shape = 21, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("white", "black")) +
  labs(x = "Distance", y = "Moran's I", fill = "P-Value") +
  theme_classic()

duration_correlogram_plot

save(duration_correlogram_plot, file = "Figures/resubmission/duration_correlogram_plot.Rdata")


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

## save results
resids_sp <- resid(pglmm_duration_sp)
rdf_sp <- mutate(mdf_phylo_spp, residuals = resids)

fit_sp = correlog(x = rdf_sp$lon, y = rdf_sp$lat, z = rdf_sp$residuals,
                  increment = 25000, latlon = F, resamp = 100)

plot(fit_sp)

## get spatial correlation model coef
re <- ranef(pglmm_duration_sp) %>% knitr::kable()
fe <- fixef(pglmm_duration_sp) %>% knitr::kable()

kableExtra::save_kable(fe, file = "Tables/resubmission/pglmm_FE_duration_spatial_correlation.html")
kableExtra::save_kable(re, file = "Tables/resubmission/pglmm_RE_duration_spatial_correlation.html")
