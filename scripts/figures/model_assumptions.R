library(phyr)
library(ape)
library(tidyverse)
library(stringr)
library(grid)
library(ggplotify)
library(geoR)

# read in data
mdf <- read.csv("data/LMM_data/phenoTraits_data_withNumObsData.csv", stringsAsFactors = FALSE) %>% 
  filter(scientificName != "Diabrotica undecimpunctata") %>% 
  mutate(numObs = scale(numObs))


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

# plot_bayes(pglmm_onset, sort = TRUE)

## Offset

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


## Duration
pglmm_dur <- pglmm(duration ~ temp + prec + temp_seas +
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


# spatial autocorrelation
jitter_dist = 12500 # 12.5 km
resid_onset <- mdf_phylo_spp %>% 
  mutate(resid = resids)

gdf <- as.geodata(resid_onset, coords.col = c("lon", "lat"), data.col = "resid")
gdf <-  geoR::jitterDupCoords(gdf, max = jitter_dist)
v1 <- variog(gdf, trend = "1st")

v1
plot(v1)

a <- as.grob(function()plot(v1))

cp_onset <- cowplot::plot_grid(on_r_hist,on_qq, a, ncol = 1,
                               labels = c("A", "B", "C"))

cp_onset

## OFFSET TIME

# residuals
resids_off <- residuals(pglmm_offset)
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
jitter_dist = 12500 # 12.5 km
resid_offset <- mdf_phylo_spp %>% 
  mutate(resid = resids_off)

gdf <- as.geodata(resid_offset, coords.col = c("lon", "lat"), data.col = "resid")
gdf <-  geoR::jitterDupCoords(gdf, max = jitter_dist)
v2 <- variog(gdf, trend = "1st")

v2
plot(v2)

b <- as.grob(function()plot(v2))

cp_offset <- cowplot::plot_grid(off_r_hist,off_qq, b, ncol = 1,
                                labels = c("D", "E", "F"))

cp_offset

## Duration TIME

# residuals
resids_dur <- residuals(pglmm_duration)
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
jitter_dist = 12500 # 12.5 km
resid_duration <- mdf_phylo_spp %>% 
  mutate(resid = resids_dur)

gdf <- as.geodata(resid_duration, coords.col = c("lon", "lat"), data.col = "resid")
gdf <-  geoR::jitterDupCoords(gdf, max = jitter_dist)
v3 <- variog(gdf, trend = "1st")

v3
plot(v3)

c <- as.grob(function()plot(v3))

cp_duration <- cowplot::plot_grid(dur_r_hist,dur_qq, c, ncol = 1,
                                  labels = c("G", "H", "I"))

cp_duration


## TOtal cp

cp_total <- cowplot::plot_grid(cp_onset, cp_offset, cp_duration, ncol = 3, nrow = 1)

cp_total

ggsave(plot = cp_total, filename = "Figures/model_assumptions.png", width = 8, height = 7)

