library(phyr)
library(ape)
library(tidyverse)
library(stringr)

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


# Onset Posterior distributions

x <- pglmm_onset
n_samp <- 1000

random_samps <- lapply(x$inla.model$marginals.hyperpar, 
                       function(x) INLA::inla.rmarginal(n_samp, INLA::inla.tmarginal(function(x) sqrt(1 / x), x))) %>% 
  setNames(names(x$random.effects))

random_samps2 <- random_samps[-length(random_samps)] %>% 
  setNames(names(x$random.effects)) %>%
  dplyr::as_tibble() %>%
  tidyr::pivot_longer(cols = dplyr::everything(),
                      names_to = "var",
                      values_to = "val") %>%
  dplyr::mutate(effect_type = "Random Effects")

fixed_samps <- lapply(x$inla.model$marginals.fixed, function(x) INLA::inla.rmarginal(n_samp, x)) %>% 
  dplyr::as_tibble() %>%
  tidyr::pivot_longer(cols = dplyr::everything(),
                      names_to = "var",
                      values_to = "val") %>%
  dplyr::mutate(effect_type = "Fixed Effects")

samps <- dplyr::bind_rows(random_samps2, fixed_samps) %>%
  dplyr::mutate(effect_type = factor(effect_type, 
                                     levels = c("Random Effects", "Fixed Effects")))

ci <- samps %>%
  dplyr::group_by(var, effect_type) %>%
  dplyr::summarise(lower = quantile(val, 0.025),
                   upper = quantile(val, 0.975),
                   mean = mean(val),
                   .groups = "drop_last")
sig_vars <- ci %>%
  dplyr::mutate(sig = ifelse(effect_type == "Random Effects",
                             "CI no overlap with zero",
                             ifelse(sign(lower) == sign(upper),
                                    "CI no overlap with zero",
                                    "CI overlaps zero"))) %>%
  dplyr::select(var, sig)

samps <- samps %>%
  dplyr::left_join(sig_vars, by = "var") %>%
  dplyr::group_by(var) %>%
  dplyr::filter(abs(val - mean(val)) < (10 * sd(val))) %>% 
  dplyr::ungroup()

pal <- c("dodgerblue2", "firebrick2", "springgreen2")

re_on <- filter(samps, effect_type == "Random Effects") %>% 
  mutate(Phenometric = "Onset")
re_ci_on <- filter(ci, effect_type == "Random Effects") %>% 
  mutate(Phenometric = "Onset")

on_re_pd <- ggplot(re_on, aes(val, var, height = ..density..)) +
  ggridges::geom_density_ridges(aes(val, var, alpha = sig, fill = sig), 
                                color = "gray70") +
  geom_point(aes(x = mean, y = var), data = re_ci_on, inherit.aes = FALSE) +
  geom_errorbarh(aes(xmin = lower, xmax = upper, y = var), data = re_ci_on,
                 inherit.aes = FALSE, height = 0.1) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
  scale_alpha_manual(values = c(0.8, 0.2)) +
  scale_fill_manual(values = rev(pal)) +
  ylab("") +
  xlab("Estimate") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 16)) + 
  ggtitle("Onset")

#ggsave("NumObsResults/Figures/onset_RE.png", dpi = 400, plot = on_re_pd)

# Offset Posterior distributions

x <- pglmm_offset
n_samp <- 1000

random_samps <- lapply(x$inla.model$marginals.hyperpar, 
                       function(x) INLA::inla.rmarginal(n_samp, INLA::inla.tmarginal(function(x) sqrt(1 / x), x))) %>% 
  setNames(names(x$random.effects))

random_samps2 <- random_samps[-length(random_samps)] %>% 
  setNames(names(x$random.effects)) %>%
  dplyr::as_tibble() %>%
  tidyr::pivot_longer(cols = dplyr::everything(),
                      names_to = "var",
                      values_to = "val") %>%
  dplyr::mutate(effect_type = "Random Effects")

fixed_samps <- lapply(x$inla.model$marginals.fixed, function(x) INLA::inla.rmarginal(n_samp, x)) %>% 
  dplyr::as_tibble() %>%
  tidyr::pivot_longer(cols = dplyr::everything(),
                      names_to = "var",
                      values_to = "val") %>%
  dplyr::mutate(effect_type = "Fixed Effects")

samps <- dplyr::bind_rows(random_samps2, fixed_samps) %>%
  dplyr::mutate(effect_type = factor(effect_type, 
                                     levels = c("Random Effects", "Fixed Effects")))

ci <- samps %>%
  dplyr::group_by(var, effect_type) %>%
  dplyr::summarise(lower = quantile(val, 0.025),
                   upper = quantile(val, 0.975),
                   mean = mean(val),
                   .groups = "drop_last")
sig_vars <- ci %>%
  dplyr::mutate(sig = ifelse(effect_type == "Random Effects",
                             "CI no overlap with zero",
                             ifelse(sign(lower) == sign(upper),
                                    "CI no overlap with zero",
                                    "CI overlaps zero"))) %>%
  dplyr::select(var, sig)

samps <- samps %>%
  dplyr::left_join(sig_vars, by = "var") %>%
  dplyr::group_by(var) %>%
  dplyr::filter(abs(val - mean(val)) < (10 * sd(val))) %>% 
  dplyr::ungroup()

pal <- c("dodgerblue2", "firebrick2", "springgreen2")

re_off <- filter(samps, effect_type == "Random Effects") %>% 
  mutate(Phenometric = "Offset")
re_ci_off <- filter(ci, effect_type == "Random Effects") %>% 
  mutate(Phenometric = "Offset")

off_re_pd <- ggplot(re_off, aes(val, var, height = ..density..)) +
  ggridges::geom_density_ridges(aes(val, var, alpha = sig, fill = sig), 
                                color = "gray70") +
  geom_point(aes(x = mean, y = var), data = re_ci_off, inherit.aes = FALSE) +
  geom_errorbarh(aes(xmin = lower, xmax = upper, y = var), data = re_ci_off,
                 inherit.aes = FALSE, height = 0.1) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
  scale_alpha_manual(values = c(0.8, 0.2)) +
  scale_fill_manual(values = "firebrick2") +
  ylab("") +
  xlab("Estimate") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 16)) + 
  ggtitle("Offset")

#ggsave("NumObsResults/Figures/offset_RE.png", dpi = 400, plot = off_re_pd)

# Duration Posterior distributions
x <- pglmm_duration
n_samp <- 1000

random_samps <- lapply(x$inla.model$marginals.hyperpar, 
                       function(x) INLA::inla.rmarginal(n_samp, INLA::inla.tmarginal(function(x) sqrt(1 / x), x))) %>% 
  setNames(names(x$random.effects))

random_samps2 <- random_samps[-length(random_samps)] %>% 
  setNames(names(x$random.effects)) %>%
  dplyr::as_tibble() %>%
  tidyr::pivot_longer(cols = dplyr::everything(),
                      names_to = "var",
                      values_to = "val") %>%
  dplyr::mutate(effect_type = "Random Effects")

fixed_samps <- lapply(x$inla.model$marginals.fixed, function(x) INLA::inla.rmarginal(n_samp, x)) %>% 
  dplyr::as_tibble() %>%
  tidyr::pivot_longer(cols = dplyr::everything(),
                      names_to = "var",
                      values_to = "val") %>%
  dplyr::mutate(effect_type = "Fixed Effects")

samps <- dplyr::bind_rows(random_samps2, fixed_samps) %>%
  dplyr::mutate(effect_type = factor(effect_type, 
                                     levels = c("Random Effects", "Fixed Effects")))

ci <- samps %>%
  dplyr::group_by(var, effect_type) %>%
  dplyr::summarise(lower = quantile(val, 0.025),
                   upper = quantile(val, 0.975),
                   mean = mean(val),
                   .groups = "drop_last")
sig_vars <- ci %>%
  dplyr::mutate(sig = ifelse(effect_type == "Random Effects",
                             "CI no overlap with zero",
                             ifelse(sign(lower) == sign(upper),
                                    "CI no overlap with zero",
                                    "CI overlaps zero"))) %>%
  dplyr::select(var, sig)

samps <- samps %>%
  dplyr::left_join(sig_vars, by = "var") %>%
  dplyr::group_by(var) %>%
  dplyr::filter(abs(val - mean(val)) < (10 * sd(val))) %>% 
  dplyr::ungroup()

pal <- c("dodgerblue2", "firebrick2", "springgreen2")

re_dur <- filter(samps, effect_type == "Random Effects") %>% 
  mutate(Phenometric = "Duration")
re_ci_dur <- filter(ci, effect_type == "Random Effects") %>% 
  mutate(Phenometric = "Duration")

dur_re_pd <- ggplot(re_dur, aes(val, var, height = ..density..)) +
  ggridges::geom_density_ridges(aes(val, var, alpha = sig, fill = sig), 
                                color = "gray70") +
  geom_point(aes(x = mean, y = var), data = re_ci_dur, inherit.aes = FALSE) +
  geom_errorbarh(aes(xmin = lower, xmax = upper, y = var), data = re_ci_dur,
                 inherit.aes = FALSE, height = 0.1) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
  scale_alpha_manual(values = c(0.8, 0.2)) +
  scale_fill_manual(values = "dodgerblue2") +
  ylab("") +
  xlab("Estimate") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 16)) + 
  ggtitle("Duration")

#ggsave(filename = "NumObsResults/Figures/duration_re_pd.png", plot = dur_re_pd, dpi = 400)

cp <- cowplot::plot_grid(on_re_pd, off_re_pd, dur_re_pd, ncol = 2,
                         labels = c("A", "B", "C"))

#ggsave(plot = cp, filename = "Figures/allPhenometrics_re_pd.png", dpi = 400,
#       width = 10, height = 8)

re_comb <- rbind(re_on, re_off, re_dur)
re_ci_comb <- rbind(re_ci_on, re_ci_off, re_ci_dur)

re_comb$Phenometric <- factor(re_comb$Phenometric, levels = c("Onset", "Offset", "Duration"))

comb_re_pd <- ggplot(re_comb, aes(val, var, height = ..density..)) +
  ggridges::geom_density_ridges(aes(val, var, fill = Phenometric), 
                                alpha = 0.35, color = "gray70") +
  geom_point(aes(x = mean, y = var), data = re_ci_comb, inherit.aes = FALSE,
             alpha = 0.75) +
  geom_errorbarh(aes(xmin = lower, xmax = upper, y = var), data = re_ci_comb,
                 inherit.aes = FALSE, height = 0.1, alpha = 0.75) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
  scale_alpha_manual(values = c(0.8, 0.2)) +
  scale_fill_manual(values = rev(pal)) +
  ylab("") +
  xlab("Estimate") +
  theme_minimal() +
  theme(axis.text = element_text(size = 14),
        strip.text = element_text(size = 16)) 

ggsave(filename = "Figures/RandomEffects_PosteriorDistribution.png", 
       dpi = 400, height = 8, width = 8)

