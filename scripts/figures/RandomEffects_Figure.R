library(phyr)
library(tidyverse)
library(stringr)

load("ModelOutputs/resubmission/pglmm_onset.Rdata") 
# note that due to huge file size of pglmm_onset.Rdata, .Rdata is not saved
# in GitHub directory. Please run /scripts/analyses/spatial_PLMM.R through line 
# 172 to make your own pglmm_onset.Rdata file

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
load("ModelOutputs/resubmission/pglmm_offset.Rdata")
# note that due to huge file size of pglmm_offset.Rdata, .Rdata is not saved
# in GitHub directory. Please run /scripts/analyses/spatial_PLMM_offset.R through line 
# 309 to make your own pglmm_offset.Rdata file

x <- pglmm_offset_sp
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
load("ModelOutputs/resubmission/pglmm_duration.Rdata")
# note that due to huge file size of pglmm_duration.Rdata, .Rdata is not saved
# in GitHub directory. Please run /scripts/analyses/spatial_PLMM_duration.R through line 
# 278 to make your own pglmm_duration.Rdata file

x <- pglmm_duration_sp
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

comb_re_pd

ggsave(filename = "Figures/resubmission/RandomEffects_PosteriorDistribution.png", 
       dpi = 400, height = 8, width = 8)
