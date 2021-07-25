library(dplyr)
library(ggplot2)

# pulling R2 values from LMM models


## Onset
# Climate only models
# R2m R2c 0.355 0.769

#Climate and Traits
# R2m R2c 0.454863 0.7679138

## Offset
# climate only models
# 0.0764 0.650

#Climate and traits
# 0.185 0.613

## Duration
# climate only model
# 0.352 0.773

# climate and traits
# 0.464 0.762

#### PGLMMs
#Onset #0.746
#Offset # 0.579
#Duration # 0.746

r2s <- data.frame(model = c("Emergence", "Emergence", "Emergence", "Emergence", "Emergence", 
                            "Termination", "Termination", "Termination","Termination", "Termination",
                            "Duration", "Duration", "Duration", "Duration", "Duration"),
                  model_type = c("Climate", "Climate", "Full", "Full", "PGLMM", 
                                 "Climate", "Climate", "Full", "Full", "PGLMM", 
                                 "Climate", "Climate", "Full", "Full", "PGLMM"),
                  R2_type = c("Marginal", "Conditional", "Marginal", "Conditional", "Conditional",
                              "Marginal", "Conditional", "Marginal", "Conditional", "Conditional",
                              "Marginal", "Conditional", "Marginal", "Conditional", "Conditional"),
                  R2 = c(0.355, 0.769, 0.454863, 0.7679138, 0.746,
                         0.0764, 0.650, 0.185, 0.613, 0.579,
                         0.352, 0.773, 0.464, 0.762, 0.746))

mdf <- r2s %>% 
  mutate(model = factor(model, levels = c("Emergence", "Termination", "Duration")))

R2_plot <- ggplot() + 
  # Plot conditional
  geom_bar( mdf,
            stat = "identity",
            position = "identity",
            mapping = aes(x = model_type, y = R2, fill = R2_type)
  ) +
  geom_bar( filter(mdf, R2_type == "Marginal"),
            stat = "identity",
            position = "identity",
            mapping = aes(x = model_type, y = R2, fill = R2_type)
  ) +
  scale_y_continuous(
    name = expression(R^2),
    limits = c(0,0.85) ,
    expand = c(0,0)
  ) +
  scale_x_discrete(
    name = "Model Type"
  ) +
  theme_classic() +
  facet_wrap(~model) +
  scale_fill_manual(values = c("turquoise", "turquoise4"), 
                    guide=guide_legend(reverse=T)) +
  theme(
    strip.background = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 0.5, size = 11),
    strip.text.x = element_text(size = 12),
    text = element_text(size = 14),
    legend.title = element_blank(),
    legend.position = "bottom"
  )

R2_plot

ggsave(plot = R2_plot, filename = "Figures/resubmission/R2.png",
       width = 5, height = 4)
