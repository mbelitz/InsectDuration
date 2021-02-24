library(dplyr)
library(ggplot2)

# pulling R2 values from LMM models


## Onset
# Climate only models
# R2m R2c 0.3560002 0.7688878

#Climate and Traits
# R2m R2c 0.454863 0.7679138

## Offset
# climate only models
# 0.07578981 0.6529253

#Climate and traits
# 0.1879277 0.609447

## Duration
# climate only model
# 0.3504058 0.7737233

# climate and traits
# 0.4622136 0.7623328

#### PGLMMs
#Onset #0.748
#Offset # 0.580
#Duration # 0.749

r2s <- data.frame(model = c("Emergence", "Emergence", "Emergence", "Emergence", "Emergence", 
                            "Termination", "Termination", "Termination","Termination", "Termination",
                            "Duration", "Duration", "Duration", "Duration", "Duration"),
                  model_type = c("Climate", "Climate", "Full", "Full", "PGLMM", 
                                 "Climate", "Climate", "Full", "Full", "PGLMM", 
                                 "Climate", "Climate", "Full", "Full", "PGLMM"),
                  R2_type = c("Marginal", "Conditional", "Marginal", "Conditional", "Conditional",
                              "Marginal", "Conditional", "Marginal", "Conditional", "Conditional",
                              "Marginal", "Conditional", "Marginal", "Conditional", "Conditional"),
                  R2 = c(0.3560002, 0.7688878, 0.454863, 0.7679138, 0.748,
                         0.07578981, 0.6529253, 0.1879277, 0.609447, 0.580,
                         0.3504058, 0.7737233, 0.4622136, 0.7623328, 0.749))

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

ggsave(plot = R2_plot, filename = "Figures/R2.png",
       width = 5, height = 4)
