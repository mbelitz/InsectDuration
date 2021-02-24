library(tidyverse)
library(jpeg)

# read in total results
phenesse_results <- read.csv(file = "data/phenesse_climate_join/duration_climate_population_data_25km.csv",
                             stringsAsFactors = FALSE) 

#read in model df
# removed D. undecimpunctata from analyses b/c there is no true diapause in this species
# removed L. cyanea b/c it wasn't in phylogeny
mdf <- read.csv("data/LMM_data/phenoTraits_data_withNumObsData.csv", stringsAsFactors = FALSE) %>% 
  filter(scientificName != "Diabrotica undecimpunctata") %>% 
  filter(scientificName != "Libellula cyanea")

# count number of species used in final analysis by order
mdf_sum <- mdf %>% 
  group_by(Order) %>% 
  summarise(cells = n())

mdf_sum2 <- mdf %>% 
  distinct(scientificName, .keep_all = TRUE) %>% 
  group_by(Order) %>% 
  summarise(count = n()) 

num_spp_traits <- left_join(mdf_sum, mdf_sum2) %>% 
  mutate(Included = "Yes")


## count number of species that could have been used if 
## trait data were available

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

mig_spp <- c("Anax junius", "Pantala flavenscens", "Pantala hymenaea", "Tramea lacerata", 
             "Sympetrum corruptum", "Sympetrum vicinum", "Libellula pulchella", "Libellula vibrans",
             "Tramea lacerata", "Tramea onusta", "Pantala flavescens", "Libellula quadrimaculata",
             "Ertythrodiplax umbrata", "Epiaeschna heros", "Tramea carolina",
             "Libellula semifasciata", "Pantala hymenaea", 
             "Spoladea recurvalis", "Ponoquina ocola", "Plutella xylostella",
             "Chrysodeixis includens", "Phoebis sennae", "Abaeis nicippe",
             "Libytheana carinenta", "Agraulis vanillae", "Junonia coenia",
             "Danaus plexippus", "Vanessa virginiensis", "Vanessa cardui",
             "Vanessa atalanta", "Danaus gilippus", "Nymphalis antiopa", 
             "Polygonia interrogationis", "Lerema accius")

datadens <- phenesse_results %>% 
  group_by(scientificName, Order) %>% 
  summarise(count = n())

has_5_cells <- filter(datadens, count >= 5)

pheno_data <- phenesse_results %>% 
  filter (scientificName != "Apis mellifera") %>% 
  filter(scientificName != "Diabrotica undecimpunctata") %>% 
  filter(scientificName %in% has_5_cells$scientificName) %>% 
  filter(!scientificName %in% eusocial) %>% 
  filter(!scientificName %in% mig_spp)

noTraitsSum <- pheno_data %>% 
  group_by(Order) %>% 
  summarise(cells = n())

noTraitsSum2 <- pheno_data %>% 
  distinct(scientificName, .keep_all = TRUE) %>% 
  group_by(Order) %>% 
  summarise(count = n()) 

num_spp_Notraits <- left_join(noTraitsSum, noTraitsSum2) %>% 
  mutate(Included = "No",
         count = count - num_spp_traits$count,
         cells = cells - num_spp_traits$cells)

## combine together

total_sums <- rbind(num_spp_traits, num_spp_Notraits)

# add in jpgs
library(png)
dip <- readPNG("data/silhouttes/diptera_blackout.png", native = T)
ode <- readPNG("data/silhouttes/Zygoptera.Maxime-Dahirel.png")
hy <- readPNG("data/silhouttes/Halictus-rubicundus.PhyloPic.2a425db2.Melissa-Broussard.png")
lep <- readPNG("data/silhouttes/Saturniidae.Mali-o-Kodis-photograph-by-Jim-Vargo.png")
cic <- readPNG("data/silhouttes/cicada_blackout.png", native = T)
cole <- readPNG("data/silhouttes/Adephaga_Carabidae_Carabina_Ca.png", native = T)

library(cowplot)

p <-  total_sums %>% 
  mutate(Order = fct_reorder(Order, desc(count))) %>% 
  mutate(Included = factor(Included, levels = c("No", "Yes"))) %>%
  ggplot() + 
  geom_bar(stat = "identity", aes (y = Order, x = count, fill = Included), width = 0.5,) +
  labs(y = "", x = "Number of Species") +
  scale_fill_manual(values = c("thistle3", "orchid4"), guide=guide_legend(reverse=T)) +
  scale_x_continuous(expand = c(0,0)) +
  theme_classic()  + 
  theme(axis.text.y = element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y = element_blank(),
        legend.position="none"
  ) 
p

pimage <- axis_canvas(p, axis = 'y') + 
  draw_image(cic, y = 5.5,  scale = 1) +
  draw_image(cole, y = 4.5,   scale = 1) +
  draw_image(dip, y = 3.5,  scale = 1) +
  draw_image(hy, y = 2.5, scale = 1) +
  draw_image(ode, y = 1.5,  scale = 1) +
  draw_image(lep, y = 0.5,scale = 1)

a <- ggdraw(insert_yaxis_grob(p, pimage, position = "left"))

## Now for phenoestimates

p2 <-  total_sums %>% 
  mutate(Order = fct_reorder(Order, desc(count))) %>% 
  mutate(Included = factor(Included, levels = c("No", "Yes"))) %>%
  ggplot() + 
  geom_bar(stat = "identity", aes (y = Order, x = cells, fill = Included), width = 0.5,) +
  labs(y = "", x = "Pheno-estimates") +
  scale_fill_manual(values = c("thistle3", "orchid4"), guide=guide_legend(reverse=T)) +
  scale_x_continuous(expand = c(0,0)) +
  theme_classic()  + 
  theme(axis.text.y = element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y = element_blank(),
        legend.position="bottom"
  ) 
p2

pimage <- axis_canvas(p, axis = 'y') + 
  draw_image(cic, y = 5.5,  scale = 1) +
  draw_image(cole, y = 4.5,   scale = 1) +
  draw_image(dip, y = 3.5,  scale = 1) +
  draw_image(hy, y = 2.5, scale = 1) +
  draw_image(ode, y = 1.5,  scale = 1) +
  draw_image(lep, y = 0.5,scale = 1)

b <- ggdraw(insert_yaxis_grob(p2, pimage, position = "left"))

## combine together

cp <- cowplot::plot_grid(a,b, labels = c("A", "B"), ncol = 1, nrow = 2)
cp

ggsave(cp, filename = "Figures/studySpp.png",
       width = 4, height = 6)