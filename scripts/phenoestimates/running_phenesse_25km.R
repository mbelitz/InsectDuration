library(data.table)
library(dplyr)
library(sf)

setwd("/srv/gspeedq32/mbelitz/insect-pheno-duration/")

source("scripts/00_pkg_functions_V2.R")

#### import odonata data

odes <- fread("data/cleaned_data/cleaned_odes.csv")

ode_list <- list("Pachydiplax longipennis",
                 "Erythemis simplicicollis",
                 "Plathemis lydia",
                 "Libellula luctuosa",
                 "Calopteryx maculata",
                 "Perithemis tenera",
                 "Ischnura verticalis",
                 "Ischnura posita",
                 "Sympetrum corruptum",
                 "Sympetrum vicinum",
                 "Celithemis eponina",
                 "Libellula incesta",
                 "Libellula pulchella",
                 "Enallagma civile",
                 "Argia moesta",
                 "Argia apicalis",
                 "Ischnura ramburii",
                 "Libellula vibrans",
                 "Orthemis ferruginea",
                 "Argia fumipennis",
                 "Libellula saturata",
                 "Tramea lacerata",
                 "Argia sedula",
                 "Celithemis elisa",
                 "Argia vivida",
                 "Enallagma basidens",
                 "Ladona julia",
                 "Ischnura cervula",
                 "Enallagma signatum",
                 "Tramea onusta",
                 "Ischnura hastata",
                 "Pantala flavescens",
                 "Brachymesia gravida",
                 "Lestes rectangularis",
                 "Libellula quadrimaculata",
                 "Aeshna umbrosa",
                 "Argia translata",
                 "Libellula croceipennis",
                 "Enallagma exsulans",
                 "Rhionaeschna multicolor",
                 "Enallagma geminatum",
                 "Erythrodiplax umbrata",
                 "Libellula needhami",
                 "Sympetrum obtrusum",
                 "Sympetrum illotum",
                 "Ladona deplanata",
                 "Sympetrum semicinctum",
                 "Archilestes grandis",
                 "Dythemis velox",
                 "Enallagma aspersum",
                 "Argia immunda",
                 "Erythrodiplax berenice",
                 "Epiaeschna heros",
                 "Telebasis salva",
                 "Tramea carolina",
                 "Lestes congener",
                 "Dromogomphus spinosus",
                 "Aeshna canadensis",
                 "Libellula cyanea",
                 "Libellula semifasciata",
                 "Sympetrum ambiguum",
                 "Libellula forensis",
                 "Pantala hymenaea",
                 "Argia plana",
                 "Erythemis collocata",
                 "Enallagma carunculatum",
                 "Phanogomphus exilis",
                 "Dythemis fugax",
                 "Erythemis vesiculosa",
                 "Aeshna interrupta",
                 "Hetaerina titia",
                 "Sympetrum pallipes",
                 "Dythemis nigrescens",
                 "Erpetogomphus designatus",
                 "Celithemis fasciata",
                 "Brachymesia furcata",
                 "Cordulia shurtleffii",
                 "Dorocordulia libera")

ode_fun <- function(x){
  
  df <- dplyr::filter(odes, binomial == x)
  
  gridz <- tryCatch(plt_summary(cell_size = 25000, dat = df, n_per_cell = 10),
                    error = function(e) message("Not enough data, Grid failed"))
  
  duration <- tryCatch(run_phenesse(df = gridz$dat_to_use, minimum_obs = 10, earliest_year = 2015, 
                                    last_year = 2019,n_item = 50000, onset_perct = 0.05, 
                                    num_cores = 30, offset_perct = 0.95),
                       error = function(e) message("Not enough data, phenesse failed"))
  
  output <- tryCatch(duration %>% 
                       mutate(scientificName = x) %>% 
                       mutate(Order = "Odonata"), 
                     error = function(e) message("Not enough data, skipping step"))
  
  tryCatch(write.csv(x = output, file = paste("phenesse_outputs_25km_V2/odonata/", x, ".csv", sep = ""), row.names = FALSE),
           error = function(e) message("Error: no data can't write file"))
}


lapply(ode_list, FUN = ode_fun)

### Bee time

hy <- fread("data/cleaned_data/cleaned_hymenoptera.csv")

hy_list <- list(
  "Bombus ternarius",
  "Bombus bimaculatus",
  "Halictus ligatus",
  "Bombus melanopygus",
  "Dasymutilla occidentalis",
  "Monobia quadridens",
  "Bombus sonorus",
  "Augochlora pura",
  "Xylocopa micans",
  "Xylocopa tabaniformis",
  "Melissodes bimaculatus",
  "Agapostemon virescens",
  "Xylocopa varipuncta",
  "Bombus huntii",
  "Bombus fervidus",
  "Bombus borealis",
  "Apis mellifera",
  "Bombus impatiens",
  "Xylocopa virginica",
  "Bombus griseocollis",
  "Bombus pensylvanicus",
  "Polistes fuscatus",
  "Dolichovespula maculata",
  "Polistes dominula",
  "Polistes exclamans",
  "Sphex ichneumoneus",
  "Sphecius speciosus",
  "Vespula maculifrons",
  "Scolia dubia",
  "Polistes metricus",
  "Eumenes fraternus",
  "Vespula squamosa",
  "Vespula pensylvanica",
  "Vespa crabro",
  "Polistes dorsalis",
  "Pelecinus polyturator",
  "Sphex pensylvanicus",
  "Dolichovespula arenaria",
  "Bicyrtes quadrifasciatus"
)

hy_fun <- function(x){
  
  df <- dplyr::filter(hy, binomial == x)
  
  gridz <- tryCatch(plt_summary(cell_size = 25000, dat = df, n_per_cell = 10),
                    error = function(e) message("Not enough data, Grid failed"))
  
  duration <- tryCatch(run_phenesse(df = gridz$dat_to_use, minimum_obs = 10, earliest_year = 2015, 
                                    last_year = 2019,n_item = 50000, onset_perct = 0.05, 
                                    num_cores = 30, offset_perct = 0.95),
                       error = function(e) message("Not enough data, phenesse failed"))
  
  output <- tryCatch(duration %>% 
                       mutate(scientificName = x) %>% 
                       mutate(Order = "Hymenoptera"), 
                     error = function(e) message("Not enough data, skipping step"))
  
  tryCatch(write.csv(x = output, file = paste("phenesse_outputs_25km_V2/hymenoptera/", x, ".csv", sep = ""), row.names = FALSE),
           error = function(e) message("Error: no data can't write file"))
}


lapply(hy_list, FUN = hy_fun)

### Beetles make the world go round -- potential problems here. come back and examine data

cole <- fread("data/cleaned_data/cleaned_coleoptera.csv")

cole_list <- list(
  "Diabrotica undecimpunctata",
  "Popillia japonica",
  "Cicindela sexguttata",
  "Tetraopes tetrophthalmus",
  "Chauliognathus pensylvanicus",
  "Monochamus scutellatus",
  "Chauliognathus marginatus",
  "Ellychnia corrusca",
  "Pelidnota punctata",
  "Exomala orientalis",
  "Chrysochus auratus",
  "Megacyllene robiniae",
  "Lucanus capreolus",
  "Cicindela punctulata",
  "Photinus pyralis",
  "Cicindela repanda",
  "Typocerus velutinus",
  "Calosoma scrutator",
  "Orthosoma brunneum")

cole_fun <- function(x){
  
  df <- dplyr::filter(cole, binomial == x)
  
  gridz <- tryCatch(plt_summary(cell_size = 25000, dat = df, n_per_cell = 10),
                    error = function(e) message("Not enough data, Grid failed"))
  
  duration <- tryCatch(run_phenesse(df = gridz$dat_to_use, minimum_obs = 10, earliest_year = 2015, 
                                    last_year = 2019,n_item = 50000, onset_perct = 0.05, 
                                    num_cores = 1, offset_perct = 0.95),
                       error = function(e) message("Not enough data, phenesse failed"))
  
  output <- tryCatch(duration %>% 
                       mutate(scientificName = x) %>% 
                       mutate(Order = "Coleoptera"), 
                     error = function(e) message("Not enough data, skipping step"))
  
  tryCatch(write.csv(x = output, file = paste("phenesse_outputs_25km_V2/coleoptera/", x, ".csv", sep = ""), row.names = FALSE),
           error = function(e) message("Error: no data can't write file"))
}


lapply(cole_list, FUN = cole_fun)

### Flys or flies? - also potential problems here? check output csvs

dip <- fread("data/cleaned_data/cleaned_diptera.csv")

dip_list <- list(
  "Toxomerus marginatus",
  "Toxomerus geminatus",
  "Eristalis transversa",
  "Allograpta obliqua",
  "Eristalis tenax",
  "Delphinia picta",
  "Xenox tigrinus",
  "Plecia nearctica",
  "Toxomerus politus",
  "Asphondylia auripila",
  "Chrysopilus thoracicus",
  "Clogmia albipunctata",
  "Eurosta solidaginis",
  "Eristalis dimidiata",
  "Promachus hinei",
  "Tabanus atratus",
  "Helophilus fasciatus",
  "Bombylius major"
)

dip_fun <- function(x){
  
  df <- dplyr::filter(dip, binomial == x)
  
  gridz <- tryCatch(plt_summary(cell_size = 25000, dat = df, n_per_cell = 10),
                    error = function(e) message("Not enough data, Grid failed"))
  
  duration <- tryCatch(run_phenesse(df = gridz$dat_to_use, minimum_obs = 10, earliest_year = 2015, 
                                    last_year = 2019,n_item = 50000, onset_perct = 0.05, 
                                    num_cores = 30, offset_perct = 0.95),
                       error = function(e) message("Not enough data, phenesse failed"))
  
  output <- tryCatch(duration %>% 
                       mutate(scientificName = x) %>% 
                       mutate(Order = "Diptera"), 
                     error = function(e) message("Not enough data, skipping step"))
  
  tryCatch(write.csv(x = output, file = paste("phenesse_outputs_25km_V2/diptera/", x, ".csv", sep = ""), row.names = FALSE),
           error = function(e) message("Error: no data can't write file"))
}


lapply(dip_list, FUN = dip_fun)

#### import Lep

insect_adult1 <- fread("data/vijay_adults/Insect_adult.csv")
insect_adult2 <- fread("data/vijay_adults/Insect_adult1.csv")
insect_adult3 <- fread("data/vijay_adults/insect_adult2.csv")

tots_adults <- rbind(insect_adult1, insect_adult2, insect_adult3) %>% 
  mutate(binomial = scientificName) %>% 
  dplyr::select(decimalLongitude, decimalLatitude, everything())

### Lep time

lep_list <- list(
  "Danaus plexippus",
  "Vanessa cardui",
  "Junonia coenia",
  "Papilio glaucus",
  "Agraulis vanillae",
  "Vanessa atalanta",
  "Papilio polyxenes",
  "Limenitis arthemis",
  "Pieris rapae",
  "Hylephila phyleus",
  "Phyciodes tharos",
  "Battus philenor",
  "Pyrrharctia isabella",
  "Strymon melinus",
  "Epargyreus clarus",
  "Danaus gilippus",
  "Hyles lineata",
  "Vanessa virginiensis",
  "Atteva aurea",
  "Antheraea polyphemus",
  "Nymphalis antiopa",
  "Euptoieta claudia",
  "Actias luna",
  "Polygonia interrogationis",
  "Halysidota tessellaris",
  "Cupido comyntas",
  "Colias eurytheme",
  "Limenitis archippus",
  "Hyphantria cunea",
  "Lophocampa caryae",
  "Asterocampa celtis",
  "Atalopedes campestris",
  "Papilio troilus",
  "Libytheana carinenta",
  "Estigmene acrea",
  "Euchaetes egle",
  "Phoebis sennae",
  "Hypena scabra",
  "Papilio cresphontes",
  "Pontia protodice",
  "Hypercompe scribonia",
  "Acronicta americana",
  "Eacles imperialis",
  "Speyeria cybele",
  "Malacosoma americana",
  "Polygonia comma",
  "Anartia jatrophae",
  "Automeris io",
  "Orgyia leucostigma",
  "Heliconius charithonia",
  "Coenonympha tullia",
  "Lymantria dispar",
  "Malacosoma disstria",
  "Mythimna unipuncta",
  "Abaeis nicippe",
  "Chlosyne lacinia",
  "Lerema accius",
  "Nathalis iole",
  "Polites peckius",
  "Hemaris diffinis",
  "Papilio rutulus",
  "Poanes zabulon",
  "Dryocampa rubicunda",
  "Papilio zelicaon",
  "Manduca sexta",
  "Asterocampa clyton",
  "Spilosoma virginica",
  "Ascalapha odorata",
  "Pyrisitia lisa",
  "Phyciodes phaon",
  "Cercyonis pegala",
  "Cisseps fulvicollis",
  "Megisto cymela",
  "Colias philodice",
  "Hemiargus ceraunus",
  "Urbanus proteus",
  "Udea rubigalis",
  "Echinargus isola",
  "Lophocampa maculata",
  "Hemaris thysbe",
  "Leptotes marina",
  "Campaea perlata",
  "Galgula partita",
  "Prochoerodes lineola",
  "Noctua pronuba",
  "Burnsius communis",
  "Anartia fatima",
  "Ancyloxypha numitor",
  "Microcrambus elegans",
  "Spodoptera ornithogalli",
  "Costaconvexa centrostrigaria",
  "Brephidium exilis",
  "Nomophila nearctica",
  "Ctenucha virginica",
  "Lethe anthedon",
  "Nymphalis californica",
  "Calycopis cecrops",
  "Chlosyne nycteis",
  "Limenitis lorquini",
  "Papilio multicaudata",
  "Papilio rumiko",
  "Euphydryas chalcedona",
  "Nadata gibbosa",
  "Haematopis grataria",
  "Eurytides marcellus",
  "Celastrina echo",
  "Idia aemula",
  "Scopula limboundata",
  "Poanes melane",
  "Icaricia acmon",
  "Erynnis horatius",
  "Ceratomia undulosa",
  "Euphyes vestris",
  "Glaucopsyche lygdamus",
  "Adelpha californica",
  "Idia americalis",
  "Helicoverpa zea",
  "Darapsa myron",
  "Hypoprepia fucosa",
  "Paonias excaecata",
  "Zerene cesonia",
  "Hyalophora cecropia",
  "Ascia monuste",
  "Xylophanes tersa",
  "Alypia octomaculata",
  "Burnsius oileus",
  "Celastrina neglecta",
  "Anthanassa texana",
  "Dryas iulia",
  "Vanessa annabella",
  "Megalopyge opercularis",
  "Phyciodes cocyta",
  "Callophrys gryneus",
  "Macaria pustularia",
  "Phyciodes mylitta",
  "Siproeta stelenes",
  "Erynnis funeralis",
  "Papilio canadensis",
  "Hermeuptychia sosybius",
  "Amorpha juglandis",
  "Poanes hobomok",
  "Eutrapela clemataria",
  "Ochlodes sylvanoides",
  "Thymelicus lineola",
  "Eumorpha pandorus",
  "Orthonama obstipata",
  "Epimecis hortaria",
  "Marimatha nigrofimbria",
  "Palthis asopialis",
  "Choristoneura rosaceana",
  "Pleuroprucha insulsaria",
  "Papilio eurymedon",
  "Spoladea recurvalis",
  "Chlorochlamys chloroleucaria",
  "Lacinipolia renigera",
  "Iridopsis defectaria",
  "Dione moneta",
  "Palthis angulalis",
  "Eudryas grata",
  "Panoquina ocola",
  "Agrotis ipsilon",
  "Halysidota harrisii",
  "Crambus agitatellus",
  "Erynnis baptisiae",
  "Zale lunata",
  "Phoebis agarithe",
  "Nematocampa resistaria",
  "Plutella xylostella",
  "Urola nivalis",
  "Hypsopygia olinalis",
  "Lerodea eufala",
  "Lycaena phlaeas",
  "Papilio palamedes",
  "Mestra amymone",
  "Anaea andria",
  "Autographa precationis",
  "Pseudeustrotia carneola",
  "Anicla infecta",
  "Amphipyra pyramidoides",
  "Eusarca confusaria",
  "Paonias myops",
  "Leuconycta diphteroides",
  "Plusiodonta compressipalpis",
  "Hypena baltimoralis",
  "Sunira bicolorago",
  "Hymenia perspectalis",
  "Melipotis indomita",
  "Lambdina fiscellaria",
  "Erynnis juvenalis",
  "Rivula propinqualis",
  "Anticarsia gemmatalis",
  "Euptoieta hegesia",
  "Uresiphita reversalis",
  "Urbanus dorantes",
  "Aglais milberti",
  "Cydia latiferreana",
  "Anthocharis sara",
  "Patalene olyzonaria",
  "Atlides halesus",
  "Chrysodeixis includens",
  "Drepana arcuata",
  "Ectropis crepuscularia",
  "Argyrotaenia velutinana",
  "Ennomos magnaria",
  "Thorybes pylades",
  "Manduca rustica",
  "Strymon istapa",
  "Acrolophus popeanella",
  "Satyrium calanus",
  "Eurema daira",
  "Maliattha synochitis",
  "Cycnia tenera",
  "Parallelia bistriaris",
  "Haploa clymene",
  "Clepsis peritana",
  "Hyles gallii",
  "Elophila obliteralis",
  "Acharia stimulea",
  "Protodeltote muscosula",
  "Leptophobia aripa",
  "Icaricia icarioides",
  "Gluphisia septentrionis",
  "Leptotes cassius",
  "Phyllodesma americana",
  "Sphinx kalmiae",
  "Nephelodes minians",
  "Caenurgina erechtea",
  "Synchlora aerata",
  "Marpesia petreus",
  "Schizura unicornis",
  "Chioides albofasciatus",
  "Euclea delphinii",
  "Apatelodes torrefacta",
  "Citheronia regalis",
  "Pachysphinx modesta",
  "Cucullia convexipennis",
  "Synchlora frondaria",
  "Lascoria ambigualis",
  "Calycopis isobeon",
  "Polygrammate hebraeicum",
  "Eudryas unio",
  "Euchromius ocellea",
  "Pyralis farinalis",
  "Anatrytone logan",
  "Hyalophora euryalus",
  "Polites themistocles",
  "Spodoptera frugiperda",
  "Callophrys augustinus",
  "Erynnis tristis",
  "Ochropleura implecta",
  "Achyra rantalis",
  "Euphydryas phaeton",
  "Catocala ilia",
  "Chlosyne theona",
  "Apodemia virgulti",
  "Pieris marginalis",
  "Thyridopteryx ephemeraeformis",
  "Boloria bellona",
  "Catocala relicta",
  "Utetheisa ornatrix",
  "Acronicta oblinita",
  "Lycaena hyllus",
  "Panopoda rufimargo",
  "Dione juno",
  "Asterocampa leilia",
  "Trichodezia albovittata",
  "Orgyia definita",
  "Eumorpha fasciatus",
  "Acronicta insularis",
  "Anavitrinella pampinaria",
  "Kricogonia lyside",
  "Idaea dimidiata",
  "Celastrina lucia",
  "Peridroma saucia",
  "Idia lubricalis",
  "Protoboarmia porcelaria",
  "Raphia frater",
  "Speyeria atlantis",
  "Chlosyne janais"
)

lep_fun <- function(x){
  
  df <- dplyr::filter(tots_adults, binomial == x)
  
  gridz <- tryCatch(plt_summary(cell_size = 25000, dat = df, n_per_cell = 10),
                    error = function(e) message("Not enough data, Grid failed"))
  
  duration <- tryCatch(run_phenesse(df = gridz$dat_to_use, minimum_obs = 10, earliest_year = 2015, 
                                    last_year = 2019,n_item = 50000, onset_perct = 0.05, 
                                    num_cores = 30, offset_perct = 0.95),
                       error = function(e) message("Not enough data, phenesse failed"))
  
  output <- tryCatch(duration %>% 
                       mutate(scientificName = x) %>% 
                       mutate(Order = "Lepidoptera"), 
                     error = function(e) message("Not enough data, skipping step"))
  
  tryCatch(write.csv(x = output, file = paste("phenesse_outputs_25km_V2/lepidoptera/", x, ".csv", sep = ""), row.names = FALSE),
           error = function(e) message("Error: no data can't write file"))
}


lapply(lep_list, FUN = lep_fun)

# Cicadas, last but not least -- except in comparable species ;)

cic <- fread("data/cleaned_data/cleaned_cicadas.csv")

cicadas_list <- list("Neotibicen superbus",
                     "Megatibicen resh")

cic_fun <- function(x){
  
  df <- dplyr::filter(cic, binomial == x)
  
  gridz <- tryCatch(plt_summary(cell_size = 25000, dat = df, n_per_cell = 10),
                    error = function(e) message("Not enough data, Grid failed"))
  
  duration <- tryCatch(run_phenesse(df = gridz$dat_to_use, minimum_obs = 10, earliest_year = 2015, 
                                    last_year = 2019,n_item = 50000, onset_perct = 0.05, 
                                    num_cores = 30, offset_perct = 0.95),
                       error = function(e) message("Not enough data, phenesse failed"))
  
  output <- tryCatch(duration %>% 
                       mutate(scientificName = x) %>% 
                       mutate(Order = "Hemiptera"), 
                     error = function(e) message("Not enough data, skipping step"))
  
  tryCatch(write.csv(x = output, file = paste("phenesse_outputs_25km_V2/hemiptera/", x, ".csv", sep = ""), row.names = FALSE),
           error = function(e) message("Error: no data can't write file"))
}


lapply(cicadas_list, FUN = cic_fun)
