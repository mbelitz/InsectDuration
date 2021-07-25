# InsectDuration
Code and data used to generate results and figures of "Climate drivers of adult insect activity are conditioned by life history traits" manuscript.
If you download or clone this repository, you can fully replicate LMM and PLMM results and figures. We also release the cleaned occurrence records datasets that were used to estimate phenology metrics. Scored trait data can be downloaded in the data/traits subdirectory. 

## A brief overview of what is found in each directory and subdirectory

## data
Contains data that was used to estimate phenology metrics and data used in modeling framework.

### LMM_data
Data used in linear mixed model and phylogeneitc linear mixed modeling framework. 
- phenoTraits_data.csv have all estimated phenology metrics for species with trait data.
- phenoTraits_data_withNumObsData.csv has climate variables in original values instead of scaled values to figure out what 1 s.d. of a variable was equal to.
- phenoTraits_data_withNumObsData.csv has number of Obs term added to dataset. Note this is the dataset used in all mixed models.
- phenoTraits_data_withNumObsData_updatedSeasTraits.csv has updates seasonality trait based on Hopkins-corrected methodology.

### occurrence records
Occurrence records used to estimate phenometrics. These data were cleaned from original data sources to remove records that were not in our species list determined by needing 1000 iNaturalist observations within our 6 groups. We also removed records that did not have eventDate, decimalLongitude, or decimalLatitude information.
- cleaned_occurrences.zip is a zip file with the occurrence records used to make phenology estimates for cicadas, Hymenoptera, Coleoptera, Odonata, and Diptera. 
- adult_leps.zip is a zip file with butterfly iNaturalist records annotated as "Adult". We did only included observations annotated as "Adult" because many Lepidoptera observations are of caterpillars.

### phenesse_climate_join
This subdirectory contains a single csv with all species with enough data to estimate phenology. This includes records of species that were removed because they were migratory, do not overwinter, were eusocial, or did not have trait information.

### phylogeny
This subdirectory contains a single .tre file that is a subtree representing the insects with traits used in our analysis. 

### silhouttes
Subdirectory contianing .pngs of insect silhouttes used to generate Figure 1. Silhouttes taken from phylopic. See manuscript for attribution.

### traits
Subdirectory containing trait information for species list, including the species that did not have completed trait information.

## Figures
Subdirectory to save figures to. Contains raw .pngs of figures used in manuscript included Supporting Information figures. Figure 1 and 2 of main text remained the same, so aare found in main Figures directory. The updated figures made as part of revising the manuscript due to reviewer and editor comments can be found in the resubmission subdirectory.

## scripts
Contains R scripts used to generate phenology estimates, run analyses, and generate figures. Contains three subdirectories for these three aims.

### analyses
Code to run data analyses (linear mixed models (LMM) and phylogenetic linear mixed models (PLMM))
- spatial_PLMM.R runs script to generate onset LMM and onset PLMM.
- spatial_PLMM_offset.R runs script to generate offset LMM, offset PLMM, and offset PLMM with spatial term.
- spatial_PLMM_duration.R runs script to generate duration LMM, duration PLMM, and duration PLMM with spatial term.
- seasonalityTrait_byHopkinsLaw.R script to generate Hopkin's-corrected seasonality trait as described in Methods.

### figures
Code to generate figures used in manuscript.
- model_assumptions.R runs script to generate SI Fig. 2
- MoransI_Fig.R script to generate SI Fig. 1
- R2_Fig.R runs script to generate SI Fig. 4
- RandomEffects_Figure.r runs script to generate SI Fig. 3
- interaction_plots_duration.R script to generate Fig. 5
- interaction_plots_offset.R script to generate Fig. 4
- interaction_plots_onset.R script to generate Fig. 3
- studyArea.R holds code to generate Fig. 2
- studySpp.R has the code to generate Fig. 1

### phenoestimates
Code to generate phenology estimates. NOTE THAT THESE SCRIPTS WILL NOT RUN IF YOU SOURCE THEM. In order to run these scripts without a bug, you will need to unzip the data in the occurrence_records subdirectory and change the paths in the scripts to the correct locations. SECOND NOTE: Running this script would take a while, due to bootstrapping to get confidence intervals around the phenology estimate. 
- 00_pkg_functions.R holds code for functions used in the following scripts.
- combine_estimates_25km.R script that combines phenology estimate outputs into a single dataframe.
- running_phenesse_25km.R script that generates phenology estimates through multiple taxon groups. 

## Tables
Subdirectory to save Table result outputs. Contains manipulated .xlsx file with all LMM results in single location. PGLMM results are saved as .html files and split up by fixed effects (FE) and random effects (RE) results.


