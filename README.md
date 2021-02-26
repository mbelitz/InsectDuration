# InsectDuration
Code and data used to generate results and figures of "Climate drivers of adult insect activity are conditioned by life history traits" manuscript.
If you download or clone this repository, you can fully replicate LMM and PLMM results and figures. We also release the cleaned occurrence records datasets that were used to estimate phenology metrics.

## A brief overview of what is found in each directory and subdirectory

## data
Contains data that was used to estimate phenology metrics and data used in modeling framework.

### LMM_data
Data used in linear mixed model and phylogeneitc linear mixed modeling framework. 
- phenoTraits_data.csv have all estimated phenology metrics for species with trait data.
- phenoTraits_data_withNumObsData.csv has climate variables in original values instead of scaled values to figure out what 1 s.d. of a variable was equal to.
- phnoTraits_data_withNumObsData.csv has number of Obs term added to dataset. Note this is the dataset used in all mixed models.

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
Subdirectory to save figures to. Contains raw .pngs of figures used in manuscript included Supporting Information figures. 



