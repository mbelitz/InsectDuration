Metadata for phenoTraits_Data_withNumObsData_updatedSeasTrat.csv, which was used to complete all models that were reported in the Manuscript.

Each column represents a new value or trait. Each row is a new year x cell x species combination that had enough observations to generate a phenology estimate. 

Column A: number of total columns
Column B: onset -- estimate of the 5% phenometric (unit is days)
Column C: onset_low -- low 95% CI value of the onset estimate (unit is days)
Column D: onset_high -- high 95% CI value of the onset estimate (unit is days)
Column E: year -- year of estimate
Column F: offset -- estimate of 95% phenometric (unit is days)
Column G: offset_low -- low 95% CI value of the offset estimate (unit is days)
Column H: offset_high -- high 95% CI value of the offset estimate (unit is days)
Column I: duration -- difference between offset and onset estimates (unit is days)
Column J: scientificName -- species with phenoestimate
Column K: Order -- Order of species with phenoestimate
Column L: id_cells -- ID of the cell with phenoestimate
Column M: lon -- x-coordinate value of cell centroid (unit is meters)
Column N: lat -- y-coordinate value of cell centroid (unit is meters)
Column O: OS -- Order:scientitifName of species with phenoestimate
Column P: prec -- scaled annual precipitation value for cell & year combination (unit is mm)
Column Q: temp -- scaled mean annual temperature value for cell & year combination (unit is degrees C)
Column R: bio4 -- unscaled temperature seasonality value for cell & year combination (unit is degrees C)
Column S: bio15 -- unscaled precipitation seasonality value for cell & year combination (unit is mm)
Column T: pop -- scaled human population density for cell 
Column U: numObs -- scaled number of observations for cell, year, species combination
Column V: prec_seas -- scaled values of column S
Column W: temp_seas -- scaled values of column R
Column X: higher_taxon -- order of species with phenoestimate
Column Y: flights -- voltinism categorization of species (univoltine or not univoltine)
Column Z: development -- development categorization of species (holometabolous or hemimetabolous)
Column AA: immature.habitat -- habitat of larvae of species (above ground, fresh water, underground)
Column AB: diapause.stage -- lifestage that species overwinters (Larvae, Pupae, Adult, Egg)
Column AC: larval.diet -- What does the species eat as a larvae? (carnivorous, live plants, detritivorous)
Column AD: seas -- season of emergence of species, not accounting for Hopkin's Law (Fall, Spring, Summer)
Column AE: elevation -- elevation of cell (unit is meters)
Column AF: lon_cor -- longitudinal distance from most southernly cell (unit is meters)
Column AG: lat_cor -- latitudinal distance from most southernly cell (unit is meters)
Column AH: elev_cor -- difference in elevation between cell and most southernly cell (unit is meters)
Column AI: hopkins_cor -- Number of days earlier or later that onset should occur in relation to our most southernly cell, given a Hopkin's Law correction (unit in days)
Aolumn AJ: hopkins_onset -- Hopkin's corrected onset for species-cell-year combination (unit is days)
mean_onset: Average Hopkin's corrected onset values for a particular species
sd_onset: standard deviation of onset values for a particular species
seas2: Hopkin's Law-corrected season categorization based on emergence of species (Fall, Spring, Summer) 


