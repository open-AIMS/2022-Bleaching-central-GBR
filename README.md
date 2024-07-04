# Coral bleaching severity #
---
Script that fits and plots generalised linear mixed effects models to quantify bleaching severity across taxa, colony sizes, and environmental gradients during the 
2022 mild bleaching event that affected the central Great Barrier Reef.


## Description of the file structure ##

There are three folders: **Data** (contains the raw data), **Code** (contains an R script that fits the 
generalised linear mixed effects models and generates the plots), and **Outputs** (contains a copy of the plots).


### Description of the data ##
The Data folder has one file (*bleaching_data.csv*), that includes the following variables:

* Plot - unique identification for each 12x6m plot (four replicates per site)
* Site - unique identification for each site (replicate sites per reef)
* Reef- unique identification for each reef:
    * OCCH - Chicken Reef,
    * OCDA- Davies Reef,
    * OCLB- Little Broadhurst Reef,
    * PA -Palm Islands
* TL_id - colony identification at the plot level
* Taxa - taxonomic idenfitication of each colony:
    * Adgt - *Acropora* digitate spp.,
    * Pmas - *Porites* massive spp.,
    * Pver - *Pocillopora verrucosa*,
    * Acor - *Acropora* corymbose spp.,
    * Lobo_sp - *Lobophyllia* spp.,
    * Menc- *Montipora* encrusting spp.,
    * Gspp- *Goniastrea* spp.,
    * Spis - *Stylophora pistillata*,
    * Shys- *Seriatopora hystrix*,
    * Abra- *Acropora* branching spp.,
    * Pdam - *Pocillopora damicornis*,
    * Plat - *Platygyra* spp.,
    * Atab - *Acropora* tabular spp.
* Bleaching - bleaching score:
    * 1 - no bleaching,
    * 2 - 1-50% of the colony area was bleached,
    * 3- 51-95% of the colony area was bleached,
    * 4 - 96-100% of the colony area was bleached,
    * 5 - partial mortality,
    * 6- whole colony mortality.
* Bleaching_res - bleaching response (bleached = 1, not bleached = 0)
* Bleaching_minor - whether the colony had minor bleaching (1-50% of the colony) =1, or no bleaching =0. NAs indicate more severe bleaching levels.
* Bleaching_major - whether the colony had major bleaching (51-95% of the colony) = 1, or minor to no bleaching = 0. NAs indicate more severe bleaching levels.
* Bleaching_fully - whether the colony was fully bleached (>96% of the colony, partial mortality, or fully dead) =1, or major to no bleaching = 0.
* lat - site latitude
* long - site longitude
* Zone - depth zone (S = shallow, 3-8 m; D = deep, 9-15 m)
* ubedmean = horizontal water velocity at the bed (m/s). Extracted from Callaghan et al 2023-  https://espace.library.uq.edu.au/view/UQ:8246441 .
* Date_2022 = date of the initial survey in January 2022
* Date_bleaching = date of the bleaching survey in March - April 2022
* DHW - degree heating weeks. Extracted from NOAA's daily 5km satellite dervied estimates: https://coralreefwatch.noaa.gov/product/5km/index_5km_dhw.php.
* rugosity_plane: rugosity at the site, extracted from the reconstructed 3D models
* fd_fugosity_plane: fractal dimension at the site, extracted from the reconstructed 3D models
* slope: the slope of the site (in degrees), extracted from the reconstructed 3D models
* Cross_shelf: cross shelf position (inshore or offshore)
* log_area : colony area in m2 (log10)
