# Quantifying bleaching severity across taxa, colony sizes, and environmental gradients during the 2022 mild bleaching event that affected the central Great Barrier Reef
---
This repository contains code and data needed to reproduce the article:

Álvarez-Noriega M<sup>1</sup>, Aston E<sup>1</sup>, Becker M<sup>1</sup>, Fabricius KE<sup>1</sup>, Figueira W<sup>2</sup>, Gordon S<sup>1</sup>, Krensel R<sup>2</sup>, Lechene M<sup>1,3</sup>, Remmers T<sup>1,3</sup>, Toor M<sup>1</sup>,  Ferrari R<sup>1</sup> *Challenging paradigms around the role of colony size, taxa, and environment on bleaching susceptibility.* Global Change Biology. <br> <br> 
Institutional Affiliations: <br>
 <sup>1</sup> Australian Institute of Marine Science, Townsville, QLD 4810, Australia <br>
 <sup>2</sup>  School of Life and Environmental Sciences, University of Sydney, Sydney, NSW 2006, Australia <br>
 <sup>3</sup>  College of Science and Engineering, James Cook University, 1 James Cook Dr, Douglas 
<br>



### Description of the file structure ##

There are three folders: **Data** (contains the raw data), **Code** (contains an R script that fits the 
cumulative link mixed models -CLMMs- and generates the plots), and **Outputs** (contains a copy of the plots).


### Description of the data ##
The Data folder has one file (*bleaching_data.csv*), that includes the following variables:
* Reef- unique identification for each reef:
    * OCCH - Chicken Reef,
    * OCDA- Davies Reef,
    * OCLB- Little Broadhurst Reef,
    * PA -Palm Islands
* Exposure- prevailing wind exposure:
    * BA- back (leeward, protected)
    * FR- front (windward, exposed)
    * FL- flank (intermediate exposure)
    * LA- lagoon
* Zone - depth zone (S = shallow, 3-8 m; D = deep, 9-15 m)
* Site - unique identification for each site (replicate sites per reef). Labelled as: Reef_Exposure/ReplicateNumber/Zone
* Plot - unique identification for each 12x6m plot (four replicates per site). Labelled as: Site_PlotReplicate
* TL_id - colony identification number (unique for each colony within plot)
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
* lat - site latitude
* long - site longitude
* ubedmean = horizontal water velocity at the bed (m/s). Extracted from Callaghan et al 2023-  https://espace.library.uq.edu.au/view/UQ:8246441 .
* Date_2022 = date of the initial survey in January 2022
* Date_bleaching = date of the bleaching survey in March - April 2022
* DHW - degree heating weeks (°C-weeks). Extracted from NOAA's daily 5km satellite dervied estimates: https://coralreefwatch.noaa.gov/product/5km/index_5km_dhw.php.
* rugosity_plane: rugosity at the site, extracted from the reconstructed 3D models
* fd_fugosity_plane: fractal dimension at the site, extracted from the reconstructed 3D models
* slope: the slope of the site (in degrees), extracted from the reconstructed 3D models
* Cross_shelf: cross shelf position (inshore or offshore)
* log_area : colony area in (m<sup>2</sup>) (log<sub>10</sub>)

### Description of the code ##
The Code folder has the script (*Analyses_CLMM.R*) that fits the CLMMs and generates the plots in the manuscript.

### Description of the outputs ##
The Outputs folder includes the figures generated by the script *Analyses_CLMM.R*.
   * **Figure2.png**: Proportion of observed colonies in each bleaching category depending on taxa (panel a), and cross shelf gradient, depth, and wave exposure (panels b and c). We only had lagoon sites in the shallow areas of the offshore reefs.
   * **Figure3.png**: Fitted model predictions of the probability of a colony presenting different levels of bleaching severity as a function of fractal dimension (panel a), degree heating weeks (DHWs) (panel b), slope of the reef (panel c), and water velocity at bed in inshore (panel d) and offshore reefs (panel e) in shallow sites. The coloured surfaces indicate the bleaching score predicted by the CLMM model for the given variables, and the points underneath the plot show the raw data (n = 3559). Each colour represents a bleaching score, reflecting different levels of bleaching severity. 
   * **Figure4.png**: Fitted model predictions of the probability of a colony presenting different levels of bleaching severity as a function of colony size for the different taxa in shallow sites. The coloured surfaces indicate the bleaching score predicted by the CLMM model for the given colony size, and the points underneath the plot show the raw data (n = 3559). Each colour represents a bleaching score, reflecting different levels of bleaching severity. Colony sizes are shown in log-scale. If assuming a circular shape, 10-3m2 would be approximately equivalent to a colony with a 4 cm in diameter, 10-2m2 to 11 cm in diameter, 10-1m2 to 35 cm in diameter, 1 m2 to 113 cm in diameter, and 10 m2 to 356 cm in diameter. 
   * **Figure5.png**: Fitted model predictions of the probability of a colony presenting different levels of bleaching severity as a function of depth and taxa. The colours indicate the bleaching score predicted by the CLMM model. Only taxa with at least 40 colonies per cross shelf gradient and depth were included. Only reefs with shallow and deep sites were included. Points underneath the plot show the raw data (n = 2056).

