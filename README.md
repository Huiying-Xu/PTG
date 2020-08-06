# Leaf traits data from Gongga Mountain, China
This is the git repository for plant trait and climate data used in the paper: Xu et al. (under review in Tree Physiology). The repository provides plant traits and climate data from 18 sites along a 3000-m elevation gradient in the Gongga Mountain region, China. 
## Project information
This project reports on leaf traits of plants and climate data at 18 sites from 1143 to 4361 m in Gongga Mountain region (29° 22' to 29° 55' N and 101° 1' to 102° 9' E). The trait data were collected during the active growing season in 2018 and 2019, as part of a Chinese research project. The vegetation type changes from deciduous broad-leaved forest dominated by Betulaceae, Urticaceae, Caprifoliaceae and Rosaceae, to evergreen needle-leaved forest and deciduous shrubland dominated by Pinaceae and/or Rosaceae and Ericaceae with increasing elevation.
## Data
- The data are provided in the form of three tables connected by ‘Site ID’. The tables are (1) site information, (2) trait values by species at each site, and (3) climate data for each site. Each file is a csv file
- he dataset is stored on **Leaf traits data in Gongga Mountain, China**: [Zenodo](link here)
### *Site Table*
- Site ID
- Longitude and latitude (decimal degrees)
- Elevation (m)
- Collection month and year
- Vegetation type
### *Leaf trait Table*
- Site ID
- Species name
- Leaf phenology
- leaf mass per unit area (Ma, gC m–2): standard protocols (Cornelissen et al., 2003)
- the maximum capacity of carboxylation at 25 ˚C (Vcmax25, μmolC m–2 s–1): leaf gas-exchange measurements (LI-6400) and one-point method (De Kauwe et al., 2016)
- the ratio of leaf-internal to ambient CO2 partial pressure (χ, unitless): estimated by carbon isotopic values (δ13C) using the method of Ubierna and Farquhar (2014)
- leaf nitrogen content per unit area (Narea, gN m–2): standard protocols (Cornelissen et al., 2003)
### *Climate Table*
The Simple Process-Led Algorithms for Simulating Habitats (SPLASH) model (Davis et al., 2017) was used to estimate the bioclimatic variables. The file contains the following fields:
- Site ID
- the growing-season mean value of temperature (Tg), where growing season is defined as the period when the daily temperature is above 0 ˚C
- mean daytime temperature in July (TdJ)
- growing-season mean vapour pressure deficit (D0), where growing season is defined as the period when the daily temperature is above 0 ˚C
- growing-season mean photosynthetically active radiation (R0), where growing season is defined as the period when the daily temperature is above 0 ˚C
- the ratio of growing-season length to the number of days in the year (f), where growing season is defined as the period when the daily temperature is above 0 ˚C
- the leaf-area-index weighted R0 (RLAI), where leaf area index was the mean value for July and August 2018 and 2019 
- the ratio of annual actual evapotranspiration to annual potential evapotranspiration (αp)
### *Raw climate data source*
- Monthly maximum and minimum temperature, fraction of sunshine hours, water vapor pressure and precipitation were interpolated to these 18 sites with data from [meteorological stations around Gongga Mountain region](http://data.cma.cn/data/cdcdetail/dataCode/SURF_CLI_CHN_MUL_MON.html) using ANUSPLIN (Hutchinson and Xu, 2004)
- Leaf area index from July and August in 2018 and 2019 was derived from the MODIS leaf area index product [MCD15A3H](https://modis.gsfc.nasa.gov/)
- All the climate data were the average of values from January 2015 to December 2017
## Code
- The data was processed and analysed using R. All code is stored on [github](https://github.com/Huiying-Xu/PTG/tree/master/photosynthesis)
- To download data, the code in `download.R` can be used. 
- The key code to produce predicted trait values in Xu et al. (under review in Tree Physiology) is in `predict.R`
## Reference
Cornelissen, J. H. C., Lavorel, S., Garnier, E. S., Buchmann, N., and Gurvich: A handbook of protocols for standardised and easy measurement of plant functional traits worldwide, Aust. J. Bot., 51, 335-380, 2003.
Ubierna, N., and Farquhar, G. D.: Advances in measurements and models of photosynthetic carbon isotope discrimination in C3 plants, Plant, Cell & Environment, 37, 1494-1498, 10.1111/pce.12346, 2014.
De Kauwe, M. G., Lin, Y. S., Wright, I. J., Medlyn, B. E., Crous, K. Y., Ellsworth, D. S., Maire, V., Prentice, I. C., Atkin, O. K., and Rogers, A.: A test of the 'one-point method' for estimating maximum carboxylation capacity from field-measured, light-saturated photosynthesis, New Phytol., 210, 1130-1144, 2016.
Hutchinson, M. F., and Xu, T.: Anusplin version 4.2 user guide, Centre for Resource and Environmental Studies, The Australian National University, Canberra, 54, 2004.
Davis, T. W., Prentice, I. C., Stocker, B. D., Thomas, R. T., Whitley, R. J., Wang, H., Evans, B. J., Gallego-Sala, A. V., Sykes, M. T., and Cramer, W.: Simple process-led algorithms for simulating habitats (SPLASH v.1.0): robust indices of radiation, evapotranspiration and plant-available moisture, Geosci. Model Dev., 10, 1-25, 2017.
