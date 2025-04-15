# Light Niches in Lizard Island
Data and code used for analysis of light niches in Lizard Island.

## File descriptions
### data
Empty as raw data (DEMs and 903 hobo readings) were too heavy. Output contains aggregated environmental data and geometrical variables extracted from DEMs with the R scripts provided. For access to raw data, please contact me.

### output
- **env_int.csv**: light readings and environmental data (PAR and tide) summarized per unit
- **env_long**: environmental data in long format (all the readings are maintained here)
- **geom_new**: table of complexity values
- **df_int_tot**: environmental (proportion of light integrals reaching substratum) and complexity data in long format (multiscale)
- **df_tot**: environmental and complexity data per unit (multiscale)

### R
- **read_environment.R**: loops to read the hobos, and add par and tides data.
- **geometry.R**: get complexity variables (multiscale) 
- **local_variables.R**: code for merging data and model fitting an outputs
- **figures.R**: code for figures not produced during model fitting
  
