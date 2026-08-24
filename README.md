# Coral reef habitat complexity decreases and diversifies local availability of light

Data and code used for analysis of light niches in Lizard Island - Brambilla et al. 2026 (Limnology and Oceanography). Figure also provided.

## File descriptions
### data
Empty as raw data (DEMs and 903 hobo readings) are too heavy. Output contains aggregated environmental data and geometrical variables extracted from DEMs with the R scripts provided. For access to raw data, please contact me.

### output
- **env_int.csv**: light readings and environmental data (PAR and tide) summarized per unit
- **env_long**: environmental data in long format (all the readings are maintained here, code cleans it)
- **geom_new**: table of complexity values at different scales
- **df_int_tot**: environmental (proportion of light integrals reaching substratum) and complexity data in long format (multiscale)
- **df_tot**: environmental and complexity data per unit (multiscale)

### R
- **read_environment.R**: loops to read the hobos recordings, and add par and tides data.
- **geometry.R**: get complexity variables (multiscale) 
- **local_variables.R**: code for merging data and model fitting + create merged outputs
- **figures.R**: code for figures not produced during model fitting

### figures
- **Figure 1**: Study sites and sampling scheme example
- **Figure 2**: Available light and complexity variables distributions
- **Figure 3**: Proportion of available light distributions across the structural complexity space measured for patch scales of S=0.5m and models results across scales
