# Light Niches in Lizard Island
Data and code used for analysis of light niches in Lizard

## File descriptions
### Data
Empty as raw data (DEMs and 900+hobo readings) were too heavy. See outputs for variables extracted (env and geom)

### output
- **env_int.csv**: light readings and environmental data (PAR and tide) summarized per unit
- **env_long**: environmental data in long format (all the readings are maintained here,)
- **geom_new**: table of complexity values with the new clipped R and D
- **df_int_tot**: environmental and complexity data in long format (3 scales)
- **df_tot**: environmental and complexity data per unit (with attenuated variables, 3 scales)

### R
- **read_environment.R**: loops to read the hobos, and merge in par and tides.
- **geometry.R**: get complexity variables (new Rclipped, 3 scales) 
- **fit_models.R**: handle data to make tables for model fitting
- **play.R**:bits of code for quick graphs and checks I used throughout the analysis. Very random...

