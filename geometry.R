# script to get all the explanatory variables and prepare tables for the analysis

#### libraries #####

library(rgl))
library(Rvcg)
library(raster)
library(rgdal)
library(stringr)
library(sp)
library(ggplot2)
library(tidyverse)

#### loop for light  ####

# TB17 is the crs reference
TB17dem <- raster("D:/Dropbox/My Dropbox/NC-RR_environment_data/dems/TB17.tif")

# read dems names
# rec are "RRYY.tif" format: RR=record site code, YY = year
recs <- dir("D:/Dropbox/My Dropbox/NC-RR_environment_data/dems")

# set the areas for sensitivity analysis. dim is the side lenght of the square to consider

# read deployment schemes
deployment_unit18 <- read.csv("D:/Dropbox/My Dropbox/NC-RR_environment_data/environment/HB18_units_clean.csv", h=T)
depl_scheme18 <- filter(deployment_unit18, RR.COL == "RR")
str(depl_scheme18)
deployment_unit19 <- read.csv("D:/Dropbox/My Dropbox/NC-RR_environment_data/environment/HB_UnitsScheme19.csv", h=T)
str(deployment_unit19)
store <- data.frame()

for (rec in recs) {
  RRname <- str_sub(rec, 1, 2) #get site name
  year <- str_sub(rec, 3, 4) #get year
  data <-
    raster(paste0("D:/Dropbox/My Dropbox/NC-RR_environment_data/dems/", rec))
  data <- projectRaster(data, crs = crs(TB17dem))
  shp <-
    readOGR(
      paste0(
        "D:/Dropbox/My Dropbox/NC-RR_environment_data/shps/",
        RRname,
        year,
        ".shp"
      )
    )
  shp <- spTransform(shp, crs(data))
  coords <- shp@coords
  storeb <- data.frame()
  ## start image
  # pdf(paste0("D:/Dropbox/My Dropbox/NC-RR_environment_data/outputs/", RRname, year,"_", dim, "cut.pdf"),
  #     width = 10, height = 10)
  #
  # plot(data, asp=1, col=terrain.colors(50))
  #
  ## loop for every hobo
  for (i in 1:nrow(coords)) {
    lon <- coords[i, 1]  #get longitude
    lat <- coords[i, 2]  #get latitude
    unit <- shp@data[i, 1]  #get unit name
    dpHB <- raster::extract(data,cbind(lon,lat)) #get depth of the hobo
    
    for (L in c(0.25, 0.5, 0.75)) {
      # list all the sizes you want to compute around coord pair
      scl <-
        L / c(1, 2, 4, 8, 16, 32, 64) # Scales, they cover 2 orders magnitude
      L0 <- min(scl) # Min grain
      x0 <- lon - L / 2
      y0 <- lat - L / 2
      
      temp <- crop(data, extent(x0, x0 + L, y0, y0 + L))
      temp2 <- temp
      temp2@data@values[temp2@data@values < dpHB] <- dpHB
      H <-
        max(temp2[!is.na(temp2)]) - min(temp2[!is.na(temp2)]) #heigh range
      temp3 <- as(temp2, 'SpatialGridDataFrame')
      Rplane <- surfaceArea(temp3) / (sum(temp2@data@values>=dpHB)*res(temp)[1]*res(temp)[2]) #Surface area
      filter <- temp < dpHB
      temp4 <- mask(temp, filter, maskvalue=1)
      temp5 <- as(temp4, 'SpatialGridDataFrame')
      Rclip <- surfaceArea(temp5)/(sum(temp@data@values>=dpHB)*res(temp)[1]*res(temp)[2])
      Dplane <-
        3 - log10(H / (sqrt(2) * L0 * sqrt(Rplane ^ 2 - 1))) / log10(L / L0) #Fractal dimension
      Dclip <-
        3 - log10(H / (sqrt(2) * L0 * sqrt(Rclip ^ 2 - 1))) / log10(L / L0)
      storeb <-
        rbind(storeb,
              data.frame(
                rec = RRname,
                year,
                lon,
                lat,
                Dclip,
                Rclip,
                Dplane,
                Rplane,
                H,
                L,
                unit,
                dpHB,
                new_unit = NA
              ))
    } # closes dim within unit loop
    
  } # closes coord pairs (units) within rec

  # add coherent unit ID
  
  # in 2017 we have only nembers (1 to 50) in the shps, so add H or H0 before that
  # in 2018 we have units numbers (U01 to U30). Gather hobo IDs from the deployment scheme
  # in 2019 we have light hobo IDs - nothing needed
  if (year == 17) {
    storeb$new_unit[str_length(storeb$unit) == 1] <-
      paste0("H0", storeb$unit[str_length(storeb$unit) == 1])
    storeb$new_unit[str_length(storeb$unit) == 2] <-
      paste0("H", storeb$unit[str_length(storeb$unit) == 2])
  }
  if (year == 18) {
    prov <- merge(
      storeb,
      depl_scheme18,
      by.x = c("rec", "unit"),
      by.y = c("site", "unit.ID")
    )
    storeb$new_unit <-
      as.character(prov$light[match(storeb$unit, prov$unit)])
  }
  if (year == 19) {
    storeb$new_unit <- storeb$unit
  }
  
  store <- rbind(store, storeb)
  write.csv(store, "output/geom_new.csv", row.names = FALSE)
} # closes rec within dim loop

write.csv(
  store,
  paste0(
    "D:/Dropbox/My Dropbox/NC-RR_environment_data/outputs/geometric_traits_",
    dim,
    "cut.csv"
  ),
  row.names = FALSE
)


str(store)
unique(store$new_unit)

# END

# geometric variables distribution

ggplot(store, aes(Dplane))+
  geom_histogram()+
  facet_grid(~rec)+
  theme_minimal()

ggplot(store, aes(Rplane))+
  geom_histogram()+
  facet_wrap(~rec)+
  theme_minimal()

ggplot(store, aes(Dclip))+
  geom_histogram()+
  facet_wrap(~rec)+
  theme_minimal()

ggplot(store, aes(Rclip))+
  geom_histogram()+
  facet_wrap(~rec)+
  theme_minimal()