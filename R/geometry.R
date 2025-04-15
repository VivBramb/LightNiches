# script to get all the explanatory variables and prepare tables for the analysis

#### libraries #####

library(rgl)
library(Rvcg)
library(raster)
library(rgdal)
library(stringr)
library(sp)
library(ggplot2)
library(tidyverse)

#### loop for light  ####

# TB17 is the crs reference
TB17dem <- raster("~/Dropbox/My Dropbox/NC-RR_environment_data/dems/TB17.tif")

# read dems names
# rec are "RRYY.tif" format: RR=record site code, YY = year
recs <- dir("~/Dropbox/My Dropbox/NC-RR_environment_data/dems")
recs <- recs[recs != "L119_old.tif"]
# set the areas for sensitivity analysis. dim is the side lenght of the square to consider

# read deployment schemes
deployment_unit18 <- read.csv("~/Dropbox/My Dropbox/NC-RR_environment_data/environment/HB18_units_clean.csv", h=T)
depl_scheme18 <- filter(deployment_unit18, RR.COL == "RR")
str(depl_scheme18)
deployment_unit19 <- read.csv("~/Dropbox/My Dropbox/NC-RR_environment_data/environment/HB_UnitsScheme19.csv", h=T)
str(deployment_unit19)
store <- data.frame()

for (rec in recs) {
  RRname <- str_sub(rec, 1, 2) #get site name
  year <- str_sub(rec, 3, 4) #get year
  data <-
    raster(paste0("~/Dropbox/My Dropbox/NC-RR_environment_data/dems/", rec))
  data <- projectRaster(data, crs = crs(TB17dem))
  shp <-
    readOGR(
      paste0(
        "~/Dropbox/My Dropbox/NC-RR_environment_data/shps/",
        RRname,
        year,
        ".shp"
      )
    )
  shp <- spTransform(shp, crs(data))
  coords <- shp@coords
  storeb <- data.frame()
  ## start image
  # pdf(paste0("~/Dropbox/My Dropbox/NC-RR_environment_data/outputs/", RRname, year,"_", dim, "cut.pdf"),
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
    
    for (L in c(0.25, 0.5, 0.75, 1)) {
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
      R <- surfaceArea(temp3) / (L*L)
      Rplane <- surfaceArea(temp3) / (sum(temp2@data@values>=dpHB)*res(temp)[1]*res(temp)[2]) #Surface area
      filter <- temp < dpHB
      temp4 <- mask(temp, filter, maskvalue=1)
      temp5 <- as(temp4, 'SpatialGridDataFrame')
      Rclip <- surfaceArea(temp5)/(sum(temp@data@values>=dpHB)*res(temp)[1]*res(temp)[2])
      Dplane <-
        3 - log10(H / (sqrt(2) * L0 * sqrt(Rplane ^ 2 - 1))) / log10(L / L0) #Fractal dimension
      Dclip <-
        3 - log10(H / (sqrt(2) * L0 * sqrt(Rclip ^ 2 - 1))) / log10(L / L0)
      D <-
        3 - log10(H / (sqrt(2) * L0 * sqrt(R ^ 2 - 1))) / log10(L / L0)
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
                R,
                D,
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
  write.csv(store, "output/store_2024.csv", row.names = FALSE)
} # closes rec within dim loop

write.csv(
  store,
  paste0(
    "~/outputs/geometric_traits_",
    dim,
    "cutanduncut.csv"
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

##### full site metric #####


store05 <- data.frame()
store1 <- data.frame()
store2 <- data.frame()
store4 <- data.frame()
store8 <- data.frame()


coords <- read.table("data/RR_Coords.txt", head = TRUE)
Ccoord <- coords[,c("x", "y")]
c2 <- SpatialPoints(Ccoord)

#View(files_scheme)
# read coords of the centers and transpose to right projection
c2@proj4string <- CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
c <- spTransform(c2, CRS("+proj=tmerc +lat_0=-14.6989736557007 +lon_0=145.448257446289 
                         +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"))
coords_tr<- as.data.frame(c)
Ccoords_tr<- as.data.frame(cbind(site = coords$Site,
                                 x= coords_tr$x, 
                                 y= coords_tr$y))
Ccoords_tr$x <- as.numeric(Ccoords_tr$x)
Ccoords_tr$y <- as.numeric(Ccoords_tr$y)
store <- data.frame()

L <- 2
scl <- L / c(1, 2, 4, 8, 16, 32, 64, 128) # Scales, they cover 2 orders magnitude
L0 <- min(scl) # Min grain


recs <- recs[c(1:6, 8:24)]


# get site variables from 2018
recs_2019 <- recs[str_detect(recs, "19")]

for (rec in recs_2019) {
  
  # Load geotif for reef record
  RRname <- str_sub(rec, 1, 2) #get site name
  year <- str_sub(rec, 3, 4) #get year
  data <-
    raster(paste0("~/Dropbox/My Dropbox/NC-RR_environment_data/dems/", rec))
  if(crs(data)@projargs != crs(TB17dem)@projargs){
    data <- projectRaster(data, crs = crs(TB17dem))
  }
  
  # Get lower corner of 8x8m bounding box;
  xb <- Ccoords_tr$x[which(Ccoords_tr$site == RRname)] - 4
  yb <- Ccoords_tr$y[which(Ccoords_tr$site == RRname)] - 4
  
  # Iterate through  1x1 m quadrats (reps = 64) in reef record
  # Scope (extent), scales of variation, and resolution (grain)
  # L <- 1 # Scope, 1 by 1 m reef patches
  # scl <- L / c(1, 2, 4, 8, 16, 32, 64) # Scales, aim for 2 orders magnitude
  # L0 <- min(scl) # Grain, resolution of processing ~ 16 mm
  # 
  # rep <- 1
  # grid1 <- seq(0,7, by = 1)
  # for (i in grid1) {
  #   for (j in grid1) {
  #     x0 <- xb + i
  #     y0 <- yb + j
  #     temp <- crop(data, extent(x0, x0 + L, y0, y0 + L))
  #     H <- max(temp[!is.na(temp)]) - min(temp[!is.na(temp)]) #heigh range
  #     temp2 <- as(temp, 'SpatialGridDataFrame')
  #     R <- surfaceArea(temp2) / (sum(!is.na(temp@data@values))*res(temp)[1]*res(temp)[2])
  #     D <- 3 - log10(H / (sqrt(2) * L0 * sqrt(R^2 - 1))) / log10(L / L0)
  #     store1 <- rbind(store1, data.frame(rec=rec, rep=rep, x0 = x0, y0 = y0, H=H, R=R, D=D))
  #     write.csv(store1, "output/grid1_2024.csv", row.names=FALSE)
  #     rep <- rep + 1
  #   }
  # }
  # 
  # 
  # # Iterate through  2x2 m quadrats (reps = 16) in reef record
  # # Scope (extent), scales of variation, and resolution (grain)
  # L <- 2 # Scope, 2 by 2 m reef patches
  # scl <- L / c(1, 2, 4, 8, 16, 32, 64) # Scales, aim for 2 orders magnitude
  # L0 <- min(scl) # Grain, resolution of processing

  rep <- 1
  grid2 <- seq(0,6, by = 2)
  for (i in grid2) {
    for (j in grid2) {
      x0 <- xb + i
      y0 <- yb + j
      temp <- crop(data, extent(x0, x0 + L, y0, y0 + L))
      H <- max(temp[!is.na(temp)]) - min(temp[!is.na(temp)]) #heigh range
      temp2 <- as(temp, 'SpatialGridDataFrame')
      R <- surfaceArea(temp2) / (sum(!is.na(temp@data@values))*res(temp)[1]*res(temp)[2])
      D <- 3 - log10(H / (sqrt(2) * L0 * sqrt(R^2 - 1))) / log10(L / L0)
      store2 <- rbind(store2, data.frame(rec=rec, rep=rep, x0 = x0, y0 = y0, H=H, R=R, D=D))
      write.csv(store2, "output/grid2_2024.csv", row.names=FALSE)
      rep <- rep + 1
    }
  }


    # Iterate through  4x4 m quadrats (reps = 4) in reef record
    # # Scope (extent), scales of variation, and resolution (grain)
    # L <- 4 # Scope, 4 by 4 m reef patches
    # scl <- L / c(1, 2, 4, 8, 16, 32, 64) # Scales, aim for 2 orders magnitude
    # L0 <- min(scl) # Grain, resolution of processing
    # 
    # rep <- 1
    # for (i in c(0,4)) {
    #   for (j in c(0,4)) {
    #     x0 <- xb + i
    #     y0 <- yb + j
    #     temp <- crop(data, extent(x0, x0 + L, y0, y0 + L))
    #     H <- max(temp[!is.na(temp)]) - min(temp[!is.na(temp)]) #heigh range
    #     temp2 <- as(temp, 'SpatialGridDataFrame')
    #     R <- surfaceArea(temp2) / (sum(!is.na(temp@data@values))*res(temp)[1]*res(temp)[2])
    #     D <- 3 - log10(H / (sqrt(2) * L0 * sqrt(R^2 - 1))) / log10(L / L0)
    #     store4 <- rbind(store4, data.frame(rec=rec, rep=rep, x0 = x0, y0 = y0, H=H, R=R, D=D))
    #     write.csv(store4, "output/grid4_2024.csv", row.names=FALSE)
    #     rep <- rep + 1
    #   }
    # }

    # Iterate through  8x8m quadrats (reps = 1) in reef record
    # Scope (extent), scales of variation, and resolution (grain)

  
  # # get variables
  # temp <- crop(data, extent(xb, xb + 8, yb, yb + 8))
  # H <-    max(temp[!is.na(temp)]) - min(temp[!is.na(temp)]) #height range
  # temp2 <- as(temp, 'SpatialGridDataFrame')
  # Asurf<- surfaceArea(temp2)
  # R <- Asurf / (sum(!is.na(temp@data@values))*res(temp)[1]*res(temp)[2]) #Surface area
  # D <- 3 - log10(H / (sqrt(2) * L0 * sqrt(R ^ 2 - 1))) / log10(L / L0) #Fractal dimension
  # R64 <- Asurf / 64 #Surface area
  # D64 <- 3 - log10(H / (sqrt(2) * L0 * sqrt(R64 ^ 2 - 1))) / log10(L / L0) #Fractal dimension
  # 
  # #store and save file
  # store <- rbind(store, data.frame(rec = RRname, year, x0 = xb+4, y0 = yb+4, D, R, H, R64, D64, Asurf))
  # write.csv(store, "output/geom_sites_all.csv", row.names = FALSE)
}
