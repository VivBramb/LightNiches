# script to get all the explanatory variables and prepare tables for the analysis

#### libraries #####

library(rgl)
library(parallel)
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
      Rplane <- surfaceArea(temp3) / (sum(temp2@data@values>=dpHB)*res(temp)[1]*res(temp)[2]) #Surface area
      filter <- temp < dpHB
      temp4 <- mask(temp, filter, maskvalue=1)
      plot(temp5)
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

} # closes dim loop

str(store)
unique(store$new_unit)

# END

# geometric variables distribution

ggplot(store, aes(D))+
  geom_histogram()+
  facet_grid(~rec)+
  theme_minimal()

ggplot(store, aes(R))+
  geom_histogram()+
  facet_wrap(~rec)+
  theme_minimal()

ggplot(store, aes(C))+
  geom_histogram()+
  facet_wrap(~rec)+
  theme_minimal()

ggplot(store, aes(dp))+
  geom_histogram()+
  facet_wrap(~rec)+
  theme_minimal()

ggplot(store, aes(a))+
  geom_histogram()+
  facet_wrap(~rec)+
  theme_minimal()

ggplot(store, aes(b))+
  geom_histogram()+
  facet_wrap(~rec)+
  theme_minimal()


#### loop for temperature ####

# TB17 is the crs reference
TB17dem <- raster("D:/Dropbox/My Dropbox/NC-RR_environment_data/dems/TB17.tif")

# read dems names
# rec are "RRYY.tif" format: RR=record site code, YY = year
recs <- dir("D:/Dropbox/My Dropbox/NC-RR_environment_data/dems")

# set the areas for sensitivity analysis. dim is the side lenght of the square to consider
dims <- c(0.25,0.5,0.75,1)

# read deployment schemes
deployment_unit18 <- read.csv("D:/Dropbox/My Dropbox/NC-RR_environment_data/environment/HB18_units_clean.csv", h=T)
depl_scheme18 <- filter(deployment_unit18, RR.COL == "RR")
str(depl_scheme18)
depl_scheme19 <- read.csv("D:/Dropbox/My Dropbox/NC-RR_environment_data/environment/HB_UnitsScheme19.csv", h=T)
str(depl_scheme19)

for (dim in dims) {
  
  # get windows for D
  wins <- c(dim, dim/2, dim/4, dim/8, dim/16, dim/32)
  
  store <- data.frame() # store the concatenation of all the RR
  storeb <- data.frame() # store a RR at a time
  
  for (rec in recs) {
    RRname <- str_sub(rec,1,2) #get site name
    year <- str_sub(rec,3,4) #get year
    data <- raster(paste0("D:/Dropbox/My Dropbox/NC-RR_environment_data/dems/", rec))
    data <- projectRaster(data, crs = crs(TB17dem))
    shp <- readOGR(paste0("D:/Dropbox/My Dropbox/NC-RR_environment_data/shps/", RRname, year, ".shp"))
    shp <- spTransform(shp, crs(data))
    
    ### leaving this here for when cover is needed
    # cover <- readOGR(paste0("F:/environment_paper/data/shps/", RRname, "_cover.shp"))
    # cover <- spTransform(cover,crs(data))
    # cov_sf <- st_as_sf(cover, coords = c("long", "lat"))
    # cov_pol <- as(cov_sf,"Spatial")
    # cv <- gIntersection(soft_cover ,as(extent(data), "SpatialPolygons"))
    # hard_cover <- cov_pol[(cov_pol$type == "coral" ||
    #                       cov_pol$type == "hard" ||
    #                         cov_pol$type == "hard?" ||
    #                       cov_pol$type == "NA"),]
    # soft_cover <- cov_pol[(cov_pol$type == "soft"||
    #                       cov_pol$type == "soft?"),]
    
    ## read hobos coordinates
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
      
      lon <- coords[i,1]  #get longitude
      lat <- coords[i,2]  #get latitude
      unit <- shp@data[i,1]  #get unit name
      
      x0 <- lon - dim/2 #get area bottom right longitude
      y0 <- lat - dim/2 #get area bottom right latitude
      dpHB <- raster::extract(data,cbind(lon,lat)) #get depth of the hobo
      
      # enter the loop only if the square around the hobo is within the map and
      # not too close to the margins
      if((x0 > extent(data)[1] & (x0+dim) < extent(data)[2]) &
         (y0 > extent(data)[3] & (y0+dim) < extent(data)[4])) {
        temp <- data.frame()
        
        # get D with the variation method
        for (i in wins) {
          ss <- seq(i, dim, i) - i
          x <- rep(ss, dim/i)
          y <- rep(ss, each=dim/i)
          # If you can't get multicore working, change "mcmapply" to "mapply" below:
          temp <- rbind(temp, data.frame(i=i, x=x, y=y, e=mcmapply(r_func, x, y)))
        }
        
        temp2 <- aggregate(e ~ i, temp, mean)
        D <- 3 - coef(lm(log10(e) ~ log10(i), temp2))[2]
        
        # get curvature
        temp3 <- crop(data, extent(x0, x0 + dim, y0, y0 + dim))
        curv <- curvature(temp3, type = "mcnab")
        C <- sum(curv@data@values[!is.na(curv@data@values)])
        
        # get easting and northing
        t3 <- rasterToPoints(temp3)
        colnames(t3) <- c("x","y","z")
        t3 <- as.data.frame(t3)
        bfplane <- lm(data =t3, z~x+y)
        a <- as.numeric(coef(bfplane)[2]) # easting
        b <- as.numeric(coef(bfplane)[3]) # northing
        dp <- dpHB # depth
        hr <- max(temp3[!is.na(temp3)]) - min(temp3[!is.na(temp3)]) #heigh range
        
        # cover for when will be needed
        #PCs <- PC(temp3, soft_cover)
        # SCs <- SC(temp3, soft_cover)
        #PCh <- PC(temp3, hard_cover)
        # SCh <- SC(temp3, hard_cover)
        
        # get R
        temp3 <- as(temp3, 'SpatialGridDataFrame')
        R <- surfaceArea(temp3) / dim^2
        
        #plot rectangle and values
        # rect(x0, y0, x0 + dim, y0 + dim, border="black", lwd=2)
        # text(mean(c(x0, x0 + dim)), mean(c(y0, y0 + dim)), 
        #      paste(unit, "\n", "D=", round(D, 2), "\nR=", round(R, 2), "\nhr=", round(hr,2)),
        #      col="black", cex=0.5)
        storeb <- rbind(storeb,data.frame(rec=RRname, year, lon, lat, D, R, dp, C,a,b, hr, dim,
                                          unit, new_unit = NA))
        
      } # closes if extent overlap loop
    } # closes hohos within rec within dim loop
    
    # add coherent unit ID 
    
    # in 2017 we have only nembers (1 to 50) in the shps, so add H or H0 before that
    # in 2018 we have units numbers (U01 to U30). Gather hobo IDs from the deployment scheme
    # in 2019 we have light hobo IDs - nothing needed
    if (year == 17){
      storeb$new_unit[str_length(storeb$unit) == 1] <- paste0("H0", storeb$unit[str_length(storeb$unit) == 1])
      storeb$new_unit[str_length(storeb$unit) == 2] <- paste0("H", storeb$unit[str_length(storeb$unit) == 2])
    }
    if (year ==18){
      prov <- merge(storeb,depl_scheme18, by.x = c("rec","unit"),
                    by.y = c("site","unit.ID"))
      storeb$new_unit <- as.character(prov$temp[match(storeb$unit, prov$unit)])
    }
    if (year ==19){
      storeb$new_unit <- as.character(depl_scheme19$temp[match(storeb$unit, depl_scheme19$light)])
    }
    
    # dev.off()
    store <- rbind(store, storeb)
    
  } # closes rec within dim loop
  
  write.csv(store, paste0("D:/Dropbox/My Dropbox/NC-RR_environment_data/outputs/geometric_traits_",dim,"_uncut.csv"), row.names=FALSE)
  
} # closes dim loop

str(store)
unique(store$new_unit)

# END

# geometric variables distribution

ggplot(store, aes(D))+
  geom_histogram()+
  facet_grid(~rec)+
  theme_minimal()

ggplot(store, aes(R))+
  geom_histogram()+
  facet_wrap(~rec)+
  theme_minimal()

ggplot(store, aes(C))+
  geom_histogram()+
  facet_wrap(~rec)+
  theme_minimal()

ggplot(store, aes(dp))+
  geom_histogram()+
  facet_wrap(~rec)+
  theme_minimal()

ggplot(store, aes(a))+
  geom_histogram()+
  facet_wrap(~rec)+
  theme_minimal()

ggplot(store, aes(b))+
  geom_histogram()+
  facet_wrap(~rec)+
  theme_minimal()


© 2022 GitHub, Inc.
Terms
Privacy
Security
Status
Docs
Contact GitHub
Pricing
API
Training
Blog
About
