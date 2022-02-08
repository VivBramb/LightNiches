###############################################
#####   Read out response variables       #####
###############################################

# Script to read light and prepare tables for the analysis
# It includes producing merged df with par and tide info
# 
# files <-list.files("D:/Dropbox/My Dropbox/NC-RR_environment_data/dems/data/shps", full.names = TRUE)
# 
# sapply(files,FUN=function(eachPath){
#   file.rename(from=eachPath,to=sub(pattern="CB",replacement="CB17",eachPath))
# })


#### libraries #####

library(stringr)
library(tidyverse)
library(lubridate)
library(ggplot2)
library("TideHarmonics")
library("RColorBrewer")
library(zoo)

#### functions ####
ld <- function(lux, depth) {
  lux_log <- log(lux)
  k <- 4/50
  lux_log <- lux_log - depth * k
  return(exp(lux_log))
}

## load deployment dates for each sampling event for all the variable
deployment <- read.csv("D:/Dropbox/My Dropbox/NC-RR_environment_data/environment/deployment_scheme.csv", h=T)
deployment$datetime0 <- as.POSIXct(as.character(deployment$from),format = "%d/%m/%Y %H:%M")
deployment$datetime1 <- as.POSIXct(as.character(deployment$to),format = "%d/%m/%Y %H:%M")
str(deployment)

# prepare dataframe for storing light reading form the hobos dataloggers (h)
store_h_all <- data.frame()

# get folder names
folders <- dir("D:/Dropbox/My Dropbox/NC-RR_environment_data/environment/light")

## read unit composition from 2018 (when more than 1 logger were deployed at each location)
# each unit since 2018 had a light and a temp logger and since they were flooding a lot
# the pair changed from a deployment to the other

deployment_unit18 <- read.csv("D:/Dropbox/My Dropbox/NC-RR_environment_data/environment/HB18_units_clean.csv", h=T)
str(deployment_unit18)
# filter units used at the RR scale
depl_scheme18 <- filter(deployment_unit18, RR.COL == "RR")
str(depl_scheme18) 

# read hobos looping in all the folders
for (folder in folders) {
  
  site <- str_sub(folder,1,2)
  date <- str_sub(folder,3,4)
  rep <- str_sub(folder,5,5)
  hbs <- dir(paste0("D:/Dropbox/My Dropbox/NC-RR_environment_data/environment/light/", folder,"/csv/"))
  
  if (date == 17) {
    hbs_cut <- unique(str_sub(hbs,1,4))
    for (hb in hbs_cut) {
      HB <- read.csv(paste0("D:/Dropbox/My Dropbox/NC-RR_environment_data/environment/light/", folder,"/csv/", hb,".csv"), skip = 1)
      HB<- HB[,1:4]
      HB$hoboID <- as.factor(hb)
      HB$site <- site
      HB$date <- date
      HB$rep <- rep
      HB[,2] <- as.POSIXct(as.character(HB[,2]), 
                           format = "%m/%d/%y %I:%M:%S %p")
      HB$unit <- paste0("H", str_sub(hb,3,4))
      HB$folder <- folder
      colnames(HB) <- c("x","datetime","temp","light","hoboID","site","year","rep","unit", "folder")
      HB_depl <- HB %>%
        filter(datetime >= deployment$datetime0[deployment$file_name %in% folder & deployment$variable %in% "light"], 
               datetime < deployment$datetime1[deployment$file_name %in% folder & deployment$variable %in% "light"])
      store_h_all <- rbind(store_h_all, HB_depl)
    }
  }
  if (date == 18) {
    site_depl_scheme18 <-
      filter(depl_scheme18, site == str_sub(folder, 1, 2))
    hbs_cut <- as.character(site_depl_scheme18$light)
    for (hb in hbs_cut) {
      if (any(hb == unique(str_sub(hbs, 1, 3))) == TRUE) {
        HB <- read.csv(paste0("D:/Dropbox/My Dropbox/NC-RR_environment_data/environment/light/",
                              folder, "/csv/",hb,".csv"), skip = 1)
        HB <- HB[, 1:4]
        HB$hoboID <- as.factor(hb)
        HB$site <- site
        HB$date <- date
        HB$rep <- rep
        HB[, 2] <- as.POSIXct(as.character(HB[, 2]), format = "%m/%d/%y %I:%M:%S %p")
        HB$unit <-
          site_depl_scheme18[which(site_depl_scheme18$light == hb), 6]
        HB$folder <- folder
        colnames(HB) <- c("x", "datetime", "temp", "light", "hoboID", "site",
                          "year", "rep", "unit", "folder")
        HB_depl <- HB %>%
          filter(datetime >= deployment$datetime0[deployment$file_name %in% folder &
                                                    deployment$variable %in% "light"],
                 datetime < deployment$datetime1[deployment$file_name %in% folder &
                                                   deployment$variable %in% "light"])
        store_h_all <- rbind(store_h_all, HB_depl)
      }
    }
  }
  
  if (date == 19) {
    hbs_cut <- unique(str_sub(hbs,1,3))
    for (hb in hbs_cut) {
      HB <- read.csv(paste0("D:/Dropbox/My Dropbox/NC-RR_environment_data/environment/light/", folder,"/csv/", hb,".csv"), skip = 1)
      HB<- HB[,1:4]
      HB$hoboID <- as.factor(hb)
      HB$site <- site
      HB$date <- date
      HB$rep <- rep
      HB[,2] <- as.POSIXct(as.character(HB[,2]), format = "%m/%d/%y %I:%M:%S %p")
      HB$unit <- hb
      HB$folder <- folder
      colnames(HB) <- c("x","datetime","temp","light","hoboID","site","year","rep","unit","folder")
      HB_depl <- HB %>%
        filter(datetime >= deployment$datetime0[deployment$file_name %in% folder & deployment$variable %in% "light"]) %>%
        filter(datetime < deployment$datetime1[deployment$file_name %in% folder & deployment$variable %in% "light"])
      store_h_all <- rbind(store_h_all, HB_depl)
    }
  }
}


#### plot light ####

ggplot(data = store_h_all, aes(x = datetime, y = light, col = site, alpha = .1)) +
  geom_point(shape = 20)+
  geom_smooth()+
  theme_minimal()+
  facet_wrap(.~folder, scales = "free_x")

ggplot(data = store_h_all, aes(x = datetime, y = light, col = hoboID, alpha = .1)) +
  geom_point(shape = 20)+
  geom_smooth()+
  theme_minimal()+
  theme(legend.title = element_blank(), legend.position = "none")+
  facet_wrap(.~folder, scales = "free_x")

#### check individual datasets - sanity checks ####

HS19a_subs <- store_h_all %>% filter(folder == "HS19a")
RS19a_subs <- store_h_all %>% filter(folder == "RS19a")
RS17a_subs <- store_h_all %>% filter(folder == "RS17a")
HS18b_subs <- store_h_all %>% filter(folder == "HS18b")
TB19a_subs <- store_h_all %>% filter(folder == "TB19a")
ggplot(data = TB19a_subs, aes(x = datetime, y = light)) +
  geom_point(shape = 20)+
  geom_smooth(method = "loess")+
  theme_minimal()+
  theme(legend.title = element_blank(), legend.position = "none")+
  facet_wrap(.~hoboID, scales = "free_x")
summary(TB19a_subs)
length(unique(HS18b_subs$hoboID))

ggplot(data = HS19a_subs, aes(x = datetime, y = light, col = hoboID, alpha = .1)) +
  geom_point(shape = 20)+
  geom_smooth()+
  theme_minimal()+
  theme(legend.title = element_blank(), legend.position = "none")
length(unique(HS19a_subs$hoboID))

ggplot(data = RS17a_subs[RS17a_subs$hoboID == "RS37",], aes(x = datetime, y = light, col = hoboID, alpha = .1)) +
  geom_point(shape = 20)+
  geom_smooth()+
  theme_minimal()+
  theme(legend.title = element_blank(), legend.position = "none")+
  facet_wrap(.~hoboID, scales = "free_x")
length(unique(HS19a_subs$hoboID))

## sanity checks 
# Hobos per folder

count_h <- store_h_all %>%
  group_by(folder) %>%
  summarise(n_distinct(hoboID))

# some folders miss units with respect to raw data - check

store_h_all %>%
  filter(folder == "CB18a") %>%
  distinct(hoboID)
#H32 in light folder but used for temp, so not read

store_h_all %>%
  filter(folder == "HB18b") %>%
  distinct(hoboID)
#H41 abd H51 are not in the deployment scheme file, so not read

store_h_all %>%
  filter(folder == "N319a") %>%
  distinct(hoboID)
# H67 set wrong sampling frequence - it finishes before deployment is done

store_h_all %>%
  filter(folder == "RS17a") %>%
  distinct(hoboID)
# H01 and H39 didn't start (still had 2018 data from last deplyment)
# H33 set wrong sampling frequence- it finishes before deployment is done

store_h_all %>%
  filter(folder == "L119a") %>%
  distinct(hoboID)
# J12 didn't start(still had 2018 data from last deplyment)
# H51 has wrong day - didn't start.

#### filter out hobo handling  errors ####
## I set the initial date/time of recording wrong when launching them. 
## they are all from 2017 so I wouldn't be surprised...
## take them out
store_h_all <- as.data.frame(store_h_all%>%
                               filter(hoboID != "TM25",
                                      hoboID != "SE53",
                                      hoboID != "HS41",
                                      hoboID != "RS38",
                                      hoboID != "RS40",
                                      hoboID != "TB53",
                                      hoboID != "RS27",
                                      hoboID != "RS37",
                                      hoboID != "TB29",
                                      hoboID != "SE19")) 


#### summarize readings for light  ####

light <- store_h_all
str(light)
light$mol_light <- 0.0185*light$light

# light12h <- light %>% 
#   filter(lubridate::hour(datetime) <= 18 &
#            lubridate::hour(datetime) > 6)
# light12h$time <- hms::as_hms(light12h$datetime) 
# light12h$date <- as_date(light12h$datetime) 
# 
# ggplot(data = light12h, aes (x = time, y = mol_light, col = unit)) +
#   geom_point() +
#   geom_smooth()+
#   facet_wrap(~ folder, scales = "free_x")+
#   theme(legend.position = "none")+
#   theme_minimal()

light$time <- hms::as_hms(light$datetime) 
light$date <- as_date(light$datetime) 
light$minute <- minute(light$time)
light <- as.data.frame(light %>%
                            filter(minute %in% seq(0,55,5)))

unique(light$datetime)

light_b <- as.data.frame (light %>%
  #dplyr::arrange(site,datetime) %>%  #arrange() orders the rows of a data frame by the values of selected columns
  group_by(folder) %>% #group_by() takes an existing tbl and converts it into a grouped tbl where operations are performed "by group"
  mutate(light_mol_3p = zoo::rollmean(mol_light, k = 3, fill = 0),
                light_mol_5p = zoo::rollmean(mol_light, k =5, fill = 0),
                light_3p = zoo::rollmean(light, k = 3, fill = 0),
                light_5p = zoo::rollmean(light, k =5, fill = 0)))


ggplot(data = light_b) +
  geom_point(aes(x = datetime, y = light_3p, col = hoboID, alpha = .1), shape = 20)+
  theme_minimal()+
  theme(legend.title = element_blank(), legend.position = "none")+
  facet_wrap(.~folder, scales = "free_x")

ggplot(data = light_b) +
  geom_point(aes(x = datetime, y = light_5p, col = hoboID, alpha = .1), shape = 20)+
  theme_minimal()+
  theme(legend.title = element_blank(), legend.position = "none")+
  facet_wrap(.~folder, scales = "free_x")


light <- data.frame(light_b %>%
                         group_by(folder,hoboID,unit) %>%
                         mutate(integral_mol = sum(mol_light),
                                integral_lux = sum(light),
                                integral_mol_3p = sum(light_mol_3p),
                                integral_lux_3p = sum(light_3p),
                                integral_mol_5p = sum(light_mol_5p),
                                integral_lux_5p = sum(light_5p)))

tail(sort(light$integral_lux),10L)

length(unique(light$integral_lux[light$integral_lux>100000000]))
unique(light$integral_lux[light$integral_lux>100000000])
unique(light$folder[light$integral_lux>100000000])
nrow(light)
nrow(light[light$light>40000,])


L118a_subs <- light_b %>% filter(folder == "L118a")
names(L118a_subs)
ggplot(data = L118a_subs, aes(x = datetime, y = light_mol_5p)) +
  geom_line()+
  #geom_smooth(method = "loess")+
  theme_minimal()+
  theme(legend.title = element_blank(), legend.position = "none")+
  facet_wrap(.~hoboID, scales = "free_x")

summary(lm(light$integral_lux~light$integral_lux_3p))
plot(light$integral_lux~light$integral_lux_3p)
plot(light$integral_lux~light$integral_lux_5p)
summary(lm(light$integral_lux~light$integral_lux_5p))

#### read PAR ####
# some values were missing, 
# so I assumed linear increment/decrement of the value 
# at timestamps the observations wasn't made

years <- c("2017","2018","2019")
store_PAR <- data.frame()

for(year in years) {
  par <- read.csv(paste0("D:/Dropbox/My Dropbox/NC-RR_environment_data/environment/PAR/PAR_Lizard", 
                         year,".csv"),h = T)
  
  if (year == "2019") {
    par$datetime_pos <- gsub("T", " ", as.character(par$date.time))
    par$datetime_pos <- gsub("+1000", "", as.character(par$datetime_pos))
    par$datetime_pos <- as.POSIXct(as.character(par$datetime_pos), 
                                   format = "%Y-%m-%d %H:%M:%S")
  }
  else {
    par$datetime_pos <- as.POSIXct(as.character(par$datetime), 
                                   format = "%d/%m/%Y %H:%M")
  }
  
  new_df <- data.frame(datetime_pos = seq(min(par$datetime_pos),
                                          max(par$datetime_pos), "min"))
  new_df$par <- 0
  new_df$par <- par$value[match(new_df$datetime_pos,
                                par$datetime_pos)]
  
  ## loop to fill in  missing values
  for (i in seq(from = 1,
                to = nrow(new_df) - 10,
                by = 10)) {
    if (is.na(new_df$par[i]) == FALSE) {
      for (j in 1:9) {
        new_df$par[i + j] <-
          new_df$par[i] + 0.1 * j * (new_df$par[i + 10] - new_df$par[i])
      }
    }
    
  }
  
  
  parNAind <- which(is.na(new_df$par))
  lenght_count <- rle(is.na(new_df$par))
  lenghts <- lenght_count$lengths
  values <- lenght_count$values
  lenght_start <- c(1, 1 + which(diff(is.na(new_df$par)) != 0))
  fd_par_fix <- data.frame(cbind(lenghts, values, lenght_start))
  
  ## loop to fill in  missing values when more than 9 minutes are missing
  # this is a patch, could be coded better into the previous loop
  # but does the job
  
  for (i in 1:length(lenght_start)) {
    if (fd_par_fix$values[i] == 1) {
      for (j in fd_par_fix$lenght_start[i]:(fd_par_fix$lenght_start[i] + fd_par_fix$lenghts[i] -
                                            1)) {
        delta <-
          (new_df$par[lenght_start[i + 1]] - new_df$par[lenght_start[i] - 1]) / (fd_par_fix$lenghts[i])
        new_df$par[j] <- new_df$par[j - 1] + delta
      }
    }
  }
  
  new_df$par
  store_PAR <- rbind(store_PAR, new_df)
}

store_PAR$year <- as.numeric(format(store_PAR$datetime_pos,'%Y'))

ggplot(store_PAR, aes(x = datetime_pos, y = par)) +
  geom_line()+
  theme_minimal()+
  facet_wrap(~as.factor(year), scales = "free_x")

summary(store_PAR)
View(store_PAR)

#### compute tides ####

years <- c("2017","2018","2019")
store_tides <- data.frame()

for( year in years){
  
  tides <- read.csv(paste0("D:/Dropbox/My Dropbox/NC-RR_environment_data/environment/tides/tides_Lizard",year,".csv"), h = T)
  #str(tides)
  tides$datetime_char <- as.character(tides$datetime)
  tides$datetime_pos <- as.POSIXct(tides$datetime_char, format = "%d/%m/%Y %H:%M")
  # plot(tides$datetime_pos[!is.na(tides$cm)],tides$cm[!is.na(tides$cm)],
  #      type = "l", col = tides$status[!is.na(tides$cm)])
  # 
  # ggplot(tides[!is.na(tides$cm),], aes(datetime_pos,cm, col = status))+
  #   #geom_point(size = 2)+
  #   geom_line(inherit.aes = FALSE,aes(datetime_pos,cm))+
  #   theme_minimal()
  
  tides$t <- c(1:nrow(tides))
  
  index <- c("min", "max","i", "o")
  valuesC <- c("#CA0020", "#0571B0", "#92C5DE",  "#F4A582")
  tides$col <- valuesC[match(tides$status, index)]
  
  #### commented out: other harmonics fitted to check the one suggested for Lizard was best fit
  # m1_hc7 <- ftide(tides$cm[!is.na(tides$cm)], tides$datetime_pos[!is.na(tides$cm)], hc7)
  # plot_int <- tides$datetime_pos
  # plot(m1_hc7,plot_int[1], plot_int[length(plot_int)])
  # #plot(m1_hc7,predict_int[1], predict_int[length(predict_int)], split = TRUE)
  # points(tides[!is.na(tides$cm),]$cm ~ tides[!is.na(tides$cm),]$datetime_pos,
  #        pch = 19, col = tides[!is.na(tides$cm),]$col)
  # plot(resid(m1_hc7),pch = 19, col = tides[!is.na(tides$cm),]$col)
  # sum(resid(m1_hc7))/nrow(tides[!is.na(tides$cm),])
  # 
  # 
  # m1_hc4 <- ftide(tides$cm[!is.na(tides$cm)], tides$datetime_pos[!is.na(tides$cm)], hc4)
  # plot(m1_hc4,plot_int[1], plot_int[length(plot_int)])
  # #plot(m1_hc4,predict_int[1], predict_int[length(predict_int)], split = TRUE)
  # points(tides[!is.na(tides$cm),]$cm ~ tides[!is.na(tides$cm),]$datetime_pos,
  #        pch = 19, col = tides[!is.na(tides$cm),]$col)
  # plot(resid(m1_hc4),pch = 19, col = tides[!is.na(tides$cm),]$col)
  # sum(resid(m1_hc4))/nrow(tides[!is.na(tides$cm),])
  plot_int <- tides$datetime_pos
  m1_hc37 <- ftide(tides$cm[!is.na(tides$cm)], tides$datetime_pos[!is.na(tides$cm)], hc37)
  plot(m1_hc37,plot_int[1], plot_int[length(plot_int)])
  #plot(m1_hc4,predict_int[1], predict_int[length(predict_int)], split = TRUE)
  points(tides[!is.na(tides$cm),]$cm ~ tides[!is.na(tides$cm),]$datetime_pos,
         pch = 19, col = tides[!is.na(tides$cm),]$col)
  plot(resid(m1_hc37),pch = 19, col = tides[!is.na(tides$cm),]$col)
  sum(resid(m1_hc37))/nrow(tides[!is.na(tides$cm),])
  
  # tide_minmax <- 
  #   tides %>% filter(status == "min" | status == "max")
  # 
  # str(tide_minmax)
  # m2_hc7 <- ftide(tide_minmax$cm, tide_minmax$datetime_pos, hc7)
  # plot(m2_hc7,plot_int[1], plot_int[length(plot_int)])
  # points(tides[!is.na(tides$cm),]$cm ~ tides[!is.na(tides$cm),]$datetime_pos,
  #        pch = 19, col = tide_minmax$col)
  # plot(resid(m2_hc7),pch = 19, col = tide_minmax$col)
  # sum(resid(m2_hc7))/nrow(tide_minmax)
  # 
  # 
  # m2_hc4 <- ftide(tide_minmax$cm, tide_minmax$datetime_pos, hc4)
  # plot(m2_hc4,plot_int[1], plot_int[length(plot_int)])
  # points(tides[!is.na(tides$cm),]$cm ~ tides[!is.na(tides$cm),]$datetime_pos,
  #        pch = 19, col = tide_minmax$col)
  # plot(resid(m2_hc4),pch = 19, col = tide_minmax$col)
  # sum(resid(m2_hc7))/nrow(tide_minmax)
  # 
  # m2_hc37 <- ftide(tide_minmax$cm, tide_minmax$datetime_pos, hc37)
  # plot(m2_hc37,plot_int[1], plot_int[length(plot_int)])
  # points(tides[!is.na(tides$cm),]$cm ~ tides[!is.na(tides$cm),]$datetime_pos,
  #        pch = 19, col = tide_minmax$col)
  # plot(resid(m2_hc37),pch = 19, col = tide_minmax$col)
  # sum(resid(m2_hc37))/nrow(tide_minmax)
  # }
  
  fulltide <- predict(m1_hc37,  plot_int[1], plot_int[length(plot_int)], by = 1/60)
  str(fulltide)
  plot(fulltide, type = "l")
  points(tides[!is.na(tides$cm),]$cm ~ tides[!is.na(tides$cm),]$t,
         pch = 19, col = tides[!is.na(tides$cm),]$col)
  tides$fulltide <- NA
  tides$fulltide [1:length(fulltide)] <- fulltide
  tides$year <- year
  store_tides <- rbind(store_tides,tides)
}

store_tides$fulltide[store_tides$year == "2017"] <-store_tides$fulltide[store_tides$year == "2017"]*100

ggplot(store_tides, aes(datetime_pos,fulltide))+
  geom_point()+
  facet_wrap(.~year, scale = "free_x", nrow=3)+
  theme_minimal()
summary(store_tides)

#### merge dataframes ####
## merge light, tides and PAR
## add differential light as a variable


tides <- store_tides
par <- store_PAR

str(tides)
str(par)
str(light)

colnames(light)[2] <- "datetime_pos"
env <- merge(light, tides, by = "datetime_pos")
env <- merge(env, par, by = "datetime_pos", all.x = TRUE)
str(env)
dim(env)

envclean <- env[,c(1,3:24,33,35)]
colnames(envclean) <-  c("datetime_pos", "temp", "light_lux", "hoboID", "site",
                            "year", "rep", "unit", "folder", "light_mol", 
                            "time", "date", "minute", "light_mol_3p", 
                            "light_mol_5p", "light_3p", "light_5p", 
                            "integral_mol", "integral_lux", "integral_mol_3p", 
                            "integral_lux_3p", "integral_mol_5p","integral_lux_5p", 
                            "tide_cm", "surface_par" )
summary(envclean)

write.csv(envclean, "output/env_long.csv", row.names=FALSE)

env_sum <- unique(envclean[,c("hoboID", "site",
                              "year", "rep", "unit", "folder",
                              "integral_mol", "integral_lux", "integral_mol_3p", 
                              "integral_lux_3p", "integral_mol_5p","integral_lux_5p")])

write.csv(env_sum, "output/env_int.csv", row.names=FALSE)







