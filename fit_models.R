##### complexity and attenuation curves #####


#### libraries #####

library(ggplot2)
library(tidyverse)
library(viridis)
library(mgcv)


#### functions ####
ld <- function(par, depth) {
  par_log <- log(par)
  k <- 4/50
  par_log <- par_log - depth * k
  return(exp(par_log))
  }

#### read environmental data ####

env <-  read.csv("D:/Dropbox/LightNiches/output/env_long.csv", h=T)
env$mergekey <- paste0(env$site,env$year,env$unit) # create unique ID for each HB
str(env)
dim(env)
env$datetime_pos <- as.POSIXct(env$datetime_pos)
env$time2 <- hms::as_hms(env$datetime_pos) 
env$minute <- minute(env$time2)
env <- as.data.frame(env %>%
                       filter(minute %in% seq(0,55,5)))
dim(env)
##### fit light ######

dims <- c(0.25,0.5,0.75)
# dims <- 0.5
# initialize storing lists and dfs

df_int_tot <- data.frame() 
gam_list<-list()
for (dim in dims) {
  
  geom <- read.csv(paste0("D:/Dropbox/My Dropbox/NC-RR_environment_data/outputs_medium old/geometric_traits_", dim, "cut.csv"), h = T)
  geomT <- read.csv(paste0("D:/Dropbox/My Dropbox/NC-RR_environment_data/outputs_medium old/geometric_traits_", dim, "_uncut.csv"), h = T)
  
  geom2 <- read.csv(paste0("D:/Dropbox/My Dropbox/NC-RR_environment_data/outputs_medium old/geometric_traits_", dim, "_TBSE18_cut.csv"), h = T)
  geom2T <- read.csv(paste0("D:/Dropbox/My Dropbox/NC-RR_environment_data/outputs_medium old/geometric_traits_", dim, "_TBSE18uncut.csv"), h = T)
  
  geom <- rbind(geom,geom2)
  geomT <- rbind(geomT,geom2T)
  geom$mergekey <- paste0(geom$rec, geom$year, geom$new_unit) #unique ID for hobo deployment
  geom$Du <- geomT$D
  geom$Ru <- geomT$R
  geom$Hu <- geomT$hr
  geom$Cu <- geomT$C
  geom$au <- geomT$a
  geom$bu <- geomT$b
  ## merge datasets
  df <- merge(env,geom, by = "mergekey")
  
  # str(df)
  summary(df)
  
  ##some surface par are <0.
  nrow(df[df$surface_par<0,])
  
  df$surface_par[df$surface_par<0] <- 0
  
  # correct light values as if there was always maximum registered par at the surface
  df$water_col <- df$tide_cm/100 - df$dp
  df$pard[df$surface_par !=0] <- ld(df$surface_par[df$surface_par !=0], df$water_col[df$surface_par !=0])
  df$pard [df$surface_par ==0] <- 0
  #summary(df$pard)
  df$light_pard <- 0
  df$light_pard[df$pard !=0] <- df$light_mol_5p[df$pard !=0]/df$pard[df$pard !=0]
  summary(df$light_pard)

  df_long <- as.data.frame(df %>%
                        group_by(mergekey,lat,lon) %>%
                          filter(rep == "a") %>%
                        mutate(integral_ld = sum(light_pard),
                               integral_pard = sum(pard),
                               par_integral = sum(surface_par)))
  
# get n
  df_int <- unique(df_long[,c("integral_ld","integral_pard", "hr", "site","year.x",
                         "integral_mol","integral_lux","integral_mol_3p","integral_lux_3p", 
                         "integral_mol_5p", "integral_lux_5p", "par_integral","folder", 
                         "D",  "Du", "R","dp","Ru", "Hu", "C","lon","lat",
                         "hoboID", "a","b","mergekey", "Cu","au","bu")])
  length(unique(df$mergekey)) # 903//936!
  length(unique(df_int$mergekey)) # 903
  dim(df_int) # 905  25
  sum(duplicated(df_int$mergekey)) # there are 2 duplicates: two same IDs in the same record
  df_int[duplicated(df_int$mergekey),]
  df_int[524:527,] #RS17H44
  df_int[409:412,] # N317H34
  
  #leave one out for now and then go to the shp and the paper location annotation
  df_int <- df_int[!duplicated(df_int$mergekey),]
  
  dim(df_int) # 903 26
  
  # scale values and set factors
  # tide is in cm, so divide it to get m
  
  df_int$year <- as.factor(df_int$year)
  df_int$Rl <- as.numeric(log10(df_int$R))
  df_int$Hl <- as.numeric(log10(df_int$hr))
  df_int$aa <- abs(df_int$a)
  df_int$Rul <- log(df_int$Ru)
  df_int$Hul <- log(df_int$Hu)
  df_int$dim <- dim
  
  ## get D_ theory ##
  df_int$Dt <- 3 - log10(df_int$Hu / (sqrt(2) * (dim/32) * sqrt((df_int$Ru)^2 - 1))) / log10(dim /(dim/32))
  
  ## store_df ##
  df_int$rec<- as.factor(df_int$site)
  
  slc1 <- gam(integral_ld~s(Dt,Rl)+ s(lon, lat, bs = "gp", m=2),
              data = df_int, method="REML")
  slc2 <- gam(integral_ld~Dt*Rl+ s(lon, lat, bs = "gp", m=2),
              data = df_int, method="REML")
  
  lab <-paste0("slc1",dim)
  gam_list[[lab]] <- slc1
  #summary(slc1)
  #png(paste0("D:/Dropbox/LightNiches/output/gamViz ",lab, ".png"))
  plot(slc1, residuals = TRUE, pages = 1, 
       pch = 19, cex = .5, shade = TRUE, seWithMean = TRUE)
  #dev.off()
  lab <-paste0("slc2",dim)
  gam_list[[lab]] <- slc2
  summary(slc2)
  
  df_int_tot<- df_long
  vis.gam(slc1, view=c("Dt", "Rl"), type="response", plot.type="contour", main = "GAM / s(D,R)",
          color = "gray", ylab = expression(paste("log(R "[dim], ")")), xlab = expression(paste("D"[50])))
  points(df_int$Rl~df_int$Dt, pch = 19, cex = .6, 
         col = rgb(red=0, green=0, blue=0, alpha=0.3))
  
  summary(slc1)
  vis.gam(slc2, view=c("Dt", "Rl"), type="response", plot.type="contour", main = "LM / D*R",
          color = "gray", ylab = expression(paste("log(R "[dim], ")")), xlab = expression(paste("D"[50])))
  points(df_int$Rl~df_int$Dt, pch = 19, cex = .6, 
         col = rgb(red=0, green=0, blue=0, alpha=0.3))
  summary(slc2)
}


###########
source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")

ylab.light = expression(paste("light integral, μMol/mm"^"    2"))
df_intALL <- df_int
df_intALL$folder = "zall"
df_intALL$rec = "zall"
df_int2 <- rbind(cbind(df_intALL, rec2 = df_int$rec),cbind(df_int, rec2 = df_int$rec))

gl <- ggplot(data = df_int2[df_int2$integral_ld <300,], 
             aes(x = folder, y = integral_ld, fill = rec)) +
  geom_point(aes(y = integral_ld, color = rec2),
             position = position_jitter(width = .15), alpha = 0.6, size = 1) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0), scale = "width") +
  expand_limits(x = 1) +
  #coord_flip() + # flip or not?
  scale_color_viridis("site",labels=c("zall" = "all",
                                      "CB17a" = "CB","CB18a" = "CB","CB19a" = "CB",
                                      "GT17a" = "GT","GT19a" = "GT",
                                      "HS17a" = "HS","HS19a" = "HS",
                                      "L117a" = "L1","L119a" = "L1",
                                      "N317a" = "N3","N318a" = "N3","N319a" = "N3",
                                      "RS17a" = "RS","RS18a" = "RS","RS19a" = "RS",
                                      "TB17a" = "TB","TB18a" = "TB","TB19a" = "TB",
                                      "TM17a" = "TM","TM18a" = "TM","TM19a" = "TM",
                                      "SE17a" = "SE","SE18a" = "SE","SE19a" = "SE"), 
                      discrete = TRUE,guide = FALSE) +
  scale_x_discrete(labels=c("zall" = "all",
                            "CB17a" = "CB","CB18a" = "CB","CB19a" = "CB",
                            "GT17a" = "GT","GT19a" = "GT",
                            "HS17a" = "HS","HS19a" = "HS",
                            "L117a" = "L1","L119a" = "L1",
                            "N317a" = "N3","N318a" = "N3","N319a" = "N3",
                            "RS17a" = "RS","RS18a" = "RS","RS19a" = "RS",
                            "TB17a" = "TB","TB18a" = "TB","TB19a" = "TB",
                            "TM17a" = "TM","TM18a" = "TM","TM19a" = "TM",
                            "SE17a" = "SE","SE18a" = "SE","SE19a" = "SE")) +
  scale_fill_viridis("site",labels=c("zall" = "all",
                                     "CB17a" = "CB","CB18a" = "CB","CB19a" = "CB",
                                     "GT17a" = "GT","GT19a" = "GT",
                                     "HS17a" = "HS","HS19a" = "HS",
                                     "L117a" = "L1","L119a" = "L1",
                                     "N317a" = "N3","N318a" = "N3","N319a" = "N3",
                                     "RS17a" = "RS","RS18a" = "RS","RS19a" = "RS",
                                     "TB17a" = "TB","TB18a" = "TB","TB19a" = "TB",
                                     "TM17a" = "TM","TM18a" = "TM","TM19a" = "TM",
                                     "SE17a" = "SE","SE18a" = "SE","SE19a" = "SE"), 
                     discrete = TRUE) +
  theme_bw() +
  labs(x ="\n sampling events", y = ylab.light)+
  guides(fill = FALSE,
         color = FALSE, 
         shape = FALSE,
         size = FALSE,
         alpha = FALSE) +
  theme(panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text =  element_text(size = 12),
        text = element_text(size = 12),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16))

gl

slc1 <- gam(int_scaled~s(Dt,Rl)+ b + aa +s(lon, lat, bs = "gp", m=2),
            data = df_int, method="REML")
lab <-paste0("slc1",dim)
gam_list[[lab]] <- slc1
png(paste0("D:/Dropbox/LightNiches/output/gamViz ",lab, ".png"))
plot(slc1, residuals = TRUE, pages = 1, pch = 19, cex = .5, shade = TRUE, seWithMean = TRUE, shift = coef(sr)[1])
dev.off()

pred <- predict(slc1, df_int)
preddf <-cbind(pred = pred, obs = df_int$int_scaled, mod = lab)
write.csv(preddf, 
          paste0("D:/Dropbox/LightNiches/output/PrevObs",lab,".csv"), 
          row.names=FALSE)
sum <- summary(slc1)
s.table <- cbind(par = row.names(sum$s.table), sum$s.table)

write.csv(s.table, 
          paste0("D:/Dropbox/LightNiches/output/model_summary_s ",lab,"light.csv"),
          row.names=FALSE)
write.csv(cbind(sum$r.sq,sum$dev.expl,sum$n), 
          paste0("D:/Dropbox/LightNiches/output/model_summary_p ",lab,"light.csv"),
          row.names=FALSE)

vis.gam(slc1, view=c("Dt", "Rl"), type="response", plot.type="contour", main = "50cm-side\n",
        color = "gray", ylab = expression(paste("log(R "[50], ")")), xlab = expression(paste("D"[50])))
points(df_int_tot[df_int_tot$dim == 0.5,]$Rl~df_int_tot[df_int_tot$dim == 0.5,]$Dt, pch = 19, cex = .6, 
       col = rgb(red=0, green=0, blue=0, alpha=0.3))


plot(int_scaled ~ Rl, df_int[df_int$int_scaled < 10000,])
plot(int_scaled ~ log(Ru), df_int[df_int$int_scaled < 10000,])
