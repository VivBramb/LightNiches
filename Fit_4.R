### Interaction



#### libraries #####
library(lme4)
library(mgcv)
library(ggplot2)
library(mgcViz)
library(tidyverse)
library(plot3D)
library(viridis)
library(sjPlot)
library(stargazer)
library(GGally)
library(ggpubr)
library(itsadug)
##
library(MuMIn)


#### functions ####
ld <- function(lux, depth) {
  lux_log <- log(lux)
  k <- 4/50
  lux_log <- lux_log - depth * k
  return(exp(lux_log))
}

#### read environmental data ####

env <-  read.csv("D:/Dropbox/My Dropbox/NC-RR_environment_data/outputs_medium old/env12hclean2.csv", h=T)
env$mergekey <- paste0(env$site,env$yr,env$unit) # create unique ID for each HB
#str(env)
env$datetime_pos <- as.POSIXct(env$datetime_pos)
env <- as.data.frame(env%>%
                       filter(rep == "a")) #only get 1 rep per site
env$time2 <- hms::as_hms(env$datetime_pos) 
env$minute <- minute(env$time2)
env <- as.data.frame(env %>%
                       filter(minute %in% seq(0,55,5)))

##### fit light ######

dims <- c(0.25,0.5,0.75)
# dims <- 0.5
# initialize storing lists and dfs
gam_list<- list()
lmer_list <- list()
df_int_tot <- data.frame() 
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
  df <-merge(env,geom, by = "mergekey")
  
  # str(df)
  summary(df)
  
  ##some surface par are <0.
  nrow(df[df$par<0,])
  
  ##coerce to 0
  df$par[df$par<0] <- 0
  
  ## get maximum surface par at each time of measurement
  df <- data.frame(df%>%
                     group_by(time) %>%
                     mutate(par_max_time = max(par)))
  
  # correct light values as if there was always maximum registered par at the surface
  df$water_col <- df$tide/100 - df$dp
  df$pard <- ld(df$par, df$water_col)
  df$light_pard <- df$light/df$pard
  ###this needs fixing!!
  
  sum(df$pard == 0)

  df$light_pard[is.infinite(df$light_pard)] <- 0
  df$light_pard[is.na(df$light_pard)] <- 0
  
  ## set the initial date/time of recording wrong when launching them. 
  ## they are all from 2017 so I wouldn't be surprise.
  # take them out
  
  df <- as.data.frame(df%>%
                        filter(mergekey != "TM17H25",
                               mergekey != "SE17H53",
                               mergekey != "HS17H41",
                               mergekey != "RS17H38",
                               mergekey != "RS17H40",
                               mergekey != "TB17H53",
                               mergekey != "RS17H27",
                               mergekey != "RS17H37",
                               mergekey != "TB18H29",
                               mergekey != "SE18H19")) 
  
  # get integrals of par at the surface, 
  # par at the bottom,  mean tide height 
  # and water column at 12:00
  df <- data.frame(df%>%
                     group_by(mergekey,lat,lon) %>%
                     mutate(par_integral = sum(par),
                            mean_tide = mean(tide),
                            int_scaled = sum(light_pard ),
                            light_noscale = sum(mol_light),
                            watercol12 = tide[time == "12:00"]/100 - dp))
  df <- data.frame(df%>%
                     group_by(folder) %>%
                     mutate(int_scale2 = as.numeric(scale(light_noscale))))
  df$dp_tide <- df$dp + df$tide/100
  # get n
  df_int <- unique(df[,c("light_noscale","int_scaled","int_scale2", "hr", "mean_tide", "rec","year",
                         "par_integral","folder", "D",  "Du", "R","dp","Ru", "Hu", "C","lon","lat","watercol12",
                         "hoboID", "a","b","mergekey", "Cu","au","bu")])
  length(unique(df$mergekey)) # 903
  length(unique(df_int$mergekey)) # 903
  dim(df_int) # 905  25
  sum(duplicated(df_int$mergekey)) # there are 2 duplicates: two same IDs in the same record
  
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
  df_int$tdp <- (df_int$mean_tide/100) - df_int$dp
  df_int$aa <- abs(df_int$a)
  df_int$Rul <- log(df_int$Ru)
  df_int$Hul <- log(df_int$Hu)
  df_int$dim <- dim
  
  ## get D_ theory ##
  df_int$Dt <- 3 - log10(df_int$Hu / (sqrt(2) * (dim/32) * sqrt((df_int$Ru)^2 - 1))) / log10(dim /(dim/32))
  
  ## store_df ##
  df_int$rec<- as.factor(df_int$rec)
  
  
  ### D*R ####
  sc1 <- gam(light_noscale~s(Dt,Rl)+ s(watercol12)+
              s(par_integral)+
              s(b) + s(aa) +s(lon, lat, bs = "gp", m=2),
             data = df_int, method="REML")
  lab <-paste0("sc1",dim)
  gam_list[[lab]] <- sc1
  png(paste0("D:/Dropbox/LightNiches/output/gamViz",lab, ".png"))
  plot(sc1, residuals = TRUE, pages = 1, pch = 19,
       cex = .5, shade = TRUE, seWithMean = TRUE, shift = coef(sc1)[1])
  dev.off()
  
  pred <- predict(sc1, df_int)
  preddf <-cbind(pred = pred, obs = df_int$light_noscale, mod = lab)
  write.csv(preddf, 
            paste0("D:/Dropbox/LightNiches/output/PrevObs",lab,".csv"), 
            row.names=FALSE)
  sum <- summary(sc1)
  s.table <- cbind(par = row.names(sum$s.table), sum$s.table)
  
  write.csv(s.table, 
            paste0("D:/Dropbox/LightNiches/output/model_summary_s ",lab,"light.csv"),
            row.names=FALSE)
  write.csv(cbind(sum$r.sq,sum$dev.expl,sum$n), 
            paste0("D:/Dropbox/LightNiches/output/model_summary_p ",lab,"light.csv"),
            row.names=FALSE)
  
  slc1 <- gam(light_noscale~s(Dt,Rl)+ watercol12+
                par_integral+
                b + aa +s(lon, lat, bs = "gp", m=2),
              data = df_int, method="REML")
  lab <-paste0("slc1",dim)
  gam_list[[lab]] <- slc1
  png(paste0("D:/Dropbox/LightNiches/output/gamViz ",lab, ".png"))
  plot(slc1, residuals = TRUE, pages = 1, pch = 19, cex = .5, shade = TRUE, seWithMean = TRUE, shift = coef(sr)[1])
  dev.off()
  
  pred <- predict(slc1, df_int)
  preddf <-cbind(pred = pred, obs = df_int$light_noscale, mod = lab)
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
  
  sr1 <- gam(light_noscale~s(Dt,Rl)+ s(watercol12)+s(par_integral)+
              s(b) + s(aa) + s(rec, bs = "re"), data = df_int, method="REML")
  lab <-paste0("sr1",dim)
  gam_list[[lab]] <- sr1
  png(paste0("D:/Dropbox/LightNiches/output/gamViz ",lab, ".png"))
  plot(sr1, residuals = TRUE, pages = 1, pch = 19, cex = .5, shade = TRUE, seWithMean = TRUE, shift = coef(sr)[1])
  dev.off()
  
  pred <- predict(sr1, df_int)
  preddf <-cbind(pred = pred, obs = df_int$light_noscale, mod = lab)
  write.csv(preddf, 
            paste0("D:/Dropbox/LightNiches/output/PrevObs",lab,".csv"), 
            row.names=FALSE)
  sum <- summary(sr1)
  s.table <- cbind(par = row.names(sum$s.table), sum$s.table)
  
  write.csv(s.table, 
            paste0("D:/Dropbox/LightNiches/output/model_summary_s ",lab,"light.csv"),
            row.names=FALSE)
  write.csv(cbind(sum$r.sq,sum$dev.expl,sum$n), 
            paste0("D:/Dropbox/LightNiches/output/model_summary_p ",lab,"light.csv"),
            row.names=FALSE)
  
  lc1 <- gam(light_noscale~Dt*Rl+ watercol12+
               par_integral+
              b + aa +s(lon, lat, bs = "gp", m=2), 
             data = df_int, method="REML")
  lab <-paste0("lc1",dim)
  gam_list[[lab]] <- lc1
  png(paste0("D:/Dropbox/LightNiches/output/gamViz ",lab, ".png"))
  plot(lc1, residuals = TRUE, pages = 1, pch = 19, cex = .5, shade = TRUE, seWithMean = TRUE, shift = coef(lc)[1])
  lc1$var.summary
  lc1$pterms
  
  
  dev.off()
  
  pred <- predict(lc1, df_int)
  preddf <-cbind(pred = pred, obs = df_int$light_noscale, mod = lab)
  write.csv(preddf, 
            paste0("D:/Dropbox/LightNiches/output/PrevObs",lab,".csv"), 
            row.names=FALSE)
  sum <- summary(lc1)
  s.table <- cbind(par = row.names(sum$s.table), sum$s.table)
  
  write.csv(s.table, 
            paste0("D:/Dropbox/LightNiches/output/model_summary_s ",lab,"light.csv"),
            row.names=FALSE)
  write.csv(cbind(sum$r.sq,sum$dev.expl,sum$n), 
            paste0("D:/Dropbox/LightNiches/output/model_summary_p ",lab,"light.csv"),
            row.names=FALSE)
  
  
  lr1 <- gam(light_noscale~Dt*Rl+ watercol12+
              scale(par_integral)+
              b + aa + s(rec, bs = "re"), data = df_int, method="REML")
  lab <-paste0("lr1",dim)
  gam_list[[lab]] <- lr1
  png(paste0("D:/Dropbox/LightNiches/output/gamViz ",lab, ".png"))
  plot(lr1, residuals = TRUE, pages = 1, pch = 19, cex = .5, shade = TRUE, seWithMean = TRUE, shift = coef(lr)[1])
  lr1$var.summary
  lr1$pterms
  
  
  dev.off()
  
  pred <- predict(lr1, df_int)
  preddf <-cbind(pred = pred, obs = df_int$light_noscale, mod = lab)
  write.csv(preddf, 
            paste0("D:/Dropbox/LightNiches/output/PrevObs",lab,".csv"), 
            row.names=FALSE)
  sum <- summary(lr1)
  s.table <- cbind(par = row.names(sum$s.table), sum$s.table)
  
  write.csv(s.table, 
            paste0("D:/Dropbox/LightNiches/output/model_summary_s ",lab,"light.csv"),
            row.names=FALSE)
  write.csv(cbind(sum$r.sq,sum$dev.expl,sum$n), 
            paste0("D:/Dropbox/LightNiches/output/model_summary_p ",lab,"light.csv"),
            row.names=FALSE)
  
  lm1 <- lmer(light_noscale ~ Dt*Rl + watercol12 +
               scale(par_integral) +
               b + aa + (1|rec), data = df_int)
  lab <-paste0("lm1",dim)
  lmer_list[[lab]] <- lm1
  png(paste0("D:/Dropbox/LightNiches/output/",lab, " coeff.png"))
  print(plot_model(lm1))
  dev.off()
  stargazer(lm1, type = "text",
            digits = 3,
            star.cutoffs = c(0.05, 0.01, 0.001),
            digit.separator = "",
            out = paste0("D:/Dropbox/LightNiches/output/",lab, " coeff.html"))
  
  df_int_tot<-rbind(df_int_tot,df_int)
  
  
  #### H*R #####
  
  sc2 <- gam(light_noscale~s(Hl,Rl)+ s(watercol12)+
               s(par_integral)+
               s(b) + s(aa) +s(lon, lat, bs = "gp", m=2),
             data = df_int, method="REML")
  lab <-paste0("sc2",dim)
  gam_list[[lab]] <- sc2
  png(paste0("D:/Dropbox/LightNiches/output/gamViz",lab, ".png"))
  plot(sc2, residuals = TRUE, pages = 1, pch = 19,
       cex = .5, shade = TRUE, seWithMean = TRUE, shift = coef(sc2)[1])
  dev.off()
  
  pred <- predict(sc2, df_int)
  preddf <-cbind(pred = pred, obs = df_int$light_noscale, mod = lab)
  write.csv(preddf, 
            paste0("D:/Dropbox/LightNiches/output/PrevObs",lab,".csv"), 
            row.names=FALSE)
  sum <- summary(sc2)
  s.table <- cbind(par = row.names(sum$s.table), sum$s.table)
  
  write.csv(s.table, 
            paste0("D:/Dropbox/LightNiches/output/model_summary_s ",lab,"light.csv"),
            row.names=FALSE)
  write.csv(cbind(sum$r.sq,sum$dev.expl,sum$n), 
            paste0("D:/Dropbox/LightNiches/output/model_summary_p ",lab,"light.csv"),
            row.names=FALSE)
  
  slc2 <- gam(light_noscale~s(Hl,Rl)+ watercol12+
                par_integral+
                b + aa +s(lon, lat, bs = "gp", m=2),
              data = df_int, method="REML")
  lab <-paste0("slc2",dim)
  gam_list[[lab]] <- slc2
  png(paste0("D:/Dropbox/LightNiches/output/gamViz ",lab, ".png"))
  plot(slc2, residuals = TRUE, pages = 1, pch = 19, cex = .5, shade = TRUE, seWithMean = TRUE, shift = coef(sr)[1])
  dev.off()
  
  pred <- predict(slc2, df_int)
  preddf <-cbind(pred = pred, obs = df_int$light_noscale, mod = lab)
  write.csv(preddf, 
            paste0("D:/Dropbox/LightNiches/output/PrevObs",lab,".csv"), 
            row.names=FALSE)
  sum <- summary(slc2)
  s.table <- cbind(par = row.names(sum$s.table), sum$s.table)
  
  write.csv(s.table, 
            paste0("D:/Dropbox/LightNiches/output/model_summary_s ",lab,"light.csv"),
            row.names=FALSE)
  write.csv(cbind(sum$r.sq,sum$dev.expl,sum$n), 
            paste0("D:/Dropbox/LightNiches/output/model_summary_p ",lab,"light.csv"),
            row.names=FALSE)
  
  
  sr2 <- gam(light_noscale~s(Hl,Rl)+ s(watercol12)+s(par_integral)+
               s(b) + s(aa) + s(rec, bs = "re"), data = df_int, method="REML")
  lab <-paste0("sr2",dim)
  gam_list[[lab]] <- sr2
  png(paste0("D:/Dropbox/LightNiches/output/gamViz ",lab, ".png"))
  plot(sr2, residuals = TRUE, pages = 1, pch = 19, cex = .5, shade = TRUE, seWithMean = TRUE, shift = coef(sr)[1])
  dev.off()
  
  pred <- predict(sr2, df_int)
  preddf <-cbind(pred = pred, obs = df_int$light_noscale, mod = lab)
  write.csv(preddf, 
            paste0("D:/Dropbox/LightNiches/output/PrevObs",lab,".csv"), 
            row.names=FALSE)
  sum <- summary(sr2)
  s.table <- cbind(par = row.names(sum$s.table), sum$s.table)
  
  write.csv(s.table, 
            paste0("D:/Dropbox/LightNiches/output/model_summary_s ",lab,"light.csv"),
            row.names=FALSE)
  write.csv(cbind(sum$r.sq,sum$dev.expl,sum$n), 
            paste0("D:/Dropbox/LightNiches/output/model_summary_p ",lab,"light.csv"),
            row.names=FALSE)
  
  lc2 <- gam(light_noscale~Hl*Rl+ watercol12+
               par_integral+
               b + aa +s(lon, lat, bs = "gp", m=2), 
             data = df_int, method="REML")
  lab <-paste0("lc2",dim)
  gam_list[[lab]] <- lc2
  png(paste0("D:/Dropbox/LightNiches/output/gamViz ",lab, ".png"))
  plot(lc2, residuals = TRUE, pages = 1, pch = 19, cex = .5, shade = TRUE, seWithMean = TRUE, shift = coef(lc)[1])
  lc2$var.summary
  lc2$pterms
  
  
  dev.off()
  
  pred <- predict(lc2, df_int)
  preddf <-cbind(pred = pred, obs = df_int$light_noscale, mod = lab)
  write.csv(preddf, 
            paste0("D:/Dropbox/LightNiches/output/PrevObs",lab,".csv"), 
            row.names=FALSE)
  sum <- summary(lc2)
  s.table <- cbind(par = row.names(sum$s.table), sum$s.table)
  
  write.csv(s.table, 
            paste0("D:/Dropbox/LightNiches/output/model_summary_s ",lab,"light.csv"),
            row.names=FALSE)
  write.csv(cbind(sum$r.sq,sum$dev.expl,sum$n), 
            paste0("D:/Dropbox/LightNiches/output/model_summary_p ",lab,"light.csv"),
            row.names=FALSE)
  
  
  lr2 <- gam(light_noscale~Hl*Rl+ watercol12+
               scale(par_integral)+
               b + aa + s(rec, bs = "re"), data = df_int, method="REML")
  lab <-paste0("lr2",dim)
  gam_list[[lab]] <- lr2
  png(paste0("D:/Dropbox/LightNiches/output/gamViz ",lab, ".png"))
  plot(lr2, residuals = TRUE, pages = 1, pch = 19, cex = .5, shade = TRUE, seWithMean = TRUE, shift = coef(lr)[1])
  lr2$var.summary
  lr2$pterms
  
  
  dev.off()
  
  pred <- predict(lr2, df_int)
  preddf <-cbind(pred = pred, obs = df_int$light_noscale, mod = lab)
  write.csv(preddf, 
            paste0("D:/Dropbox/LightNiches/output/PrevObs",lab,".csv"), 
            row.names=FALSE)
  sum <- summary(lr2)
  s.table <- cbind(par = row.names(sum$s.table), sum$s.table)
  
  write.csv(s.table, 
            paste0("D:/Dropbox/LightNiches/output/model_summary_s ",lab,"light.csv"),
            row.names=FALSE)
  write.csv(cbind(sum$r.sq,sum$dev.expl,sum$n), 
            paste0("D:/Dropbox/LightNiches/output/model_summary_p ",lab,"light.csv"),
            row.names=FALSE)
  
  lm2 <- lmer(light_noscale ~ Hl*Rl + watercol12 +
                scale(par_integral) +
                b + aa + (1|rec), data = df_int)
  lab <-paste0("lm2",dim)
  lmer_list[[lab]] <- lm2
  png(paste0("D:/Dropbox/LightNiches/output/",lab, " coeff.png"))
  print(plot_model(lm2))
  dev.off()
  stargazer(lm2, type = "text",
            digits = 3,
            star.cutoffs = c(0.05, 0.01, 0.001),
            digit.separator = "",
            out = paste0("D:/Dropbox/LightNiches/output/",lab, " coeff.html"))
  
  df_int_tot<-rbind(df_int_tot,df_int)
  
  #### H*D #####
  
  sc3 <- gam(light_noscale~s(Hl,Dt)+ s(watercol12)+
               s(par_integral)+
               s(b) + s(aa) +s(lon, lat, bs = "gp", m=2),
             data = df_int, method="REML")
  lab <-paste0("sc3",dim)
  gam_list[[lab]] <- sc3
  png(paste0("D:/Dropbox/LightNiches/output/gamViz",lab, ".png"))
  plot(sc3, residuals = TRUE, pages = 1, pch = 19,
       cex = .5, shade = TRUE, seWithMean = TRUE, shift = coef(sc3)[1])
  dev.off()
  
  pred <- predict(sc3, df_int)
  preddf <-cbind(pred = pred, obs = df_int$light_noscale, mod = lab)
  write.csv(preddf, 
            paste0("D:/Dropbox/LightNiches/output/PrevObs",lab,".csv"), 
            row.names=FALSE)
  sum <- summary(sc3)
  s.table <- cbind(par = row.names(sum$s.table), sum$s.table)
  
  write.csv(s.table, 
            paste0("D:/Dropbox/LightNiches/output/model_summary_s ",lab,"light.csv"),
            row.names=FALSE)
  write.csv(cbind(sum$r.sq,sum$dev.expl,sum$n), 
            paste0("D:/Dropbox/LightNiches/output/model_summary_p ",lab,"light.csv"),
            row.names=FALSE)
  
  slc3 <- gam(light_noscale~s(Hl,Dt)+ watercol12+
                par_integral+
                b + aa +s(lon, lat, bs = "gp", m=2),
              data = df_int, method="REML")
  lab <-paste0("slc3",dim)
  gam_list[[lab]] <- slc3
  png(paste0("D:/Dropbox/LightNiches/output/gamViz ",lab, ".png"))
  plot(slc3, residuals = TRUE, pages = 1, pch = 19, cex = .5, shade = TRUE, seWithMean = TRUE, shift = coef(sr)[1])
  dev.off()
  
  pred <- predict(slc3, df_int)
  preddf <-cbind(pred = pred, obs = df_int$light_noscale, mod = lab)
  write.csv(preddf, 
            paste0("D:/Dropbox/LightNiches/output/PrevObs",lab,".csv"), 
            row.names=FALSE)
  sum <- summary(slc3)
  s.table <- cbind(par = row.names(sum$s.table), sum$s.table)
  
  write.csv(s.table, 
            paste0("D:/Dropbox/LightNiches/output/model_summary_s ",lab,"light.csv"),
            row.names=FALSE)
  write.csv(cbind(sum$r.sq,sum$dev.expl,sum$n), 
            paste0("D:/Dropbox/LightNiches/output/model_summary_p ",lab,"light.csv"),
            row.names=FALSE)
  
  
  
  sr3 <- gam(light_noscale~s(Hl,Dt)+ s(watercol12)+s(par_integral)+
               s(b) + s(aa) + s(rec, bs = "re"), data = df_int, method="REML")
  lab <-paste0("sr3",dim)
  gam_list[[lab]] <- sr3
  png(paste0("D:/Dropbox/LightNiches/output/gamViz ",lab, ".png"))
  plot(sr3, residuals = TRUE, pages = 1, pch = 19, cex = .5, shade = TRUE, seWithMean = TRUE, shift = coef(sr)[1])
  dev.off()
  
  pred <- predict(sr3, df_int)
  preddf <-cbind(pred = pred, obs = df_int$light_noscale, mod = lab)
  write.csv(preddf, 
            paste0("D:/Dropbox/LightNiches/output/PrevObs",lab,".csv"), 
            row.names=FALSE)
  sum <- summary(sr3)
  s.table <- cbind(par = row.names(sum$s.table), sum$s.table)
  
  write.csv(s.table, 
            paste0("D:/Dropbox/LightNiches/output/model_summary_s ",lab,"light.csv"),
            row.names=FALSE)
  write.csv(cbind(sum$r.sq,sum$dev.expl,sum$n), 
            paste0("D:/Dropbox/LightNiches/output/model_summary_p ",lab,"light.csv"),
            row.names=FALSE)
  
  lc3 <- gam(light_noscale~Hl*Dt+ watercol12+
               par_integral+
               b + aa +s(lon, lat, bs = "gp", m=2), 
             data = df_int, method="REML")
  lab <-paste0("lc3",dim)
  gam_list[[lab]] <- lc3
  png(paste0("D:/Dropbox/LightNiches/output/gamViz ",lab, ".png"))
  plot(lc3, residuals = TRUE, pages = 1, pch = 19, cex = .5, shade = TRUE, seWithMean = TRUE, shift = coef(lc)[1])
  lc3$var.summary
  lc3$pterms
  
  
  dev.off()
  
  pred <- predict(lc3, df_int)
  preddf <-cbind(pred = pred, obs = df_int$light_noscale, mod = lab)
  write.csv(preddf, 
            paste0("D:/Dropbox/LightNiches/output/PrevObs",lab,".csv"), 
            row.names=FALSE)
  sum <- summary(lc3)
  s.table <- cbind(par = row.names(sum$s.table), sum$s.table)
  
  write.csv(s.table, 
            paste0("D:/Dropbox/LightNiches/output/model_summary_s ",lab,"light.csv"),
            row.names=FALSE)
  write.csv(cbind(sum$r.sq,sum$dev.expl,sum$n), 
            paste0("D:/Dropbox/LightNiches/output/model_summary_p ",lab,"light.csv"),
            row.names=FALSE)
  
  
  lr3 <- gam(light_noscale~Hl*Dt+ watercol12+
               scale(par_integral)+
               b + aa + s(rec, bs = "re"), data = df_int, method="REML")
  lab <-paste0("lr3",dim)
  gam_list[[lab]] <- lr3
  png(paste0("D:/Dropbox/LightNiches/output/gamViz ",lab, ".png"))
  plot(lr3, residuals = TRUE, pages = 1, pch = 19, cex = .5, shade = TRUE, seWithMean = TRUE, shift = coef(lr)[1])
  lr3$var.summary
  lr3$pterms
  
  
  dev.off()
  
  pred <- predict(lr3, df_int)
  preddf <-cbind(pred = pred, obs = df_int$light_noscale, mod = lab)
  write.csv(preddf, 
            paste0("D:/Dropbox/LightNiches/output/PrevObs",lab,".csv"), 
            row.names=FALSE)
  sum <- summary(lr3)
  s.table <- cbind(par = row.names(sum$s.table), sum$s.table)
  
  write.csv(s.table, 
            paste0("D:/Dropbox/LightNiches/output/model_summary_s ",lab,"light.csv"),
            row.names=FALSE)
  write.csv(cbind(sum$r.sq,sum$dev.expl,sum$n), 
            paste0("D:/Dropbox/LightNiches/output/model_summary_p ",lab,"light.csv"),
            row.names=FALSE)
  
  lm3 <- lmer(light_noscale ~ Hl*Dt + watercol12 +
                scale(par_integral) +
                b + aa + (1|rec), data = df_int)
  lab <-paste0("lm3",dim)
  lmer_list[[lab]] <- lm3
  png(paste0("D:/Dropbox/LightNiches/output/",lab, " coeff.png"))
  print(plot_model(lm3))
  dev.off()
  stargazer(lm3, type = "text",
            digits = 3,
            star.cutoffs = c(0.05, 0.01, 0.001),
            digit.separator = "",
            out = paste0("D:/Dropbox/LightNiches/output/",lab, " coeff.html"))
  
  df_int_tot<-rbind(df_int_tot,df_int)
} 
#closes dim loop

## gam models
#### look up results #####
lc150 <- gam_list[["lc10.5"]]
lc125 <- gam_list[["lc10.25"]]
lc175 <- gam_list[["lc10.75"]]
lr150 <- gam_list[["lr10.5"]]
lr125 <- gam_list[["lr10.25"]]
lr175 <- gam_list[["lr10.75"]]
sc150 <- gam_list[["sc10.5"]]
sc125 <- gam_list[["sc10.25"]]
sc175 <- gam_list[["sc10.75"]]
sr150 <- gam_list[["sr10.5"]]
sr125 <- gam_list[["sr10.25"]]
sr175 <- gam_list[["sr10.75"]]

slc150 <- gam_list[["slc10.5"]]
slc125 <- gam_list[["slc10.25"]]
slc175 <- gam_list[["slc10.75"]]
slc250 <- gam_list[["slc20.5"]]
slc225 <- gam_list[["slc20.25"]]
slc275 <- gam_list[["slc20.75"]]
slc350 <- gam_list[["slc30.5"]]
slc325 <- gam_list[["slc30.25"]]
slc375 <- gam_list[["slc30.75"]]

sr150 <- gam_list[["sr10.5"]]
sr125 <- gam_list[["sr10.25"]]
sr175 <- gam_list[["sr10.75"]]

lm125 <- lmer_list[["lm10.25"]]
lm150 <- lmer_list[["lm10.5"]]
lm175 <- lmer_list[["lm10.75"]]

lc250 <- gam_list[["lc20.5"]]
lc225 <- gam_list[["lc20.25"]]
lc275 <- gam_list[["lc20.75"]]
lr250 <- gam_list[["lr20.5"]]
lr225 <- gam_list[["lr20.25"]]
lr275 <- gam_list[["lr20.75"]]
sc250 <- gam_list[["sc20.5"]]
sc225 <- gam_list[["sc20.25"]]
sc275 <- gam_list[["sc20.75"]]
sr250 <- gam_list[["sr20.5"]]
sr225 <- gam_list[["sr20.25"]]
sr275 <- gam_list[["sr20.75"]]

lm225 <- lmer_list[["lm20.25"]]
lm250 <- lmer_list[["lm20.5"]]
lm275 <- lmer_list[["lm20.75"]]

lc350 <- gam_list[["lc30.5"]]
lc325 <- gam_list[["lc30.25"]]
lc375 <- gam_list[["lc30.75"]]
lr350 <- gam_list[["lr30.5"]]
lr325 <- gam_list[["lr30.25"]]
lr375 <- gam_list[["lr30.75"]]
sc350 <- gam_list[["sc30.5"]]
sc325 <- gam_list[["sc30.25"]]
sc375 <- gam_list[["sc30.75"]]
sr350 <- gam_list[["sr30.5"]]
sr325 <- gam_list[["sr30.25"]]
sr375 <- gam_list[["sr30.75"]]

lm325 <- lmer_list[["lm30.25"]]
lm350 <- lmer_list[["lm30.5"]]
lm375 <- lmer_list[["lm30.75"]]

## coeff tables ####

##lc####
lc_coefs <- tab_model(lc125, lc150, lc175,lc225, lc250, lc275, lc325, lc350, lc375, show.stat = TRUE, collapse.se = TRUE, show.ci = FALSE,
                      # pred.labels = c("intercept", "fractal \ndimension", "surface\nrugosity", "water column",
                      #                 "surface \nPAR integral", "northing", "sloping", "D*R", "s(lon,lat)"),
                      dv.labels = c("25-cm side predictors", "50-cm side predictors","75-cm side predictors",
                                    "25-cm side predictors", "50-cm side predictors","75-cm side predictors",
                                    "25-cm side predictors", "50-cm side predictors","75-cm side predictors"),
                      string.pred = "coefficients",
                      string.est = "estimates\n(SE)",
                      string.p = "p-value",
                      string.stat = "statistic", 
                      title = "linear predictors + coords")

lc_coefs1 <- tab_model(lc125, lc150, lc175, show.stat = TRUE, collapse.se = TRUE, show.ci = FALSE,
                       # pred.labels = c("intercept", "fractal \ndimension", "surface\nrugosity", "water column",
                       #                 "surface \nPAR integral", "northing", "sloping", "D*R", "s(lon,lat)"),
                       dv.labels = c("25-cm side predictors", "50-cm side predictors","75-cm side predictors"),
                       string.pred = "coefficients",
                       string.est = "estimates\n(SE)",
                       string.p = "p-value",
                       string.stat = "statistic")
lc_coefs2 <- tab_model(lc225, lc250, lc275, show.stat = TRUE, collapse.se = TRUE, show.ci = FALSE,
                       # pred.labels = c("intercept", "fractal \ndimension", "surface\nrugosity", "water column",
                       #                 "surface \nPAR integral", "northing", "sloping", "D*R", "s(lon,lat)"),
                       dv.labels = c("25-cm side predictors", "50-cm side predictors","75-cm side predictors"),
                       string.pred = "coefficients",
                       string.est = "estimates\n(SE)",
                       string.p = "p-value",
                       string.stat = "statistic")
lc_coefs3 <- tab_model(lc325, lc350, lc375, show.stat = TRUE, collapse.se = TRUE, show.ci = FALSE,
                       # pred.labels = c("intercept", "fractal \ndimension", "surface\nrugosity", "water column",
                                       # "surface \nPAR integral", "northing", "sloping", "D*R", "s(lon,lat)"),
                       dv.labels = c("25-cm side predictors", "50-cm side predictors","75-cm side predictors"),
                       string.pred = "coefficients",
                       string.est = "estimates\n(SE)",
                       string.p = "p-value",
                       string.stat = "statistic")

##sc####

sc_coefs <- tab_model(sc125, sc150, sc175,
                      sc225, sc250, sc275,
                      sc325, sc350, sc375,
                      show.stat = TRUE, show.df = TRUE, show.ci = FALSE,
                      # pred.labels = c("intercept", "fractal \ndimension", "surface\nrugosity", "mean water\ncolumn",
                      #                 "water column\nat noon", "height\nrange","surface \nPAR integral", "northing", "sloping", "s(lon,lat)"),
                      dv.labels = c("25-cm side predictors", "50-cm side predictors","75-cm side predictors",
                                    "25-cm side predictors", "50-cm side predictors","75-cm side predictors",
                                    "25-cm side predictors", "50-cm side predictors","75-cm side predictors"),
                      string.pred = "coefficients",
                      string.est = "estimates",
                      string.p = "p-value",
                      string.df = "edf",
                      string.stat = "F", 
                      title = "smoothed predictors + coords")
##slc####
slc_coefs <- tab_model(slc125, slc150, slc175,
                      slc225, slc250, slc275, 
                      slc325, slc350, slc375,
                      show.stat = TRUE, collapse.se = TRUE, show.ci = FALSE,
                      # pred.labels = c("intercept", "fractal \ndimension", "surface\nrugosity", "water column",
                      #                 "surface \nPAR integral", "northing", "sloping", "D*R", "s(lon,lat)"),
                      dv.labels = c("25-cm side predictors", "50-cm side predictors","75-cm side predictors",
                                    "25-cm side predictors", "50-cm side predictors","75-cm side predictors",
                                    "25-cm side predictors", "50-cm side predictors","75-cm side predictors"),
                      #string.pred = "coefficients",
                      #string.est = "estimates\n(SE)",
                      #string.p = "p-value",
                      string.stat = "statistic", 
                      title = "smooted and linear predictors + coords")

###lm ####
lm_coefs <- tab_model(lm125, lm150, lm175,
                      lm225, lm250, lm275,
                      lm325, lm350, lm375,
                      show.stat = TRUE, show.ci = FALSE, collapse.se = TRUE, 
                      # pred.labels = c("intercept", "fractal \ndimension", "surface\nrugosity", "mean water\ncolumn",
                      #                 "water column\nat noon", "height\nrange","surface \nPAR integral", "northing", "sloping"),
                      dv.labels = c("25-cm side predictors", "50-cm side predictors","75-cm side predictors",
                                    "25-cm side predictors", "50-cm side predictors","75-cm side predictors",
                                    "25-cm side predictors", "50-cm side predictors","75-cm side predictors"),
                      string.pred = "coefficients",
                      string.est = "estimates",
                      string.p = "p-value",
                      string.stat = "t test", 
                      title = "linear predictor + site (lm)")
####lr ####
lr_coefs <- tab_model(lr125, lr150, lr175,
                      lr225, lr250, lr275,
                      lr325, lr350, lr375,
                      show.stat = TRUE, collapse.se = TRUE, show.ci = FALSE,
                      # pred.labels = c("intercept", "fractal \ndimension", "surface\nrugosity", "mean water\ncolumn",
                      #                 "water column\nat noon", "height\nrange","surface \nPAR integral", "northing", "sloping", "s(lon,lat)"),
                      dv.labels = c("25-cm side predictors", "50-cm side predictors","75-cm side predictors",
                                    "25-cm side predictors", "50-cm side predictors","75-cm side predictors",
                                    "25-cm side predictors", "50-cm side predictors","75-cm side predictors"),
                      string.pred = "coefficients",
                      string.est = "estimates\n(SE)",
                      string.p = "p-value",
                      string.stat = "statistic", 
                      title = "linear predictor + site")
###sr ####
sr_coefs <- tab_model(sr125, sr150, sr175,
                      sr225, sr250, sr275,
                      sr325, sr350, sr375,
                      show.stat = TRUE, show.df = TRUE, show.ci = FALSE,
                      # pred.labels = c("intercept", "fractal \ndimension", "surface\nrugosity", "mean water\ncolumn",
                      #                 "water column\nat noon", "height\nrange","surface \nPAR integral", "northing", "sloping", "s(lon,lat)"),
                      dv.labels = c("25-cm side predictors", "50-cm side predictors","75-cm side predictors",
                                    "25-cm side predictors", "50-cm side predictors","75-cm side predictors",
                                    "25-cm side predictors", "50-cm side predictors","75-cm side predictors"),
                      string.pred = "coefficients",
                      string.est = "estimates",
                      string.p = "p-value",
                      string.df = "edf",
                      string.stat = "F", 
                      title = "smoothed predictor + site")


##### plots ####
####getViz ####

## slc ####
slc50v <- getViz(slc150)
png("D:/Dropbox/LightNiches/output/viz gam slc150.png",
    width = 600, height = 400)
print(plot(slc50v, allTerms = T, select = c(1:8),residuals = TRUE)+ 
        l_fitLine(colour = "black") +l_ciPoly() +
        l_ciLine(colour = "black", linetype = 2) + 
        l_points(shape = 19, size = 1, alpha = 0.3) + theme_classic(), pages = 1) # Calls print.plotGam()
dev.off()
slc25v <- getViz(slc125)
png("D:/Dropbox/LightNiches/output/viz gam slc125.png",
    width = 600, height = 400)
print(plot(slc25v, allTerms = T, select = c(1:8),residuals = TRUE)+ 
        l_fitLine(colour = "black") +l_ciPoly() +
        l_ciLine( colour = "black", linetype = 2) + 
        l_points(shape = 19, size = 1, alpha = 0.3) + theme_classic(), pages = 1) # Calls print.plotGam()
dev.off()
png("D:/Dropbox/LightNiches/output/viz gam slc175.png",
    width = 600, height = 400)
slc75v <- getViz(slc175)
print(plot(slc75v, allTerms = T, select = c(1:8),residuals = TRUE)+ 
        l_fitLine(colour = "black") +l_ciPoly() +
        l_ciLine( colour = "black", linetype = 2) + 
        l_points(shape = 19, size = 1, alpha = 0.3) + theme_classic(), pages = 1) # Calls print.plotGam()

dev.off()
png("D:/Dropbox/LightNiches/output/viz gam lc75.png",
    width = 600, height = 400)
lc75v <- getViz(lc75)
print(plot(lc75v, allTerms = TRUE, select = c(2:9))+ 
        l_fitLine(colour = "black") +l_ciPoly() +
        l_ciLine( colour = "black", linetype = 2) + 
        l_points(shape = 19, size = 1, alpha = 0.3) + theme_classic(), pages = 1) # Calls print.plotGam()
dev.off()
png("D:/Dropbox/LightNiches/output/viz gam lc55.png",
    width = 600, height = 400)
lc50v <- getViz(lc50)
print(plot(lc50v, allTerms = TRUE, select = c(2:9))+ 
        l_fitLine(colour = "black") +l_ciPoly() +
        l_ciLine( colour = "black", linetype = 2) + 
        l_points(shape = 19, size = 1, alpha = 0.3) + theme_classic(), pages = 1) # Calls print.plotGam()
dev.off()
png("D:/Dropbox/LightNiches/output/viz gam lc25.png",
    width = 600, height = 400)
lc25v <- getViz(lc25)
print(plot(lc25v, allTerms = TRUE, select = c(2:9))+ 
        l_fitLine(colour = "black") +l_ciPoly() +
        l_ciLine( colour = "black", linetype = 2) + 
        l_points(shape = 19, size = 1, alpha = 0.3) + theme_classic(), pages = 1) # Calls print.plotGam()



DR50 <- plot( sm(slc50v, 1))
DR75 <- plot( sm(slc75v, 1))
DR25 <- plot( sm(slc25v, 1))

DR25
DR75 
DR50

R50 <- plot( sm(slc50v,2)) + l_ciPoly() +
  l_fitLine(colour = "black") +
  geom_smooth(method = "lm", color = "#35978F")+
  #l_ciLine(colour = "black", linetype = 2) + 
  l_points(shape = 19, size = 1, alpha = 0.3) + theme_classic() + ylab("s(R)") +xlab("Surface rugosity (R)")

mw50 <- plot( sm(sc50v, 3)) +  l_ciPoly() +
  l_fitLine(colour = "black") +
  geom_smooth(method = "lm", color = "#35978F")+
  #l_ciLine(colour = "black", linetype = 2) + 
  l_points(shape = 19, size = 1, alpha = 0.3) + theme_classic() + ylab("s(mw)") +xlab("water column height (mw)")

w1250 <- plot( sm(sc50v, 4)) + l_ciPoly() +
  l_fitLine(colour = "black") + 
  geom_smooth(method = "lm", color = "#35978F")+
  #l_ciLine(colour = "black", linetype = 2) + 
  l_points(shape = 19, size = 1, alpha = 0.3) + theme_classic() + ylab("s(w12)") +xlab("water column at noon (w12)")

H50 <- plot( sm(sc50v, 5)) + l_ciPoly() +
  l_fitLine(colour = "black") + 
  geom_smooth(method = "lm", color = "#35978F")+
  #l_ciLine(colour = "black", linetype = 2) + 
  l_points(shape = 19, size = 1, alpha = 0.3) + theme_classic() + ylab("s(H)") +xlab("height range (H)")

PAR50 <- plot( sm(sc50v, 6)) + l_ciPoly() +
  l_fitLine(colour = "black") + 
  geom_smooth(method = "lm", color = "#35978F")+
  #l_ciLine(colour = "black", linetype = 2) + 
  l_points(shape = 19, size = 1, alpha = 0.3) + theme_classic() + ylab("s(PAR)") +xlab("surface PAR integral (PAR)")

n50 <- plot( sm(sc50v, 7)) + l_ciPoly() +
  l_fitLine(colour = "black") + 
  geom_smooth(method = "lm", color = "#35978F")+
  #l_ciLine(colour = "black", linetype = 2) + 
  l_points(shape = 19, size = 1, alpha = 0.3) + theme_classic() + ylab("s(b)") +xlab("northing (b)")

s50 <- plot( sm(sc50v, 8)) + l_ciPoly() +
  l_fitLine(colour = "black") + 
  geom_smooth(method = "lm", color = "#35978F")+
  #l_ciLine(colour = "black", linetype = 2) + 
  l_points(shape = 19, size = 1, alpha = 0.3) + theme_classic() + ylab("s(a)") +xlab("sloping(a)")
png("D:/Dropbox/LightNiches/output/vizgam sc 50 chapt.png",
    width = 550, height = 700)
gridPrint(D50,n50,R50,s50,H50,PAR50,mw50,w1250, ncol = 2)
dev.off()

check.gamViz(sc50v,
             a.qq = list(method = "tnorm", 
                         a.cipoly = list(fill = "light blue")), 
             a.respoi = list(size = 0.5), 
             a.hist = list(bins = 10))

#### light in the geometric space #####
png("D:/Dropbox/LightNiches/output/light in the geometric space slc.png",
    width = 700, height = 600)
par(mfrow=c(3,3), cex = 1.1, mai=c(0.9,1,0.4,0.4))

vis.gam(slc125, view=c("Dt", "Rl"), type="response", plot.type="contour", main = "25cm-side\n",
        color = "gray", ylab = expression(paste("log(R "[25], ")")), xlab = expression(paste("D"[25])))
points(df_int_tot[df_int_tot$dim == 0.25,]$Rl~df_int_tot[df_int_tot$dim == 0.25,]$Dt, pch = 19, cex = .6, 
       col = rgb(red=0, green=0, blue=0, alpha=0.3))
vis.gam(slc150, view=c("Dt", "Rl"), type="response", plot.type="contour", main = "50cm-side\n",
        color = "gray", ylab = expression(paste("log(R "[50], ")")), xlab = expression(paste("D"[50])))
points(df_int_tot[df_int_tot$dim == 0.5,]$Rl~df_int_tot[df_int_tot$dim == 0.5,]$Dt, pch = 19, cex = .6, 
       col = rgb(red=0, green=0, blue=0, alpha=0.3))
vis.gam(slc175, view=c("Dt", "Rl"), type="response", plot.type="contour", main = "75cm-side\n",
        color = "gray", ylab = expression(paste("log(R "[75], ")")), xlab = expression(paste("D"[75])))
points(df_int_tot[df_int_tot$dim == 0.75,]$Rl~df_int_tot[df_int_tot$dim == 0.75,]$Dt, pch = 19, cex = .6, 
       col = rgb(red=0, green=0, blue=0, alpha=0.3))

vis.gam(slc225, view=c("Hl", "Rl"), type="response", plot.type="contour", main = "",
        color = "gray", ylab = expression(paste("log(H"[25], ")")), xlab = expression(paste("log(R"[25], ")")))
points(df_int_tot[df_int_tot$dim == 0.25,]$Rl~df_int_tot[df_int_tot$dim == 0.25,]$Hl, pch = 19, cex = .6, 
       col = rgb(red=0, green=0, blue=0, alpha=0.3))
vis.gam(slc250, view=c("Hl", "Rl"), type="response", plot.type="contour", main = "",
        color = "gray", ylab = expression(paste("log(H"[50], ")")), xlab = expression(paste("log(R"[50], ")")))
points(df_int_tot[df_int_tot$dim == 0.5,]$Rl~df_int_tot[df_int_tot$dim == 0.5,]$Hl, pch = 19, cex = .6, 
       col = rgb(red=0, green=0, blue=0, alpha=0.3))
vis.gam(slc275, view=c("Hl", "Rl"), type="response", plot.type="contour", main = "",
        color = "gray", ylab = expression(paste("log(H"[75], ")")), xlab = expression(paste("log(R"[75], ")")))
points(df_int_tot[df_int_tot$dim == 0.75,]$Rl~df_int_tot[df_int_tot$dim == 0.75,]$Hl, pch = 19, cex = .6, 
       col = rgb(red=0, green=0, blue=0, alpha=0.3))

vis.gam(slc325, view=c("Dt", "Hl"), type="response", plot.type="contour",main = "",
        color = "gray", xlab = expression(paste("D"[25], "")), ylab = expression(paste("log(H "[25], ")")))
points(df_int_tot[df_int_tot$dim == 0.25,]$Hl~df_int_tot[df_int_tot$dim == 0.25,]$Dt, pch = 19, cex = .6, 
       col = rgb(red=0, green=0, blue=0, alpha=0.3))
vis.gam(slc350, view=c("Dt", "Hl"), type="response", plot.type="contour", main = "",
        color = "gray", xlab = expression(paste("D"[50], "")), ylab = expression(paste("log(H "[50], ")")))
points(df_int_tot[df_int_tot$dim == 0.5,]$Hl~df_int_tot[df_int_tot$dim == 0.5,]$Dt, pch = 19, cex = .6, 
       col = rgb(red=0, green=0, blue=0, alpha=0.3))
vis.gam(slc375, view=c("Dt", "Hl"), type="response", plot.type="contour", main = "",
        color = "gray", xlab = expression(paste("D"[75], "")), ylab = expression(paste("log(H "[75], ")")))
points(df_int_tot[df_int_tot$dim == 0.75,]$Hl~df_int_tot[df_int_tot$dim == 0.75,]$Dt, pch = 19, cex = .6, 
       col = rgb(red=0, green=0, blue=0, alpha=0.3))

dev.off()


png("D:/Dropbox/LightNiches/output/light in the geometric space sc.png",
    width = 700, height = 600)
par(mfrow=c(3,3), cex = 1.1, mai=c(0.9,1,0.4,0.4))

vis.gam(sc25, view=c("Hl", "Rl"), type="response", plot.type="contour", main = "25cm-side\n",
        color = "gray", ylab = expression(paste("log(R "[25], ")")), xlab = expression(paste("log(H"[25], ")")))
points(df_int_tot[df_int_tot$dim == 0.25,]$Rl~df_int_tot[df_int_tot$dim == 0.25,]$Hl, pch = 19, cex = .6, 
       col = rgb(red=0, green=0, blue=0, alpha=0.3))
vis.gam(sc50, view=c("Dt", "Rl"), type="response", plot.type="contour", main = "50cm-side\n",
        color = "gray", ylab = expression(paste("log(R "[50], ")")), xlab = expression(paste("log(H"[50], ")")))
points(df_int_tot[df_int_tot$dim == 0.5,]$Rl~df_int_tot[df_int_tot$dim == 0.5,]$Dt, pch = 19, cex = .6, 
       col = rgb(red=0, green=0, blue=0, alpha=0.3))
vis.gam(sc75, view=c("Hl", "Rl"), type="response", plot.type="contour", main = "75cm-side\n",
        color = "gray", ylab = expression(paste("log(R "[75], ")")), xlab = expression(paste("log(H"[75], ")")))
points(df_int_tot[df_int_tot$dim == 0.75,]$Rl~df_int_tot[df_int_tot$dim == 0.75,]$Hl, pch = 19, cex = .6, 
       col = rgb(red=0, green=0, blue=0, alpha=0.3))

vis.gam(sc25, view=c("Dt", "Rl"), type="response", plot.type="contour", main = "",
        color = "gray", ylab = expression(paste("D "[25], "")), xlab = expression(paste("log(R"[25], ")")))
points(df_int_tot[df_int_tot$dim == 0.25,]$Rl~df_int_tot[df_int_tot$dim == 0.25,]$Dt, pch = 19, cex = .6, 
       col = rgb(red=0, green=0, blue=0, alpha=0.3))
vis.gam(sc50, view=c("Dt", "Rl"), type="response", plot.type="contour", main = "",
        color = "gray", ylab = expression(paste("D "[50], "")), xlab = expression(paste("log(R"[50], ")")))
points(df_int_tot[df_int_tot$dim == 0.5,]$Rl~df_int_tot[df_int_tot$dim == 0.5,]$Dt, pch = 19, cex = .6, 
       col = rgb(red=0, green=0, blue=0, alpha=0.3))
vis.gam(sc75, view=c("Dt", "Rl"), type="response", plot.type="contour", main = "",
        color = "gray", ylab = expression(paste("D "[75], "")), xlab = expression(paste("log(R"[75], ")")))
points(df_int_tot[df_int_tot$dim == 0.75,]$Rl~df_int_tot[df_int_tot$dim == 0.75,]$Dt, pch = 19, cex = .6, 
       col = rgb(red=0, green=0, blue=0, alpha=0.3))

vis.gam(sc25, view=c("Hl", "Dt"), type="response", plot.type="contour",main = "",
        color = "gray", xlab = expression(paste("D"[25], "")), ylab = expression(paste("log(H "[25], ")")))
points(df_int_tot[df_int_tot$dim == 0.25,]$Dt~df_int_tot[df_int_tot$dim == 0.25,]$Hl, pch = 19, cex = .6, 
       col = rgb(red=0, green=0, blue=0, alpha=0.3))
vis.gam(sc50, view=c("Hl", "Dt"), type="response", plot.type="contour", main = "",
        color = "gray", xlab = expression(paste("D"[50], "")), ylab = expression(paste("log(H "[50], ")")))
points(df_int_tot[df_int_tot$dim == 0.5,]$Dt~df_int_tot[df_int_tot$dim == 0.5,]$Hl, pch = 19, cex = .6, 
       col = rgb(red=0, green=0, blue=0, alpha=0.3))
vis.gam(sc75, view=c("Hl", "Dt"), type="response", plot.type="contour", main = "",
        color = "gray", xlab = expression(paste("D"[75], "")), ylab = expression(paste("log(H "[75], ")")))
points(df_int_tot[df_int_tot$dim == 0.75,]$Dt~df_int_tot[df_int_tot$dim == 0.75,]$Hl, pch = 19, cex = .6, 
       col = rgb(red=0, green=0, blue=0, alpha=0.3))

dev.off()


png("D:/Dropbox/LightNiches/output/light in the geometric space linear predictors lc.png",
    width = 700, height = 600)
par(mfrow=c(3,3), cex = 1.1, mai=c(0.9,1,0.4,0.4))

vis.gam(lc25, view=c("Hl", "Rl"), type="response", plot.type="contour", main = "25cm-side\n",
        color = "gray", ylab = expression(paste("log(R "[25], ")")), xlab = expression(paste("log(H"[25], ")")))
points(df_int_tot[df_int_tot$dim == 0.25,]$Rl~df_int_tot[df_int_tot$dim == 0.25,]$Hl, pch = 19, cex = .6, 
       col = rgb(red=0, green=0, blue=0, alpha=0.3))
vis.gam(lc50, view=c("Hl", "Rl"), type="response", plot.type="contour", main = "50cm-side\n",
        color = "gray", ylab = expression(paste("log(R "[50], ")")), xlab = expression(paste("log(H"[50], ")")))
points(df_int_tot[df_int_tot$dim == 0.5,]$Rl~df_int_tot[df_int_tot$dim == 0.5,]$Hl, pch = 19, cex = .6, 
       col = rgb(red=0, green=0, blue=0, alpha=0.3))
vis.gam(lc75, view=c("Hl", "Rl"), type="response", plot.type="contour", main = "75cm-side\n",
        color = "gray", ylab = expression(paste("log(R "[75], ")")), xlab = expression(paste("log(H"[75], ")")))
points(df_int_tot[df_int_tot$dim == 0.75,]$Rl~df_int_tot[df_int_tot$dim == 0.75,]$Hl, pch = 19, cex = .6, 
       col = rgb(red=0, green=0, blue=0, alpha=0.3))

vis.gam(lc25, view=c("Dt", "Rl"), type="response", plot.type="contour", main = "",
        color = "gray", ylab = expression(paste("D "[25], "")), xlab = expression(paste("log(R"[25], ")")))
points(df_int_tot[df_int_tot$dim == 0.25,]$Rl~df_int_tot[df_int_tot$dim == 0.25,]$Dt, pch = 19, cex = .6, 
       col = rgb(red=0, green=0, blue=0, alpha=0.3))
vis.gam(lc50, view=c("Dt", "Rl"), type="response", plot.type="contour", main = "",
        color = "gray", ylab = expression(paste("D "[50], "")), xlab = expression(paste("log(R"[50], ")")))
points(df_int_tot[df_int_tot$dim == 0.5,]$Rl~df_int_tot[df_int_tot$dim == 0.5,]$Dt, pch = 19, cex = .6, 
       col = rgb(red=0, green=0, blue=0, alpha=0.3))
vis.gam(lc75, view=c("Dt", "Rl"), type="response", plot.type="contour", main = "",
        color = "gray", ylab = expression(paste("D "[75], "")), xlab = expression(paste("log(R"[75], ")")))
points(df_int_tot[df_int_tot$dim == 0.75,]$Rl~df_int_tot[df_int_tot$dim == 0.75,]$Dt, pch = 19, cex = .6, 
       col = rgb(red=0, green=0, blue=0, alpha=0.3))

vis.gam(lc25, view=c("Hl", "Dt"), type="response", plot.type="contour",main = "",
        color = "gray", xlab = expression(paste("D"[25], "")), ylab = expression(paste("log(H "[25], ")")))
points(df_int_tot[df_int_tot$dim == 0.25,]$Dt~df_int_tot[df_int_tot$dim == 0.25,]$Hl, pch = 19, cex = .6, 
       col = rgb(red=0, green=0, blue=0, alpha=0.3))
vis.gam(lc50, view=c("Hl", "Dt"), type="response", plot.type="contour", main = "",
        color = "gray", xlab = expression(paste("D"[50], "")), ylab = expression(paste("log(H "[50], ")")))
points(df_int_tot[df_int_tot$dim == 0.5,]$Dt~df_int_tot[df_int_tot$dim == 0.5,]$Hl, pch = 19, cex = .6, 
       col = rgb(red=0, green=0, blue=0, alpha=0.3))
vis.gam(lc75, view=c("Hl", "Dt"), type="response", plot.type="contour", main = "",
        color = "gray", xlab = expression(paste("D"[75], "")), ylab = expression(paste("log(H "[75], ")")))
points(df_int_tot[df_int_tot$dim == 0.75,]$Dt~df_int_tot[df_int_tot$dim == 0.75,]$Hl, pch = 19, cex = .6, 
       col = rgb(red=0, green=0, blue=0, alpha=0.3))

dev.off()

png("D:/Dropbox/LightNiches/output/light in the geometric space linear predictors sr.png",
    width = 700, height = 500)
par(mfrow=c(3,3), cex = 1.1, mai=c(0.9,1,0.4,0.4))

vis.gam(sr25, view=c("Hl", "Rl"), type="response", plot.type="contour", main = "25cm-side\n",
        color = "gray", ylab = expression(paste("log(R "[25], ")")), xlab = expression(paste("log(H"[25], ")")))
points(df_int_tot[df_int_tot$dim == 0.25,]$Rl~df_int_tot[df_int_tot$dim == 0.25,]$Hl, pch = 19, cex = .6, 
       col = rgb(red=0, green=0, blue=0, alpha=0.3))
vis.gam(sr50, view=c("Hl", "Rl"), type="response", plot.type="contour", main = "50cm-side\n",
        color = "gray", ylab = expression(paste("log(R "[50], ")")), xlab = expression(paste("log(H"[50], ")")))
points(df_int_tot[df_int_tot$dim == 0.5,]$Rl~df_int_tot[df_int_tot$dim == 0.5,]$Hl, pch = 19, cex = .6, 
       col = rgb(red=0, green=0, blue=0, alpha=0.3))
vis.gam(sr75, view=c("Hl", "Rl"), type="response", plot.type="contour", main = "75cm-side\n",
        color = "gray", ylab = expression(paste("log(R "[75], ")")), xlab = expression(paste("log(H"[75], ")")))
points(df_int_tot[df_int_tot$dim == 0.75,]$Rl~df_int_tot[df_int_tot$dim == 0.75,]$Hl, pch = 19, cex = .6, 
       col = rgb(red=0, green=0, blue=0, alpha=0.3))

vis.gam(sr25, view=c("Dt", "Rl"), type="response", plot.type="contour", main = "",
        color = "gray", ylab = expression(paste("D "[25], "")), xlab = expression(paste("log(R"[25], ")")))
points(df_int_tot[df_int_tot$dim == 0.25,]$Rl~df_int_tot[df_int_tot$dim == 0.25,]$Dt, pch = 19, cex = .6, 
       col = rgb(red=0, green=0, blue=0, alpha=0.3))
vis.gam(sr50, view=c("Dt", "Rl"), type="response", plot.type="contour", main = "",
        color = "gray", ylab = expression(paste("D "[50], "")), xlab = expression(paste("log(R"[50], ")")))
points(df_int_tot[df_int_tot$dim == 0.5,]$Rl~df_int_tot[df_int_tot$dim == 0.5,]$Dt, pch = 19, cex = .6, 
       col = rgb(red=0, green=0, blue=0, alpha=0.3))
vis.gam(sr75, view=c("Dt", "Rl"), type="response", plot.type="contour", main = "",
        color = "gray", ylab = expression(paste("D "[75], "")), xlab = expression(paste("log(R"[75], ")")))
points(df_int_tot[df_int_tot$dim == 0.75,]$Rl~df_int_tot[df_int_tot$dim == 0.75,]$Dt, pch = 19, cex = .6, 
       col = rgb(red=0, green=0, blue=0, alpha=0.3))

vis.gam(sr25, view=c("Hl", "Dt"), type="response", plot.type="contour",main = "",
        color = "gray", xlab = expression(paste("D"[25], "")), ylab = expression(paste("log(H "[25], ")")))
points(df_int_tot[df_int_tot$dim == 0.25,]$Dt~df_int_tot[df_int_tot$dim == 0.25,]$Hl, pch = 19, cex = .6, 
       col = rgb(red=0, green=0, blue=0, alpha=0.3))
vis.gam(sr50, view=c("Hl", "Dt"), type="response", plot.type="contour", main = "",
        color = "gray", xlab = expression(paste("D"[50], "")), ylab = expression(paste("log(H "[50], ")")))
points(df_int_tot[df_int_tot$dim == 0.5,]$Dt~df_int_tot[df_int_tot$dim == 0.5,]$Hl, pch = 19, cex = .6, 
       col = rgb(red=0, green=0, blue=0, alpha=0.3))
vis.gam(sr75, view=c("Hl", "Dt"), type="response", plot.type="contour", main = "",
        color = "gray", xlab = expression(paste("D"[75], "")), ylab = expression(paste("log(H "[75], ")")))
points(df_int_tot[df_int_tot$dim == 0.75,]$Dt~df_int_tot[df_int_tot$dim == 0.75,]$Hl, pch = 19, cex = .6, 
       col = rgb(red=0, green=0, blue=0, alpha=0.3))

dev.off()

png("D:/Dropbox/LightNiches/output/light in the geometric space linear predictors lr.png",
    width = 800, height = 800)
par(mfrow=c(3,3), cex = 1.1, mai=c(0.9,1,0.4,0.4))

vis.gam(lr25, view=c("Hl", "Rl"), type="response", plot.type="contour", main = "25cm-side\n",
        color = "gray", ylab = expression(paste("log(R "[25], ")")), xlab = expression(paste("log(H"[25], ")")))
points(df_int_tot[df_int_tot$dim == 0.25,]$Rl~df_int_tot[df_int_tot$dim == 0.25,]$Hl, pch = 19, cex = .6, 
       col = rgb(red=0, green=0, blue=0, alpha=0.3))
vis.gam(lr50, view=c("Hl", "Rl"), type="response", plot.type="contour", main = "50cm-side\n",
        color = "gray", ylab = expression(paste("log(R "[50], ")")), xlab = expression(paste("log(H"[50], ")")))
points(df_int_tot[df_int_tot$dim == 0.5,]$Rl~df_int_tot[df_int_tot$dim == 0.5,]$Hl, pch = 19, cex = .6, 
       col = rgb(red=0, green=0, blue=0, alpha=0.3))
vis.gam(lr75, view=c("Hl", "Rl"), type="response", plot.type="contour", main = "75cm-side\n",
        color = "gray", ylab = expression(paste("log(R "[75], ")")), xlab = expression(paste("log(H"[75], ")")))
points(df_int_tot[df_int_tot$dim == 0.75,]$Rl~df_int_tot[df_int_tot$dim == 0.75,]$Hl, pch = 19, cex = .6, 
       col = rgb(red=0, green=0, blue=0, alpha=0.3))

vis.gam(lr25, view=c("Dt", "Rl"), type="response", plot.type="contour", main = "",
        color = "gray", ylab = expression(paste("D "[25], "")), xlab = expression(paste("log(R"[25], ")")))
points(df_int_tot[df_int_tot$dim == 0.25,]$Rl~df_int_tot[df_int_tot$dim == 0.25,]$Dt, pch = 19, cex = .6, 
       col = rgb(red=0, green=0, blue=0, alpha=0.3))
vis.gam(lr50, view=c("Dt", "Rl"), type="response", plot.type="contour", main = "",
        color = "gray", ylab = expression(paste("D "[50], "")), xlab = expression(paste("log(R"[50], ")")))
points(df_int_tot[df_int_tot$dim == 0.5,]$Rl~df_int_tot[df_int_tot$dim == 0.5,]$Dt, pch = 19, cex = .6, 
       col = rgb(red=0, green=0, blue=0, alpha=0.3))
vis.gam(lr75, view=c("Dt", "Rl"), type="response", plot.type="contour", main = "",
        color = "gray", ylab = expression(paste("D "[75], "")), xlab = expression(paste("log(R"[75], ")")))
points(df_int_tot[df_int_tot$dim == 0.75,]$Rl~df_int_tot[df_int_tot$dim == 0.75,]$Dt, pch = 19, cex = .6, 
       col = rgb(red=0, green=0, blue=0, alpha=0.3))

vis.gam(lr25, view=c("Hl", "Dt"), type="response", plot.type="contour",main = "",
        color = "gray", xlab = expression(paste("D"[25], "")), ylab = expression(paste("log(H "[25], ")")))
points(df_int_tot[df_int_tot$dim == 0.25,]$Dt~df_int_tot[df_int_tot$dim == 0.25,]$Hl, pch = 19, cex = .6, 
       col = rgb(red=0, green=0, blue=0, alpha=0.3))
vis.gam(lr50, view=c("Hl", "Dt"), type="response", plot.type="contour", main = "",
        color = "gray", xlab = expression(paste("D"[50], "")), ylab = expression(paste("log(H "[50], ")")))
points(df_int_tot[df_int_tot$dim == 0.5,]$Dt~df_int_tot[df_int_tot$dim == 0.5,]$Hl, pch = 19, cex = .6, 
       col = rgb(red=0, green=0, blue=0, alpha=0.3))
vis.gam(lr75, view=c("Hl", "Dt"), type="response", plot.type="contour", main = "",
        color = "gray", xlab = expression(paste("D"[75], "")), ylab = expression(paste("log(H "[75], ")")))
points(df_int_tot[df_int_tot$dim == 0.75,]$Dt~df_int_tot[df_int_tot$dim == 0.75,]$Hl, pch = 19, cex = .6, 
       col = rgb(red=0, green=0, blue=0, alpha=0.3))

dev.off()


png("D:/Dropbox/LightNiches/output/light in the geometric space sc 05.png",
    width = 800, height = 300)
par(mfrow=c(1,3), cex = 1.1, mai=c(1,1,0,0.2))


vis.gam(sc50, view=c("Hl", "Rl"), type="response", plot.type="contour",
        color = "gray", ylab = expression(paste("log(R "[50], ")")), xlab = expression(paste("log(H"[50], ")")))
points(df_int_tot[df_int_tot$dim == 0.5,]$Rl~df_int_tot[df_int_tot$dim == 0.5,]$Hl, pch = 19, cex = .6, 
       col = rgb(red=0, green=0, blue=0, alpha=0.3))

vis.gam(sc50, view=c("Dt", "Rl"), type="response", plot.type="contour",
        color = "gray", ylab = expression(paste("D "[50], "")), xlab = expression(paste("log(R"[50], ")")))
points(df_int_tot[df_int_tot$dim == 0.5,]$Rl~df_int_tot[df_int_tot$dim == 0.5,]$Dt, pch = 19, cex = .6, 
       col = rgb(red=0, green=0, blue=0, alpha=0.3))

vis.gam(sc50, view=c("Hl", "Dt"), type="response", plot.type="contour",
        color = "gray", xlab = expression(paste("D"[50], "")), ylab = expression(paste("log(H "[50], ")")))
points(df_int_tot[df_int_tot$dim == 0.5,]$Dt~df_int_tot[df_int_tot$dim == 0.5,]$Hl, pch = 19, cex = .6, 
       col = rgb(red=0, green=0, blue=0, alpha=0.3))


dev.off()



## gam viz shift by intercept

par(mfrow=c(3,3), cex = 1.1, mai=c(0.9,1,0.4,0.4))

#png(paste0("D:/Dropbox/LightNiches/output/gamViz dim = 25.png"))
plot.gam(sc50, residuals = TRUE, pages = 1, pch = 19, cex = .5, shade = TRUE, seWithMean = TRUE, shift = coef(sc50)[1])

#dev.off()
#### predict vs observed ######
df25 <-df_int_tot[df_int_tot$dim == 0.25,]
df50 <-df_int_tot[df_int_tot$dim == 0.50,]
df75 <-df_int_tot[df_int_tot$dim == 0.75,]

### slc ####
df25$pa125slc <- predict.gam(slc125,df25)
df50$pa150slc <- predict.gam(slc150,df50)
df75$pa175slc <- predict.gam(slc175,df75)
df25$pa225slc <- predict.gam(slc225,df25)
df50$pa250slc <- predict.gam(slc250,df50)
df75$pa275slc <- predict.gam(slc275,df75)
df25$pa325slc <- predict.gam(slc325,df25)
df50$pa350slc <- predict.gam(slc350,df50)
df75$pa375slc <- predict.gam(slc375,df75)
pvoslc <- as.data.frame(
  rbind(
    cbind(obs = df25$light_noscale,
          pred1 = df25$pa125slc, 
          pred2 = df25$pa225slc,
          pred3 = df25$pa325slc, 
          dim = "25cm-side predictors",
          rec = as.character(df25$rec)),
    cbind(obs = df50$light_noscale, 
          pred1 = df50$pa150slc, 
          pred2 = df50$pa250slc,
          pred3 = df50$pa350slc, 
          dim = "50cm-side predictors",
          rec = as.character(df50$rec)),
    cbind(obs = df75$light_noscale, 
          pred1 = df75$pa175slc, 
          pred2 = df75$pa275slc,
          pred3 = df75$pa375slc, 
          dim = "75cm-side predictors",
          rec = as.character(df75$rec))))


pvoslc$dim <- as.factor(pvoslc$dim)
pvoslc$obs <- as.numeric(pvoslc$obs)
pvoslc$pred1 <- as.numeric(pvoslc$pred1)
pvoslc$pred2 <- as.numeric(pvoslc$pred2)
pvoslc$pred3 <- as.numeric(pvoslc$pred3)
pvoslc_p <- ggplot(pvoslc, aes(obs, pred1, col = rec))+
  geom_point(cex = 3, alpha = 0.4) +
  xlab("\nObserved light integral")+
  ylab("Predicted light integrals,\nlinear predictors\n") +
  scale_colour_viridis_d(name = "site")+
  geom_abline(linetype = 2, cex = 1, slope = 1, intercept = 0)+
  geom_smooth(method = "lm", se = FALSE, col = "black")+
  theme_minimal()+
  facet_wrap (~dim)+
  theme(text = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(vjust = 0.5),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 12))



df25$pf25 <- predict.gam(lc25,df25, exclude = "s(lon,lat)")
df50$pf50 <- predict.gam(lc50,df50, exclude = "s(lon,lat)")
df75$pf75 <- predict.gam(lc75,df75, exclude = "s(lon,lat)")

df25$pa25lc <- predict.gam(lc25,df25)
df50$pa50lc <- predict.gam(lc50,df50)
df75$pa75lc <- predict.gam(lc75,df75)

df25$pa25sc <- predict.gam(sc25,df25)
df50$pa50sc <- predict.gam(sc50,df50)
df75$pa75sc <- predict.gam(sc75,df75)

df25$pa125slc <- predict.gam(slc125,df25)
df50$pa150slc <- predict.gam(slc150,df50)
df75$pa175slc <- predict.gam(slc175,df75)
df25$pa225slc <- predict.gam(slc225,df25)
df50$pa250slc <- predict.gam(slc250,df50)
df75$pa275slc <- predict.gam(slc275,df75)
df25$pa325slc <- predict.gam(slc325,df25)
df50$pa350slc <- predict.gam(slc350,df50)
df75$pa375slc <- predict.gam(slc375,df75)
pvo <- as.data.frame(rbind(cbind(obs = df25$light_noscale,predlc = df25$pa25lc, predsr =  df25$pa25sr,
                                 predsc = df25$pa25sc, predlr =  df25$pa25lr, predlm =  df25$pa25lm,
                                 dim = "25cm-side predictors",rec = as.character(df25$rec)),
                           cbind(obs = df50$light_noscale, predlc = df50$pa50lc, predsr =  df50$pa50sr,
                                 predsc = df50$pa50sc, predlr =  df50$pa50lr, predlm =  df50$pa50lm,
                                 dim = "50cm-side predictors",rec = as.character(df50$rec)),
                           cbind(obs = df75$light_noscale,predlc = df75$pa75lc, predsr =  df75$pa75sr,
                                 predsc = df75$pa75sc, predlr =  df75$pa75lr, predlm =  df75$pa75lm,
                                 dim = "75cm-side predictors",rec = as.character(df75$rec))))

df25$pa25lr <- predict.gam(lr25,df25)
df50$pa50lr <- predict.gam(lr50,df50)
df75$pa75lr <- predict.gam(lr75,df75)

df25$pa25sr <- predict.gam(sr25,df25)
df50$pa50sr <- predict.gam(sr50,df50)
df75$pa75sr <- predict.gam(sr75,df75)

df25$pa25lm <- predict(lm25,df25)
df50$pa50lm <- predict(lm50,df50)
df75$pa75lm <- predict(lm75,df75)


pvo <- as.data.frame(rbind(cbind(obs = df25$light_noscale,predlc = df25$pa25lc, predsr =  df25$pa25sr,
                                 predsc = df25$pa25sc, predlr =  df25$pa25lr, predlm =  df25$pa25lm,
                                 dim = "25cm-side predictors",rec = as.character(df25$rec)),
                           cbind(obs = df50$light_noscale, predlc = df50$pa50lc, predsr =  df50$pa50sr,
                                 predsc = df50$pa50sc, predlr =  df50$pa50lr, predlm =  df50$pa50lm,
                                 dim = "50cm-side predictors",rec = as.character(df50$rec)),
                           cbind(obs = df75$light_noscale,predlc = df75$pa75lc, predsr =  df75$pa75sr,
                                 predsc = df75$pa75sc, predlr =  df75$pa75lr, predlm =  df75$pa75lm,
                                 dim = "75cm-side predictors",rec = as.character(df75$rec))))
pvo$dim <- as.factor(pvo$dim)
pvo$obs <- as.numeric(pvo$obs)
pvo$predlc <- as.numeric(pvo$predlc)
pvo$predsc <- as.numeric(pvo$predsc)
pvo$predlr <- as.numeric(pvo$predlr)
pvo$predsr <- as.numeric(pvo$predsr)
pvo$predlm <- as.numeric(pvo$predlm)
pvo$rec <- as.factor(pvo$rec)
str(pvo)
pvoglc <- ggplot(pvo, aes(obs, predlc, col = rec))+
  geom_point(cex = 3, alpha = 0.4) +
  xlab("\nObserved light integral")+
  ylab("Predicted light integrals,\nlinear predictors\n") +
  scale_colour_viridis_d(name = "site")+
  geom_abline(linetype = 2, cex = 1, slope = 1, intercept = 0)+
  geom_smooth(method = "lm", se = FALSE, col = "black")+
  theme_minimal()+
  facet_wrap (~dim)+
  theme(text = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(vjust = 0.5),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 12))
pvogsc <- ggplot(pvo, aes(obs, predsc, col = rec))+
  geom_point(cex = 3, alpha = 0.4) +
  xlab("\nObserved light integral")+
  ylab("Predicted light integrals,\nsmoothed predictors") +
  scale_colour_viridis_d(name = "site")+
  geom_abline(linetype = 2, cex = 1, slope = 1, intercept = 0)+
  geom_smooth(method = "lm", se = FALSE, col = "black")+
  theme_minimal()+
  facet_wrap (~dim)+
  theme(text = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(vjust = 0.5),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 12))
png("D:/Dropbox/LightNiches/output/pred vs obs full.png",
    width = 700, height = 700)
ggarrange(pvogsc,pvoglc, ncol = 1 )
dev.off()


srsc <- ggplot(pvo, aes(predsr, predsc, col = rec))+
  geom_point(cex = 2, alpha = .4) +
  ggtitle("a. Models with smooth basis for geometric predictors")+
  xlab("\nPredicted light integrals,site as random effect")+
  ylab("Predicted light integrals,\ncoordinates with \nGaussian process basis") +
  scale_colour_viridis_d(name = "site")+
  geom_abline(linetype = 2, cex = 1, slope = 1, intercept = 0)+
  theme_minimal()+
  facet_wrap (~dim)+
  theme(text = element_text(size = 10),
        plot.title = element_text(size=14, face="bold.italic"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(vjust = 0.5),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 12))

lrlc <- ggplot(pvo, aes(predlr, predlc, col = rec))+
  geom_point(cex = 2, alpha = .4) +
  ggtitle("b. Models with linear geometric predictors")+
  xlab("\nPredicted light integrals, site as random effect")+
  ylab("Predicted light integrals,\ncoordinates with \nGaussian process basis") +
  scale_colour_viridis_d(name = "site")+
  geom_abline(linetype = 2, cex = 1, slope = 1, intercept = 0)+
  theme_minimal()+
  facet_wrap (~dim)+
  theme(text = element_text(size = 10),
        plot.title = element_text(size=14, face="bold.italic"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(vjust = 0.5),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 12))


pvplm <- ggplot(pvo, aes(predlc, predlm, col = rec))+
  geom_point(cex = 1) +
  xlab("\nObserved light integral")+
  ylab("Predicted light integral\n") +
  scale_colour_viridis_d(name = "site")+
  geom_abline(linetype = 2, cex = 1, slope = 1, intercept = 0)+
  geom_smooth(method = "lm", se = FALSE, col = "black")+
  theme_minimal()+
  facet_wrap (~dim)+
  theme(text = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(vjust = 0.5),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 12))

lrsr <- ggplot(pvo, aes(predlr, predsr, col = rec))+
  geom_point(cex = 2, alpha = .4) +
  ggtitle("c. Models with site as random effect")+
  xlab("\nPredicted light integrals, linear predictors")+
  ylab("Predicted light integrals,\nsmooth basis for \ngeometric predictors") +
  scale_colour_viridis_d(name = "site")+
  geom_abline(linetype = 2, cex = 1, slope = 1, intercept = 0)+
  #geom_smooth(method = "lm", se = FALSE, col = "black")+
  theme_minimal()+
  facet_wrap (~dim)+
  theme(text = element_text(size = 10),
        plot.title = element_text(size=14, face="bold.italic"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(vjust = 0.5),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 12))

lcsc <- ggplot(pvo, aes(predlc, predsc, col = rec))+
  geom_point(cex = 2, alpha = .4) +
  ggtitle("d. Models with coordinates with Gaussian process basis")+
  xlab("\nPredicted light integrals, linear predictors")+
  ylab("Predicted light integrals,\nsmooth basis for \ngeometric predictors") +
  scale_colour_viridis_d(name = "site")+
  geom_abline(linetype = 2, cex = 1, slope = 1, intercept = 0)+
  #geom_smooth(method = "lm", se = FALSE, col = "black")+
  theme_minimal()+
  facet_wrap (~dim)+
  theme(text = element_text(size = 10),
        plot.title = element_text(size=14, face="bold.italic"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(vjust = 0.5),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 12))


plmVpsc <- ggplot(pvo, aes(predlm, predsc, col = rec))+
  geom_point(cex = 3, alpha = .4) +
  xlab("\nObserved light integral")+
  ylab("Predicted light integral\n") +
  scale_colour_viridis_d(name = "site")+
  geom_abline(linetype = 2, cex = 1, slope = 1, intercept = 0)+
  geom_smooth(method = "lm", se = FALSE, col = "black")+
  theme_minimal()+
  facet_wrap (~dim)+
  theme(text = element_text(size = 9),
        axis.title.x = element_text(size = 9),
        axis.title.y = element_text(size = 9),
        axis.text = element_text(size = 9),
        axis.text.x = element_text(vjust = 0.5),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 9),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 9))

png("D:/Dropbox/LightNiches/output/preds vs preds.png",
    width = 700, height = 900)
ggarrange(srsc, lrlc, lrsr, lcsc, ncol = 1 )
dev.off()

###plot for chapter

sc <- gam(light_noscale~s(Dt,Rl)+ s(tdp)+ s(watercol12)+
            s(par_integral)+
            s(b) + s(aa) +s(lon, lat, bs = "gp", m=2), data = df_int, method="REML")
lab <-paste0("sc",dim)
gam_list[[lab]] <- sc
png(paste0("D:/Dropbox/LightNiches/output/gamViz ",lab, ".png"))
plot(sc, residuals = TRUE, pages = 1, pch = 19, cex = .5, shade = TRUE, seWithMean = TRUE, shift = coef(sc)[1])
dev.off()

pred <- predict(sc, df_int)
preddf <-cbind(pred = pred, obs = df_int$light_noscale, mod = lab, rec = as.character(df_int$rec))


str(preddf)
preddf <- as.data.frame(preddf)
preddf$pred <- as.numeric(preddf$pred)
preddf$obs <- as.numeric(preddf$obs)
preddf$rec = as.factor(preddf$rec)

sc05po <- ggplot(preddf, aes(obs, pred, col = rec))+
  geom_point(cex = 3, alpha = 0.4) +
  xlab("\nObserved light integral")+
  ylab("Predicted light integrals\n") +
  scale_colour_viridis_d(name = "site")+
  geom_abline(linetype = 2, cex = 1, slope = 1, intercept = 0)+
  geom_smooth(method = "lm", se = FALSE, col = "black")+
  theme_minimal()+
  theme(text = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(vjust = 0.5),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 12))
png("D:/Dropbox/LightNiches/output/preds vs obs 0.5.png",
    width = 600, height = 400)
plot(sc05po)
dev.off()


### D theory vs D variation method #####

Dcomp25 <- ggplot(df25, aes(Du, Dt))+
  geom_point(alpha = .4,cex = 2) +
  xlab(expression(paste("D"[25], " (variation method)")))+
  ylab(expression(paste("D"[25], " (theory)"))) +
  geom_smooth(method = "lm", se = FALSE, col = "black")+
  stat_cor(method="pearson",  digits = 2) +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text =  element_text(size = 12),
        text = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12))

Dcomp50 <- ggplot(df50, aes(Du, Dt))+
  geom_point(alpha = .4,cex = 2) +
  xlab(expression(paste("D"[50], " (variation method)")))+
  ylab(expression(paste("D"[50], " (theory)"))) +
  geom_smooth(method = "lm", se = FALSE, col = "black")+
  stat_cor(method="pearson") +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text =  element_text(size = 12),
        text = element_text(size = 14),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12))


Dcomp75 <- ggplot(df75, aes(Du, Dt))+
  geom_point(alpha = .4,cex = 2) +
  xlab(expression(paste("D"[75], " (variation method)")))+
  ylab(expression(paste("D"[75], " (theory)"))) +
  geom_smooth(method = "lm", se = FALSE, col = "black")+
  stat_cor(method="pearson") +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text =  element_text(size = 12),
        text = element_text(size = 14),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12))


png("D:/Dropbox/My Dropbox/NC-RR_environment_data/outputs_fit2/D vs Dtheory.png",
    width = 250, height = 500, units = "px")
print(ggarrange(Dcomp25,Dcomp50,Dcomp75, ncol = 1))
dev.off()


#### pair plot #####

data50 <- df50[,c("light_noscale", "par_integral", "Dt", "Rl", "Hl",
                  "tdp", "a", "b", "watercol12","rec")]
colnames(data50) <- c("reef\nlight", "surface\nlight", "D", "log(R)", "log(H)",
                      "mean\n water\n column", "a", "b", "water \ncolumn\n at noon","site")
ggcorr(data50, palette = "RdYlGn", name = "rho", 
       label = FALSE, label_color = "black")
columns <- 1:ncol(data50)
pairs50 <- ggpairs(data50, title = "",
                   axisLabels = "show", 
                   columnLabels = colnames(data50[, columns])) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

png("D:/Dropbox/LightNiches/output/pairs 50.png", res = 300,
    width = 30, height = 20, units = "cm")
print(pairs50)
dev.off()


data25 <- df25[,c("light_noscale", "par_integral", "Dt", "Rl", "Hl",
                  "tdp", "a", "b", "watercol12","rec")]
colnames(data25) <- c("reef\nlight", "surface\nlight", "D", "log(R)", "log(H)",
                      "mean\n water\n column", "a", "b", "water \ncolumn\n at noon","site")
columns <- 1:ncol(data25)
pairs25 <- ggpairs(data25, title = "",  
                   axisLabels = "show", 
                   columnLabels = colnames(data25[, columns]))  +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

png("D:/Dropbox/LightNiches/output/pairs 25.png", res = 300,
    width = 30, height = 20, units = "cm")
print(pairs25)
dev.off()



data75 <- df75[,c("light_noscale", "par_integral", "Dt", "Rl", "Hl",
                  "tdp", "a", "b", "watercol12","rec")]
colnames(data75) <- c("reef\nlight", "surface\nlight", "D", "log(R)", "log(H)",
                      "mean\n water\n column", "a", "b", "water \ncolumn\n at noon","site")
ggcorr(data75, palette = "RdYlGn", name = "rho", 
       label = FALSE, label_color = "black")
columns <- 1:ncol(data75)
pairs75 <- ggpairs(data75, title = "",  
                   axisLabels = "show", 
                   columnLabels = colnames(data75[, columns])) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

png("D:/Dropbox/LightNiches/output/pairs 75.png", res = 300,
    width = 30, height = 20, units = "cm")
print(pairs75)
dev.off()

#### geometry among scales  #####

scaleD <- data.frame(key = c(df25$mergekey, df50$mergekey,df75$mergekey),
                     rec = as.factor(c(as.character(df25$rec), as.character(df50$rec),as.character(df75$rec))),
                     val = c(df25$Dt, df50$Dt,df75$Dt),
                     var = as.factor(rep("Dt",2708)),
                     dim = as.factor(c(rep("dim25",903),rep("dim50",903),rep("dim75",902))))
scaleR <- data.frame(key = c(df25$mergekey, df50$mergekey,df75$mergekey),
                     rec = as.factor(c(as.character(df25$rec), as.character(df50$rec),as.character(df75$rec))),
                     val = c(df25$Rl, df50$Rl,df75$Rl),
                     var = as.factor(rep("Rl",2708)),
                     dim = as.factor(c(rep("dim25",903),rep("dim50",903),rep("dim75",902))))
scaleH <- data.frame(key = c(df25$mergekey, df50$mergekey,df75$mergekey),
                     rec = as.factor(c(as.character(df25$rec), as.character(df50$rec),as.character(df75$rec))),
                     val = c(df25$Hl, df50$Hl,df75$Hl),
                     var = as.factor(rep("Hl",2708)),
                     dim = as.factor(c(rep("dim25",903),rep("dim50",903),rep("dim75",902))))
scale <- rbind(scaleD, scaleR, scaleH)
str(scale)

ggplot(scale, aes(key, val, col = key))+
  geom_point()+
  facet_wrap(~var)+
  theme_minimal()+
  theme(legend.text = element_blank())


