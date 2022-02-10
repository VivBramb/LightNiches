##### complexity and attenuation curves #####


#### libraries #####

library(ggplot2)
library(tidyverse)
library(viridis)
library(mgcv)
library(lme4)

#### functions ####
ld <- function(par, depth) {
  par_log <- log(par)
  k <- 4/50
  par_log <- par_log - depth * k
  return(exp(par_log))
}

ld1 <- function(par, depth) {
  par_log <- log(par)
  k <- 40/50
  par_log <- par_log - depth * k
  return(exp(par_log))
}

ld2 <- function(par, depth) {
  par_log <- log(par)
  k <- 4/500
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

# initialize storing lists and dfs

df_int_tot <- data.frame() 
df_tot <- data.frame()
gam_list<-list()

geom_all <- read.csv(paste0("output/geom_new.csv"), h = T)
dims <- unique(geom_all$L)

for (dim in dims) {

  geom <- geom_all[geom_all$L == dim,]
  geom$mergekey <- paste0(geom$rec, geom$year, geom$new_unit) #unique ID for hobo deployment
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
  
  
  df$pard1[df$surface_par !=0] <- ld1(df$surface_par[df$surface_par !=0], df$water_col[df$surface_par !=0])
  df$pard1 [df$surface_par ==0] <- 0
  #summary(df$pard)
  df$light_pard1 <- 0
  df$light_pard1[df$pard1 !=0] <- df$light_mol_5p[df$pard1 !=0]/df$pard1[df$pard1 !=0]
  summary(df$light_pard1)
  
  
  df$pard2[df$surface_par !=0] <- ld2(df$surface_par[df$surface_par !=0], df$water_col[df$surface_par !=0])
  df$pard2 [df$surface_par ==0] <- 0
  #summary(df$pard)
  df$light_pard2 <- 0
  df$light_pard2[df$pard2 !=0] <- df$light_mol_5p[df$pard2 !=0]/df$pard2[df$pard2 !=0]
  summary(df$light_pard2)
  
  

  df_long <- as.data.frame(df %>%
                        group_by(mergekey,lat,lon) %>%
                          filter(rep == "a") %>%
                        mutate(integral_props = sum(light_pard),
                               integral_pard = sum(pard),
                               integral_props1 = sum(light_pard1),
                               integral_pard1 = sum(pard1),
                               integral_props2 = sum(light_pard2),
                               integral_pard2 = sum(pard2),
                               bottom_integral = sum(light_mol_5p)))
  df_long$prop_integrals = df_long$bottom_integral/df_long$integral_pard
  df_long$prop_integrals1 = df_long$bottom_integral/df_long$integral_pard1
  df_long$prop_integrals2 = df_long$bottom_integral/df_long$integral_pard2
  
  df_tot <- rbind(df_long,df_tot)
# get n
  df_int <- unique(df_long[,c("integral_props","prop_integrals",
                              "integral_props1","prop_integrals1",
                              "integral_props2","prop_integrals2",
                              "H","site","year.x",
                         "integral_pard","integral_mol_5p", "folder",
                         "Dclip",  "Dplane", "Rclip","dpHB","Rclip", 
                         "Rplane","lon","lat",
                         "hoboID", "mergekey","L")])
  
  # plot(data = df_int, integral_props~prop_integrals)
  # length(unique(df$mergekey))
  # length(unique(df_int$mergekey))
  # dim(df_int)
  # sum(duplicated(df_int$mergekey)) # there are 2 duplicates: two same IDs in the same record
  # df_int[duplicated(df_int$mergekey),]
  # df_int[841:915,]
  # df_int[524:527,] #RS17H44
  # df_int[409:412,] # N317H34
  # 
  # leave one out for now and then go to the shp and the paper location annotation
  df_int <- df_int[!duplicated(df_int$mergekey),]
  df_int <- df_int[df_int$folder != "RS19a",]
  dim(df_int) 
  #removeNAs
  df_int <- df_int[!is.na(df_int$Rclip)& !is.na(df_int$Rplane),]
  
  # scale values and set factors
  # tide is in cm, so divide it to get m
  
  df_int$year <- as.factor(df_int$year)
  df_int$Rcl <- as.numeric(log10(df_int$Rclip))
  df_int$Rpl <- as.numeric(log10(df_int$Rplane))
  df_int$Hl <- as.numeric(log10(df_int$H))
  df_int$dim <- dim
  df_int$rec<- as.factor(df_int$site)
  df_int_tot <- rbind(df_int_tot, df_int)
  
  ss <- gamm(prop_integrals~s(Dclip,Rcl)+ s(lon, lat, bs = "gp", m=2), family = "binomial",
              data = df_int, method="REML")
  ls <- gamm(prop_integrals~Dclip*Rcl+ s(lon, lat, bs = "gp", m=2), family = "binomial",
              data = df_int, method="REML")
  sl <- gamm(prop_integrals~s(Dclip,Rcl), random = list(rec=~1), family = "binomial",
               data = df_int, method="REML")
  ll <- gamm(prop_integrals~Dclip*Rcl, random = list(rec=~1), family = "binomial",
              data = df_int, method="REML")
  
  lab <-paste0("ss",dim)
  gam_list[[lab]] <- ss
  lab <-paste0("ls",dim)
  gam_list[[lab]] <- ls
  lab <-paste0("sl",dim)
  gam_list[[lab]] <- sl
  lab <-paste0("ll",dim)
  gam_list[[lab]] <- ll

}


par(mfrow=c(1,3))
vis.gam(gam_list["ss0.25"]$ss0.25$gam,view=c("Dclip", "Rcl"), type="response", plot.type="contour", main = "GAM25 / s(D,R) + s(lon,lat)",
  color = "gray", ylab = expression(paste("log(R "[dim], ")")), xlab = expression(paste("D"[50])))
points(df_int_tot[df_int_tot$L == 0.25,]$Rcl~df_int_tot[df_int_tot$L == 0.25,]$Dclip, pch = 19, cex = .6,
  col = rgb(red=0, green=0, blue=0, alpha=0.7))

vis.gam(gam_list["ss0.5"]$ss0.5$gam,view=c("Dclip", "Rcl"), type="response", plot.type="contour", main = "GAM50 / s(D,R) + s(lon,lat)",
        color = "gray", ylab = expression(paste("log(R "[dim], ")")), xlab = expression(paste("D"[50])))
points(df_int_tot[df_int_tot$L == 0.5,]$Rcl~df_int_tot[df_int_tot$L == 0.5,]$Dclip, pch = 19, cex = .6,
       col = rgb(red=0, green=0, blue=0, alpha=0.3))

vis.gam(gam_list["ss0.75"]$ss0.75$gam,view=c("Dclip", "Rcl"), type="response", plot.type="contour", main = "GAM75 / s(D,R) + s(lon,lat)",
        color = "gray", ylab = expression(paste("log(R "[dim], ")")), xlab = expression(paste("D"[50])))
points(df_int_tot[df_int_tot$L == 0.75,]$Rcl~df_int_tot[df_int_tot$L == 0.75,]$Dclip, pch = 19, cex = .6,
       col = rgb(red=0, green=0, blue=0, alpha=0.7))


vis.gam(gam_list["ls0.25"]$ls0.25$gam,view=c("Dclip", "Rcl"), type="response", plot.type="contour", main = "LM/GAM25 / D*R  + s(lon,lat)",
        color = "gray", ylab = expression(paste("log(R "[dim], ")")), xlab = expression(paste("D"[50])))
points(df_int_tot[df_int_tot$L == 0.25,]$Rcl~df_int_tot[df_int_tot$L == 0.25,]$Dclip, pch = 19, cex = .6,
       col = rgb(red=0, green=0, blue=0, alpha=0.3))

vis.gam(gam_list["ls0.5"]$ls0.5$gam,view=c("Dclip", "Rcl"), type="response", plot.type="contour", main = "LM/GAM50 /  D*R  + s(lon,lat)",
        color = "gray", ylab = expression(paste("log(R "[dim], ")")), xlab = expression(paste("D"[50])))
points(df_int_tot[df_int_tot$L == 0.5,]$Rcl~df_int_tot[df_int_tot$L == 0.5,]$Dclip, pch = 19, cex = .6,
       col = rgb(red=0, green=0, blue=0, alpha=0.3))

vis.gam(gam_list["ls0.75"]$ls0.75$gam,view=c("Dclip", "Rcl"), type="response", plot.type="contour", main = "LM/GAM75 /  D*R + s(lon,lat)",
        color = "gray", ylab = expression(paste("log(R "[dim], ")")), xlab = expression(paste("D"[50])))
points(df_int_tot[df_int_tot$L == 0.75,]$Rcl~df_int_tot[df_int_tot$L == 0.75,]$Dclip, pch = 19, cex = .6,
       col = rgb(red=0, green=0, blue=0, alpha=0.3))



vis.gam(gam_list["sl0.25"]$sl0.25$gam,view=c("Dclip", "Rcl"), type="response", plot.type="contour", main = "GAM/LM25 /s(D,R) + 1|site",
        color = "gray", ylab = expression(paste("log(R "[dim], ")")), xlab = expression(paste("D"[50])))
points(df_int_tot[df_int_tot$L == 0.25,]$Rcl~df_int_tot[df_int_tot$L == 0.25,]$Dclip, pch = 19, cex = .6,
       col = rgb(red=0, green=0, blue=0, alpha=0.3))

vis.gam(gam_list["sl0.5"]$sl0.5$gam,view=c("Dclip", "Rcl"), type="response", plot.type="contour", main = "GAM/LM50 /  s(D,R) + 1|site",
        color = "gray", ylab = expression(paste("log(R "[dim], ")")), xlab = expression(paste("D"[50])))
points(df_int_tot[df_int_tot$L == 0.5,]$Rcl~df_int_tot[df_int_tot$L == 0.5,]$Dclip, pch = 19, cex = .6,
       col = rgb(red=0, green=0, blue=0, alpha=0.3))

vis.gam(gam_list["sl0.75"]$sl0.75$gam,view=c("Dclip", "Rcl"), type="response", plot.type="contour", main = "LM/GAM75 /  s(D,R)  + 1|site",
        color = "gray", ylab = expression(paste("log(R "[dim], ")")), xlab = expression(paste("D"[50])))
points(df_int_tot[df_int_tot$L == 0.75,]$Rcl~df_int_tot[df_int_tot$L == 0.75,]$Dclip, pch = 19, cex = .6,
       col = rgb(red=0, green=0, blue=0, alpha=0.3))



vis.gam(gam_list["ll0.25"]$ll0.25$gam,view=c("Dclip", "Rcl"), type="response", plot.type="contour", main = "LM25 / D*R + 1|site",
        color = "gray", ylab = expression(paste("log(R "[dim], ")")), xlab = expression(paste("D"[50])))
points(df_int_tot[df_int_tot$L == 0.25,]$Rcl~df_int_tot[df_int_tot$L == 0.25,]$Dclip, pch = 19, cex = .6,
       col = rgb(red=0, green=0, blue=0, alpha=0.3))

vis.gam(gam_list["ll0.5"]$ll0.5$gam,view=c("Dclip", "Rcl"), type="response", plot.type="contour", main =  "LM50 / D*R + 1|site",
        color = "gray", ylab = expression(paste("log(R "[dim], ")")), xlab = expression(paste("D"[50])))
points(df_int_tot[df_int_tot$L == 0.5,]$Rcl~df_int_tot[df_int_tot$L == 0.5,]$Dclip, pch = 19, cex = .6,
       col = rgb(red=0, green=0, blue=0, alpha=0.3))

vis.gam(gam_list["ll0.75"]$ll0.75$gam,view=c("Dclip", "Rcl"), type="response", plot.type="contour", main =  "LM75 / D*R + 1|site",
        color = "gray", ylab = expression(paste("log(R "[dim], ")")), xlab = expression(paste("D"[50])))
points(df_int_tot[df_int_tot$L == 0.75,]$Rcl~df_int_tot[df_int_tot$L == 0.75,]$Dclip, pch = 19, cex = .6,
       col = rgb(red=0, green=0, blue=0, alpha=0.3))


1###########
source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")

ylab.light = expression(paste("Proportion of the integrals"))
df_intALL <- df_int
df_intALL$folder = "zall"
df_intALL$rec = "zall"
df_int2 <- rbind(cbind(df_intALL, rec2 = df_int$rec),cbind(df_int, rec2 = df_int$rec))

gl <- ggplot(data = df_int2, 
             aes(x = folder, y = prop_integrals, fill = rec)) +
  geom_point(aes(y = prop_integrals, color = rec2),
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

dim(df)

######

df_tot %>% 
  filter(folder == "GT17a" ) %>%
  ggplot() +
  geom_line(aes(x = datetime_pos, y = light_mol), col = "grey80")+
  geom_line(aes(x = datetime_pos, y = surface_par), col = "red")+
  geom_line(aes(x = datetime_pos, y = light_mol_5p), col = "black")+
  #geom_smooth(method = "loess")+
  theme_minimal()+
  theme(legend.title = element_blank(), legend.position = "none")+
  facet_wrap(.~hoboID,nrow = 4, scales = "free_x")

df_tot %>% 
  filter(folder == "GT19a" ) %>%
  ggplot() +
  geom_line(aes(x = datetime_pos, y = light_mol), col = "grey80")+
  geom_line(aes(x = datetime_pos, y = surface_par), col = "red")+
  geom_line(aes(x = datetime_pos, y = light_mol_5p), col = "black")+
  #geom_smooth(method = "loess")+
  theme_minimal()+
  theme(legend.title = element_blank(), legend.position = "none")+
  facet_wrap(.~hoboID,nrow = 4, scales = "free_x")

example <- as.data.frame(df_tot %>%
                           filter((folder == "GT17a" & hoboID == "GT03")| (folder == "GT19a"& hoboID == "H02")))
example%>%
  ggplot() +
  geom_line(aes(x = datetime_pos, y = light_mol), col = "grey80")+
  geom_line(aes(x = datetime_pos, y = surface_par), col = "red")+
  geom_line(aes(x = datetime_pos, y = light_mol_5p), col = "black")+
  #geom_smooth(method = "loess")+
  theme_minimal()+
  theme(legend.title = element_blank(), legend.position = "none")+
  facet_wrap(.~hoboID,nrow = 4, scales = "free_x")

example%>%
  ggplot() +
  geom_line(aes(x = datetime_pos, y = light_mol), col = "grey80")+
  geom_line(aes(x = datetime_pos, y = surface_par), alpha = .2, col = "red")+
  geom_line(aes(x = datetime_pos, y = pard), col = "red")+
  geom_line(aes(x = datetime_pos, y = light_mol_5p), col = "black")+
  #geom_smooth(method = "loess")+
  theme_minimal()+
  theme(legend.title = element_blank(), legend.position = "none")+
  facet_wrap(.~hoboID,nrow = 4, scales = "free_x")

df_tot %>% 
  filter(folder == "GT17a" ) %>%
  ggplot() +
  ylim(0,2000)+
  geom_line(aes(x = datetime_pos, y = surface_par), col = "red")+
  geom_point(aes(x = datetime_pos, y = light_mol_5p), alpha = .3, col = "grey70")+
  geom_point(aes(x = datetime_pos, y = pard2), alpha = .3, col = "black")+
  #geom_smooth(method = "loess")+
  theme_minimal()+
  theme(legend.title = element_blank(), legend.position = "none")

#######

library(GGally)
for (dim in dims) {}
data <- df_int_tot[df_int_tot$L==dim,c("prop_integrals","integral_mol_5p", "H", "Rcl", "Dclip", "dpHB",
                  "year","site")]
data$H<- log(data$H)
colnames(data) <- c("prop of\nintegrals", "bottom\nintegrals", "log(Height\nrange)", "log(Rugosity)",
                      "Fractal\n dimension", "HOBO\ndepth","year","site")
ggcorr(data, palette = "RdYlGn", name = "rho", 
       label = FALSE, label_color = "black")
columns <- 1:ncol(data)
pairs <- ggpairs(data, title = "",
                   axisLabels = "show", 
                   columnLabels = colnames(data[, columns])) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
pairs


data <- df_int_tot[df_int_tot$L==dim,c("prop_integrals","integral_mol_5p", "H", "Rcl", "Dclip", "dpHB",
                                       "year","site")]
colnames(data) <- c("prop of\nintegrals", "bottom\nintegrals", "Height\nrange", "Rugosity",
                    "Fractal\n dimension", "HOBO\ndepth","year","site")
ggcorr(data, palette = "RdYlGn", name = "rho", 
       label = FALSE, label_color = "black")
columns <- 1:ncol(data)
pairs <- ggpairs(data, title = "",
                 axisLabels = "show", 
                 columnLabels = colnames(data[, columns])) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
pairs



png("D:/Dropbox/My Dropbox/NC-RR_environment_data/outputs_fit4/pairs 50.png", res = 300,
    width = 30, height = 20, units = "cm")
print(pairs50)
dev.off()




