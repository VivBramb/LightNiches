library(ggplot2)
library(tidyverse)
library(viridis)
library(mgcv)
library(lme4)
library(hms)
library(lubridate)
library(stringr)
library(MuMIn)
library(brms)
library(report)

#### functions ####
ld <- function(par, depth) {
  par_log <- log(par)
  k <- 4/50
  par_log <- par_log - depth * k
  return(exp(par_log))
}

#### read environmental data ####
env <-  read.csv("output/env_long.csv", h=T)
env<- env[env$folder != "L119b",]

env$mergekey <- paste0(env$site,env$year) # create unique ID for merging with geometry df
str(env)
dim(env)
env$datetime_pos <- as.POSIXct(env$datetime_pos)
unique(env$folder)
summary(env)

## all sampling have same n?
env_check <- env %>%
  group_by(folder) %>%
  mutate(n_distinct(datetime_pos))
unique(env_check$`n_distinct(datetime_pos)`)
# yes!

## all sam
env_check <- env %>%
  group_by(folder) %>%
  mutate(n_distinct(datetime_pos))%>%
  count(datetime_pos)
#View(env_check)

#View(env[env$datetime_pos == "2019-11-23 00:00:00",]) 
#fix it
env_fix <- env[env$datetime_pos == "2019-11-23 00:00:00" & env$tide_cm == 33.0398,]
env <- rbind(env[env$datetime_pos != "2019-11-23 00:00:00",], env_fix)
env_check <- env %>%
  group_by(folder) %>%
  mutate(n_distinct(datetime_pos))%>%
  count(datetime_pos)
#View(env_check)
#yes!


# ## read geometry
# geom2 <- read.csv(paste0("output/grid2_corrected.csv"), h = T)
# geom2$mergekey <- as.factor(str_sub(geom2$rec, start =1, end = 4))
# geom2<- as.data.frame(geom2 %>%
#                         group_by(mergekey)%>%
#                         mutate(Dmean = mean(D),
#                                Rmean = mean(R),
#                                Hmean = mean(H))%>%
#                         distinct(mergekey, Dmean, Rmean, Hmean))
# 
# 
# geom1 <- read.csv("output/geom_sites_8m.csv")
# geom1$mergekey <- paste0(geom1$rec,geom1$year)
# 
# geom <- merge(geom1,geom2)
# 
# df_tot <- merge(env,geom)
# names(df_tot)


### merge with depths
geom <- read.csv("output/store_2024.csv", h = T)
geom <- geom [,-c(5:8)]
geom <- as.data.frame(geom %>% pivot_wider(names_from = L, values_from = c(R,D,H)))

store2 <- read.csv("output/grid2_2024.csv", h = T)
geom_site  <- as.data.frame(store2 %>% 
  group_by(rec) %>%
  mutate(Rm_log = mean(log10(R)),
         Hm = mean(H)) %>%
  distinct(rec,Rm_log,Hm))

geom_site$rec <- str_sub(geom_site$rec, start =1, end = 2)

geom <- merge(geom,geom_site)

dim(geom)
# get 75cm squares
unique(env$unit)
unique(geom$new_unit)
env$mergekeyHB <- paste0(env$mergekey, env$unit)
geom$mergekeyHB <- paste0(geom$rec, geom$year, geom$new_unit)

#unique(env$mergekeyHB)
#unique(geom$mergekeyHB)

df <- merge(geom,env,by = "mergekeyHB")

summary(df)
names(df)
dim(df)


# correct light values as if there was always maximum registered par at the surface
nrow(df[df$surface_par<0,])
df$surface_par[df$surface_par<0] <- 0

df$water_col <- df$tide_cm/100 - df$dpHB
df$pard[df$surface_par !=0] <- ld(df$surface_par[df$surface_par !=0], df$water_col[df$surface_par !=0])
df$pard [df$surface_par ==0] <- 0
#summary(df$pard)
df$light_pard <- 0
df$light_pard[df$pard !=0] <- df$light_mol_5p[df$pard !=0]/df$pard[df$pard !=0]
summary(df$light_pard)


# df$pard1[df$surface_par !=0] <- ld1(df$surface_par[df$surface_par !=0], df$water_col[df$surface_par !=0])
# df$pard1 [df$surface_par ==0] <- 0
# #summary(df$pard)
# df$light_pard1 <- 0
# df$light_pard1[df$pard1 !=0] <- df$light_mol_5p[df$pard1 !=0]/df$pard1[df$pard1 !=0]
# summary(df$light_pard1)
# 
# 
# df$pard2[df$surface_par !=0] <- ld2(df$surface_par[df$surface_par !=0], df$water_col[df$surface_par !=0])
# df$pard2 [df$surface_par ==0] <- 0
# #summary(df$pard)
# df$light_pard2 <- 0
# df$light_pard2[df$pard2 !=0] <- df$light_mol_5p[df$pard2 !=0]/df$pard2[df$pard2 !=0]
# summary(df$light_pard2)
# 


df_long <- as.data.frame(df %>%
                           group_by(mergekeyHB,lat,lon) %>%
                           filter(rep == "a") %>%
                           mutate(water_col_mean = mean(water_col),
                                  integral_props = sum(light_pard),
                                  integral_pard = sum(pard),
                                  # integral_props1 = sum(light_pard1),
                                  # integral_pard1 = sum(pard1),
                                  # integral_props2 = sum(light_pard2),
                                  # integral_pard2 = sum(pard2),
                                  bottom_integral = sum(light_mol_5p)))
df_long$prop_integrals = df_long$bottom_integral/df_long$integral_pard
# df_long$prop_integrals1 = df_long$bottom_integral/df_long$integral_pard1
# df_long$prop_integrals2 = df_long$bottom_integral/df_long$integral_pard2

# get n

names(df_long)
df_int <- unique(df_long[,c( # "integral_props",
  "prop_integrals",
  # "integral_props1","prop_integrals1",
  # "integral_props2","prop_integrals2",
  "R_0.25","R_0.5","R_0.75","R_1",
  "D_0.25","D_0.5","D_0.75","D_1",
  "H_0.25","H_0.5","H_0.75","H_1",
  "Hm","Rm_log",#"Dmean",
  "water_col_mean",
  "site","year.x",
  "integral_pard","integral_mol_5p", "folder",
  "lon","lat","dpHB",
  "mergekeyHB")])

length(unique(df$mergekey))
length(unique(df_int$mergekey))
dim(df_int)
sum(duplicated(df_int$mergekey)) # there are 2 duplicates: two same IDs in the same record
# df_int[duplicated(df_int$mergekey),]
# df_int[841:915,]
# df_int[524:527,] #RS17H44
# df_int[409:412,] # N317H34

# leave one out 
df_int <- df_int[!duplicated(df_int$mergekeyHB),]
dim(df_int) 
summary(df_int)

table(df_int$folder)


#write.csv(df_int, "output/df_int_tot_040322.csv", row.names = F)
#write.csv(df_tot, "output/df_tot_040322.csv", row.names = F)
#write.csv(df_tot, "output/df_tot_0524.csv", row.names = F)
# 
# df_int1 <- read.csv("output/df_int_tot_23.csv")
# summary(df_int1)
# table(df_int1$folder)
# df_tot <- read.csv("output/df_tot_030322.csv")

#### model fitting ####

#full
# 
# set.seed(1234)
# brmsHRl_0.25_full <- brms::brm(prop_integrals ~ H_0.25*Rl_0.25 + (1|water_col_mean) + (1|folder), #not scaled
#                          data = df_int, family = Beta(), chains = 6,
#                          iter = 10000, warmup = 2000, thin = 2)
# brmsHRl_0.5_full <- brms::brm(prop_integrals ~ H_0.5*Rl_0.5+ (1|water_col_mean), #not scaled
#                                data = df_int, family = Beta(), chains = 4,
#                                iter = 10000, warmup = 2000, thin = 2)
# brmsHRl_0.75_full <- brms::brm(prop_integrals ~ H_0.75*Rl_0.75+ (1|water_col_mean), #not scaled
#                                data = df_int, family = Beta(), chains = 4,
#                                iter = 10000, warmup = 2000, thin = 2)
# brmsHRl_1_full <- brms::brm(prop_integrals ~ H_1*Rl_1+ (1|water_col_mean), #not scaled
#                                data = df_int, family = Beta(), chains = 4,
#                                iter = 10000, warmup = 2000, thin = 2)



### H and R ####
# 
# set.seed(1234)
# brmsHR_0.25 <- brms::brm(prop_integrals ~ H_0.25*R_0.25+ (1|folder), #not scaled
#                         data = df_int, family = Beta(), chains = 4,
#                         iter = 4000, warmup = 1000)
# brmsR_H_0.25 <- brms::brm(prop_integrals ~ H_0.25+R_0.25+ (1|folder), #not scaled
#                         data = df_int, family = Beta(), chains = 4,
#                         iter = 4000, warmup = 1000)
# brmsR_0.25 <- brms::brm(prop_integrals ~ R_0.25+ (1|folder), #not scaled
#                         data = df_int, family = Beta(), chains = 4,
#                         iter = 4000, warmup = 1000)
# brmsH_0.25 <- brms::brm(prop_integrals ~ H_0.25 + (1|folder), #not scaled
#                    data = df_int, family = Beta(), chains = 4,
#                    iter = 4000, warmup = 1000)
# 
# brmsHR_0.5 <- brms::brm(prop_integrals ~ H_0.5*R_0.5+ (1|folder), #not scaled
#                     data = df_int, family = Beta(), chains = 4,
#                     iter = 4000, warmup = 1000)
# brmsR_H_0.5 <- brms::brm(prop_integrals ~ H_0.5+R_0.5+ (1|folder), #not scaled
#                      data = df_int, family = Beta(), chains = 4,
#                      iter = 4000, warmup = 1000)
# brmsR_0.5 <- brms::brm(prop_integrals ~ R_0.5+ (1|folder), #not scaled
#                    data = df_int, family = Beta(), chains = 4,
#                    iter = 4000, warmup = 1000)
# brmsH_0.5 <- brms::brm(prop_integrals ~ H_0.5 + (1|folder), #not scaled
#                    data = df_int, family = Beta(), chains = 4,
#                    iter = 4000, warmup = 1000)
# 
# brmsHR_0.75 <- brms::brm(prop_integrals ~ H_0.75*R_0.75+ (1|folder), #not scaled
#                     data = df_int, family = Beta(), chains = 4,
#                     iter = 4000, warmup = 1000)
# brmsR_H_0.75 <- brms::brm(prop_integrals ~ H_0.75+R_0.75+ (1|folder), #not scaled
#                      data = df_int, family = Beta(), chains = 4,
#                      iter = 4000, warmup = 1000)
# brmsR_0.75 <- brms::brm(prop_integrals ~ R_0.75+ (1|folder), #not scaled
#                    data = df_int, family = Beta(), chains = 4,
#                    iter = 4000, warmup = 1000)
# brmsH_0.75 <- brms::brm(prop_integrals ~ H_0.75 + (1|folder), #not scaled
#                    data = df_int, family = Beta(), chains = 4,
#                    iter = 4000, warmup = 1000)
# 
# brmsHR_1 <- brms::brm(prop_integrals ~ H_1*R_1+ (1|folder), #not scaled
#                          data = df_int, family = Beta(), chains = 4,
#                          iter = 4000, warmup = 1000)
# brmsR_H_1 <- brms::brm(prop_integrals ~ H_1+R_1+ (1|folder), #not scaled
#                           data = df_int, family = Beta(), chains = 4,
#                           iter = 4000, warmup = 1000)
# brmsR_1 <- brms::brm(prop_integrals ~ R_1+ (1|folder), #not scaled
#                         data = df_int, family = Beta(), chains = 4,
#                         iter = 4000, warmup = 1000)
# brmsH_1 <- brms::brm(prop_integrals ~ H_1 + (1|folder), #not scaled
#                         data = df_int, family = Beta(), chains = 4,
#                         iter = 4000, warmup = 1000)
# 
# 
# # compute and compare  WAIC and LOOic
# model.results <- data.frame(rbind(summary(brmsHR_0.25)$fixed,
#                                   summary(brmsR_H_0.25)$fixed,
#                                   summary(brmsR_0.25)$fixed,
#                                   summary(brmsH_0.25)$fixed,
#                                   summary(brmsHR_0.5)$fixed,
#                                   summary(brmsR_H_0.5)$fixed,
#                                   summary(brmsR_0.5)$fixed,
#                                   summary(brmsH_0.5)$fixed,
#                                   summary(brmsHR_0.75)$fixed,
#                                   summary(brmsR_H_0.75)$fixed,
#                                   summary(brmsR_0.75)$fixed,
#                                   summary(brmsH_0.75)$fixed,
#                                   summary(brmsHR_1)$fixed,
#                                   summary(brmsR_H_1)$fixed,
#                                   summary(brmsR_1)$fixed,
#                                   summary(brmsH_1)$fixed))
# 
# model.results$L <- c(rep(0.25,11),rep(0.5,11),rep(0.75,11),rep(1,11))
# model.results$parameter <- rep(c("intercept","H","R","H*R",
#                                  "intercept","H","R",
#                                  "intercept","R","intercept","H"),4) 
# model.results$param_order <- rep(c(1,4,3,2,
#                                  1,3,2,
#                                  1,2,1,2),4) 
# model.results$model <- rep((c(rep("1. full interaction",4),rep("2. full additive",3),
#                            rep("3. only R",2),rep("4. only H",2))),4)
# 
themeRN <-theme(text = element_text(size = 12),
                axis.title.x = element_text(size = 14),
                axis.title.y = element_text(size = 14),
                axis.text = element_text(size = 12),
                axis.text.x = element_text(vjust = 0.5),
                legend.title = element_text(size = 14),
                legend.text = element_text(size = 14),
                panel.grid = element_blank())
# ggplot()+
#   geom_linerange(data= model.results,
#                  mapping=aes(y=reorder(parameter,desc(param_order)), xmin=l.95..CI, xmax=u.95..CI), size = 3, col = "grey80")+
#   geom_point(data= model.results,
#              mapping=aes(y=reorder(parameter,desc(param_order)), x=Estimate, size = 4))+
#   xlab("Estimates with 95% confidence intervals") +
#   ylab("Response variable") +
#   scale_colour_manual(values = c("grey60", "grey20"))+
#   facet_grid(model~L, scales = "free")+
#   # Big bold line at y=0
#   geom_vline(xintercept=0,size=1, alpha=0.3, linetype="dashed")+
#   theme(panel.grid = element_blank(),strip.text.x = element_text(size = 14), panel.grid.major.y = element_line(colour="grey90", size=0.5),legend.position="none", 
#         axis.text = element_text(size = 12), axis.title = element_text(size = 15)) +
#   themeRN
# 
# ggsave("output/local_variable_fits.pdf", width =20, height = 15, units ="cm",  dpi = 330)


#### H and log R #####

df_int$Rl_0.25 <- log10(df_int$R_0.25)
df_int$Rl_0.5 <- log10(df_int$R_0.5)
df_int$Rl_0.75 <- log10(df_int$R_0.75)
df_int$Rl_1 <- log10(df_int$R_1)

df_int <- df_int[df_int$R_0.25 >=1|df_int$R_0.5 >=1|df_int$R_0.75 >=1|df_int$R_1 >=1, ]
set.seed(1234)
brmsHRl_0.25 <- brms::brm(prop_integrals ~ H_0.25*Rl_0.25+ (1|folder), #not scaled
                         data = df_int, family = Beta(), chains = 4,
                         iter = 4000, warmup = 1000)
brmsRl_H_0.25 <- brms::brm(prop_integrals ~ H_0.25+Rl_0.25+ (1|folder), #not scaled
                          data = df_int, family = Beta(), chains = 4,
                          iter = 4000, warmup = 1000)
brmsRl_0.25 <- brms::brm(prop_integrals ~ Rl_0.25+ (1|folder), #not scaled
                        data = df_int, family = Beta(), chains = 4,
                        iter = 4000, warmup = 1000)
brmsH_0.25 <- brms::brm(prop_integrals ~ H_0.25 + (1|folder), #not scaled
                   data = df_int, family = Beta(), chains = 4,
                   iter = 4000, warmup = 1000)

brmsHRl_0.5 <- brms::brm(prop_integrals ~ H_0.5*Rl_0.5+ (1|folder), #not scaled
                        data = df_int, family = Beta(), chains = 4,
                        iter = 4000, warmup = 1000)
brmsRl_H_0.5 <- brms::brm(prop_integrals ~ H_0.5+Rl_0.5+ (1|folder), #not scaled
                         data = df_int, family = Beta(), chains = 4,
                         iter = 4000, warmup = 1000)
brmsRl_0.5 <- brms::brm(prop_integrals ~ Rl_0.5+ (1|folder), #not scaled
                       data = df_int, family = Beta(), chains = 4,
                       iter = 4000, warmup = 1000)
brmsH_0.5 <- brms::brm(prop_integrals ~ H_0.5 + (1|folder), #not scaled
                   data = df_int, family = Beta(), chains = 4,
                   iter = 4000, warmup = 1000)

brmsHRl_0.75 <- brms::brm(prop_integrals ~ H_0.75*Rl_0.75+ (1|folder), #not scaled
                         data = df_int, family = Beta(), chains = 4,
                         iter = 4000, warmup = 1000)
brmsRl_H_0.75 <- brms::brm(prop_integrals ~ H_0.75+Rl_0.75+ (1|folder), #not scaled
                          data = df_int, family = Beta(), chains = 4,
                          iter = 4000, warmup = 1000)
brmsRl_0.75 <- brms::brm(prop_integrals ~ Rl_0.75+ (1|folder), #not scaled
                        data = df_int, family = Beta(), chains = 4,
                        iter = 4000, warmup = 1000)
brmsH_0.75 <- brms::brm(prop_integrals ~ H_0.75 + (1|folder), #not scaled
                   data = df_int, family = Beta(), chains = 4,
                   iter = 4000, warmup = 1000)

brmsHRl_1 <- brms::brm(prop_integrals ~ H_1*Rl_1+ (1|folder), #not scaled
                      data = df_int, family = Beta(), chains = 4,
                      iter = 4000, warmup = 1000)
brmsRl_H_1 <- brms::brm(prop_integrals ~ H_1+Rl_1+ (1|folder), #not scaled
                       data = df_int, family = Beta(), chains = 4,
                       iter = 4000, warmup = 1000)
brmsRl_1 <- brms::brm(prop_integrals ~ Rl_1+ (1|folder), #not scaled
                     data = df_int, family = Beta(), chains = 4,
                     iter = 4000, warmup = 1000)
brmsH_1 <- brms::brm(prop_integrals ~ H_1 + (1|folder), #not scaled
                        data = df_int, family = Beta(), chains = 4,
                        iter = 4000, warmup = 1000)

# site
brmsHRl_s <- brms::brm(prop_integrals ~ Hm*Rm_log+ (1|folder), #not scaled
                       data = df_int, family = Beta(), chains = 4,
                       iter = 4000, warmup = 1000)
brmsRl_H_s <- brms::brm(prop_integrals ~ Hm+Rm_log+ (1|folder), #not scaled
                        data = df_int, family = Beta(), chains = 4,
                        iter = 4000, warmup = 1000)
brmsRl_s <- brms::brm(prop_integrals ~ Rm_log+ (1|folder), #not scaled
                      data = df_int, family = Beta(), chains = 4,
                      iter = 4000, warmup = 1000)
brmsH_s <- brms::brm(prop_integrals ~ Hm + (1|folder), #not scaled
                     data = df_int, family = Beta(), chains = 4,
                     iter = 4000, warmup = 1000)


#### WAIC and LOO


loo(brmsHRl_s,brmsRl_H_s, brmsRl_s,brmsH_s)
loo(brmsHRl_0.25,brmsRl_H_0.25, brmsRl_0.25,brmsH_0.25)
loo(brmsHRl_0.75,brmsRl_H_0.75, brmsRl_0.75,brmsH_0.75)
loo(brmsHRl_0.5,brmsRl_H_0.5, brmsRl_0.5,brmsH_0.5)



### SM 3 effect size plot ####

model.results.Rl <- data.frame(rbind(summary(brmsHRl_0.25)$fixed,
                                  summary(brmsRl_H_0.25)$fixed,
                                  summary(brmsRl_0.25)$fixed,
                                  summary(brmsH_0.25)$fixed,
                                  summary(brmsHRl_0.5)$fixed,
                                  summary(brmsRl_H_0.5)$fixed,
                                  summary(brmsRl_0.5)$fixed,
                                  summary(brmsH_0.5)$fixed,
                                  summary(brmsHRl_0.75)$fixed,
                                  summary(brmsRl_H_0.75)$fixed,
                                  summary(brmsRl_0.75)$fixed,
                                  summary(brmsH_0.75)$fixed,
                                  summary(brmsHRl_1)$fixed,
                                  summary(brmsRl_H_1)$fixed,
                                  summary(brmsRl_1)$fixed,
                                  summary(brmsH_1)$fixed,
                                  summary(brmsHRl_s)$fixed,
                                  summary(brmsRl_H_s)$fixed,
                                  summary(brmsRl_s)$fixed,
                                  summary(brmsH_s)$fixed))

#model.results.Rl <- data.frame(summary(brmsHRl_0.)$fixed)
model.results.Rl$S <- c(rep("S = 0.25m",11),rep("S = 0.5m",11),rep("S = 0.75m",11),rep("S = 1m",11), rep("site-level",11))
model.results.Rl$parameter <- rep(c("intercept","H","R_log","H*R_log",
                                 "intercept","H","R_log",
                                 "intercept","R_log","intercept","H"),5) 
model.results.Rl$param_order <- rep(c(1,4,3,2,
                                   1,3,2,
                                   1,2,1,2),5) 
model.results.Rl$model <- rep((c(rep("1. full interaction",4),rep("2. full additive",3),
                              rep("3. only R",2),rep("4. only H",2))),5)

model.results.Rl$lower <- exp(model.results.Rl$l.95..CI)/(1+exp(model.results.Rl$l.95..CI))
model.results.Rl$upper <- exp(model.results.Rl$u.95..CI)/(1+exp(model.results.Rl$u.95..CI))
model.results.Rl$estimate <- exp(model.results.Rl$Estimate)/(1+exp(model.results.Rl$Estimate))
exp(0.3)/ (1 + exp(0.3))


ggplot()+
  geom_linerange(data= model.results.Rl,
                 mapping=aes(y=reorder(parameter,desc(param_order)), xmin=lower, xmax=upper), size = 3, col = "grey80")+
  geom_point(data= model.results.Rl,
             mapping=aes(y=reorder(parameter,desc(param_order)), x=estimate, size = 4))+
  xlab("estimates with 95% credible intervals") +
  ylab("parameter") +
  scale_colour_manual(values = c("grey60", "grey20"))+
  facet_grid(model~S, scales = "free")+
  # Big bold line at y=0
  geom_vline(xintercept=0.5,size=1, alpha=0.3, linetype="dashed")+
  theme(panel.grid = element_blank(),strip.text.x = element_text(size = 14), 
        panel.grid.major.y = element_line(colour="grey90", size=0.5),
        legend.position="none", 
        axis.text = element_text(size = 12), axis.title = element_text(size = 15),
        panel.background = element_rect(fill = "white", colour = "black"))


ggsave("output/local_variable_fits_Rl_transformed.png", 
       width =20, height = 15, units ="cm",  dpi = 330)


ggplot()+
  geom_linerange(data= model.results.Rl,
                 mapping=aes(y=reorder(parameter,desc(param_order)), xmin=l.95..CI, xmax=u.95..CI), size = 3, col = "grey80")+
  geom_point(data= model.results.Rl,
             mapping=aes(y=reorder(parameter,desc(param_order)), x=Estimate, size = 4))+
  xlab("Estimates with 95% confidence intervals") +
  ylab("Parameter") +
  scale_colour_manual(values = c("grey60", "grey20"))+
  facet_grid(model~L, scales = "free")+
  # Big bold line at y=0
  geom_vline(xintercept=0,size=1, alpha=0.3, linetype="dashed")+
  theme(panel.grid = element_blank(),strip.text.x = element_text(size = 14), 
      panel.grid.major.y = element_line(colour="grey90", size=0.5),
      legend.position="none", 
      axis.text = element_text(size = 12), axis.title = element_text(size = 15),
      panel.background = element_rect(fill = "white", colour = "black"))

ggsave("output/local_variable_fits_Rl_nontrasformed.pdf", width =20, height = 15, units ="cm",  dpi = 330)


## prediction plots 0.5  ####

plot(conditional_effects(brmsHRl_0.5, effects = "H_0.5:Rl_0.5", surface = TRUE, 
                    stype = "raster", method = "posterior_predict", rug = TRUE,
                    surface_args = (interpolate = TRUE),
                    points = TRUE))

predict.df <- cbind(H_0.5 = as.numeric(brmsHRl_0.5$data$H_0.5), 
                    Rl_0.5 = as.numeric(brmsHRl_0.5$data$Rl_0.5),
                    folder = as.factor(rep("NEW", nrow(brmsHRl_0.5$data))))
predict <- predict(object = brmsHRl_0.5, 
                                        newdata = predict.df,
                                        allow_new_levels = TRUE)

#str(predict)
#class(predict)
#View(predict)
  
predict.df$predict <- predict
#brmsHRl_0.5$data
library(sgeostat)

grid.lines <- 100
names(brmsHRl_0.5$data)

H.pred.full <- seq(min(df_int$H_0.5), max(df_int$H_0.5), length=grid.lines)
R.pred.full <- seq(min(df_int$Rl_0.5), max(df_int$Rl_0.5), length=grid.lines)

mat <- matrix(1, length(R.pred.full), length(H.pred.full))
maty <- matrix(H.pred.full, length(R.pred.full), length(H.pred.full), byrow = TRUE)
matx <- matrix(R.pred.full, length(R.pred.full), length(H.pred.full))

mat.s <- expand.grid(Rl_0.5=R.pred.full, H_0.5=H.pred.full)
pred.df.full <- data.frame(Rl_0.5 = mat.s[,1], H_0.5 = mat.s[,2], 
                           #sdp = rep(mean(mod1@frame$sdp), nrow(mat.s)),
                           folder = as.factor(rep("NEW", nrow(mat.s))))

pred.s <- predict(object = brmsHRl_0.5, newdata = pred.df.full, allow_new_levels=TRUE)

pred.df.full$pred <- pred.s[,1]
#write.csv(pred.df,"data/output/pred2.df.csv")

pred.mat.full = xtabs(pred ~ Rl_0.5  + H_0.5, data = pred.df.full)
# layout(matrix(c(rep(1, 12), rep(2, 12)), nrow=5, byrow=TRUE))
pred.mat.store <- pred.mat.full  
pred.mat.full <- pred.mat.store

pred.mat.full[!in.chull(matx, maty, 
                        df_int$Rl_0.5[chull(df_int$Rl_0.5, df_int$H_0.5)], 
                        df_int$H_0.5[chull(df_int$Rl_0.5, df_int$H_0.5)])] <- NA


png("output/pred_full_0.5_rev.png",  width =20, height = 16, units ="cm", res = 300)

image(R.pred.full,H.pred.full, z = pred.mat.full, col=mako(n=100, begin=0.3), #col=hcl.colors(100, "blues",  alpha=0.9), # axes=FALSE, 
      xlim=c(min(df_int$Rl_0.5), max(df_int$Rl_0.5)), 
      ylim=c(min(df_int$H_0.5), max(df_int$H_0.5)), 
      xlab = expression("log"[10]*"(surface rugosity), dimensionless (R_log)"), 
      ylab = "height range, m (H)" , 
      cex.lab=1.7,
      xaxt = "n", yaxt = "n")

axis(1, at = c(0, 0.25, 0.5, 0.75), cex.axis = 1.5)
axis(2, at = c(0, 0.5, 1), cex.axis = 1.5)
#axis(2, at=seq(1, 2.6, 0.full))
#axis(1, at=seq(0.75, 1.87, 0.full))
points(df_int$Rl_0.5, df_int$H_0.5, col=rgb(1,1,1,0.5), pch=16) 
contour(R.pred.full,H.pred.full, pred.mat.full, method = "edge", levels=c(0.3, 0.4, 0.5),
        add=TRUE, labcex = 1.4, lwd = 1.8, vfont =c("sans serif", "bold")
        )
dev.off()


### site level predictions ####

H.pred.full <- seq(min(df_int$Hm), max(df_int$Hm), length=grid.lines)
R.pred.full <- seq(min(df_int$Rm_log), max(df_int$Rm_log), length=grid.lines)

mat <- matrix(1, length(R.pred.full), length(H.pred.full))
maty <- matrix(H.pred.full, length(R.pred.full), length(H.pred.full), byrow = TRUE)
matx <- matrix(R.pred.full, length(R.pred.full), length(H.pred.full))

mat.s <- expand.grid(Rm_log=R.pred.full, Hm=H.pred.full)
pred.df.full <- data.frame(Rm_log = mat.s[,1], Hm= mat.s[,2], 
                           #sdp = rep(mean(mod1@frame$sdp), nrow(mat.s)),
                           folder = as.factor(rep("NEW", nrow(mat.s))))

pred.s <- predict(object = brmsHRl_s, newdata = pred.df.full, allow_new_levels=TRUE)

pred.df.full$pred <- pred.s[,1]
#write.csv(pred.df,"data/output/pred2.df.csv")

pred.mat.full = xtabs(pred ~ Rm_log  + Hm, data = pred.df.full)
# layout(matrix(c(rep(1, 12), rep(2, 12)), nrow=5, byrow=TRUE))
pred.mat.store <- pred.mat.full  
pred.mat.full <- pred.mat.store

pred.mat.full[!in.chull(matx, maty, 
                        df_int$Rm_log[chull(df_int$Rm_log, df_int$Hm)], 
                        df_int$Hm[chull(df_int$Rm_log, df_int$Hm)])] <- NA


png("output/pred_full_site.png",  width = 12, height = 10, units ="cm", res = 300)

image(R.pred.full,H.pred.full, z = pred.mat.full, col=mako(n=100, begin=0.3), #col=hcl.colors(100, "blues",  alpha=0.9), # axes=FALSE, 
       xlim=c(min(df_int$Rm_log), max(df_int$Rm_log)), 
       ylim=c(min(df_int$Hm), max(df_int$Hm)), 
       xlab = expression("log"[10]*"(Surface Rugosity) (R_log)"), 
       ylab = "Height range, m (H)" , 
       cex.lab=1.1)
#axis(2, at=seq(1, 2.6, 0.full))
#axis(1, at=seq(0.75, 1.87, 0.full))
points(df_int$Rm_log, df_int$Hm, col=rgb(1,1,1,0.25), pch=16) 
contour(R.pred.full,H.pred.full, pred.mat.full, method = "edge", levels=c(0.3, 0.4, 0.25),
        add=TRUE, labcex = 1.4, lwd = 1.8, vfont =c("sans serif", "bold")
)
dev.off()


### prediction plots .25 #####
grid.lines <- 100
names(brmsHRl_0.25$data)

H.pred.full <- seq(min(df_int$H_0.25), max(df_int$H_0.25), length=grid.lines)
R.pred.full <- seq(min(df_int$Rl_0.25), max(df_int$Rl_0.25), length=grid.lines)

mat <- matrix(1, length(R.pred.full), length(H.pred.full))
maty <- matrix(H.pred.full, length(R.pred.full), length(H.pred.full), byrow = TRUE)
matx <- matrix(R.pred.full, length(R.pred.full), length(H.pred.full))

mat.s <- expand.grid(Rl_0.25=R.pred.full, H_0.25=H.pred.full)
pred.df.full <- data.frame(Rl_0.25 = mat.s[,1], H_0.25 = mat.s[,2], 
                           #sdp = rep(mean(mod1@frame$sdp), nrow(mat.s)),
                           folder = as.factor(rep("NEW", nrow(mat.s))))

pred.s <- predict(object = brmsHRl_0.25, newdata = pred.df.full, allow_new_levels=TRUE)

pred.df.full$pred <- pred.s[,1]
#write.csv(pred.df,"data/output/pred2.df.csv")

pred.mat.full = xtabs(pred ~ Rl_0.25  + H_0.25, data = pred.df.full)
# layout(matrix(c(rep(1, 12), rep(2, 12)), nrow=5, byrow=TRUE))
pred.mat.store <- pred.mat.full  
pred.mat.full <- pred.mat.store

pred.mat.full[!in.chull(matx, maty, 
                        df_int$Rl_0.25[chull(df_int$Rl_0.25, df_int$H_0.25)], 
                        df_int$H_0.25[chull(df_int$Rl_0.25, df_int$H_0.25)])] <- NA


png("output/pred_full_0.25.png",    width = 12, height = 10, units ="cm", res = 300)

image(R.pred.full,H.pred.full, z = pred.mat.full, col=mako(n=100, begin=0.3), #col=hcl.colors(100, "blues",  alpha=0.9), # axes=FALSE, 
      xlim=c(min(df_int$Rl_0.25), max(df_int$Rl_0.25)), 
      ylim=c(min(df_int$H_0.25), max(df_int$H_0.25)), 
      xlab = expression("log"[10]*"(Surface Rugosity) (R_log)"), 
      ylab = "Height range, m (H)" , 
      cex.lab=1.1)

#axis(2, at=seq(1, 2.6, 0.full))
#axis(1, at=seq(0.75, 1.87, 0.full))
points(df_int$Rl_0.25, df_int$H_0.25, col=rgb(1,1,1,0.25), pch=16) 
contour(R.pred.full,H.pred.full, pred.mat.full, method = "edge", levels=c(0.3, 0.4, 0.25),
        add=TRUE, labcex = 1.4, lwd = 1.8, vfont =c("sans serif", "bold")
)
dev.off()

#### prediction plots 1m ####

grid.lines <- 100
names(brmsHRl_1$data)

H.pred.full <- seq(min(df_int$H_1), max(df_int$H_1), length=grid.lines)
R.pred.full <- seq(min(df_int$Rl_1), max(df_int$Rl_1), length=grid.lines)

mat <- matrix(1, length(R.pred.full), length(H.pred.full))
maty <- matrix(H.pred.full, length(R.pred.full), length(H.pred.full), byrow = TRUE)
matx <- matrix(R.pred.full, length(R.pred.full), length(H.pred.full))

mat.s <- expand.grid(Rl_1=R.pred.full, H_1=H.pred.full)
pred.df.full <- data.frame(Rl_1 = mat.s[,1], H_1 = mat.s[,2], 
                           #sdp = rep(mean(mod1@frame$sdp), nrow(mat.s)),
                           folder = as.factor(rep("NEW", nrow(mat.s))))

pred.s <- predict(object = brmsHRl_1, newdata = pred.df.full, allow_new_levels=TRUE)

pred.df.full$pred <- pred.s[,1]
#write.csv(pred.df,"data/output/pred2.df.csv")

pred.mat.full = xtabs(pred ~ Rl_1  + H_1, data = pred.df.full)
# layout(matrix(c(rep(1, 12), rep(2, 12)), nrow=5, byrow=TRUE))
pred.mat.store <- pred.mat.full  
pred.mat.full <- pred.mat.store

pred.mat.full[!in.chull(matx, maty, 
                        df_int$Rl_1[chull(df_int$Rl_1, df_int$H_1)], 
                        df_int$H_1[chull(df_int$Rl_1, df_int$H_1)])] <- NA


png("output/pred_full_1.png",   width = 12, height = 10, units ="cm", res = 300)

image(R.pred.full,H.pred.full, z = pred.mat.full, col=mako(n=100, begin=0.3), #col=hcl.colors(100, "blues",  alpha=0.9), # axes=FALSE, 
      xlim=c(min(df_int$Rl_1), max(df_int$Rl_1)), 
      ylim=c(min(df_int$H_1), max(df_int$H_1)), 
      xlab = expression("log"[10]*"(Surface Rugosity) (R_log)"), 
      ylab = "Height range, m (H)" , 
      cex.lab=1.1)

#axis(2, at=seq(1, 2.6, 0.full))
#axis(1, at=seq(0.75, 1.87, 0.full))
points(df_int$Rl_1, df_int$H_1, col=rgb(1,1,1,1), pch=16) 
contour(R.pred.full,H.pred.full, pred.mat.full, method = "edge", levels=c(0.3, 0.4, 1),
        add=TRUE, labcex = 1.4, lwd = 1.8, vfont =c("sans serif", "bold")
)
dev.off()

###prediction plots .75#####
H.pred.full <- seq(min(df_int$H_0.75), max(df_int$H_0.75), length=grid.lines)
R.pred.full <- seq(min(df_int$Rl_0.75), max(df_int$Rl_0.75), length=grid.lines)

mat <- matrix(1, length(R.pred.full), length(H.pred.full))
maty <- matrix(H.pred.full, length(R.pred.full), length(H.pred.full), byrow = TRUE)
matx <- matrix(R.pred.full, length(R.pred.full), length(H.pred.full))

mat.s <- expand.grid(Rl_0.75=R.pred.full, H_0.75=H.pred.full)
pred.df.full <- data.frame(Rl_0.75 = mat.s[,1], H_0.75 = mat.s[,2], 
                           #sdp = rep(mean(mod1@frame$sdp), nrow(mat.s)),
                           folder = as.factor(rep("NEW", nrow(mat.s))))

pred.s <- predict(object = brmsHRl_0.75, newdata = pred.df.full, allow_new_levels=TRUE)

pred.df.full$pred <- pred.s[,1]
#write.csv(pred.df,"data/output/pred2.df.csv")

pred.mat.full = xtabs(pred ~ Rl_0.75  + H_0.75, data = pred.df.full)
# layout(matrix(c(rep(1, 12), rep(2, 12)), nrow=5, byrow=TRUE))
pred.mat.store <- pred.mat.full  
pred.mat.full <- pred.mat.store

pred.mat.full[!in.chull(matx, maty, 
                        df_int$Rl_0.75[chull(df_int$Rl_0.75, df_int$H_0.75)], 
                        df_int$H_0.75[chull(df_int$Rl_0.75, df_int$H_0.75)])] <- NA


png("output/pred_full_0.75.png",  width = 12, height = 10, units ="cm", res = 300)

image(R.pred.full,H.pred.full, z = pred.mat.full, col=mako(n=100, begin=0.3), #col=hcl.colors(100, "blues",  alpha=0.9), # axes=FALSE, 
      xlim=c(min(df_int$Rl_0.75), max(df_int$Rl_0.75)), 
      ylim=c(min(df_int$H_0.75), max(df_int$H_0.75)), 
      xlab = expression("log"[10]*"(Surface Rugosity) (R_log)"), 
      ylab = "Height range, m (H)" , 
      cex.lab=1.1)

#axis(2, at=seq(1, 2.6, 0.full))
#axis(1, at=seq(0.75, 1.87, 0.full))
points(df_int$Rl_0.75, df_int$H_0.75, col=rgb(1,1,1,0.75), pch=16) 
contour(R.pred.full,H.pred.full, pred.mat.full, method = "edge", levels=c(0.3, 0.4, 0.75),
        add=TRUE, labcex = 1.4, lwd = 1.8, vfont =c("sans serif", "bold")
)
dev.off()




## coeff plot for inset #####

model.results.05 <- as.data.frame(summary(brmsHRl_0.5)$fixed)
model.results.05$lower <- exp(model.results.05$'l-95% CI')/(1+exp(model.results.05$'l-95% CI'))
model.results.05$upper <- exp(model.results.05$'u-95% CI')/(1+exp(model.results.05$'u-95% CI'))
model.results.05$estimate <- exp(model.results.05$Estimate)/(1+exp(model.results.05$Estimate))
model.results.05$param_order <- c(1,4,3,2)

model.results.05$parameter <- c("intercept","H","R_log","H*R_log")
                                    
ggplot()+
  geom_linerange(data= model.results.05,
                 mapping=aes(y=reorder(parameter,desc(param_order)), xmin=lower, xmax=upper), size = 3, col = "grey80")+
  geom_point(data= model.results.05,
             mapping=aes(y=reorder(parameter,desc(param_order)), x=estimate, size = 4))+
  geom_vline(xintercept=0.5,size=1, alpha=0.3, linetype="dashed")+
  theme(panel.grid = element_blank(),strip.text.x = element_text(size = 14), 
        panel.grid.major.y = element_line(colour="grey90", size=0.5),
        legend.position="none", 
        axis.text = element_text(size = 12), axis.title = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"))

ggsave("output/local_variable_fits_full0.5.pdf", width =5, height = 3, units ="cm",  dpi = 330)

p <- posterior_linpred(brmsHRl_0.5, transform = TRUE)
str(p)

# fig 4b - posterior distribution of the slopes of mod_rug
posterior_slope <- as.array(brmsHRl_0.5)
dim(posterior_slope)
dimnames(posterior_slope) 
p <-plot(posterior_slope, pars = "b_Rl_0.5")

library(bayesplot)
d <- mcmc_areas(
  posterior_slope, 
  pars = c("b_H_0.5","b_Rl_0.5","b_H_0.5:Rl_0.5"),
  prob = .95, # 90% intervals
  prob_outer = 1, # 99%
  point_est = "mean"
) 

d


### fig 5 ####

model.results.full <- model.results.Rl[model.results.Rl$model == "1. full interaction",]

ggplot()+
  geom_linerange(data= model.results.full,
                 mapping=aes(y=reorder(parameter,desc(param_order)), xmin=lower, xmax=upper), size = 3, col = "grey80")+
  geom_point(data= model.results.full,
             mapping=aes(y=reorder(parameter,desc(param_order)), x=estimate, size = 4))+
  scale_x_continuous(labels = scales::label_number(accuracy = 0.1), breaks = c(0.2, 0.5, 0.8),
) +
  xlab("estimates with 95% credible interval") +
  scale_colour_manual(values = c("grey60", "grey20"))+
  facet_grid(~S)+
  # Big bold line at y=0
  geom_vline(xintercept=0.5,size=1, alpha=0.3, linetype="dashed")+
  theme(panel.grid = element_blank(),strip.text.x = element_text(size = 14), 
        panel.grid.major.y = element_line(colour="grey90", size=0.5),
        legend.position="none", 
        axis.text = element_text(size = 12), axis.title.x = element_text(size = 15),
        axis.title.y = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"))


ggsave("output/local_variable_fits_Rl_transformed.png", 
       width =15, height = 6, units ="cm",  dpi = 330)


# fig 4c - marginal effect of status of the nubbin in mod_stat

s <- plot(marginal_effects(mod_stat))

fig4c <- s[[1]] +
  #geom_count(inherit.aes = FALSE,data = df.tot, aes(x = srug, y = s.pres), fill = "black", alpha = 0.2) +
  #theme_minimal()+ ylim(c(0,1))+
  #scale_color_grey() +
  #scale_fill_grey() +
  labs(x = "\n adult coral status", y = "presence probability \n") +
  #geom_point(aes(x = srug, y = tot.ab, color = alive.dead), inherit.aes = FALSE, size = 2, alpha = 0.1,
  #           position = position_jitter(width = 0.01)) +
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 13), 
        legend.margin = margin(3,3,3,3), legend.title = element_text(size = 13),
        legend.text = element_text(size = 12), legend.position =c(0.8,0.6), 
        panel.background = element_rect(fill = "transparent", color = NA), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent", color = NA), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent", color = NA),
        axis.ticks = element_blank()) # get rid of legend panel)) # get rid of legend panel


#### get totals

summary(df_int$prop_integrals)
df_int %>% 
  group_by(folder) %>% 
  mutate(min = min(prop_integrals),
                                       max= max(prop_integrals),
                                       min_pard = min(integral_pard),
                                       max_pard = max(integral_pard)) %>%
  select(folder,min,max, min_pard, max_pard)

df_int %>% 
  group_by(folder) %>% 
  mutate(min = min(prop_integrals),
         max= max(prop_integrals),
         min_pard = min(integral_pard)*5*60,
         max_pard = max(integral_pard)*5*60,
         min_mol = min(integral_mol_5p)*5*60,
         max_mol = max(integral_mol_5p)*5*60,
         range_mol = max_mol - min_mol,
         range_prop = max-min) %>%
  select(folder,min,max, min_pard, max_pard, min_mol, max_mol, range_prop, range_mol) %>% 
  unique() %>% View()

df_int %>% 
  group_by(folder) %>% 
  mutate(min = min(prop_integrals),
         max= max(prop_integrals),
         min_pard = min(integral_pard)*5*60,
         max_pard = max(integral_pard)*5*60,
         min_mol = min(integral_mol_5p)*5*60,
         max_mol = max(integral_mol_5p)*5*60,
         range_mol = max_mol - min_mol,
         range_prop = max-min) %>%
  select(folder,min,max, min_pard, max_pard, min_mol, max_mol, range_prop, range_mol) %>% 
  unique() %>% summary()

df_int %>% 
  group_by(folder) %>% 
  mutate(min = min(prop_integrals),
         max= max(prop_integrals),
         min_pard = min(integral_pard)*5*60,
         max_pard = max(integral_pard)*5*60,
         min_mol = min(integral_mol_5p)*5*60,
         max_mol = max(integral_mol_5p)*5*60,
         range_mol = max_mol - min_mol,
         range_prop = max-min) %>%
  select(folder,min,max, min_pard, max_pard, min_mol, max_mol, range_prop, range_mol) %>% 
  unique() %>% View()


sort(unique(df_long$datetime_pos))
