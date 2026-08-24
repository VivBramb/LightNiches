library(raster)
library(rgdal)
library(RColorBrewer)
source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")
library(ggplot2)
library(dplyr)
library(viridis) 
library(cowplot)

## figure 1 ####

### here a map of CB that year 
TB17dem <- raster("~/Dropbox/My Dropbox/NC-RR_environment_data/dems/TB17.tif")
data <- raster("~/Dropbox/My Dropbox/NC-RR_environment_data/dems/CB18.tif")
data <- projectRaster(data, crs = crs(TB17dem))
shp <- readOGR("~/Dropbox/My Dropbox/NC-RR_environment_data/shps/CB18.shp")
shp <- spTransform(shp, crs(data))
coords <- shp@coords

brewer.pal(n = 8, name = "BrBG")
# palette
ncols <- 100
newcol <- colorRampPalette(c("#8C510A", "#BF812D",
                             "#DFC27D", "#F6E8C3", "#C7EAE5", "#80CDC1", "#35978F", "#01665E"))
cols <- newcol(ncols)  #apply the function to get 100 colours

# scales example
shp[shp$Unit == "U18",]@coords

x <- -683.8083
y <- 2649.94

newcolcrop <- colorRampPalette(c("#DFC27D", "#F6E8C3", "#C7EAE5", "#80CDC1", "#35978F", "#01665E"))
colscrop <- newcolcrop(100)  #apply the function to get 100 colours

png("~/Dropbox/LightNiches/output/area example.png",  
    width = 500, height = 500, units = "px")
crop <- crop(data, extent(x-0.75, x+0.75, y-0.75, y+0.75))
plot(crop, asp=1, col=colscrop, axes = FALSE, bty="n", box=FALSE, legend = FALSE)
rect(x-0.5,y-0.5,x+0.5, y+0.5, border="black", lwd=4, lty = 6)
rect(x-0.25,y-0.25,x+0.25, y+0.25, border="black", lwd=6)
rect(x-0.375,y-0.375,x+0.375, y+0.375, border="black", lwd=4, lty = 2)
rect(x-0.125,y-0.125,x+0.125, y+0.125, border="black", lwd=4, lty = 3)

text(x, y+0.16, cex = 1.4,"S = 0.25m")
text(x, y+0.285, cex = 1.4,"S = 0.5m")
text(x, y+0.41, cex = 1.4,"S = 0.75m")
text(x, y+0.535, cex = 1.4,"S = 1m")
points(x,y, lwd=1, pch = 16)

dev.off()


# plot dem and rectangles

png("~/Dropbox/LightNiches/output/CB map.png",  
    width = 500, height = 500, units = "px")
plot(data, asp=1, col=cols,axes = FALSE, box=FALSE)
dim = .5
for (i in 1:nrow(coords)) {
  
  lon <- coords[i,1]
  lat <- coords[i,2]
  unit <- shp@data[i,1]
  x0 <- lon - dim/2
  y0 <- lat - dim/2
  text(lon, lat, paste(" ",unit,"\n \n"), col="black", cex = 1)
  rect(x- dim/2, y- dim/2, x + dim/2, y + dim/2, border="black", lwd=4)
  rect(x0, y0, x0 + dim, y0 + dim, border="black", lwd=1)
}

dev.off()

#### fig 2 #####

df_int <- as.data.frame(df_int %>% group_by(site) %>% 
                          mutate(Rm = mean(Rmean), Rr = rank(Rm), Rl_0.5 = log10(R_0.5)) %>% ungroup())

gl <- ggplot(data = df_int, 
             aes(x = reorder(site, Rm_log), y = prop_integrals, fill = site)) +
=======
                          mutate(Rm = mean(Rl_0.5), Rr = rank(Rm)) %>% ungroup())

gl <- ggplot(data = df_int, 
             aes(x = reorder(site, Rm), y = prop_integrals, fill = site)) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .9) +
  geom_point(aes(y = prop_integrals, color = site), 
             position = position_jitter(width = .15), size = 1, alpha = 0.7) +
  expand_limits(x = 1) +
  #coord_flip() + # flip or not?
  scale_color_viridis("site",labels = unique(df_int$site), discrete = TRUE) +
  scale_fill_viridis("site",labels = unique(df_int$site), discrete = TRUE) +
  #labs(tag = "a")+
  theme_bw() +
  labs(x ="site", y = "proportion of light \n available at substratum")+
  theme(panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(size = 11),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12))

ggplot(data = df_int, 
       aes(x = reorder(site, Rm_log), y = prop_integrals, fill = site)) +
  geom_flat_violin(position = position_nudge(x = .2), alpha = .9) +
  geom_point(aes(color = site),
             position = position_jitter(width = .15),
             size = 1, alpha = 0.7) +
  geom_flat_violin(
    data = df_int,
    aes(x = length(unique(df_int$site)) + 1,
        y = prop_integrals),
    inherit.aes = FALSE,
    fill = "black",
    color = "black",
    alpha = 1,
    width = 0.8
  ) +
  expand_limits(x = length(unique(df_int$site)) + 1) +
  scale_color_viridis("site", labels = unique(df_int$site), discrete = TRUE) +
  scale_fill_viridis("site", labels = unique(df_int$site), discrete = TRUE) +
  theme_bw() +
  labs(x ="site", y = "proportion of light \n available at substratum") +
  theme(panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(size = 11),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12))
        text = element_text(size = 15),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16))
gl
ggsave("~/Dropbox/LightNiches/output/lightreadings_ordered.png",
       width = 5, height = 4, dpi = 150, units = "in", device='png')

#review 
df_fig <- rbind(df_int, df_int)
df_fig$site[1:903] <- "all"
df_fig$Rm[1:903] <- 5

f2a <- ggplot(data = df_fig, 
       aes(x = reorder(site, Rm_log), y = prop_integrals, fill = site)) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .9) +
  geom_point(aes(y = prop_integrals, color = site), 
             position = position_jitter(width = .15), size = 1, alpha = 0.4) +
  #coord_flip() + # flip or not?
  scale_color_viridis("site",labels = unique(df_fig$site), discrete = TRUE) +
  scale_fill_viridis("site",labels = unique(df_fig$site), discrete = TRUE) +
  #labs(tag = "a")+
  theme_bw() +
  labs(x ="site", y = "proportion of light \n available at substratum")+
  theme(panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(size = 11),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12))

f2a <- ggplot(data = df_fig, 
              aes(x = reorder(site, Rm_log), y = prop_integrals, fill = site)) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .9) +
  geom_point(aes(y = prop_integrals), 
             position = position_jitter(width = .15), size = 0.8, alpha = 0.4) +
  #coord_flip() + # flip or not?
  scale_color_viridis("site",labels = unique(df_fig$site), discrete = TRUE) +
  scale_fill_viridis("site",labels = unique(df_fig$site), discrete = TRUE) +
  #labs(tag = "a")+
  theme_bw() +
  labs(x ="site", y = "proportion of light \n available at substratum")+
  theme(panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(size = 11),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.key = element_blank()
        )
f2a

f2a <- ggplot(data = df_int, 
              aes(x = reorder(site, Rm_log), y = prop_integrals)) +
  geom_flat_violin(
    position = position_nudge(x = .2),
    fill = "grey80",
    color = NA,
    alpha = .9
  ) +
  geom_point(position = position_jitter(width = .15),
             size = 0.8, alpha = 0.4) +
  geom_flat_violin(
    data = df_int,
    aes(x = length(unique(df_int$site)) + 2,
        y = prop_integrals),
    inherit.aes = FALSE,
    fill = "grey30",
    color = NA,
    width = 0.8
  ) +
  expand_limits(x = length(unique(df_int$site)) + 2) +
  scale_x_discrete(
    limits = c(levels(reorder(df_int$site, df_int$Rm)),"",  "all sites")
  ) +
  theme_bw() +
  labs(x ="site                                 ",
       y = "proportion of light \n available at substratum\n") +
  theme(panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(size = 12),
        axis.title.x = element_text(size = 13),
        axis.ticks.x=element_blank(),
        axis.title.y = element_text(size = 13),
        legend.position = "none")   # remove legend since colors no longer map




f2a
ggsave("~/Dropbox/LightNiches/output/Fig2a_reviewed.png",
       width = 18, height = 10, dpi = 150, units = "cm", device='png', bg = "white")


## Fig 2b ##

ggplot(data = df_int) + 
  geom_point( aes(Rl_0.5,H_0.5, color = site), size = 2, alpha = 0.5) +

ggplot(data = df_int) + 
  geom_point( aes(Rl_0.5,H_0.5, color = site), size = 3, alpha = 0.5) +
  ylab("height range, m (H)") +
  xlab(expression("log"[10]*"(surface rugosity), dimensionless (R_log)")) +
  #guides(fill = guide_legend(order = 0.1))+
  theme_bw()+ 
  #labs(tag = "b")+
  scale_color_viridis("site",labels = unique(df_int$site), discrete = TRUE) +
  scale_fill_viridis("site",labels = unique(df_int$site), discrete = TRUE) +
  theme(panel.border = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(size = 15),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16))

# review
f2b <- ggplot(data = df_int) + 
  geom_point(aes(Rl_0.5, H_0.5,
                 color = prop_integrals,
                 size  = prop_integrals),
             alpha = 0.25) +
  ylab("height range, m (H)") +
  xlab(expression("log"[10]*"(surface rugosity), dimensionless (R_log)")) +
  scale_color_viridis_c(name = " proportion \n of light \n available \n at substratum",
                        option = "mako") +
  scale_size(name = " proportion \n of light \n available \n at substratum") +
  guides(color = guide_legend(),
         size  = guide_legend()) +
  theme_bw()+ 
  theme(panel.border = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(size = 12),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.key = element_blank())
f2b
ggsave("~/Dropbox/LightNiches/output/fig2b_reviewed.png",
       width = 15, height = 10, dpi = 300, units = "cm", device='png')



f2b_narrow <- ggdraw() +
  draw_plot(f2b, x = 0.08, width =0.84)

plot_grid(
  f2a,
  NULL,
  f2b_narrow,
  ncol = 1,
  rel_heights = c(1.2, 0.1, 1.2),
  labels = c("a", "", "b")
)

ggsave("~/Dropbox/LightNiches/output/fig2accepted.png",
       width = 16, height = 17, 
       dpi = 300, units = "cm", device='png',bg = "white")

=======
ggsave("~/Dropbox/LightNiches/output/complexity_ordered.png",
       width = 5, height = 4, dpi = 150, units = "in", device='png')


gl
themeRN

R <- ggplot(data = df_int, 
            aes(x = Rl_0.5, y = prop_integrals, fill = site)) +
  #geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .9) +
  geom_point(aes(y = prop_integrals, color = site), size = 1, alpha = 0.7) +
  expand_limits(x = 1) +
  #coord_flip() + # flip or not?
  scale_color_viridis("site",labels = unique(df_int$site), discrete = TRUE) +
  scale_fill_viridis("site",labels = unique(df_int$site), discrete = TRUE) +
  #labs(tag = "a")+
  theme_bw() +
  labs(x =expression("log"[10]*"(surface rugosity)"), y = "proportion of light \navailable at substratum")+
  theme(panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(size = 15),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16))

R  

ggsave("output/SMrugosityVSlight.png",
       width = 5, height = 4, dpi = 150, units = "in", device='png')

H <- ggplot(data = df_int, 
            aes(x = H_0.5, y = prop_integrals, fill = site)) +
  #geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .9) +
  geom_point(aes(y = prop_integrals, color = site), 
             position = position_jitter(width = .15), size = 1, alpha = 0.7) +
  #geom_smooth()+
  expand_limits(x = 1) +
  #coord_flip() + # flip or not?
  scale_color_viridis("site",labels = unique(df_int$site), discrete = TRUE) +
  scale_fill_viridis("site",labels = unique(df_int$site), discrete = TRUE) +
  #labs(tag = "b")+
  theme_bw() +
  labs(x ="height range, m", y = "proportion of light \navailable at substratum")+
  theme(panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(size = 15),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16))

H  
ggsave("output/SMheightVSlight.png",
       width = 5, height = 4, dpi = 150, units = "in", device='png')

R <- ggplot(data = df_int, 
            aes(x = l_0.5, y = prop_integrals, fill = site)) +
  #geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .9) +
  geom_point(aes(y = prop_integrals, color = site), 
             position = position_jitter(width = .15), size = 1, alpha = 0.7) +
  expand_limits(x = 1) +
  #coord_flip() + # flip or not?
  scale_color_viridis("site",labels = unique(df_int$site), discrete = TRUE) +
  scale_fill_viridis("site",labels = unique(df_int$site), discrete = TRUE) +
  labs(tag = "a")+
  theme_bw() +
  labs(x ="site", y = "proportion of light available at substratum")+
  theme(panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(size = 15),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16))

R  

### SM 1 ####

example <- as.data.frame(df_long %>%
                           filter(folder == "GT17a" & hoboID == "GT03" ))
example <- as.data.frame(df_long %>%
                           filter(folder == "CB17a" & hoboID == "CB13"& datetime_pos > "2017-11-16 07:00:00"))
example%>%
  ggplot() +
  #geom_line(aes(x = datetime_pos, y = light_mol), col = "grey80")+
  geom_area(aes(x = datetime_pos, y = surface_par), alpha = .2, fill = "black")+
  geom_area(aes(x = datetime_pos, y = pard), fill = "black")+
  geom_area(aes(x = datetime_pos, y = light_mol_5p), fill = "#3490A8FF")+
  #geom_smooth(method = "loess")+
  ylab(expression(paste("light (μMol/s/m" ^"2"*")")))+
  xlab("Time") +
  theme_minimal()+
  theme(panel.border = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        text = element_text(size = 15),
        #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16))


example%>%
  ggplot() +
  #geom_line(aes(x = datetime_pos, y = light_mol), col = "grey80")+
  geom_area(aes(x = datetime_pos, y = surface_par), alpha = .2, fill = "black")+
  geom_area(aes(x = datetime_pos, y = pard), fill = "black")+
  geom_area(aes(x = datetime_pos, y = light_mol_5p), fill = "#3490A8FF")+
  geom_line(aes(x = datetime_pos, y = light_mol_5p * 1600/pard), color = "red", size = 2) +
  #geom_area(aes(x = datetime_pos, y = light_mol_5p * 1600/pard), color = "red",fill = "red", alpha = 0.3) +
  scale_y_continuous(name = expression(paste("light (μMol/s/m" ^"2"*")")),
                     sec.axis = sec_axis(~ . / 1600, name = "proportion of light \navailable at substratum\n")
  ) +
  #geom_smooth(method = "loess")+
  xlab("time") +
  theme_minimal()+
  theme(panel.border = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        text = element_text(size = 15),
        #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16))

ggsave("~/Dropbox/LightNiches/output/response_variable.png",
       width = 6, height = 3.5, dpi = 150, units = "in", device='png')

example%>%
  ggplot() +
  #geom_line(aes(x = datetime_pos, y = light_mol), col = "grey80")+
  geom_line(aes(x = datetime_pos, y = surface_par), alpha = .2, fill = "black", size = 2)+
  geom_line(aes(x = datetime_pos, y = pard), fill = "black", size = 2)+
  geom_line(aes(x = datetime_pos, y = light_mol_5p * 1600/pard), color = "red", size = 2) +
  geom_area(aes(x = datetime_pos, y = light_mol_5p * 1600/pard), color = "red",fill = "red", alpha = 0.3) +
  geom_line(aes(x = datetime_pos, y = light_mol_5p), color = "#3490A8FF", size = 2)+
  scale_y_continuous(name = expression(paste("light (μMol/s/m" ^"2"*")")),
                     sec.axis = sec_axis(~ . / 1600, name = "Proportion of light \navailable at substratum\n")
  ) +
  #geom_smooth(method = "loess")+
  xlab("Time") +
  theme_minimal()+
  theme(panel.border = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        text = element_text(size = 15),
        #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16))


### get mean pairwise distance ####

mean_distance <- df_int%>%
  group_by(folder) %>%
  summarise(
    mean_pairwise_dist = mean(
      as.vector(dist(cbind(lon, lat)))),
    median_pairwise_dist = median(
      as.vector(dist(cbind(lon, lat)))),
    min_pairwise_dist = min(
      as.vector(dist(cbind(lon, lat)))),
    max_pairwise_dist = max(
      as.vector(dist(cbind(lon, lat))))
  )
mean_distance$min_mean <- mean(mean_distance$min_pairwise_dist)
mean_distance$min_median <- median(mean_distance$min_pairwise_dist)

###

df_int%>%
  group_by(folder) %>%
  mutate(min_int = min(integral_mol_5p)) %>% View()

#### SM2 ####

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

png("output/pred_full_0.5_SM2rev.png",  width =10, height = 10, units ="cm", res = 300)

image(R.pred.full,H.pred.full, z = pred.mat.full, col=mako(n=100, begin=0.3), #col=hcl.colors(100, "blues",  alpha=0.9), # axes=FALSE, 
      xlim=c(min(df_int$Rl_0.5), max(df_int$Rl_0.5)), 
      ylim=c(min(df_int$H_0.5), max(df_int$H_0.5)), 
      xlab = expression("log"[10]*"(surface rugosity), dimensionless"), 
      ylab = "height range, m" , 
      xaxt = "n", yaxt = "n")

axis(1, at = c(0, 0.25, 0.5, 0.75))
axis(2, at = c(0, 0.5, 1))
#axis(2, at=seq(1, 2.6, 0.full))
#axis(1, at=seq(0.75, 1.87, 0.full))
points(df_int$Rl_0.5, df_int$H_0.5, col=rgb(1,1,1,0.2), pch=16) 
contour(R.pred.full,H.pred.full, pred.mat.full, method = "edge", levels=c(0.3, 0.4, 0.5),
        add=TRUE, lwd = 1.4, vfont =c("sans serif", "bold")
)
dev.off()

# S = 0.25

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

png("output/pred_full_0.25_SM2rev.png",  width =10, height = 10, units ="cm", res = 300)

image(R.pred.full,H.pred.full, z = pred.mat.full, col=mako(n=100, begin=0.3), #col=hcl.colors(100, "blues",  alpha=0.9), # axes=FALSE, 
      xlim=c(min(df_int$Rl_0.25), max(df_int$Rl_0.25)), 
      ylim=c(min(df_int$H_0.25), max(df_int$H_0.25)), 
      xlab = expression("log"[10]*"(surface rugosity), dimensionless"), 
      ylab = "height range, m" , 
      xaxt = "n", yaxt = "n")

axis(1, at = c(0, 0.25, 0.5, 0.75))
axis(2, at = c(0, 0.5, 1))
#axis(2, at=seq(1, 2.6, 0.full))
#axis(1, at=seq(0.75, 1.87, 0.full))
points(df_int$Rl_0.25, df_int$H_0.25, col=rgb(1,1,1,0.2), pch=16) 
contour(R.pred.full,H.pred.full, pred.mat.full, method = "edge", levels=c(0.3, 0.4, 0.5),
        add=TRUE, labcex = 1.1, lwd = 1.4, vfont =c("sans serif", "bold")
)
dev.off()

# S = 0.75

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

png("output/pred_full_0.75_SM2rev.png", width =10, height = 10, units ="cm", res = 300)

image(R.pred.full,H.pred.full, z = pred.mat.full, col=mako(n=100, begin=0.3), #col=hcl.colors(100, "blues",  alpha=0.9), # axes=FALSE, 
      xlim=c(min(df_int$Rl_0.75), max(df_int$Rl_0.75)), 
      ylim=c(min(df_int$H_0.75), max(df_int$H_0.75)), 
      xlab = expression("log"[10]*"(surface rugosity), dimensionless"), 
      ylab = "height range, m" , 
      xaxt = "n", yaxt = "n")

axis(1, at = c(0, 0.25, 0.5, 0.75))
axis(2, at = c(0, 0.5, 1))
#axis(2, at=seq(1, 2.6, 0.full))
#axis(1, at=seq(0.75, 1.87, 0.full))
points(df_int$Rl_0.75, df_int$H_0.75, col=rgb(1,1,1,0.2), pch=16) 
contour(R.pred.full,H.pred.full, pred.mat.full, method = "edge", levels=c(0.3, 0.4, 0.5),
        add=TRUE, labcex = 1.1, lwd = 1.4, vfont =c("sans serif", "bold")
)
dev.off()


# S = 1

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

png("output/pred_full_1_SM2rev.png", width =10, height = 10, units ="cm", res = 300)

image(R.pred.full,H.pred.full, z = pred.mat.full, col=mako(n=100, begin=0.3), #col=hcl.colors(100, "blues",  alpha=0.9), # axes=FALSE, 
      xlim=c(min(df_int$Rl_1), max(df_int$Rl_1)), 
      ylim=c(min(df_int$H_1), max(df_int$H_1)), 
      xlab = expression("log"[10]*"(surface rugosity), dimensionless"), 
      ylab = "height range, m" , 
      xaxt = "n", yaxt = "n")

axis(1, at = c(0, 0.25, 0.5, 0.75))
axis(2, at = c(0, 0.5, 1))
#axis(2, at=seq(1, 2.6, 0.full))
#axis(1, at=seq(0.75, 1.87, 0.full))
points(df_int$Rl_1, df_int$H_1, col=rgb(1,1,1,0.2), pch=16) 
contour(R.pred.full,H.pred.full, pred.mat.full, method = "edge", levels=c(0.3, 0.4, 0.5),
        add=TRUE, labcex = 1.1, lwd = 1.4, vfont =c("sans serif", "bold")
)
dev.off()

# patch = full site

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


png("output/pred_full_site_SM2rev.png",  width =10, height = 10, units ="cm", res = 300)

image(R.pred.full,H.pred.full, z = pred.mat.full, col=mako(n=100, begin=0.3), #col=hcl.colors(100, "blues",  alpha=0.9), # axes=FALSE, 
      xlim=c(min(df_int$Rm_log), max(df_int$Rm_log)), 
      ylim=c(min(df_int$Hm), max(df_int$Hm)), 
      xlab = expression("log"[10]*"(surface rugosity), dimensionless"), 
      ylab = "height range, m" , 
      xaxt = "n", yaxt = "n")

axis(1, at = c(0, 0.25, 0.5, 0.75))
axis(2, at = c(0, 0.5, 1))
#axis(2, at=seq(1, 2.6, 0.full))
#axis(1, at=seq(0.75, 1.87, 0.full))
points(df_int$Rm_log, df_int$Hm, col=rgb(1,1,1,0.2), pch=16) 
contour(R.pred.full,H.pred.full, pred.mat.full, method = "edge", levels=c(0.3, 0.4, 0.5),
        add=TRUE, labcex = 1.1, lwd = 1.4, vfont =c("sans serif", "bold")
)
dev.off()

