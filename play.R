df_int <- df_int_tot [df_int_tot$dim == .5,]

df_int %>%
  group_by(rec) %>%
  sum(R_sum = sum(R))

df_a <- stats::aggregate(R ~ rec + year, df_int, mean)
df_l <- stats::aggregate(light_noscale~ rec +year, df_int, var)

df_al <-merge(df_a, df_l)
   
plot(light_noscale ~ log10(R), df_al)
   
class(df_int$R)
           
names(df_int)
df_int
names(df_int)
df_means <- as.data.frame(df_int %>%
                            group_by(rec,year) %>%
                            filter(int_scaled < 10000)%>%
                            mutate(R_mean = mean(Rl),
                                   Ru_mean = mean(Ru),
                                   light_var = var(int_scaled)))

df <- as.data.frame(df %>%
                            group_by(site,yr) %>% filter(int_scaled <10000) %>%
                            mutate(Ru_mean = mean(Ru),
                                   light_var = var(light_pard)))

names(df)

unique(df_means$light_var)

plot(int_scaled~log(Ru_mean),df_means)
boxplot(int_scaled~log(Ru_mean),df_means)
