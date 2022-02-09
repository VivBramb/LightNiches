#### diagnostic plots #######


df_int <- df_int_tot [df_int_tot$dim == .25,]

df_int %>%
  group_by(rec) %>%
  mutate(R_sum = sum(R))

df_a <- stats::aggregate(R ~ rec + year.x, df_int, mean)
df_l <- stats::aggregate(light_lux~ rec +year.x, df_int, var)
df_al <-merge(df_a, df_l)
plot(light_lux ~ log10(R), df_al)

df_ld <- stats::aggregate(integral_ld~ rec +year.x, df_int, var)   
df_ald <-merge(df_a, df_ld)
plot(integral_ld ~ log10(R), df_ald)


           
names(df_int)
df_int
names(df_int)
df_means <- as.data.frame(df_int %>%
                            group_by(rec,year.x) %>%
                            mutate(R_mean = mean(Rclip),
                                   R_var = var(Rclip),
                                   prop_var = var(prop_integrals),
                                   light_var = var(integral_mol_5p)))

df <- as.data.frame(df %>%
                            group_by(site,yr) %>% filter(int_scaled <10000) %>%
                            mutate(Ru_mean = mean(Ru),
                                   light_var = var(light_pard)))

names(df)

unique(df_means$light_var)

plot(prop_var~log(R_mean),df_means)
abline(lm(prop_var~log(R_mean),df_means))
plot(prop_var~log(R_var),df_means)
abline(lm(prop_var~log(R_var),df_means))

plot(light_var~log(R_mean),df_means)
abline(lm(light_var~log(R_mean),df_means))
plot(light_var~log(R_var),df_means)
abline(lm(light_var~log(R_var),df_means))



boxplot(int_scaled~log(Ru_mean),df_means)
