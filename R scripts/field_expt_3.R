library(tidyverse)
library(readxl)
library(vegan)
library(ggord)
source('theme_javier.R')
theme_set(theme_javier())

field_3 <- 
  read_excel('data/Field 3.xlsx') 

# PERMANOVA---------
permanova <-
  adonis(
    sqrt(sqrt(field_3[,-c(1:4)])) ~ Treatment_or_control * Line_position    * Sampling_time  ,
    method = 'bray',
    permutations = 999,
    data = field_3
  )

permanova$aov.tab %>% 
  data.frame() %>%
  rownames_to_column(var = "term") %>% 
  write_csv('tables/field_3/permanova_table.csv',na = "")

# 04 MDS ordination -------------
mds <- metaMDS( sqrt(sqrt(field_3[,-c(1:4)])) , distance = 'bray')

## MDS biplot 
mds_plot <-
ggord(
mds,
grp_in = paste(field_3$Treatment_or_control,field_3$Sampling_time, sep = '-'),
poly = F,
alpha = 1,
ellipse = T,
arrow = 0,
repel = T,
text = .01,
vec_ext = .7
) +
theme_javier() +
scale_color_viridis_d(end = .8) +
annotate(geom = 'text', label = paste0("Stress = ",round(mds$stress,2)), x = .4, y = .4)
mds_plot

# save plot biplot -----
ggsave(
  plot = mds_plot,
  filename = "figures/field_3/mds_plot.svg",
  width = 8,
  height = 8)

# B&W version---
mds_plot_bw <- 
ggord(
  mds,
  grp_in = paste(field_3$Treatment_or_control,field_3$Sampling_time, sep = '-'),
  poly = F,
  alpha = 1,
  ellipse = T,
  arrow = 0,
  repel = T,
  text = .01,
  vec_ext = .7
) +
  theme_javier() +
  scale_color_grey(start =  0 ,end = 0) +
  scale_shape_manual(values = c(0,1,2,5)) +
  annotate(geom = 'text', label = paste0("Stress = ",round(mds$stress,2)), x = .4, y = .4)
mds_plot_bw

# save plot biplot -----
ggsave(
  plot = mds_plot_bw,
  filename = "figures/field_3/mds_plot_b&w.svg",
  width = 8,
  height = 8)


# Total cover 
field_3 <- 
  field_3 %>% 
  mutate(total_cover = Bryozoans_pooled	+Hydroids	+ Other_fouling,
         Sampling_time = fct_relevel(Sampling_time, "Pre")) 
  
bar_plot_n <-
  ggplot(
    field_3,
    aes(
      y = total_cover,
      x = Treatment_or_control,
      group =  Sampling_time,
      fill = Sampling_time
    )
  ) +
  stat_summary(
    fun = mean,
    geom = "bar",
    position = position_dodge(width = .9),
    color = 1
  ) +
  stat_summary(
    fun.data = mean_se,
    position = position_dodge(width = .9),
    geom = "errorbar",
    width = 0.2
  ) +
  theme_javier(base_family = 12) +
  scale_fill_manual(values = c('gray40', 'transparent'), name = 'Sampling time') +
  labs(x = "Treatment", y = "Mean total cover ± S.E.")
bar_plot_n
ggsave(
  plot = bar_plot_n,
  filename = "figures/field_3/bar_plot_total_cover.svg",
  width = 5.4,
  height = 2.6)

anova_n <- aov(total_cover~Treatment_or_control * Line_position    * Sampling_time, data = field_3)
anova(anova_n)

tidy(anova_n) %>%
  rename(F = "statistic") %>%
  write_csv('tables/field_3/anova_table_total_cover.csv', na = "")
