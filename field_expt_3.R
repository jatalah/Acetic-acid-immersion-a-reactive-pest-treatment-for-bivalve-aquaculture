library(tidyverse)
library(readxl)
library(vegan)
library(ggord)
source('theme_javier.R')
theme_set(theme_javier())

field_3 <- read_excel('data/Field 3.xlsx')

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
mds_plot

# save plot biplot -----
ggsave(
  plot = mds_plot_bw,
  filename = "figures/field_3/mds_plot_b&w.svg",
  width = 8,
  height = 8)
