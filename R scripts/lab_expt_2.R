library(tidyverse)
library(readxl)
library(ggfortify)
library(broom)
library(lsmeans)
library(caret)
source('theme_javier.R')
theme_set(theme_javier())

# Lab experiment 2 --------
lab_2_long <-
  read_excel('data/Laboratory 2.xlsx') %>%
  mutate(Treat = paste(Time_in_bath , `time_air-drying`, sep = " | "),
         Treatment = fct_relevel(if_else(Acetic_acid ==0, "Control", Treat), "Control"), #merge controls
         C_vs_T = fct_relevel(if_else(Acetic_acid ==0, "Control", "Treated"), "Control")) %>%
  gather(taxa, migrated,  Perna_migrated:Mytilus_migrated) %>%
  mutate(taxa = str_remove(taxa, "_migrated"),
         acetic_acid = Acetic_acid,
         Acetic_acid = fct_recode(
           factor(Acetic_acid),
           Controls = "0",
           `2 %` = "2",
           `4 %` = "4"
         ),
         Treat1 = paste(Treat, Acetic_acid, sep = "_"),
         Treat1 = fct_collapse(Treat1, Controls = c("0 | 60_Controls", "60 | 0_Controls")))

lab2 <- 
lab_2_long %>% 
  spread(taxa, migrated)

janitor::tabyl(lab_2, Time_in_bath , `time_air-drying`)
janitor::tabyl(lab_2, Treat, Acetic_acid)
janitor::tabyl(lab_2_long, Treat1)


# bar plot ------------
lab2_plot_migrated <- 
  ggplot(lab_2_long,
       aes(x = Treatment,
           y = migrated,
           fill = taxa)) +
  facet_grid(cols = vars(Acetic_acid), scales = "free_x", space = "free_x") +
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
  scale_fill_grey(name = "Taxa") +
  labs(x = "Treatment (bath | air) ", y = 'Migrated (%)')

lab2_plot_migrated

# save plot -----
ggsave(
  plot = lab2_plot_migrated,
  filename = "figures/lab_2/lab2_plot_migrated.svg",
  width = 8,
  height = 2.5
)

# ANOVAs migrated ------------

# Transform data -----
t_param_lab2 <-
  lab2 %>%
  dplyr::select(Mytilus,Perna) %>% 
  data.frame %>% 
  preProcess(method = "BoxCox") # ((x ^ lambda) - mean(x ^ lambda)) / sd(x ^ lambda)


lab2_t <-
  predict(t_param_lab2, data.frame(lab2 %>% dplyr::select(Mytilus:Perna))) %>% 
  bind_cols(lab2 %>% select (Acetic_acid :Treat, Treat1), .)


Anova_lab2 <- 
  lab2_t %>%
  gather(taxa, migrated,  Perna:Mytilus) %>%
  mutate(Treatment = fct_relevel(if_else(Acetic_acid=="Controls", "Control", Treat),"Control")) %>% 
  group_by(taxa) %>%
  nest() %>%
  mutate(
    anova_controls = map(.x = data, ~ tidy(anova(lm(migrated ~Treat , data = .x)))),
    anova_all = map(.x = data, ~ lm(migrated ~Treat1 , data = .x)),
    anova_treated = map(.x = data, ~ tidy(anova(lm(migrated ~Acetic_acid*`time_air-drying`*Time_in_bath , 
                                                   data = filter(.x, Treat!="Control"))))))

Anova_lab2 %>%
  select(anova_controls) %>% 
  unnest(cols = anova_controls) %>% 
  write_csv('tables/lab_2/ANOVA_lab2_control_vs_treated_migration.csv')


Anova_lab2 %>%
  select(anova_treated) %>% 
  unnest(cols = anova_treated) %>% 
  write_csv('tables/lab_2/ANOVA_lab2_treated_migration.csv')

## All pairwise---
Anova_lab2 %>% 
mutate(pairwise_all =  map(.x = anova_all,
                        ~lsmeans(.x,
                                 pairwise ~Treat1,
                                 adjust = "none") %>% 
                          .$contrast %>% 
                          data.frame())) %>% 
  select(pairwise_all) %>% 
  unnest(pairwise_all) %>% 
  write_csv('tables/lab_2/pairwise_all.csv', na = "")
