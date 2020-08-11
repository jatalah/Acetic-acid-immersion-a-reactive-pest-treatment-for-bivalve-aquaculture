library(tidyverse)
library(readxl)
library(ggfortify)
library(broom)
library(caret)
library(lsmeans)
source('theme_javier.R')
theme_set(theme_javier())
library(janitor)

# Read and prepare data-------------
field_2 <- 
  read_excel('data/Field 2.xlsx') %>%
  mutate(Treat = paste(Time_in_bath , `Time_air-drying`, sep = " | "),
         Treatment = fct_relevel(if_else(Time_in_bath==0, "Control", Treat),"Control"),
         Sampling_time = fct_recode(factor(Sampling_time),  `2 mo` = "2",`48 h`= "48"),
         Sampling_time = fct_relevel(Sampling_time, "48 h"),
         Treat1 = paste(Treat, Sampling_time, sep = "_"),
         Treat1 = fct_collapse(Treat1, 
                               "Control_2 mo" = c("0 | 0_2 mo", "0 | 45_2 mo","0 | 60_2 mo"),
                               "Control_48 h" = c("0 | 0_48 h", "0 | 45_48 h","0 | 60_48 h")))


tabyl(field_2, `Time_air-drying`)
tabyl(field_2, `Total_ treatment_time`)
tabyl(field_2, Time_in_bath)
tabyl(field_2, Treat1)


f2_long <- 
  field_2 %>% 
  gather(key, value, Number_GSM:Weight_fouling) %>%  
  drop_na(value) %>% 
  mutate(key = fct_recode(key,
                          `GSM (no. ind.)` = "Number_GSM",
                          `GSM biomass (g)` = "Weight_GSM",
                          `Fouling biomass (g)` = "Weight_fouling"),
         key = fct_relevel(key, "GSM (no. ind.)", "GSM biomass (g)" ))

                          
bar_plot_f2 <- 
  ggplot(f2_long, aes(
  x = Treatment,
  y = value ,
  group = factor(Sampling_time),
  fill = factor(Sampling_time)
)) +
  facet_grid(rows = vars(key), cols = vars(`Total_ treatment_time`), scales = "free", space = 'free_x') +
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
  # scale_y_log10() +
  scale_fill_grey(name = 'Treatment \n time') +
  labs(x = "Time in bath | time in air (s)") +
  theme(axis.title.y = element_blank())
bar_plot_f2

# save plot bm plot -----
ggsave(
  plot = bar_plot_f2,
  filename = "figures/field_2/bar_plot_f2.svg",
  width = 6,
  height = 5
)

# transform response variable using Box-Cox---------------
t_param2 <-
  field_2 %>%
  dplyr::select(Number_GSM:Weight_fouling) %>% 
  data.frame %>% 
  preProcess(method = "BoxCox") # ((x ^ lambda) - mean(x ^ lambda)) / sd(x ^ lambda)


field_2_t <-
  predict(t_param2, data.frame(field_2 %>% dplyr::select(Number_GSM:Weight_fouling))) %>% 
  bind_cols(field_2 %>% select (`Total_ treatment_time` :Replicate,Treat, Treatment, Treat1), .)

f2_long_t <- 
  field_2_t %>% 
  gather(key, value, Number_GSM:Weight_fouling) %>% 
  drop_na(value) %>% 
  mutate(key = fct_recode(key,
                          `GSM (no. ind.)` = "Number_GSM",
                          `GSM biomass (g)` = "Weight_GSM",
                          `Fouling biomass (g)` = "Weight_fouling"),
         key = fct_relevel(key, "GSM (no. ind.)", "GSM biomass (g)" ))


# Compare controls -------
f2_control_anovas <-
  f2_long %>%
  filter(Treatment=="Control") %>% 
  group_by(key) %>%
  nest() %>%
  mutate(
    lm_models = map(.x = data, ~ lm(value ~ Treat+Sampling_time, data = .x)),
    lm_anova = map(lm_models, anova),
    lm_summary = map(.x = lm_anova, tidy),
    lm_pairwise_treat = map(.x = lm_models,
                            ~lsmeans(.x,
                                     pairwise ~Treat,
                                     adjust = "none") %>% 
                              .$contrast %>% 
                              data.frame())
  )

f2_control_anovas %>% 
  select(lm_summary) %>% 
  unnest(lm_summary) %>% 
  write_csv('tables/field_2/ANOVA_tables_controls.csv')

# run anovas---------
f2_anovas <-
  f2_long_t %>%
  group_by(key) %>%
  nest() %>%
  mutate(
    lm_models = map(.x = data, ~ lm(value ~ Treatment* Sampling_time, data = .x)),
    lm_anova = map(lm_models, anova),
    lm_summary = map(.x = lm_anova, tidy),
    lm_pairwise_treat = map(.x = lm_models,
                      ~lsmeans(.x,
                               pairwise ~Treatment,
                               adjust = "none") %>% 
                        .$contrast %>% 
                        data.frame())
  )

f2_long_t %>%
  group_by(key) %>%
  nest() %>%
  mutate(
    lm_models = map(.x = data, ~ lm(value ~ Treat, data = .x)),
    pairwise_all = map(.x = lm_models,
                            ~lsmeans(.x,
                                     pairwise ~Treat,
                                     adjust = "none") %>% 
                              .$contrast %>% 
                              data.frame())
  ) %>% 
  select(pairwise_all) %>% 
  unnest(pairwise_all) %>% 
  write_csv('tables/field_2/pairwise_all.csv', na = "")


# separate all treatments -----------------
f2_long_t %>%
  group_by(key) %>%
  nest() %>%
  mutate(
    lm_models = map(.x = data, ~ lm(value ~ Treat1, data = .x)),
    pairwise_all = map(.x = lm_models,
                       ~lsmeans(.x,
                                pairwise ~Treat1,
                                adjust = "none") %>% 
                         .$contrast %>% 
                         data.frame())
  ) %>% 
  select(pairwise_all) %>% 
  unnest(pairwise_all) %>% 
  write_csv('tables/field_2/pairwise_all.csv', na = "")



f2_long_t %>%
  group_by(key) %>%
  nest() %>%
  mutate(
    lm_models = map(.x = data, ~ lm(value ~ Treatment, data = .x)),
    pairwise_all = map(.x = lm_models,
                       ~lsmeans(.x,
                                pairwise ~Treatment,
                                adjust = "none") %>% 
                         .$contrast %>% 
                         data.frame())
  ) %>% 
  select(pairwise_all) %>% 
  unnest(pairwise_all) %>% 
  write_csv('tables/field_2/pairwise_all_pooled_controls.csv', na = "")


# save ANOVA tables--------
f2_anovas %>% 
  select(lm_summary) %>% 
  unnest(lm_summary) %>% 
  write_csv('tables/field_2/ANOVA_tables.csv')

# save pairwise tables-----
f2_anovas %>% 
  select(lm_pairwise_treat) %>% 
  unnest(lm_pairwise_treat) %>% 
  write_csv('tables/field_2/pairwise_treatment.csv')

