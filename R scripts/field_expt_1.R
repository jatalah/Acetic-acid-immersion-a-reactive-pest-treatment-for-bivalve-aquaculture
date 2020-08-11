library(tidyverse)
library(readxl)
library(ggfortify)
library(broom)
library(lsmeans)
library(caret)
source('theme_javier.R')
theme_set(theme_javier())


# Field experiment 1 data preparation----------------
field_1 <-
  read_excel('data/Field 1.1.xlsx') %>%
  mutate(
    Treated = if_else(Acetic_acid == 0, "Control", "Acetic acid"),
    Sampling_time = fct_recode(
      factor(Sampling_time),
      `48 h` = "48",
      `1 mo` = "1",
      `3 mo` = "3"
    ),
    Time_in_bath = fct_recode(
      factor(Time_in_bath),
      `No bath` = "0",
      `30 sec` = "30",
      `60 sec` = "60"
    ),
    Sampling_time = fct_relevel(Sampling_time, "48 h"),
    acetic_acid = Acetic_acid,
    Acetic_acid = fct_recode(
      factor(Acetic_acid),
      Control = "0",
      `1 %` = "1",
      `2 %` = "2",
      `4 %` = "4"
    ),
    Treated = if_else(Acetic_acid=="Control", "Control","Treated"),
    Number_BM = Number_BM_4mm + Number_BM_2mm + Number_BM_1mm
  ) %>% 
  drop_na(Number_GSM) %>% 
  filter(Number_GSM>10) %>% 
  slice(-c(5,18))# remove outlier 

# transform response variable using Box-Cox---------------
t_param <-
  field_1 %>%
  dplyr::select(Number_GSM:Weight_GSM,
                Number_BM,
                Weight_BM_4mm,
                Weight_fouling) %>%
  data.frame %>%
  preProcess(method = "BoxCox")# ((x ^ lambda) - mean(x ^ lambda)) / sd(x ^ lambda)



field_1_t <-
  predict(t_param, data.frame(
    field_1 %>% dplyr::select(
      Number_GSM:Weight_GSM,
      Number_BM,
      Weight_BM_4mm,
      Weight_fouling
    )
  )) %>%
  bind_cols(field_1 %>% dplyr::select (acetic_acid,Acetic_acid:Replicate, Treated), .)


# Treated vs controls plots--------
f1_long <- 
  field_1 %>% 
  gather(key, value, c(Number_GSM:Weight_GSM,Number_BM, Weight_BM_4mm, Weight_fouling))


f1_long %>% 
  ggplot(aes(value)) +
  geom_histogram() +
  facet_wrap(~key, scales = 'free')


f1_long_t <- 
  field_1_t %>% 
  gather(key, value, c(Number_GSM:Weight_GSM,Number_BM, Weight_BM_4mm, Weight_fouling)) %>% 
  drop_na(value)

ggplot(f1_long, aes(
           x = Treated,
           y = value ,
           group = factor(Sampling_time),
           fill = factor(Sampling_time)
         )) +
  facet_wrap(~key, scales = 'free') +
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
  scale_fill_grey(name = "Sampling time") 

t_vs_c_lms <-
  f1_long_t %>%
  group_by(key) %>%
  nest() %>%
  mutate(
    lm_models = map(.x = data, ~ lm(value ~ Treated * Sampling_time, data = .x)),
    lm_anova = map(lm_models, ~car::Anova(.x, type = 2)),
    lm_summary = map(.x = lm_anova, tidy),
    lm_pairwise = map(.x = lm_models,
                      ~lsmeans(.x,
                          pairwise ~Treated*Sampling_time,
                          adjust = "none") %>% 
      .$contrast %>% 
      data.frame()
  ))

# save ANOVA tables--------
t_vs_c_lms %>% 
  dplyr::select(lm_summary) %>% 
  unnest(lm_summary) %>% 
  write_csv('tables/field_1/ANOVA_tables_treated_vs_controls.csv', na = "")

# save pairwise tables-----
t_vs_c_lms %>% 
  dplyr::select(lm_pairwise) %>% 
  unnest(lm_pairwise) %>% 
  write_csv('tables/field_1/pairwise_treated_vs_controls.csv', na = "")

# All plots -----
plot_all_f1 <- 
  ggplot(f1_long, 
       aes(x = Acetic_acid,
           y = value,
           group =Sampling_time,
           fill = Sampling_time
       )) +
  facet_grid(key ~Time_in_bath, scales = "free", space = "free_x") +
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
  scale_fill_grey(name = "Sampling time") +
  labs(x = "Acetic acid concentration", y = "Mean ± SE")


# save plot bm plot -----
ggsave(
  plot = plot_all_f1,
  filename = "figures/field_1/plot_all_f1.svg",
  width = 8,
  height = 8)
  
### Treated ANOVAs-----------
treated_lms <-
  f1_long_t %>%
  filter(Treated =="Treated") %>% 
  group_by(key) %>%
  nest() %>%
  mutate(
    lm_models = map(.x = data, ~ lm(value ~ Time_in_bath*Sampling_time*Acetic_acid, data = .x)),
    lm_anova = map(lm_models, anova),
    lm_summary = map(.x = lm_anova, tidy),
    pairwise_aa = map(.x = lm_models,
                      ~lsmeans(.x,
                               pairwise ~Acetic_acid,
                               adjust = "tukey") %>% 
                        .$contrast %>% 
                        data.frame()),
    pairwise_stime = map(.x = lm_models,
                      ~lsmeans(.x,
                               pairwise ~Sampling_time,
                               adjust = "tukey") %>% 
                        .$contrast %>% 
                        data.frame()),
    pairwise_btime = map(.x = lm_models,
                         ~lsmeans(.x,
                                  pairwise ~Time_in_bath,
                                  adjust = "tukey") %>% 
                           .$contrast %>% 
                           data.frame()),
    pairwise_aa_stime = map(.x = lm_models,
                      ~lsmeans(.x,
                               pairwise ~Acetic_acid*Sampling_time,
                               adjust = "tukey") %>% 
                        .$contrast %>% 
                        data.frame()),
    pairwise_aa_btime = map(.x = lm_models,
                           ~lsmeans(.x,
                                    pairwise ~Acetic_acid*Sampling_time,
                                    adjust = "tukey") %>% 
                             .$contrast %>% 
                             data.frame())
    )


# save ANOVA tables--------
treated_lms %>% 
  dplyr::select(lm_summary) %>% 
  unnest(lm_summary) %>% 
  write_csv('tables/field_1/ANOVA_tables_treated.csv', na = "")

# save pairwise tables-----
treated_lms %>% 
  dplyr::select(pairwise_aa) %>% 
  unnest(pairwise_aa) %>% 
  write_csv('tables/field_1/pairwise_treated_acetic_acid.csv', na = "")


# save pairwise tables-----
treated_lms %>% 
  dplyr::select(pairwise_stime) %>% 
  unnest(pairwise_stime) %>% 
  write_csv('tables/field_1/pairwise_treated_sampling_time.csv', na = "")


# save pairwise tables-----
treated_lms %>% 
  dplyr::select(pairwise_stime) %>% 
  unnest(pairwise_stime) %>% 
  write_csv('tables/field_1/pairwise_treated_bath_time.csv', na = "")

# save pairwise tables-----
treated_lms %>% 
  select(pairwise_aa_stime) %>% 
  unnest(pairwise_aa_stime) %>% 
  write_csv('tables/field_1/pairwise_treated_AA_s_time.csv', na = "")

# save pairwise tables-----
treated_lms %>% 
  select(pairwise_aa_btime) %>% 
  unnest(pairwise_aa_btime) %>% 
  write_csv('tables/field_1/pairwise_treated_AA_bath_time.csv', na = "")


# Pairwise all------------
f1_long_t %>%
  group_by(key) %>%
  nest() %>%
  mutate(
    lm_models = map(.x = data, ~ lm(value ~ Time_in_bath*Sampling_time*Acetic_acid, data = .x)),
    pairwise_all = map(.x = lm_models,
                       ~lsmeans(.x,
                                pairwise ~Time_in_bath*Sampling_time*Acetic_acid,
                                adjust = "none") %>% 
                         .$contrast %>% 
                         data.frame())
  ) %>% 
  select(pairwise_all) %>% 
  unnest(pairwise_all) %>% 
  write_csv('tables/field_1/pairwise_all.csv', na = "")



# individual analysis for plot blue mussels-----
field1_bm_plot <- 
  ggplot(bm_count,
       aes(
         x = Acetic_acid,
         y = count,
         group =Sampling_time,
         fill = Sampling_time
       )) +
  facet_grid(cols = vars(Time_in_bath), scales = "free_x", space = "free_x") +
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
  scale_y_log10() +
  scale_fill_grey(name = "Sampling time") +
  labs(x = "Acetic acid concentration", y = "Mean no. of blue mussel ± SE")
field1_bm_plot

 # save plot bm plot -----
ggsave(
  plot = field1_bm_plot,
  filename = "figures/field_1/field1_bm_plot.svg",
  width = 8,
  height = 2.5
)

# # ANOVA treated ----
# bm_count_t <- filter(bm_count, Treated=="Yes")
# 
# m_f1 <- lm(log(count+1)~Time_in_bath*Sampling_time*Acetic_acid, bm_count_t)
# autoplot(m_f1)
# anova(m_f1) 
# write_csv(tidy(m_f1), 'tables/field_1/ANOVA_field1_treated.csv') 
# 
# # post-hoc test ----------
# # interaction term
# posthoc_m_f1_interaction <-
#   lsmeans(m_f1,
#           pairwise ~ Sampling_time:Acetic_acid ,
#           adjust = "tukey") %>% 
#   write_csv('tables/field_1/posthoc_treated_vs_controls_interaction.csv')
# 
# posthoc_m_f1_interaction$contrasts %>% 
#   data.frame() %>% 
#   write_csv('tables/field_1/field1_posthoc_treated_interaction.csv')
# 
# # Acetic acid------
# posthoc_m_f1_AA <-
#   lsmeans(m_f1,
#           pairwise ~ Acetic_acid,
#           adjust = "none") %>%
#   .$contrast %>%
#   data.frame() %>%
#   write_csv('tables/field_1/field_1_pairwise_AA.csv')

# weight biofouling-----
plot_biof_ww <- 
  ggplot(field_1,
       aes(
         x = Acetic_acid,
         y = Weight_fouling ,
         group = factor(Sampling_time),
         fill = factor(Sampling_time)
       )) +
  facet_grid(cols = vars(Time_in_bath),
             scales = "free_x",
             space = "free_x") +
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
  scale_y_log10() +
  scale_fill_grey(name = "Sampling time") +
  labs(x = "Acetic Acid concentration", y = "Mean biofouling weight (g ± SE)")
plot_biof_ww

# save biofouling plot -----
ggsave(
  plot = plot_biof_ww,
  filename = "figures/field_1/plot_biof_ww.svg",
  width = 8,
  height = 2.5
)


# # ANOVA biofoling treated vs untreated ----
# m_f1 <- lm(Weight_fouling~Time_in_bath*Sampling_time*Acetic_acid, bm_count_t)
# autoplot(m_f1)
# anova(m_f1) 
# write_csv(tidy(m_f1), 'tables/field_1/ANOVA_field1_treated.csv') 

# # post-hoc test ----------
# # interaction term
# posthoc_m_f1_interaction <-
#   lsmeans(m_f1,
#           pairwise ~ Sampling_time:Acetic_acid ,
#           adjust = "tukey") %>% 
#   write_csv('tables/field_1/posthoc_treated_vs_controls_interaction.csv')
# 
# posthoc_m_f1_interaction$contrasts %>% 
#   data.frame() %>% 
#   write_csv('tables/field_1/field1_posthoc_treated_interaction.csv')
# 
# 
# 
# # ANOVA biofouling  ----
# biofoiling_t <- filter(bm_count, Treated=="Yes")
# 
# m_f1 <- lm(log(count+1)~Time_in_bath*Sampling_time*Acetic_acid, bm_count_t)
# autoplot(m_f1)
# anova(m_f1) 
# write_csv(tidy(m_f1), 'tables/field_1/ANOVA_field1_treated.csv') 
# 
# # post-hoc test ----------
# # interaction term
# posthoc_m_f1_interaction <-
#   lsmeans(m_f1,
#           pairwise ~ Sampling_time:Acetic_acid ,
#           adjust = "tukey") %>% 
#   write_csv('tables/field_1/posthoc_treated_vs_controls_interaction.csv')
# 
# posthoc_m_f1_interaction$contrasts %>% 
#   data.frame() %>% 
#   write_csv('tables/field_1/field1_posthoc_treated_interaction.csv')


# number and weight GSM -----
field_1 %>%
  ggplot(aes(Number_GSM, Weight_GSM, color = Sampling_time)) +
  geom_point() +
  geom_smooth(method = lm)


## GSM plots
field_1 %>%
  gather(key, value, Number_GSM:Weight_GSM) %>%
  mutate(key= fct_recode(key, `No. individulas` = "Number_GSM", `Weight (g)` = "Weight_GSM")) %>% 
  ggplot(aes(
    x = Acetic_acid,
    y = value ,
    group = Sampling_time,
    fill = Sampling_time
  )) +
  facet_grid(key ~ Time_in_bath,
             scales = "free_x",
             space = "free_x") +
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
  scale_y_log10() +
  scale_fill_grey(name = "Sampling time") +
  labs(x = "Acetic acid concentration", y = "Mean (± SE)")
