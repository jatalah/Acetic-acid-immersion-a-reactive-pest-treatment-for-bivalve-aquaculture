library(tidyverse)
library(readxl)
library(ggfortify)
library(broom)
library(lsmeans)
library(emmeans)
library(ggpmisc)
source('theme_javier.R')
theme_set(theme_javier())

# Lab experiment 1----------------
lab_1 <- 
  read_excel('data/Laboratory 1.xlsx') %>% 
  gather(size, survival, `Survival_2-5mm _spat`:`Survival_12-15mm _spat`) %>% 
  mutate(size = str_remove(size,"Survival_"),
         survival_raw = survival,
         survival=(survival/20)*100,
         size = str_remove(size," _spat"),
         size = fct_relevel(size, "2-5mm","8-10mm"),
         size = fct_recode(size, "2 - 5 mm" = "2-5mm","8 - 10 mm"= "8-10mm", "12 - 15 mm"="12-15mm"),
         aa = Acetic_acid,
         Acetic_acid = fct_recode(
           factor(Acetic_acid),
           Control = "0",
           `1 %` = "1",
           `2 %` = "2",
           `4 %` = "4",
           `8 %` = "8"
         ),
         Time_in_bath = fct_recode(
           factor(Time_in_bath),
           `10 sec` = "10",
           `30 sec` = "30",
           `60 sec` = "60"
         ))
# Bar plot survival-----------
lab1_plot_surv <- 
  ggplot(lab_1,
       aes(
         x = Acetic_acid,
         y = survival,
         group = Time_in_bath,
         fill = Time_in_bath
         )) +
  facet_wrap(~ size, scales = 'free_x') +
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
  scale_fill_grey(name = "Bath time") +
  labs(x= "Acetic acid", y = 'Survival (%)')
print(lab1_plot_surv)

# save plot bm plot -----
ggsave(
  plot = lab1_plot_surv,
  filename = "figures/lab_1/lab1_plot_surv.svg",
  width = 8,
  height = 2.5
)

### Pool Time in bath
lab1_plot_surv_pooled <- 
  ggplot(lab_1,
       aes(
         x = Acetic_acid,
         y = survival
       )) +
  # facet_grid(Time_in_bath~ size, scales = 'free_x') +
  facet_wrap(~ size, scales = 'free_x') +
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
  scale_fill_grey(name = "Bath time") +
  labs(x= "Acetic acid", y = 'Survival (%)') +
  geom_hline(yintercept = 75, lty = 2)

ggsave(
  plot = lab1_plot_surv_pooled,
  filename = "figures/lab_1/lab1_plot_surv_pooled.svg",
  width = 8,
  height = 2.5
)

# ANOVA survival------------
m1 <- lm(survival ~ Time_in_bath * size * Acetic_acid, data  = lab_1)
anova(m1)

autoplot(m1)
write_csv(tidy(anova(m1)), 'tables/lab_1/ANOVA_lab1_surv.csv', na = "") %>% print()

# Post-hoc tests---------
posthoc_m1_interaction <-
  lsmeans(m1, pairwise ~ size * Acetic_acid, adjust = "tukey") %>%
  .$contrast %>%
  data.frame() %>%
  write_csv('tables/lab_1/lab1_survival_posthoc_interaction.csv')

posthoc_m1_time <- 
  lsmeans(m1, pairwise ~ Time_in_bath, adjust = "tukey") %>% 
  .$contrasts %>% 
  data.frame() %>% 
  write_csv('tables/lab_1/lab1_survival_posthoc_time.csv')

# Regression approach --------
lab1_regression_plot <- 
  ggplot(lab_1,
       aes(
         x = aa,
         y = survival,
         color = Time_in_bath,
         group = Time_in_bath
       )) +
  facet_wrap(~ size) +
  geom_point(position = position_jitter(0.2), alpha = .7) +
  stat_smooth(method = "lm", alpha = .1) +
  scale_y_continuous(limits = c(0,110)) +
  scale_color_discrete(name = "Time in bath") +
  labs(x = 'Acetic acid concentration (%)', y = "Viability (%)") +
  stat_poly_eq(aes(label =  stat(rr.label)),
               formula = y~x, parse = TRUE, label.x = "left",label.y ="bottom")
lab1_regression_plo


# save plot bm plot -----
ggsave(
  plot = lab1_regression_plot,
  filename = "figures/lab_1/lab1_regression_plot.svg",
  width = 8,
  height = 2.5
)

m1_reg <- lm(survival ~ size * aa, data  = lab_1)
summary(m1_reg)
tidy(summary(m1_reg)) %>% 
  write_csv('tables/lab_1/lab1_regression_coef_table.csv')
tidy(anova(m1_reg)) %>% 
  write_csv('tables/lab_1/lab1_regression_summary_table.csv')

# mean +-se plots-----
mean_se_plot <- 
  ggplot(lab_1,
       aes(
         x = aa,
         y = survival,
         color = Time_in_bath,
         group = Time_in_bath,
         shape = Time_in_bath
       )) +
  
  stat_summary(fun = mean, geom = 'line', aes(lty = Time_in_bath), color = 'gray20') +
  stat_summary(fun.data = "mean_se",
               size = .5, color = 'gray20') +
  theme_javier() +
  scale_linetype(name = "Time in bath") +
  scale_shape_manual(values = c(1,0,16), name = "Time in bath") +
  labs(x = 'Acetic acid concentration (%)', y = "Viability (%)") + 
  facet_wrap(~ size)

ggsave(
  plot = mean_se_plot,
  filename = "figures/lab_1/lab1_mean_se_plot.svg",
  width = 12,
  height = 4
)



# Time in bath handling effect---
ggplot(lab_1,
       aes(x = Time_in_bath,
           y = survival)) +
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
  scale_fill_grey(name = "Bath time") +
  labs(x = "Acetic acid", y = 'Survival (%)')


lab_1 %>%
  group_by(Time_in_bath) %>%
  summarise_at(vars(survival), .funs = list(mean = mean, se = se))
