library(brms)
library(tidyverse)
library(janitor)

library(readxl)
bird_data <- read_excel("bird_data.xlsx")

get_prior(variance ~ age*dispersed + (1|ind), data = bird_data,
          family = Beta(link = "logit"))


m1_prior <- brm(variance ~ age*dispersed + (1|ind), data = bird_data,
          family = Beta(link = "logit"),
          prior = c(prior(normal(0,1), class = "b"),
                    prior(normal(0,1), class = "Intercept"),
                    prior(cauchy(0,1), class = "sd")),
          sample_prior = "only")


m1_prior

m1 <- brm(variance ~ age*dispersed + (1|ind), data = bird_data,
          family = Beta(link = "logit"),
          prior = c(prior(normal(0,1), class = "b"),
                    prior(normal(0,1), class = "Intercept"),
                    prior(cauchy(0,1), class = "sd")))


m1
marginal_effects(m1)
pp_check(m1, type = "boxplot")

marg_m1 <- marginal_effects(m1, effects = "age:dispersed")
marg_m1_df <- as.data.frame(marg_m1$`age:dispersed`)

new_data <- expand.grid(age = unique(bird_data$age), 
                        dispersed = unique(bird_data$dispersed))

renames <- c("A_N","F_N","H_N","A_D","F_D","H_D")

posts <- as.data.frame(fitted(m1, summary = F, newdata = new_data, re_formula = NA))
colnames(posts) <- renames

bird_data_plot <- bird_data %>% 
  mutate(age = fct_relevel(age, "H", "F"))

posts_plot <-  posts %>% as_tibble() %>% 
  gather(trt, variance) %>% 
  separate(trt, c("age","dispersed")) %>% 
  mutate(age = fct_relevel(age, "H", "F"),
         iter = rep(1:4000, 6)) 

ggplot(data = posts_plot, aes(x = age, y = variance, fill = dispersed, group = dispersed)) +
  geom_boxplot(aes(group = interaction(age, dispersed)), outlier.shape = NA) +
  geom_point(data = bird_data_plot, position = position_dodge(width = 0.75)) +
  # geom_line(data = bird_data_plot, aes(group = interaction(ind, dispersed)),
  #           position = position_dodge(width = 0.75)) +
  NULL


bird_data %>% 
  ggplot(aes(x = age, y = variance, group = ind, color = dispersed))+
  geom_line() +
  geom_point()



posts %>% as_tibble() %>% 
  mutate(diff_F = F_D - F_N,
         diff_H = H_D - H_N) %>% 
  summarize(median = median(diff_F),
            low95 = quantile(diff_F, probs = 0.025),
            high95 = quantile(diff_F, probs = 0.975),
            prob_greater_0 = sum(diff_F>0)/4000,
            prob_greater_0H = sum(diff_H>0)/4000)


posts %>% as_tibble() %>% 
  mutate(diff_F = F_D - F_N) 
