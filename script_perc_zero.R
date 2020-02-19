library(brms)
library(tidyverse)
library(janitor)
library(RCurl)
library(ggridges)


# Load data and models from github ---------------------------------------------------
#raw data
bird_data <- read.csv(text = getURL("https://raw.githubusercontent.com/jswesner/liebl_birds/master/bird_data.csv?token=AETBVHH4E7BIPRHDOWFNAUK6AFUVQ"))

#prior predictive model
m1_prior <- readRDS(url("https://github.com/jswesner/liebl_birds/blob/master/m1_perc_prior_predictive.rds?raw=true"))

#final bayesian model
m1_perc <- readRDS(url("https://github.com/jswesner/liebl_birds/blob/master/m1_perc.rds?raw=true"))

# Run models --------------------------------------------------------------

#Prior predictive (model using only the priors)
# m1_prior <- brm(perc.zero ~ age*dispersed + (1|ind), data = bird_data,
#                 family = Beta(link = "logit"),
#                 prior = c(prior(normal(0,1), class = "b"),
#                           prior(normal(0,2), class = "Intercept"),
#                           prior(cauchy(0,1), class = "sd")),
#                 sample_prior = "only")

#plot prior predictive marginal distributions
marginal_effects(m1_prior, effects = "age:dispersed")

saveRDS(m1_prior, file = 'm1_perc_prior_predictive.rds')

#data to condition on
new_data <- expand.grid(age = unique(bird_data$age), 
                        dispersed = unique(bird_data$dispersed))
#new_names of columns
renames <- new_data %>% unite(names,sep = "_") %>% pull(names)

#posterior of each treatment combo and then renaming columns
prior_fit <- as.data.frame(fitted(m1_prior, summary = F, newdata = new_data, re_formula = NA))
colnames(prior_fit) <- renames

#plot posteriors of each combo to see distribution
prior_fit %>% 
  gather() %>% 
  separate(key, c("age","dispersed")) %>% 
  ggplot(aes(x = value, fill = dispersed, y = age)) +
  geom_density_ridges()

#plot above shows some peaks at 1 and 0, but also a resonable spread across the 0-1 spectrum for each combo. 
#Will stick with these priors



#Full Bayesian model. Same as above but with data added
# m1_perc <- brm(perc.zero ~ age*dispersed + (1|ind), data = bird_data,
#           family = Beta(link = "logit"),
#           prior = c(prior(normal(0,1), class = "b"),
#                     prior(normal(0,2), class = "Intercept"),
#                     prior(cauchy(0,1), class = "sd")),
#           cores = 2)


m1_perc
#prior predictive check (can the model generate dataset that resemble the original dataset?)
pp_check(m1_perc, type = "boxplot") #Yes. Yes it can. On the right track.

saveRDS(m1_perc, file = 'm1_perc.rds')

#extract posteriors to plot and summarize
marg_m1_perc <- marginal_effects(m1_perc, effects = "age:dispersed")
marg_m1__percdf <- as.data.frame(marg_m1_perc$`age:dispersed`)

new_data <- expand.grid(age = unique(bird_data$age), 
                        dispersed = unique(bird_data$dispersed))

renames <- new_data %>% unite(names,sep = "_") %>% pull(names)

posts_perc <- as.data.frame(fitted(m1_perc, summary = F, newdata = new_data, re_formula = NA))
colnames(posts_perc) <- renames



# Make plots --------------------------------------------------------------
#raw data to plot
bird_data_plot <- bird_data %>% 
  mutate(age = fct_relevel(age, "H", "F"),
         age_full = case_when(age == "H" ~ "Hatchling",
                              age == "F" ~ "Fledgling",
                              age == "A" ~ "Adult"),
         age_full = fct_relevel(age_full, "Hatchling","Fledgling"),
         offset = case_when(dispersed == "D" ~ -0.1,
                            TRUE ~ 0.1))

#posterior data to plot
posts_plot_perc <-  as_tibble(posts_perc) %>% 
  mutate(iter = 1:nrow(posts_perc)) %>% 
  filter(iter < 1001) %>% 
  gather(trt, perc.zero, -iter) %>% 
  separate(trt, c("age","dispersed")) %>% 
  mutate(age = fct_relevel(age, "H", "F"),
         age_full = case_when(age == "H" ~ "Hatchling",
                              age == "F" ~ "Fledgling",
                              age == "A" ~ "Adult"),
         age_full = fct_relevel(age_full, "Hatchling","Fledgling"),
         offset = case_when(dispersed == "D" ~ -0.1,
                            TRUE ~ 0.1))




plot_box_lines <- ggplot(data = posts_plot_perc, aes(x = as.numeric(age_full) + offset, y = perc.zero, 
                                                     group = dispersed,
                                                     color = dispersed)) +
  geom_line(aes(group = interaction(iter,dispersed)),
            alpha = 0.1) +
  geom_boxplot(aes(group = interaction(dispersed,age)),
               outlier.shape = NA, width = 0.1) +
  
  scale_fill_brewer(type = "qual") +
  scale_color_brewer(type = "qual") +
  geom_point(data = bird_data_plot, position = position_dodge(width = 0),
             aes(shape = dispersed)) +
  ylab("Proportion methylated") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        text = element_text(size = 14)) +
  scale_x_continuous(breaks=c(1, 2, 3),
                   labels=c("Hatchling", "Fledgling", "Adult")) +
  NULL 

ggsave(plot_box_lines, file = "plot_box_lines.jpg", dpi = 600, width = 5, height = 4)


# Quantitative summaries of the posterior ----------------------------------------------
#summary stats of treatments
summary_stats <- posts_plot_perc %>% 
  group_by(age, dispersed) %>% 
  summarize(mean = mean(perc.zero),
            sd = sd(perc.zero),
            low95 = quantile(perc.zero, probs = 0.025),
            high95 = quantile(perc.zero, probs = 0.975)) %>% 
  ungroup() %>% 
  mutate_if(is.numeric, round, 2)

write.csv(summary_stats, file = "summary_stats.csv")

#summary of differences
as_tibble(posts_perc)  %>% 
  mutate(diff_F = F_D - F_N,
         diff_H = H_D - H_N,
         diff_A = A_D - A_N) %>% 
  select(diff_F, diff_H, diff_A) %>% 
  gather() %>% 
  group_by(key) %>% 
  summarize(median = median(value),
            sd = sd(value),
            low95 = quantile(value, probs = 0.025),
            high95 = quantile(value, probs = 0.975),
            prob_lessthan_0 = sum(value<0)/4000)



#boxplot
# plot_boxplot <- ggplot(data = posts_plot_perc, aes(x = age_full, y = perc.zero, 
#                                                     fill = dispersed, group = dispersed)) +
#   geom_boxplot(aes(group = interaction(age, dispersed)), outlier.shape = NA,
#                width = 0.4) +
#   scale_fill_grey(start = 0.5, end = 1) +
#   geom_point(data = bird_data_plot, position = position_dodge(width = 0.4)) +
#   ylab("Proportion methylated") +
#   theme_classic() +
#   theme(axis.title.x = element_blank(),
#         text = element_text(size = 14)) +
#   NULL
# 
# ggsave(plot_boxplot, file = "plot_boxplot.jpg", dpi = 600, width = 5, height = 4)
# 
# 
# plot_violin <- ggplot(data = posts_plot_perc, aes(x = age_full, y = perc.zero, 
#                                                    fill = dispersed, group = dispersed)) +
#   geom_violin(aes(group = interaction(age, dispersed)), 
#                width = 0.4) +
#   scale_fill_grey(start = 0.5, end = 1) +
#   geom_point(data = bird_data_plot, position = position_dodge(width = 0.4)) +
#   ylab("Proportion methylated") +
#   theme_classic() +
#   theme(axis.title.x = element_blank(),
#         text = element_text(size = 14)) +
#   NULL
# 
# 
# ggsave(plot_violin, file = "plot_violin.jpg", dpi = 600, width = 5, height = 4)
# 
# 
# plot_violin_lines <- ggplot(data = posts_plot_perc, aes(x = age_full, y = perc.zero, 
#                                                   group = dispersed,
#                                                   color = dispersed)) +
#   geom_violin(aes(group = interaction(dispersed,age)),
#               position = position_dodge(width = 0)) +
#   geom_line(aes(group = interaction(dispersed, iter)),
#             alpha = 0.01) +
#   scale_fill_brewer(type = "qual") +
#   scale_color_brewer(type = "qual") +
#   geom_point(data = bird_data_plot, position = position_dodge(width = 0),
#              aes(shape = dispersed)) +
#   ylab("Proportion methylated") +
#   theme_classic() +
#   theme(axis.title.x = element_blank(),
#         text = element_text(size = 14)) +
#   NULL
# 
# ggsave(plot_violin_lines, file = "plot_violin_lines.jpg", dpi = 600, width = 5, height = 4)




