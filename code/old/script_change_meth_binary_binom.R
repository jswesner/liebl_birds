library(brms)
library(tidyverse)
library(janitor)
library(RCurl)
library(ggridges)
library(cowplot)
library(viridis)
library(readxl)

binary_change_summary <- read_excel("data/binary change summary.xlsx") %>% 
  as_tibble() %>% clean_names() %>% 
  mutate(total = change + no_change) %>% 
  mutate(total_change = loss + gain) 

binom_brm <- brm(total_change|trials(total) ~ transition*behavior + (1|ind),
                    family = binomial(link = "logit"),
                    data = binary_change_summary,
                    prior = c(prior(normal(0,2), class = "Intercept"),
                              prior(normal(0,1), class = "b"),
                              prior(cauchy(0,1), class = "sd")),
                    sample_prior = T)

saveRDS(binom_brm, file = "models/binom_brm.rds")

binom_brm <- readRDS(file = "models/binom_brm.rds")

binom_brm
pp_check(binom_brm)

# extract draws
binom_draws <- posterior_samples(binom_brm) %>% as_tibble() %>% clean_names() %>% mutate(iter = 1:nrow(.)) %>% 
  select(!contains(c("r_ind", "lp"))) 

# posterior only
binom_post <- binom_draws %>% select(!contains("prior"))

str(binom_post)


binom_post_plot_data <- binom_post %>% 
  mutate(N_FA = b_intercept + b_behavior_n,
                 N_CF = b_intercept + b_transition_hf + b_behavior_n + b_transition_hf_behavior_n,
                 D_FA = b_intercept,
                 D_CF = b_intercept + b_transition_hf) %>% 
  pivot_longer(cols = c(N_CF, N_FA, D_CF, D_FA)) %>% 
  separate(name, c("behavior", "transition")) %>% 
  mutate(y = inv_logit_scaled(value))


binary_change_plot <- binom_post_plot_data %>% 
  ggplot(aes(x = transition, y = y, fill = behavior)) +
  geom_boxplot(aes(group = interaction(transition, behavior)), outlier.shape = NA) +
  geom_point(data = binary_change_summary %>% 
               mutate(transition = case_when(transition == "HF" ~ "CF", TRUE ~ transition)), 
             aes(y = total_change/total), position = position_dodge(width = 0.75),
             color = "grey50") +
  labs(y = "Proportion methylated",
       x = "Transition") +
  scale_color_viridis_d(option = "E") +
  scale_fill_viridis_d(option = "E") +
  theme_classic() +
  theme(text = element_text(size = 14)) +
  NULL 

ggsave(binary_change_plot, file = "plots/binary_change_plot.jpg", dpi = 600, width = 5, height = 4)






# samples of outcomes and wrangle
binom_total <- binary_change_summary %>% mutate(name = case_when(name == "no_change" ~ "nochange", TRUE ~ name),
                                                behavior = str_to_lower(behavior),
                                                transition = str_to_lower(transition),
                                                id = paste0(transition, "_", behavior, "_", name),
                                                change_in_methylation = name) 


# make plot 
# data to plot raw
binom_dat_plot <- binom_total %>% 
  mutate(transition = case_when(transition == "fa" ~ "Fledgling to Adult", TRUE ~ "Hatchling to Fledgling"),
         transition = fct_relevel(transition, "Hatchling to Fledgling"),
         change_in_methylation = case_when(change_in_methylation == "nochange" ~ "none", TRUE ~ change_in_methylation),
         change_in_methylation = fct_relevel(change_in_methylation, "none"),
         behavior = case_when(behavior == "d" ~ "disperse", TRUE ~ "philopatric"))

# plot with posteriors and raw data
transition_behavior_change <- post_preds %>% 
  mutate(transition = case_when(transition == "fa" ~ "Fledgling to Adult", TRUE ~ "Hatchling to Fledgling"),
         transition = fct_relevel(transition, "Hatchling to Fledgling"),
         change_in_methylation = case_when(change_in_methylation == "nochange" ~ "none", TRUE ~ change_in_methylation),
         change_in_methylation = fct_relevel(change_in_methylation, "none"),
         behavior = case_when(behavior == "d" ~ "disperse", TRUE ~ "philopatric")) %>% 
  ggplot(aes(x = change_in_methylation, fill = behavior, y = y)) + 
  geom_boxplot(aes(group = interaction(behavior,change_in_methylation)), outlier.shape = NA) + 
  facet_wrap(~transition) +
  geom_point(data = binom_dat_plot, aes(y = value, shape = behavior, size = behavior), color = "grey50", 
             position = position_dodge(width = 0.75)) +
  scale_fill_viridis_d(option = "E") +
  theme_classic() +
  scale_shape_manual(values = c(1,24)) +
  scale_size_manual(values = c(2.5, 1.5)) +
  theme(axis.title.x = element_blank(),
        text = element_text(size = 14)) +
  labs(y = "Number of loci")

ggsave(transition_behavior_change, file = "plots/transition_behavior_change.jpg", dpi = 600, width = 8, height = 3.5)






















# hf_n_loss <- binom_post %>% select(contains(c("intercept", "transition_hf", "behavior_n", "loss"))) %>% 
#   select(!contains(c("change", "sd_", "iter"))) %>% rowSums(na.rm = T)
# # hf_n_nochange <- binom_post %>% select(contains(c("intercept", "transition_hf", "behavior_n", "change"))) %>% 
# #   select(!contains(c("loss", "sd_", "iter"))) %>% rowSums(na.rm = T)
# hf_n_gain <-  binom_post %>% select(contains(c("intercept","transition_hf", "behavior_n"))) %>% 
#   select(!contains(c("loss", "sd_", "change", "iter"))) %>% rowSums(na.rm = T)
# hf_n_nochange <- 40360 - (exp(hf_n_loss) + exp(hf_n_gain))
# 
# fa_n_loss <- binom_post %>% select(contains(c("intercept", "behavior_n", "loss"))) %>% 
#   select(!contains(c("change", "sd_", "transition_hf", "iter"))) %>% rowSums(na.rm = T)
# # fa_n_nochange <- binom_post %>% select(contains(c("intercept", "behavior_n", "change"))) %>% 
# #   select(!contains(c("loss", "sd_", "iter", "transition_hf", "iter"))) %>% rowSums(na.rm = T)
# fa_n_gain <-  binom_post %>% select(contains(c("intercept", "behavior_n"))) %>% 
#   select(!contains(c("loss", "sd_", "change", "iter", "transition_hf"))) %>% rowSums(na.rm = T)
# 
# fa_n_nochange <- 40360 - (exp(fa_n_loss) + exp(fa_n_gain))
# 
# 
# hf_d_loss <- binom_post %>% select(contains(c("intercept", "transition_hf", "loss"))) %>% 
#   select(!contains(c("change", "sd_", "iter", "behavior_n"))) %>% rowSums(na.rm = T)
# # hf_d_nochange <- binom_post %>% select(contains(c("intercept", "transition_hf", "behavior_n", "change"))) %>% 
# #   select(!contains(c("loss", "sd_", "iter", "behavior_n"))) %>% rowSums(na.rm = T)
# hf_d_gain <-  binom_post %>% select(contains(c("intercept", "transition_hf", "behavior_n"))) %>% 
#   select(!contains(c("loss", "sd_", "change", "iter", "behavior_n"))) %>% rowSums(na.rm = T)
# hf_d_nochange <- 40360 - (exp(hf_d_loss) + exp(hf_d_gain))
# 
# fa_d_loss <- binom_post %>% select(contains(c("intercept", "behavior_n", "loss"))) %>% 
#   select(!contains(c("change", "sd_", "transition_hf", "iter", "behavior_n"))) %>% rowSums(na.rm = T)
# # fa_d_nochange <- binom_post %>% select(contains(c("intercept", "behavior_n", "change"))) %>% 
# #   select(!contains(c("loss", "sd_", "iter", "transition_hf", "iter", "behavior_n"))) %>% rowSums(na.rm = T)
# fa_d_gain <-  binom_post %>% select(contains(c("intercept", "behavior_n"))) %>% 
#   select(!contains(c("loss", "sd_", "change", "iter", "transition_hf", "behavior_n"))) %>% rowSums(na.rm = T)
# fa_d_nochange <- 40360 - (exp(fa_d_loss) + exp(fa_d_gain))
