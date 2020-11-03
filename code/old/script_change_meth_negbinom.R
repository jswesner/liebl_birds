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
  pivot_longer(cols = c(no_change, gain, loss)) %>% 
  mutate(value = as.integer(value))


binary_dat <- binary_change_summary %>% filter(name!= "no_change")

get_prior(value ~ transition*behavior*name + (1|ind), 
          data = binary_dat,
          family = negbinomial(link = "log", link_shape = "identity"))

negbinom_brm <- brm(value ~ transition*behavior*name + (1|ind),
                    family = negbinomial(link = "log", link_shape = "identity"),
                    data = binary_dat,
                    prior = c(prior(normal(8,2), class = "Intercept"),
                              prior(normal(0,1), class = "b"),
                              prior(cauchy(0,1), class = "sd"),
                              prior(gamma(0.01, 0.01), class = "shape")),
                    sample_prior = T)

saveRDS(negbinom_brm, file = "models/negbinom_brm.rds")

negbinom_brm <- readRDS("models/negbinom_brm.rds")
pp_check(negbinom_brm, type = "boxplot")

# extract draws
binom_draws <- posterior_samples(negbinom_brm) %>% as_tibble() %>% clean_names() %>% mutate(iter = 1:nrow(.)) %>% 
  select(!contains(c("r_ind", "lp"))) 

# posterior only
binom_post <- binom_draws %>% select(!contains("prior"))

str(binom_post)





# samples of outcomes and wrangle

hf_n_loss <- binom_post %>% select(contains(c("intercept", "transition_hf", "behavior_n", "loss"))) %>% 
  select(!contains(c("change", "sd_", "iter"))) %>% rowSums(na.rm = T)
# hf_n_nochange <- binom_post %>% select(contains(c("intercept", "transition_hf", "behavior_n", "change"))) %>% 
#   select(!contains(c("loss", "sd_", "iter"))) %>% rowSums(na.rm = T)
hf_n_gain <-  binom_post %>% select(contains(c("intercept","transition_hf", "behavior_n"))) %>% 
  select(!contains(c("loss", "sd_", "change", "iter"))) %>% rowSums(na.rm = T)
hf_n_nochange <- 40360 - (exp(hf_n_loss) + exp(hf_n_gain))

fa_n_loss <- binom_post %>% select(contains(c("intercept", "behavior_n", "loss"))) %>% 
  select(!contains(c("change", "sd_", "transition_hf", "iter"))) %>% rowSums(na.rm = T)
# fa_n_nochange <- binom_post %>% select(contains(c("intercept", "behavior_n", "change"))) %>% 
#   select(!contains(c("loss", "sd_", "iter", "transition_hf", "iter"))) %>% rowSums(na.rm = T)
fa_n_gain <-  binom_post %>% select(contains(c("intercept", "behavior_n"))) %>% 
  select(!contains(c("loss", "sd_", "change", "iter", "transition_hf"))) %>% rowSums(na.rm = T)

fa_n_nochange <- 40360 - (exp(fa_n_loss) + exp(fa_n_gain))


hf_d_loss <- binom_post %>% select(contains(c("intercept", "transition_hf", "loss"))) %>% 
  select(!contains(c("change", "sd_", "iter", "behavior_n"))) %>% rowSums(na.rm = T)
# hf_d_nochange <- binom_post %>% select(contains(c("intercept", "transition_hf", "behavior_n", "change"))) %>% 
#   select(!contains(c("loss", "sd_", "iter", "behavior_n"))) %>% rowSums(na.rm = T)
hf_d_gain <-  binom_post %>% select(contains(c("intercept", "transition_hf", "behavior_n"))) %>% 
  select(!contains(c("loss", "sd_", "change", "iter", "behavior_n"))) %>% rowSums(na.rm = T)
hf_d_nochange <- 40360 - (exp(hf_d_loss) + exp(hf_d_gain))

fa_d_loss <- binom_post %>% select(contains(c("intercept", "behavior_n", "loss"))) %>% 
  select(!contains(c("change", "sd_", "transition_hf", "iter", "behavior_n"))) %>% rowSums(na.rm = T)
# fa_d_nochange <- binom_post %>% select(contains(c("intercept", "behavior_n", "change"))) %>% 
#   select(!contains(c("loss", "sd_", "iter", "transition_hf", "iter", "behavior_n"))) %>% rowSums(na.rm = T)
fa_d_gain <-  binom_post %>% select(contains(c("intercept", "behavior_n"))) %>% 
  select(!contains(c("loss", "sd_", "change", "iter", "transition_hf", "behavior_n"))) %>% rowSums(na.rm = T)
fa_d_nochange <- 40360 - (exp(fa_d_loss) + exp(fa_d_gain))

nochanges <- bind_cols(hf_n_nochange = hf_n_nochange,
                       hf_d_nochange = hf_d_nochange,
                       fa_n_nochange = fa_n_nochange,
                       fa_d_nochange = fa_n_nochange) %>% 
  mutate(iter = 1:nrow(.)) %>% 
  pivot_longer(cols = -iter) %>% 
  separate(name, c("transition", "behavior", "change_in_methylation")) %>% 
  mutate(y = value,
         id = paste0(transition, "_", behavior, "_", change_in_methylation))


post_preds <- bind_cols(hf_n_loss = hf_n_loss, 
          hf_n_gain = hf_n_gain,
          fa_n_loss = fa_n_loss,
          fa_n_gain = fa_n_gain,
          hf_d_loss = hf_d_loss,
          hf_d_gain = hf_d_gain,
          fa_d_loss = fa_d_loss,
          fa_d_gain = fa_d_gain) %>% mutate(iter = 1:nrow(.)) %>% 
  sample_n(1000) %>% 
  pivot_longer(cols = -iter) %>% 
  separate(name, c("transition", "behavior", "change_in_methylation")) %>% 
  mutate(y = exp(value),
         id = paste0(transition, "_", behavior, "_", change_in_methylation)) %>% 
  bind_rows(nochanges)

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

