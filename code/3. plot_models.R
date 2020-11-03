library(brms)
library(tidyverse)
library(janitor)
library(RCurl)
library(ggridges)
library(cowplot)
library(viridis)
library(readxl)


# Age x Behavior ---------------------------------------------------

bird_data <- read_excel("data/binary total methylation dispersal.xlsx")

m1_perc <- readRDS(file = "models/m1_perc.rds")

marg_m1_perc <- conditional_effects(m1_perc, effects = "age:behavior")
marg_m1__percdf <- as.data.frame(marg_m1_perc$`age:behavior`)

new_data <- expand.grid(age = unique(bird_data$age), 
                        behavior = unique(bird_data$behavior))

renames <- new_data %>% unite(names,sep = "_") %>% pull(names)

posts_perc <- as.data.frame(fitted(m1_perc, summary = F, newdata = new_data, re_formula = NA))
colnames(posts_perc) <- renames


bird_data_plot <- bird_data %>% 
  mutate(age = fct_relevel(age, "H", "F"),
         age_full = case_when(age == "H" ~ "Hatchling",
                              age == "F" ~ "Fledgling",
                              age == "A" ~ "Adult"),
         age_full = fct_relevel(age_full, "Hatchling","Fledgling"),
         offset = case_when(behavior == "D" ~ -0.1,
                            TRUE ~ 0.1),
         behavior = case_when(behavior == "D" ~ "disperse",
                              TRUE ~ "philopatric"))

#posterior data to plot
posts_plot_perc <-  as_tibble(posts_perc) %>% 
  mutate(iter = 1:nrow(posts_perc)) %>% 
  filter(iter < 1001) %>% 
  gather(trt, methylated.perc, -iter) %>% 
  separate(trt, c("age","dispersed")) %>% 
  mutate(age = fct_relevel(age, "H", "F"),
         age_full = case_when(age == "H" ~ "Hatchling",
                              age == "F" ~ "Fledgling",
                              age == "A" ~ "Adult"),
         age_full = fct_relevel(age_full, "Hatchling","Fledgling"),
         offset = case_when(dispersed == "D" ~ -0.1,
                            TRUE ~ 0.1),
         behavior = case_when(dispersed == "D" ~ "disperse",
                              TRUE ~ "philopatric"))



#combine raw and posterior data in a single plot
plot_box_lines <- ggplot(data = posts_plot_perc, aes(x = as.numeric(age_full) + offset, y = methylated.perc, 
                                                     group = behavior)) +
  geom_line(aes(group = interaction(iter,behavior),
                color = behavior),
            alpha = 0.1) +
  geom_boxplot(aes(group = interaction(behavior,age),
                   fill = behavior),
               outlier.shape = NA, width = 0.1) +
  scale_color_viridis_d(option = "E") +
  scale_fill_viridis_d(option = "E") +
  geom_point(data = bird_data_plot, position = position_dodge(width = 0),
             aes(shape = behavior), color = "grey40", size = 1.5) +
  ylab("Proportion methylated") +
  # scale_shape_manual(values = c(21,24)) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        text = element_text(size = 14)) +
  scale_x_continuous(breaks=c(1, 2, 3),
                     labels=c("Hatchling", "Fledgling", "Adult")) +
  NULL 

plot_box_lines

ggsave(plot_box_lines, file = "plots/plot_box_lines.jpg", dpi = 600, width = 6, height = 4)


# Age.transition x Behavior ---------------------------------------------------

change_in_methylation_over_time_cont_data <- read_excel("data/change in methylation over time- cont data.xlsx") %>% 
  mutate(total_loci = total.change + no.change) %>% 
  mutate(age.transition = case_when(age.transition == "CF" ~ "HF", TRUE ~ age.transition))

change_cont_brm <- readRDS(file = "models/change_cont_brm.rds")

post_check_cont <- pp_check(change_cont_brm, type = "boxplot")
ggsave(post_check_cont, file = "plots/post_check_cont.jpg", width = 6, height = 4, dpi = 600)

post_cont <- posterior_samples(change_cont_brm) %>% as_tibble() %>% clean_names() %>% mutate(iter = 1:nrow(.)) %>% 
  mutate(N_HF = b_intercept + b_behavior_n,
         N_FA = b_intercept + b_age_transition_fa + b_behavior_n + b_age_transition_fa_behavior_n,
         D_HF = b_intercept,
         D_FA = b_intercept + b_age_transition_fa)

post_cont_plotdata <- post_cont %>% 
  select(N_HF, N_FA, D_HF, D_FA) %>% 
  mutate(iter = 1:nrow(.)) %>% 
  pivot_longer(cols = c(N_HF, N_FA, D_HF, D_FA)) %>% 
  separate(name, c("behavior", "transition")) %>% 
  mutate(y = inv_logit_scaled(value),
         transition = fct_relevel(transition, "HF"))


cont_plot <- post_cont_plotdata %>% 
  ggplot(aes(x = transition, y = y, fill = behavior)) +
  geom_boxplot(aes(group = interaction(transition, behavior)), outlier.shape = NA) +
  geom_point(data = change_in_methylation_over_time_cont_data, 
             aes(y = total.change/total_loci, x = age.transition,
                 group = behavior),
             position = position_dodge(width = 0.75),
             color = "grey50") +
  labs(y = "Proportion of loci that changed",
       x = "Transition")+
  scale_color_viridis_d(option = "E") +
  scale_fill_viridis_d(option = "E") +
  theme_classic() +
  theme(text = element_text(size = 14)) +
  NULL 

cont_plot


ggsave(cont_plot, file = "plots/cont_plot.jpg", dpi = 600, width = 5, height = 4)


# summarize posteriors

post_conts <- post_cont_plotdata %>% 
  group_by(behavior, transition) %>% 
  summarize(mean = mean(y),
            median = median(y),
            sd = sd(y),
            low95 = quantile(y, probs = 0.975),
            hihg95 = quantile(y, probs = 0.025))

write.csv(post_conts, file = "data/post_conts.csv", row.names = F)

post_cont_plotdata %>% 
  select(-value) %>% 
  pivot_wider(names_from = behavior, values_from = y) %>% 
  mutate(diff_n_d = N-D) %>% 
  group_by(transition) %>% 
  summarize(mean = mean(diff_n_d),
         sd = sd(diff_n_d),
         prob_diff = sum(diff_n_d>0)/4000)


# Sex ---------------------------------------------------------------------


meth_sex.brm <- readRDS("models/meth_sex.brm")
saveRDS(meth_sex.brm, file = "models/meth_sex.brm")


pp_check(meth_sex.brm, type = "boxplot")


meth_sex_post <- posterior_samples(meth_sex.brm) %>% as_tibble() %>% 
  mutate(F = inv_logit_scaled(b_Intercept),
         M = inv_logit_scaled(b_Intercept + b_sexM))



meth_sex_plot <- meth_sex_post %>% 
  pivot_longer(cols = c(F, M), names_to = "sex", values_to = "methylated_perc") %>% 
  ggplot(aes(x = sex, y = methylated_perc)) + 
  geom_boxplot(aes(group = sex), outlier.shape = NA) +
  geom_point(data = binary_total_methylation_sex, position = position_jitter(width = 0.03), alpha = 0.4) +
  coord_cartesian(ylim = c(0,1)) +
  theme_classic() +
  labs(y = "Percent Methylated")

meth_sex_plot

ggsave(meth_sex_plot, file = "plots/meth_sex_plot.jpg" , dpi = 600, width = 4, height = 4)


# summary stats

post_meth_sex <- posterior_samples(meth_sex.brm) %>% as_tibble() %>% clean_names()

post_meth_sex_diff <- post_meth_sex %>% 
  select(b_intercept, b_sex_m) %>% 
  mutate(females = inv_logit_scaled(b_intercept),
         males = inv_logit_scaled(b_intercept + b_sex_m)) %>% 
  select(males, females) %>% 
  mutate(diff = males-females) %>% 
  pivot_longer(cols = everything())

post_meth_sex_diff %>% 
  group_by(name) %>% 
  summarize(mean = mean(value),
            sd = sd(value),
            low = quantile(value, probs = 0.025),
            high = quantile(value, probs = 0.975))



