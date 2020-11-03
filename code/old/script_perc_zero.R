library(brms)
library(tidyverse)
library(janitor)
library(RCurl)
library(ggridges)
library(viridis)
library(readxl)


# Load data and models from github ---------------------------------------------------
#raw data
# bird_olddata <- read.csv(text = getURL("https://raw.githubusercontent.com/jswesner/liebl_birds/master/bird_data.csv?token=AETBVHH4E7BIPRHDOWFNAUK6AFUVQ"))


bird_data <- read_excel("data/binary total methylation dispersal.xlsx")

# fit model
# m1_perc <- brm(methylated.perc ~ age*behavior + (1|ind), data = bird_data,
#           family = Beta(link = "logit"),
#           prior = c(prior(normal(0,1), class = "b"),
#                     prior(normal(0,2), class = "Intercept"),
#                     prior(cauchy(0,1), class = "sd")),
#           cores = 2, 
#           sample_prior = T)


m1_perc <- readRDS("models/m1_perc.rds")

#run again with tighter prior on beta
# m1_perc_t <- update(m1_perc, prior = c(prior(normal(0, 0.1), class = "b"),
#                                        prior(normal(0, 2), class = "Intercept"),
#                                        prior(cauchy(0,1), class = "sd")))
# 
# m1_perc_t2 <- update(m1_perc, prior = c(prior(normal(0, 0.3), class = "b"),
#                                        prior(normal(0, 2), class = "Intercept"),
#                                        prior(cauchy(0,1), class = "sd")))
# 
# m1_perc_t3 <- update(m1_perc, prior = c(prior(normal(0, 0.5), class = "b"),
#                                        prior(normal(0, 2), class = "Intercept"),
#                                        prior(cauchy(0,1), class = "sd")))
# 
# m1_perc_t4 <- update(m1_perc, prior = c(prior(normal(0, 0.7), class = "b"),
#                                        prior(normal(0, 2), class = "Intercept"),
#                                        prior(cauchy(0,1), class = "sd")))
# 
# m1_perc_t5 <- update(m1_perc, prior = c(prior(normal(0, 2), class = "b"),
#                                         prior(normal(0, 2), class = "Intercept"),
#                                         prior(cauchy(0,1), class = "sd")))
# 
# m1_perc_t6 <- update(m1_perc, prior = c(prior(normal(0, 3), class = "b"),
#                                         prior(normal(0, 2), class = "Intercept"),
#                                         prior(cauchy(0,1), class = "sd")))
# 
# saveRDS(m1_perc, file = "models/m1_perc.rds")
# saveRDS(m1_perc_t1, file = "models/m1_perc_t1.rds")
# saveRDS(m1_perc_t2, file = "models/m1_perc_t2.rds")
# saveRDS(m1_perc_t3, file = "models/m1_perc_t3.rds")
# saveRDS(m1_perc_t4, file = "models/m1_perc_t4.rds")
# saveRDS(m1_perc_t5, file = "models/m1_perc_t5.rds")
# saveRDS(m1_perc_t6, file = "models/m1_perc_t6.rds")

m1_perc_t1 <- readRDS("models/m1_perc_t1.rds")
m1_perc_t2 <- readRDS("models/m1_perc_t2.rds")
m1_perc_t3 <- readRDS("models/m1_perc_t3.rds")
m1_perc_t4 <- readRDS("models/m1_perc_t4.rds")
m1_perc_t5 <- readRDS("models/m1_perc_t5.rds")
m1_perc_t6 <- readRDS("models/m1_perc_t6.rds")


#prior predictive check (can the model generate dataset that resemble the original dataset?)
post_check_perc <- pp_check(m1_perc, type = "boxplot") #Yes. Yes it can. On the right track.
ggsave(post_check_perc, file = "plots/post_check_perc.jpg", dpi = 600, width = 6, height = 4)

# saveRDS(m1_perc, file = 'models/m1_perc.rds')

#extract posteriors to plot and summarize
marg_m1_perc <- conditional_effects(m1_perc, effects = "age:behavior")
marg_m1__percdf <- as.data.frame(marg_m1_perc$`age:behavior`)

new_data <- expand.grid(age = unique(bird_data$age), 
                        behavior = unique(bird_data$behavior))

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


# Quantitative summaries of the posterior ----------------------------------------------
#summary stats of treatments
summary_stats <- posts_plot_perc %>% 
  group_by(age, dispersed) %>% 
  summarize(mean = mean(methylated.perc),
            sd = sd(methylated.perc),
            low95 = quantile(methylated.perc, probs = 0.025),
            high95 = quantile(methylated.perc, probs = 0.975)) %>% 
  ungroup() %>% 
  mutate_if(is.numeric, round, 2)

write.csv(summary_stats, file = "data/summary_stats.csv", row.names = F)

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





