library(brms)
library(tidyverse)
library(janitor)
library(RCurl)
library(ggridges)

# Load data and models from github ---------------------------------------------------
#raw data
change <- read.csv(text = getURL("https://raw.githubusercontent.com/jswesner/liebl_birds/master/change_in_meth.csv"))

#prior predictive model
m2_change_prior <- readRDS(url("https://github.com/jswesner/liebl_birds/blob/master/m2_change_prior.rds?raw=true"))

#final bayesian model
m2_change_full <- readRDS(url("https://github.com/jswesner/liebl_birds/blob/master/m2_change_full.rds?raw=true"))


# Run models --------------------------------------------------------------

# Prior predictive (model using only the priors)

# m2_change_prior <- brm(number.change ~ type.change*behavior*age.transition + (1|ind),
#                  data = change,
#                  family = Gamma(link = "log"),
#                  prior = c(prior(normal(0,1), class = "b"),
#                            prior(normal(0,6), class = "Intercept"),
#                            prior(cauchy(0,1), class = "sd")),
#                  sample_prior = "only")
# 
# saveRDS(m2_change_prior, file = "m2_change_prior.rds")
marginal_effects(m2_change_prior)





#Full Bayesian model. Same as above but with data added
# m2_change_full <- brm(number.change ~ type.change*behavior*age.transition + (1|ind), data = change,
#                       family = Gamma(link = "log"),
#                       prior = c(prior(normal(0,1), class = "b"),
#                                 prior(normal(0,6), class = "Intercept"),
#                                 prior(cauchy(0,1), class = "sd")))
# saveRDS(m2_change_full, file = "m2_change_full.rds")
print(m2_change_full)
pp_check(m2_change_full, type = "boxplot")
# saveRDS(m2_change_full, file = 'm2_change_full.rds')

#extract posteriors to plot and summarize
conditions = data.frame(behavior = c("d","n"))
marg_m2_change_full <- marginal_effects(m2_change_full, effects = "type.change:age.transition", conditions = conditions)
marg_m2_changedf <- as.data.frame(marg_m2_change_full$`type.change:age.transition`)

new_data_change <- expand.grid(age.transition = unique(change$age.transition),
                        type.change = unique(change$type.change), 
                        behavior = unique(change$behavior)) 

renames_change <- new_data_change %>% unite(names,sep = "_") %>% pull(names)

posts_change <- as.data.frame(fitted(m2_change_full,newdata = new_data_change, summary = F,
                                     re_formula = NA))
colnames(posts_change) <- renames_change



plot_change_data <- posts_change %>% 
  mutate(iter = 1:4000) %>% 
  gather(key, value, -iter) %>% 
  separate(key, c("age.transition", "type.change", "behavior"), remove = F)



raw_data <- change %>% 
  mutate(age.transition = case_when(age.transition == "HF" ~  "Hatchling to Fledgling", 
                                    TRUE ~ "Fledgling to Adult"),
         age.transition = fct_relevel(age.transition, "Hatchling to Fledgling"),
         behavior = case_when(behavior == "d" ~ "disperse",
                              TRUE ~ "natal"),
         type.change = str_replace(type.change, "no", "none"),
         type.change = str_replace(type.change, "lose", "loss"),
         type.change = fct_relevel(type.change, "none")) %>% 
  unite(Group, c(type.change, behavior), sep = "_", remove = F)


plot_model <- plot_change_data %>% 
  mutate(age.transition = case_when(age.transition == "HF" ~  "Hatchling to Fledgling", 
                                    TRUE ~ "Fledgling to Adult"),
         age.transition = fct_relevel(age.transition, "Hatchling to Fledgling"),
         behavior = case_when(behavior == "d" ~ "disperse",
                              TRUE ~ "natal"),
         type.change = str_replace(type.change, "no", "none"),
         type.change = str_replace(type.change, "lose", "loss"),
         type.change = fct_relevel(type.change, "none")) %>% 
  unite(Group, c(type.change, behavior), sep = "_", remove = F)
  

plot_change_meth <- ggplot() +
  geom_boxplot(data = plot_model, 
               aes(x = type.change, y =value, fill = behavior),
               position = position_dodge(width = 0.5), width = 0.5,
               outlier.shape = NA) +
  facet_wrap(~age.transition) +
  scale_fill_viridis_d(option = "E") +
   # scale_y_log10() +
  geom_point(data = raw_data, position = position_dodge(width = 0.5), 
             aes(x = type.change, y = number.change, shape = behavior),color = "grey40", size = 1.5) +
  ylim(0,150000) +
  ylab("Number of genes") +
  xlab("Change in methylation") +
  theme_classic() +
  theme(text = element_text(size = 14)) +
  NULL

plot_change_meth


ggsave(plot_change_meth, file = "plot_change_meth.jpg", dpi = 600, width = 7, height = 4)
# ggplot_build(plot_change_meth)$data[[1]]$fill

summary_stat_change <- plot_change_data %>% 
  as_tibble() %>% 
  unite(key, c("age.transition", "type.change")) %>% 
  select(iter, key, value, behavior) %>% 
  mutate(key = fct_relevel(key, "HF_no","HF_gain", "HF_lose", "FA_no")) %>% 
  pivot_wider(names_from = behavior, values_from = value) %>% 
  mutate(diff = d-n) %>% 
  group_by(key) %>% 
  summarize(mean = mean(diff),
            sd = sd(diff),
            low95 = quantile(diff, probs = 0.025),
            high95 = quantile(diff, probs = 0.975),
            prob_diff = sum(diff>0)/max(iter)) 

write.csv(summary_stat_change, file = "summary_stat_change.csv", row.names = F)


