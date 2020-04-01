library(brms)
library(tidyverse)
library(janitor)
library(RCurl)
library(ggridges)
library(cowplot)
library(viridis)

#import data and make loss and gain the same group
change_two <- read.csv(text = getURL("https://raw.githubusercontent.com/jswesner/liebl_birds/master/change_in_meth.csv")) %>% 
  group_by(age.transition, ind, behavior) %>% 
  mutate(type.change.two = case_when(type.change == "no" ~ "no",
                                     TRUE ~ "changed")) %>% 
  group_by(age.transition, ind, behavior, type.change.two) %>% 
  summarize(changed = sum(number.change)) %>%
  group_by(age.transition, ind, behavior) %>% 
  mutate(total = sum(changed),
         not_changed = total - changed,
         p = changed/total,
         age = fct_relevel(age.transition, "HF"),
         disperse = case_when(behavior == "n" ~ "natal", TRUE ~ "disperse")) %>% 
  filter(type.change.two != "no")


  

binom_brm <- brm(changed|trials(total) ~ age.transition*behavior + (1|ind),
                 family = binomial(link = "logit"),
                 data = change_two,
                 prior = c(prior(normal(0,1), class = "Intercept"),
                           prior(normal(0,1), class = "b"),
                           prior(cauchy(0,1), class = "sd")),
                 sample_prior = T)
                 
#prior predictive               
pp_check(binom_brm, type = "boxplot")

priors <- prior_samples(binom_brm, c("b_Intercept", "b_age.transitionHF", "b_behaviorn", "b_age.transitionHF:behaviorn"))
prior_pred <- priors %>% 
  mutate(prior_ad = b_Intercept,
         prior_an = b_Intercept + b_behaviorn,
         prior_fd = b_Intercept + b_age.transitionHF,
         prior_fn = b_Intercept + b_age.transitionHF + `b_age.transitionHF:behaviorn`,
         iter = 1:nrow(.)) %>% 
  select(prior_ad, prior_an, prior_fd, prior_fn, iter) %>% 
  gather(param, log_odds, -iter) %>% 
  mutate(p = inv_logit_scaled(log_odds))

#plot prior predictive.
prior_pred %>% 
  ggplot(aes(x = p, y = param)) +
  geom_density_ridges()



#posterior distribution
plot(marginal_effects(binom_brm), points = T)

posteriors <- posterior_samples(binom_brm)

posts <- posteriors %>% 
  mutate(posts_ad = b_Intercept,
         posts_an = b_Intercept + b_behaviorn,
         posts_fd = b_Intercept + b_age.transitionHF,
         posts_fn = b_Intercept + b_age.transitionHF + b_behaviorn + `b_age.transitionHF:behaviorn`,
         iter = 1:nrow(.)) %>% 
  select(posts_ad, posts_an, posts_fd, posts_fn, iter) %>% 
  gather(param, log_odds, -iter) %>% 
  mutate(p = inv_logit_scaled(log_odds),
         age = case_when(grepl("posts_a", param) ~ "FA", TRUE ~ "HF"),
         disperse = case_when(grepl("n", param) ~ "natal", TRUE ~ "disperse"),
         age = fct_relevel(age, "HF"))


#density plot
posts %>% 
  ggplot(aes(x = p, y = param)) + 
  geom_density_ridges2()



#prior v post
prior_post <- posts %>% mutate(model = "posterior") %>% 
  bind_rows(prior_pred %>% mutate(model = "prior prediction")) %>% 
  mutate(trt = str_sub(param, -2, -1),
         trt = case_when(trt == "ad" ~ "FA_disperse",
                         trt == "an" ~ "FA_natal",
                         trt == "fd" ~ "HF_disperse",
                         trt == "fn" ~ "HF_natal")) %>% 
  as_tibble()


plot_priorvpost_meth <- prior_post %>% 
  ggplot(aes(x = p, y = trt, fill = model)) + 
  geom_density_ridges2() +
  scale_fill_viridis_d() +
  theme_classic() +
  theme(axis.title.y = element_blank()) +
  labs(x = "Proportion of loci that changed methylation")

ggsave(plot_priorvpost_meth, file = "plot_priorvpost_meth.jpg", dpi = 600, width = 7, height = 5.5)


#make plot
plot_change_meth_two <- posts %>% 
  ggplot(aes(x = age, y = p, fill = disperse)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_point(data = change_two, aes(shape = disperse),
             position = position_dodge(width = 0.5), color = "grey50", size = 1.5) +
  scale_fill_viridis_d(option = "E") +
  ylab("Proportion of loci that changed methylation") +
  xlab("Age Transition") +
  theme_classic() +
  scale_x_discrete(label = c("Hatchling to Fledgling", "Fledgling to Adult")) +
  theme(text = element_text(size = 14)) +
  ylim(0,1) +
  NULL


ggsave(plot_change_meth_two, file = "plot_change_meth_two.jpg", dpi = 600, width = 6, height = 6)




#summary_stats
posts %>% 
  group_by(age, disperse) %>% 
  summarize(mean = mean(p),
            sd = sd(p),
            low95 = quantile(p, probs = 0.025),
            high95 = quantile(p, probs = 0.975),
            prob_p = sum(p>0)/max(iter)) 

posts %>% 
  group_by(age, disperse) %>% 
  select(iter, p, age, disperse) %>% 
  pivot_wider(names_from = disperse, values_from = p) %>% 
  mutate(diff = disperse - natal) %>% 
  group_by(age) %>% 
  summarize(mean = mean(diff),
            sd = sd(diff),
            low95 = quantile(diff, probs = 0.025),
            high95 = quantile(diff, probs = 0.975),
            prob_diff = sum(diff>0)/max(iter)) 
