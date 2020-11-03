library(lme4)
library(glmmTMB)
library(lsmeans)
library(cowplot)


# Analysis for Figure 1 ------------

#run mle of perc.zero models
test <- glmmTMB(methylated.perc ~ age*behavior + (1|ind), data = bird_data,
              family = beta_family(link = "logit"))

test_means <- lsmeans(test, specs = ~age|behavior) %>% rbind() %>% as_tibble() %>% select(-SE,-df) %>% 
  pivot_longer(cols = c(lsmean, lower.CL, upper.CL)) %>% 
  mutate(value = inv_logit_scaled(value)) %>% 
  pivot_wider(names_from = name, values_from = value) %>% 
  mutate(Prior = "Inf (Max LIke)",
         mle_bayes = "Maximum Likelihood",
         model = "7") %>% 
  rename(estimate = lsmean,
         lower = lower.CL,
         upper = upper.CL)

test_sum <- summary(test)

confint(test) %>% as_tibble() %>% slice(1:6)

coefs <- c("intercept", "ageF", "ageH", "behaviorN", "ageF:behaviorN","ageH:behaviorN")

m1_perc_1 <- fixef(m1_perc) %>% as_tibble() %>% mutate(sd = "1 (model in MS)", coef = coefs)
m1_perc_0.1 <- fixef(m1_perc_t1) %>% as_tibble() %>% mutate(sd = "0.1", coef = coefs)
m1_perc_0.3 <- fixef(m1_perc_t2) %>% as_tibble() %>% mutate(sd = "0.3", coef = coefs)
m1_perc_0.5 <- fixef(m1_perc_t3) %>% as_tibble() %>% mutate(sd = "0.5", coef = coefs)
m1_perc_0.7 <- fixef(m1_perc_t4) %>% as_tibble() %>% mutate(sd = "0.7", coef = coefs)
m1_perc_2 <- fixef(m1_perc_t5) %>% as_tibble() %>% mutate(sd = "2", coef = coefs)
m1_perc_3 <- fixef(m1_perc_t6) %>% as_tibble() %>% mutate(sd = "3", coef = coefs)

m1_perc_freq <- confint(test) %>% as_tibble() %>% slice(1:6) %>% mutate(sd = "Inf (Max Lik)", coef = coefs) %>% 
  rename(Q2.5 = `2.5 %`,
         Q97.5 = `97.5 %`)

combine_ests <- bind_rows(m1_perc_1, m1_perc_0.1, m1_perc_0.3, m1_perc_0.5, m1_perc_0.7,
                          m1_perc_2, m1_perc_3) %>% 
  bind_rows(m1_perc_freq)

prior_sens_plot_perczero <-combine_ests %>% 
  mutate(mle_bayes = case_when(sd == "Inf (Max Lik)" ~ "Maximum Likelihood", TRUE ~ "Bayesian")) %>% 
  mutate(model = case_when(sd == "1 (model in MS)" ~ "Model used in MS", TRUE ~ "Alternative Priors")) %>% 
  ggplot(aes(x = reorder(coef,-Estimate), y = Estimate, ymin = Q2.5, ymax = Q97.5)) +
  geom_pointrange(aes(color = sd, shape = mle_bayes, size = model),
                  position = position_dodge(width = 0.7)) +
  scale_size_manual(values = c(0.3, 0.7)) +
  geom_hline(yintercept = 0) +
  scale_color_viridis_d() +
  labs(x = "Parameter",
       y = "Estimate") +
  ylim(-2,2) +
  # facet_grid(.~sd, scales = "free") +
  theme_classic() +
  coord_flip() +
  guides(color = guide_legend(reverse = T),
         size = F)
  
prior_sens_plot_perczero

ggsave(prior_sens_plot_perczero, file = "plots/prior_sens_plot.jpg", dpi = 600, width = 5, height = 5)  

m1_perc_t0 <- m1_perc

perc_list <- list(m1_perc_t0,
                  m1_perc_t1,
                  m1_perc_t2,
                  m1_perc_t3,
                  m1_perc_t4,
                  m1_perc_t5,
                  m1_perc_t6)

perc_ce <- lapply(perc_list, conditional_effects)

perc_stack1 <- as.matrix(unlist(perc_ce)) %>% as.data.frame() %>% rownames_to_column()
perc_stack2 <- perc_stack1 %>% filter(str_detect(rowname, 'age.behavior.estimate|age.behavior.lower|age.behavior.upper'))
perc_stack <- perc_stack2 %>% separate(rowname, c("effect", "model"), sep = "__", remove = F) %>% 
  separate(model, c("group", "model")) %>% mutate(model = case_when(is.na(model) ~ "0", TRUE ~ model)) %>% 
  separate(effect, c("a","b","measure")) %>% 
  select(measure, group, model,V1) %>% 
  pivot_wider(names_from = "measure", values_from = "V1") %>% 
  mutate(group = case_when(group == "1" ~ "A.D",
                           group == "2" ~ "A.P", 
                           group == "3" ~ "F.D",
                           group == "4" ~ "F.P",
                           group == "5" ~ "H.D",
                           TRUE ~ "H.P")) %>% 
  separate(group, c("age", "behavior")) %>% 
  mutate(Prior = case_when(model == "0" ~ '1 (MS)',
                           model == "1" ~ '0.1',
                           model == "2" ~ '0.3',
                           model == "3" ~ '0.5',
                           model == "4" ~ '0.7',
                           model == "5" ~ '2',
                           model == "6" ~ '3'),
         mle_bayes = "Bayesian") %>% 
  rbind(test_means)


perc_sens_outcome_perczero <- perc_stack %>% 
  mutate(age = fct_relevel(age, "H","F","A"),
         behavior = case_when(behavior == "N" ~ "P", TRUE ~ behavior)) %>% 
  ggplot(aes(x = age, y = estimate, ymax = upper, ymin = lower, color = behavior)) + 
  geom_pointrange(position = position_dodge(width = 0.2),
                  aes(shape = mle_bayes)) +
  geom_line(aes(group = behavior)) +
  facet_grid(.~Prior) +
  labs(x = "Age",
       y = "Proportion Methylated") +
  theme_classic() +
  theme(legend.position = "top")

ggsave(perc_sens_outcome_perczero, file = "plots/perc_sens_outcome_perczero.jpg", dpi = 600, 
       width = 8, height = 3)



#alternate prior param plot


# combine_ests %>% 
#   mutate(mle_bayes = case_when(sd == "Inf (Max Lik)" ~ "Maximum Likelihood", TRUE ~ "Bayesian")) %>% 
#   mutate(model = case_when(sd == "1 (model in MS)" ~ "Model used in ms", TRUE ~ "Alternative Priors")) %>% 
#   ggplot(aes(x = sd, y = Estimate, ymin = Q2.5, ymax = Q97.5)) +
#   geom_pointrange(aes(color = model, shape = mle_bayes)) +
#   geom_hline(yintercept = 0) +
#   labs(x = "Prior SD") +
#   ylim(-2,2) +
#   facet_wrap(~coef, scales = "free") +
#   theme_classic()



# Analysis for Figure 2 ------------------------

change_in_methylation_over_time_cont_data <- read_excel("data/change in methylation over time- cont data.xlsx") %>% 
  mutate(total_loci = total.change + no.change)


# change_cont_brm <- brm(total.change|trials(total_loci) ~ age.transition*behavior + (1|individual),
#                        data = change_in_methylation_over_time_cont_data,
#                        family = binomial(link = "logit"),
#                        prior = c(prior(normal(0,1), class = "b"),
#                                  prior(normal(0,2), class = "Intercept"),
#                                  prior(cauchy(0,1), class = "sd")),
#                        cores = 2, 
#                        sample_prior = T)

change_cont_brm <- readRDS("models/change_cont_brm.rds")

#run mle
test_cont <- glmer(cbind(total.change, total_loci - total.change) ~ age.transition*behavior + (1|individual),
                data = change_in_methylation_over_time_cont_data,
                family = binomial(link = "logit"))

test_means_cont <- lsmeans(test_cont, specs = ~age.transition|behavior) %>% rbind() %>% as_tibble() %>% select(-SE,-df) %>% 
  pivot_longer(cols = c(lsmean, asymp.LCL, asymp.UCL)) %>% 
  mutate(value = inv_logit_scaled(value)) %>% 
  pivot_wider(names_from = name, values_from = value) %>% 
  mutate(Prior = "Inf (Max LIke)",
         mle_bayes = "Maximum Likelihood",
         model = "7") %>% 
  rename(estimate = lsmean,
         lower = asymp.LCL,
         upper = asymp.UCL)

test_sum_cont <- summary(test_cont)

confint(test_cont) %>% as_tibble() %>% slice(1:6)

coefs_cont <- c("intercept", "age.transitionFA", "behaviorNA", "age.transitionFA:behaviorN")

# change_cont_brm_t0 <- readRDS("models/change_cont_brm.rds")
# 
# change_cont_brm_t1 <- update(change_cont_brm_t0,  prior = c(prior(normal(0,0.1), class = "b"),
#                                                         prior(normal(0,2), class = "Intercept"),
#                                                         prior(cauchy(0,1), class = "sd")),
#                             chains = 1)
# 
# change_cont_brm_t2 <- update(change_cont_brm_t0,  prior = c(prior(normal(0,0.3), class = "b"),
#                                                         prior(normal(0,2), class = "Intercept"),
#                                                         prior(cauchy(0,1), class = "sd")),
#                             chains = 1)
# 
# change_cont_brm_t3 <- update(change_cont_brm_t0,  prior = c(prior(normal(0,0.5), class = "b"),
#                                                             prior(normal(0,2), class = "Intercept"),
#                                                             prior(cauchy(0,1), class = "sd")),
#                              chains = 1)
# change_cont_brm_t4 <- update(change_cont_brm_t0,  prior = c(prior(normal(0,0.7), class = "b"),
#                                                             prior(normal(0,2), class = "Intercept"),
#                                                             prior(cauchy(0,1), class = "sd")),
#                              chains = 1)
# change_cont_brm_t5 <- update(change_cont_brm_t0,  prior = c(prior(normal(0,2), class = "b"),
#                                                             prior(normal(0,2), class = "Intercept"),
#                                                             prior(cauchy(0,1), class = "sd")),
#                              chains = 1)
# change_cont_brm_t6 <- update(change_cont_brm_t0,  prior = c(prior(normal(0,3), class = "b"),
#                                                             prior(normal(0,2), class = "Intercept"),
#                                                             prior(cauchy(0,1), class = "sd")),
#                              chains = 1)
# 
# 
# saveRDS(change_cont_brm_t0, file = "models/change_cont_brm_t0.rds")
# saveRDS(change_cont_brm_t1, file = "models/change_cont_brm_t1.rds")
# saveRDS(change_cont_brm_t2, file = "models/change_cont_brm_t2.rds")
# saveRDS(change_cont_brm_t3, file = "models/change_cont_brm_t3.rds")
# saveRDS(change_cont_brm_t4, file = "models/change_cont_brm_t4.rds")
# saveRDS(change_cont_brm_t5, file = "models/change_cont_brm_t5.rds")
# saveRDS(change_cont_brm_t6, file = "models/change_cont_brm_t6.rds")


change_cont_brm_t0 <- readRDS("models/change_cont_brm_t0.rds")
change_cont_brm_t1 <- readRDS("models/change_cont_brm_t1.rds")
change_cont_brm_t2 <- readRDS("models/change_cont_brm_t2.rds")
change_cont_brm_t3 <- readRDS("models/change_cont_brm_t3.rds")
change_cont_brm_t4 <- readRDS("models/change_cont_brm_t4.rds")
change_cont_brm_t5 <- readRDS("models/change_cont_brm_t5.rds")
change_cont_brm_t6 <- readRDS("models/change_cont_brm_t6.rds")


m1_cont_1 <- fixef(change_cont_brm_t0) %>% as_tibble() %>% mutate(sd = "1 (model in MS)", coef = coefs_cont)
m1_cont_0.1 <- fixef(change_cont_brm_t1) %>% as_tibble() %>% mutate(sd = "0.1", coef = coefs_cont)
m1_cont_0.3 <- fixef(change_cont_brm_t2) %>% as_tibble() %>% mutate(sd = "0.3", coef = coefs_cont)
m1_cont_0.5 <- fixef(change_cont_brm_t3) %>% as_tibble() %>% mutate(sd = "0.5", coef = coefs_cont)
m1_cont_0.7 <- fixef(change_cont_brm_t4) %>% as_tibble() %>% mutate(sd = "0.7", coef = coefs_cont)
m1_cont_2 <- fixef(change_cont_brm_t5) %>% as_tibble() %>% mutate(sd = "2", coef = coefs_cont)
m1_cont_3 <- fixef(change_cont_brm_t6) %>% as_tibble() %>% mutate(sd = "3", coef = coefs_cont)

m1_cont_freq <- confint(test_cont) %>% as_tibble() %>% slice(2:5) %>% mutate(sd = "Inf (Max Lik)", coef = coefs_cont) %>% 
  rename(Q2.5 = `2.5 %`,
         Q97.5 = `97.5 %`) %>% 
  mutate(Estimate = fixef(test_cont))

combine_ests_cont <- bind_rows(m1_cont_1, m1_cont_0.1, m1_cont_0.3, m1_cont_0.5, m1_cont_0.7,
                          m1_cont_2, m1_cont_3) %>% 
  bind_rows(m1_cont_freq)

prior_sens_plot_cont <- combine_ests_cont %>% 
  mutate(mle_bayes = case_when(sd == "Inf (Max Lik)" ~ "Maximum Likelihood", TRUE ~ "Bayesian")) %>% 
  mutate(model = case_when(sd == "1 (model in MS)" ~ "Model used in MS", TRUE ~ "Alternative Priors")) %>% 
  ggplot(aes(x = reorder(coef,-Estimate), y = Estimate, ymin = Q2.5, ymax = Q97.5)) +
  geom_pointrange(aes(color = sd, shape = mle_bayes, size = model),
                  position = position_dodge(width = 0.7)) +
  scale_size_manual(values = c(0.3, 0.7)) +
  geom_hline(yintercept = 0) +
  scale_color_viridis_d() +
  labs(x = "Parameter",
       y = "Estimate") +
  ylim(-2,2) +
  # facet_grid(.~sd, scales = "free") +
  theme_classic() +
  coord_flip() +
  guides(color = guide_legend(reverse = T),
         size = F)

prior_sens_plot_cont

ggsave(prior_sens_plot_perczero, file = "plots/prior_sens_plot.jpg", dpi = 600, width = 5, height = 5)  




# make figure conditional effects
cont_list <- list(change_cont_brm_t0,
                  change_cont_brm_t1,
                  change_cont_brm_t2,
                  change_cont_brm_t3,
                  change_cont_brm_t4,
                  change_cont_brm_t5,
                  change_cont_brm_t6)

cont_ce <- lapply(cont_list, conditional_effects)

cont_stack <- as.matrix(unlist(cont_ce)) %>% as.data.frame() %>% rownames_to_column()
cont_stack <- cont_stack %>% filter(str_detect(rowname, 'age.transition.behavior.estimate|age.transition.behavior.lower|age.transition.behavior.upper'))
cont_stack_full <- cont_stack %>% separate(rowname, c("effect", "model"), sep = "__", remove = F) %>% 
  separate(model, c("group", "model")) %>% mutate(model = case_when(is.na(model) ~ "0", TRUE ~ model)) %>% 
  separate(effect, c("a","b","c","measure")) %>%
  select(measure, group, model,V1) %>% 
  pivot_wider(names_from = "measure", values_from = "V1") %>% 
  mutate(group = case_when(group == "1" ~ "CF.P",
                           group == "2" ~ "CF.D", 
                           group == "3" ~ "FA.P",
                           group == "4" ~ "FA.D")) %>% 
  separate(group, c("age.transition", "behavior")) %>% 
  mutate(Prior = case_when(model == "0" ~ '1 (MS)',
                           model == "1" ~ '0.1',
                           model == "2" ~ '0.3',
                           model == "3" ~ '0.4',
                           model == "4" ~ '0.7',
                           model == "5" ~ '2',
                           model == "6" ~ '3'),
         mle_bayes = "Bayesian") %>% 
  rbind(test_means_cont)


cont_sens_outcome_cont <- cont_stack_full %>% 
  mutate(behavior = case_when(behavior == "N" ~ "P", TRUE ~ behavior)) %>% 
  ggplot(aes(x = age.transition, y = estimate, ymax = upper, ymin = lower, color = behavior)) + 
  geom_pointrange(position = position_dodge(width = 0.2),
                  aes(shape = mle_bayes)) +
  geom_line(aes(group = behavior)) +
  facet_grid(.~Prior) +
  labs(x = "Transition",
       y = "Proportion of Loci Changed") +
  theme_classic() +
  theme(legend.position = "top")


ggsave(cont_sens_outcome_cont, file = "plots/cont_sens_outcome_cont.jpg", dpi = 600, 
       width = 8, height = 3)


legend_param_sens <- get_legend(prior_sens_plot_perczero)
plot_predict_sens <- plot_grid(perc_sens_outcome_perczero, cont_sens_outcome_cont, nrow = 2, labels = "auto")
plot_param_sens <- plot_grid(prior_sens_plot_perczero + guides(color = F, shape = F), 
                             prior_sens_plot_cont + guides(color = F, shape = F),
                             legend_param_sens, nrow = 1, labels = c("a","b",NA), rel_widths = c(1,1.1,0.3))

ggsave(plot_predict_sens, file = "plots/plot_predict_sens.jpg", dpi = 600, width = 8, height = 5)
ggsave(plot_param_sens, file = "plots/plot_param_sens.jpg", dpi = 600, width = 8, height = 5)
