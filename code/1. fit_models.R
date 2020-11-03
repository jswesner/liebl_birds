library(brms)
library(tidyverse)
library(janitor)
library(RCurl)
library(ggridges)
library(viridis)
library(readxl)


# Age x Behavior ---------------------------------------------------

bird_data <- read_excel("data/binary total methylation dispersal.xlsx")

# fit model
m1_perc <- brm(methylated.perc ~ age*behavior + (1|ind), data = bird_data,
          family = Beta(link = "logit"),
          prior = c(prior(normal(0,1), class = "b"),
                    prior(normal(0,2), class = "Intercept"),
                    prior(cauchy(0,1), class = "sd")),
          cores = 2,
          sample_prior = T)

saveRDS(m1_perc, file = "models/m1_perc.rds")

# Age.transition x Behavior ---------------------------------------------------

change_in_methylation_over_time_cont_data <- read_excel("data/change in methylation over time- cont data.xlsx") %>% 
  mutate(total_loci = total.change + no.change)


change_cont_brm <- brm(total.change|trials(total_loci) ~ age.transition*behavior + (1|individual),
                       data = change_in_methylation_over_time_cont_data,
                       family = binomial(link = "logit"),
                       prior = c(prior(normal(0,1), class = "b"),
                                 prior(normal(0,2), class = "Intercept"),
                                 prior(cauchy(0,1), class = "sd")),
                       cores = 2, 
                       sample_prior = T)

saveRDS(change_cont_brm, file = "models/change_cont_brm.rds")


# Sex ---------------------------------------------------------------------

binary_total_methylation_sex <- read_excel("data/binary total methylation sex.xlsx") %>% as_tibble() %>% 
  clean_names()


meth_sex.brm <- brm(methylated_perc ~ sex + (1|ind),
                    data = binary_total_methylation_sex,
                    family = Beta(link = "logit"),
                    prior = c(prior(normal(0,1), class = "b"),
                              prior(normal(0,2), class = "Intercept"),
                              prior(cauchy(0,1), class = "sd")),
                    cores = 2,
                    sample_prior = T)

saveRDS(meth_sex.brm, file = "models/meth_sex.brm.rds")
