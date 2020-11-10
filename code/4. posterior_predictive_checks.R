library(brms)
library(tidyverse)
library(janitor)
library(cowplot)
library(viridis)


m1_perc <- readRDS("models/m1_perc.rds")
change_cont_brm <- readRDS(file = "models/change_cont_brm.rds")

ppc_perc <- pp_check(m1_perc, type = "boxplot") +
  labs(title = "Age x Behavior")

ppc_change <- pp_check(change_cont_brm, type = "boxplot") +
  labs(title = "Age Transition x Behavior")


posterior_predictive_checks <- plot_grid(ncol = 1, ppc_perc, ppc_change, labels = "auto")
saveRDS(posterior_predictive_checks, file = "plots/posterior_predictive_checks.rds")
ggsave(posterior_predictive_checks, file = "plots/posterior_predictive_checks.tiff", width = 5, height = 6)
