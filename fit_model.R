rm(list = ls())

# loading required libraries
library(ape)
library(brms)
library(cmdstanr)
library(rstudioapi)
library(rstan)
library(tidybayes)
library(bayestestR)
library(tidyverse)
library(emmeans)
library(PNWColors)
library(MetBrewer)
library(extrafont)
library(scales)
library(extrafont)

# path
path_folder <- "D:/non_native_tree/"

# set the maximum allowed size
options(future.globals.maxSize = 7000*1024^2)
set_cmdstan_path(file.path(Sys.getenv("HOME"), ".cmdstan/cmdstan-2.35.0/"))
treetrait_fit_one <- brm(
  formula = distribution ~ Bark_thickness + Leaf_P_per_mass + Root_depth + Seed_dry_mass +
    Specific_leaf_area + Stem_conduit_diameter + Tree_height + Wood_density + Totals +
    (1|gr(spacc_phylo, cov = treetrait_vcv_one)),
  data = treetrait_modeldata_one,
  family = bernoulli(link = "logit"), sample_prior = TRUE,
  prior = c(
    prior(normal(0,1), "b"),
    prior(normal(0,1), "Intercept"),
    prior(normal(0,.5), "sd")),
  cores = 32, chains = 4, warmup = 2000, iter = 4000,
  backend = "cmdstanr", threads = threading(8),
  data2 = list(treetrait_vcv_one = treetrait_vcv_one),
  stan_model_args = list(stanc_options = list("O1")),
  file = paste0(path_folder, "05_trait_data/treetrait_model_one.rds"))

# Bark thickness
tree_model_one %>%
  emtrends(~ Bark_thickness,
           var = "Bark_thickness",
           at = list(Leaf_P_per_mass = 0,
                     Root_depth = 0,
                     Seed_dry_mass = 0,
                     Specific_leaf_area = 0,
                     Stem_conduit_diameter = 0,
                     Tree_height = 0,
                     Wood_density = 0,
                     Totals = 0),
           epred = TRUE) %>%
  gather_emmeans_draws() %>%
  median_hdi(.width = .9)

# Root_depth
tree_model_one %>%
  emtrends(~ Root_depth,
           var = "Root_depth",
           at = list(Bark_thickness = 0,
                     Leaf_P_per_mass = 0,
                     Seed_dry_mass = 0,
                     Specific_leaf_area = 0,
                     Stem_conduit_diameter = 0,
                     Tree_height = 0,
                     Wood_density = 0,
                     Totals = 0),
           epred = TRUE) %>%
  gather_emmeans_draws() %>%
  median_hdi(.width = .9)

# seed dry mass
tree_model_one %>%
  emtrends(~ Seed_dry_mass,
           var = "Seed_dry_mass",
           at = list(Bark_thickness = 0,
                     Leaf_P_per_mass = 0,
                     Root_depth = 0,
                     Specific_leaf_area = 0,
                     Stem_conduit_diameter = 0,
                     Tree_height = 0,
                     Wood_density = 0,
                     Totals = 0),
           epred = TRUE) %>%
  gather_emmeans_draws() %>%
  median_hdi(.width = .9)

# Stem_conduit_diameter
tree_model_one %>%
  emtrends(~ Stem_conduit_diameter,
           var = "Stem_conduit_diameter",
           at = list(Bark_thickness = 0,
                     Leaf_P_per_mass = 0,
                     Root_depth = 0,
                     Seed_dry_mass = 0,
                     Specific_leaf_area = 0,
                     Tree_height = 0,
                     Wood_density = 0,
                     Totals = 0),
           epred = TRUE) %>%
  gather_emmeans_draws() %>%
  median_hdi(.width = .9)

# Wood_density
tree_model_one %>%
  emtrends(~ Wood_density,
           var = "Wood_density",
           at = list(Bark_thickness = 0,
                     Leaf_P_per_mass = 0,
                     Root_depth = 0,
                     Seed_dry_mass = 0,
                     Specific_leaf_area = 0,
                     Stem_conduit_diameter = 0,
                     Tree_height = 0,
                     Totals = 0),
           epred = TRUE) %>%
  gather_emmeans_draws() %>%
  median_hdi(.width = .9)

# Leaf P per mass
newdata_leafppermass <- expand_grid(Bark_thickness = 0,
                                    Leaf_P_per_mass = seq(from = min(treetrait_modeldata$Leaf_P_per_mass),
                                                          to = max(treetrait_modeldata$Leaf_P_per_mass),
                                                          length.out = 100),
                                    Root_depth = 0,
                                    Seed_dry_mass = 0,
                                    Specific_leaf_area = 0,
                                    Stem_conduit_diameter = 0,
                                    Tree_height = 0,
                                    Wood_density = 0,
                                    Totals = 0)

grand_mean_leafppermass <- tree_model_one %>%
  epred_draws(newdata = newdata_leafppermass,
              re_formula = NA)

save(grand_mean_leafppermass, file = paste0(path_folder,"05_trait_data/grand_mean_leafppermass.Rdata"))

load(paste0(path_folder, "05_trait_data/grand_mean_leafppermass.Rdata"))

tree_model_one %>%
  emtrends(~ Leaf_P_per_mass,
           var = "Leaf_P_per_mass",
           at = list(Bark_thickness = 0,
                     Root_depth = 0,
                     Seed_dry_mass = 0,
                     Specific_leaf_area = 0,
                     Stem_conduit_diameter = 0,
                     Tree_height = 0,
                     Wood_density = 0,
                     Totals = 0),
           epred = TRUE) %>%
  gather_emmeans_draws() %>%
  median_hdi(.width = .9)

# Specific leaf area
newdata_specificleafarea <- expand_grid(Bark_thickness = 0,
                                        Leaf_P_per_mass = 0,
                                        Root_depth = 0,
                                        Seed_dry_mass = 0,
                                        Specific_leaf_area = seq(from = min(treetrait_modeldata$Specific_leaf_area),
                                                                 to = max(treetrait_modeldata$Specific_leaf_area),
                                                                 length.out = 100),
                                        Stem_conduit_diameter = 0,
                                        Tree_height = 0,
                                        Wood_density = 0,
                                        Totals = 0)

grand_mean_specificleafarea <- tree_model_one %>%
  epred_draws(newdata = newdata_specificleafarea,
              re_formula = NA)

save(grand_mean_specificleafarea, file = paste0(path_folder,"05_trait_data/grand_mean_specificleafarea.Rdata"))

load(paste0(path_folder, "05_trait_data/grand_mean_specificleafarea.Rdata"))

tree_model_one %>%
  emtrends(~ Specific_leaf_area,
           var = "Specific_leaf_area",
           at = list(Bark_thickness = 0,
                     Leaf_P_per_mass = 0,
                     Root_depth = 0,
                     Seed_dry_mass = 0,
                     Stem_conduit_diameter = 0,
                     Tree_height = 0,
                     Wood_density = 0,
                     Totals = 0),
           epred = TRUE) %>%
  gather_emmeans_draws() %>%
  median_hdi(.width = .9)

# tree height
newdata_treeheight <- expand_grid(Bark_thickness = 0,
                                  Leaf_P_per_mass = 0,
                                  Root_depth = 0,
                                  Seed_dry_mass = 0,
                                  Specific_leaf_area = 0,
                                  Stem_conduit_diameter = 0,
                                  Tree_height = seq(from = min(treetrait_modeldata$Tree_height),
                                                    to = max(treetrait_modeldata$Tree_height),
                                                    length.out = 100),
                                  Wood_density = 0,
                                  Totals = 0)

grand_mean_treeheight <- tree_model_one %>%
  epred_draws(newdata = newdata_treeheight,
              re_formula = NA)

save(grand_mean_treeheight, file = paste0(path_folder,"05_trait_data/grand_mean_treeheight.Rdata"))

load(paste0(path_folder, "05_trait_data/grand_mean_treeheight.Rdata"))

tree_model_one %>%
  emtrends(~ Tree_height,
           var = "Tree_height",
           at = list(Bark_thickness = 0,
                     Leaf_P_per_mass = 0,
                     Root_depth = 0,
                     Seed_dry_mass = 0,
                     Specific_leaf_area = 0,
                     Stem_conduit_diameter = 0,
                     Wood_density = 0,
                     Totals = 0),
           epred = TRUE) %>%
  gather_emmeans_draws() %>%
  median_hdi(.width = .9)

setwd(path_folder)
png(file=paste0(path_folder, "05_trait_data/probability_being_introduced_variables.png"),
    bg="transparent", width = 8, height = 8, units = "cm", res = 300)
ggplot() +
  stat_lineribbon(mapping = aes(x = Leaf_P_per_mass, y = .epred), data = grand_mean_leafppermass,
                  .width = .9, color = "#5D74A5", fill = alpha("#5D74A5", 0.1), linewidth = 0.5) +
  stat_lineribbon(mapping = aes(x = Specific_leaf_area, y = .epred), data = grand_mean_specificleafarea,
                  .width = .9, color = "#F6D4A9", fill = alpha("#F6D4A9", 0.15), linewidth = 0.5) +
  stat_lineribbon(mapping = aes(x = Tree_height, y = .epred), data = grand_mean_treeheight,
                  .width = .9, color = "#A8554E", fill = alpha("#A8554E", 0.1), linewidth = 0.5) +
  labs(x = "Standardized variables",
       y = "Probability of being introduced") +
  theme(panel.background = element_blank(),
        axis.line = element_line(),
        axis.text = element_text(family = "Times New Roman", size = 10),
        axis.title = element_text(family = "Times New Roman", size = 12)) +
  scale_y_continuous(labels = c(0, sprintf("%.2f", c(0.25, 0.5, 0.75, 1))),
                     breaks = c(0, 0.25, 0.50, 0.75, 1.00),
                     limits = c(0, 1)) +
  scale_x_continuous(labels = c(-5, -2.5, 0, 2.5),
                     breaks = c(-5, -2.5, 0, 2.5))

dev.off()

###################
# Threatened tree
###################

# set the maximum allowed size
options(future.globals.maxSize = 7000*1024^2)
set_cmdstan_path(file.path(Sys.getenv("HOME"), ".cmdstan/cmdstan-2.35.0/"))
treetrait_fit_threatened <- brm(
  formula = distribution ~ Bark_thickness + Leaf_P_per_mass + Root_depth + Seed_dry_mass +
    Specific_leaf_area + Stem_conduit_diameter + Tree_height + Wood_density + Totals +
    (1|gr(spacc_phylo, cov = treetrait_vcv_threatened)), 
  data = treetrait_modeldata_threatened, 
  family = bernoulli(link = "logit"), sample_prior = TRUE,
  prior = c(
    prior(normal(0,1), "b"),
    prior(normal(0,1), "Intercept"),
    prior(normal(0,.5), "sd")),
  cores = 32, chains = 4, warmup = 2000, iter = 4000,
  backend = "cmdstanr", threads = threading(8),
  data2 = list(treetrait_vcv_threatened = treetrait_vcv_threatened),
  stan_model_args = list(stanc_options = list("O1")),
  file = paste0(path_folder, "05_trait_data/treetrait_model_threatened.rds"))

# Bark thickness
tree_model_threatened %>%
  emtrends(~ Bark_thickness,
           var = "Bark_thickness",
           at = list(Leaf_P_per_mass = 0,
                     Root_depth = 0,
                     Seed_dry_mass = 0,
                     Specific_leaf_area = 0,
                     Stem_conduit_diameter = 0,
                     Tree_height = 0,
                     Wood_density = 0,
                     Totals = 0),
           epred = TRUE) %>%
  gather_emmeans_draws() %>%
  median_hdi(.width = .9)

# Root_depth
tree_model_threatened %>%
  emtrends(~ Root_depth,
           var = "Root_depth",
           at = list(Bark_thickness = 0,
                     Leaf_P_per_mass = 0,
                     Seed_dry_mass = 0,
                     Specific_leaf_area = 0,
                     Stem_conduit_diameter = 0,
                     Tree_height = 0,
                     Wood_density = 0,
                     Totals = 0),
           epred = TRUE) %>%
  gather_emmeans_draws() %>%
  median_hdi(.width = .9)

# Seed dry mass
tree_model_threatened %>%
  emtrends(~ Seed_dry_mass,
           var = "Seed_dry_mass",
           at = list(Bark_thickness = 0,
                     Leaf_P_per_mass = 0,
                     Root_depth = 0,
                     Specific_leaf_area = 0,
                     Stem_conduit_diameter = 0,
                     Tree_height = 0,
                     Wood_density = 0,
                     Totals = 0),
           epred = TRUE) %>%
  gather_emmeans_draws() %>%
  median_hdi(.width = .9)

# Stem_conduit_diameter
tree_model_threatened %>%
  emtrends(~ Stem_conduit_diameter,
           var = "Stem_conduit_diameter",
           at = list(Bark_thickness = 0,
                     Leaf_P_per_mass = 0,
                     Root_depth = 0,
                     Seed_dry_mass = 0,
                     Specific_leaf_area = 0,
                     Tree_height = 0,
                     Wood_density = 0,
                     Totals = 0),
           epred = TRUE) %>%
  gather_emmeans_draws() %>%
  median_hdi(.width = .9)

# Wood_density
tree_model_threatened %>%
  emtrends(~ Wood_density,
           var = "Wood_density",
           at = list(Bark_thickness = 0,
                     Leaf_P_per_mass = 0,
                     Root_depth = 0,
                     Seed_dry_mass = 0,
                     Specific_leaf_area = 0,
                     Stem_conduit_diameter = 0,
                     Tree_height = 0,
                     Totals = 0),
           epred = TRUE) %>%
  gather_emmeans_draws() %>%
  median_hdi(.width = .9)

# Leaf P per mass
newdata_leafppermass_threatened <- expand_grid(Bark_thickness = 0,
                                               Leaf_P_per_mass = seq(from = min(treetrait_modeldata_threatened$Leaf_P_per_mass),
                                                                     to = max(treetrait_modeldata_threatened$Leaf_P_per_mass),
                                                                     length.out = 100),
                                               Root_depth = 0,
                                               Seed_dry_mass = 0,
                                               Specific_leaf_area = 0,
                                               Stem_conduit_diameter = 0,
                                               Tree_height = 0,
                                               Wood_density = 0,
                                               Totals = 0)

grand_mean_leafppermass_threatened <- tree_model_threatened %>%
  epred_draws(newdata = newdata_leafppermass_threatened,
              re_formula = NA)

save(grand_mean_leafppermass_threatened, file = paste0(path_folder,"05_trait_data/grand_mean_leafppermass_threatened.Rdata"))

load(paste0(path_folder, "05_trait_data/grand_mean_leafppermass_threatened.Rdata"))

tree_model_threatened %>%
  emtrends(~ Leaf_P_per_mass,
           var = "Leaf_P_per_mass",
           at = list(Bark_thickness = 0,
                     Root_depth = 0,
                     Seed_dry_mass = 0,
                     Specific_leaf_area = 0,
                     Stem_conduit_diameter = 0,
                     Tree_height = 0,
                     Wood_density = 0,
                     Totals = 0),
           epred = TRUE) %>%
  gather_emmeans_draws() %>%
  median_hdi(.width = .9)

# Specific leaf area
newdata_specificleafarea_threatened <- expand_grid(Bark_thickness = 0,
                                                   Leaf_P_per_mass = 0,
                                                   Root_depth = 0,
                                                   Seed_dry_mass = 0,
                                                   Specific_leaf_area = seq(from = min(treetrait_modeldata_threatened$Specific_leaf_area),
                                                                            to = max(treetrait_modeldata_threatened$Specific_leaf_area),
                                                                            length.out = 100),
                                                   Stem_conduit_diameter = 0,
                                                   Tree_height = 0,
                                                   Wood_density = 0,
                                                   Totals = 0)

grand_mean_specificleafarea_threatened <- tree_model_threatened %>%
  epred_draws(newdata = newdata_specificleafarea_threatened,
              re_formula = NA)

save(grand_mean_specificleafarea_threatened, file = paste0(path_folder,"05_trait_data/grand_mean_specificleafarea_threatened.Rdata"))

load(paste0(path_folder, "05_trait_data/grand_mean_specificleafarea_threatened.Rdata"))

tree_model_threatened %>%
  emtrends(~ Specific_leaf_area,
           var = "Specific_leaf_area",
           at = list(Bark_thickness = 0,
                     Leaf_P_per_mass = 0,
                     Root_depth = 0,
                     Seed_dry_mass = 0,
                     Stem_conduit_diameter = 0,
                     Tree_height = 0,
                     Wood_density = 0,
                     Totals = 0),
           epred = TRUE) %>%
  gather_emmeans_draws() %>%
  median_hdi(.width = .9)

# Tree height
newdata_treeheight_threatened <- expand_grid(Bark_thickness = 0,
                                             Leaf_P_per_mass = 0,
                                             Root_depth = 0,
                                             Seed_dry_mass = 0,
                                             Specific_leaf_area = 0,
                                             Stem_conduit_diameter = 0,
                                             Tree_height = seq(from = min(treetrait_modeldata_threatened$Tree_height),
                                                               to = max(treetrait_modeldata_threatened$Tree_height),
                                                               length.out = 100),
                                             Wood_density = 0,
                                             Totals = 0)

grand_mean_treeheight_threatened <- tree_model_threatened %>%
  epred_draws(newdata = newdata_treeheight_threatened,
              re_formula = NA)

save(grand_mean_treeheight_threatened, file = paste0(path_folder,"05_trait_data/grand_mean_treeheight_threatened.Rdata"))

load(paste0(path_folder, "05_trait_data/grand_mean_treeheight_threatened.Rdata"))

tree_model_threatened %>%
  emtrends(~ Tree_height,
           var = "Tree_height",
           at = list(Bark_thickness = 0,
                     Leaf_P_per_mass = 0,
                     Root_depth = 0,
                     Seed_dry_mass = 0,
                     Specific_leaf_area = 0,
                     Stem_conduit_diameter = 0,
                     Wood_density = 0,
                     Totals = 0),
           epred = TRUE) %>%
  gather_emmeans_draws() %>%
  median_hdi(.width = .9)

setwd(path_folder)
png(file=paste0(path_folder, "05_trait_data/probability_being_introduced_variables_threatened.png"),
    bg="transparent", width = 8, height = 8, units = "cm", res = 300)
ggplot() +
  stat_lineribbon(mapping = aes(x = Leaf_P_per_mass, y = .epred), data = grand_mean_leafppermass_threatened,
                  .width = .9, color = "#5D74A5", fill = alpha("#5D74A5", 0.1), linewidth = 0.5) +
  stat_lineribbon(mapping = aes(x = Specific_leaf_area, y = .epred), data = grand_mean_specificleafarea_threatened,
                  .width = .9, color = "#F6D4A9", fill = alpha("#F6D4A9", 0.15), linewidth = 0.5) +
  stat_lineribbon(mapping = aes(x = Tree_height, y = .epred), data = grand_mean_treeheight_threatened,
                  .width = .9, color = "#A8554E", fill = alpha("#A8554E", 0.1), linewidth = 0.5) +
  labs(x = "Standardized variables",
       y = "Probability of being introduced") +
  theme(panel.background = element_blank(),
        axis.line = element_line(),
        axis.text = element_text(family = "Times New Roman", size = 10),
        axis.title = element_text(family = "Times New Roman", size = 12)) +
  scale_y_continuous(labels = c(0, sprintf("%.2f", c(0.25, 0.5, 0.75, 1))),
                     breaks = c(0, 0.25, 0.50, 0.75, 1.00),
                     limits = c(0, 1))
dev.off()

#######################
# non-threatened tree
#######################

# set the maximum allowed size
options(future.globals.maxSize = 7000*1024^2)
set_cmdstan_path(file.path(Sys.getenv("HOME"), ".cmdstan/cmdstan-2.35.0/"))
treetrait_fit_nonthreatened_one <- brm(
  formula = distribution ~ Bark_thickness + Leaf_P_per_mass + Root_depth + Seed_dry_mass +
    Specific_leaf_area + Stem_conduit_diameter + Tree_height + Wood_density + Totals +
    (1|gr(spacc_phylo, cov = treetrait_vcv_nonthreatened_one)),
  data = treetrait_modeldata_nonthreatened_one,
  family = bernoulli(link = "logit"), sample_prior = TRUE,
  prior = c(
    prior(normal(0,1), "b"),
    prior(normal(0,1), "Intercept"),
    prior(normal(0,.5), "sd")),
  cores = 32, chains = 4, warmup = 2000, iter = 4000,
  backend = "cmdstanr", threads = threading(8),
  data2 = list(treetrait_vcv_nonthreatened_one = treetrait_vcv_nonthreatened_one),
  stan_model_args = list(stanc_options = list("O1")),
  file = paste0(path_folder, "05_trait_data/treetrait_model_nonthreatened_one.rds"))

# Bark thickness
tree_model_nonthreatened_one %>%
  emtrends(~ Bark_thickness,
           var = "Bark_thickness",
           at = list(Leaf_P_per_mass = 0,
                     Root_depth = 0,
                     Seed_dry_mass = 0,
                     Specific_leaf_area = 0,
                     Stem_conduit_diameter = 0,
                     Tree_height = 0,
                     Wood_density = 0,
                     Totals = 0),
           epred = TRUE) %>%
  gather_emmeans_draws() %>%
  median_hdi(.width = .9)

# Root_depth
tree_model_nonthreatened_one %>%
  emtrends(~ Root_depth,
           var = "Root_depth",
           at = list(Bark_thickness = 0,
                     Leaf_P_per_mass = 0,
                     Seed_dry_mass = 0,
                     Specific_leaf_area = 0,
                     Stem_conduit_diameter = 0,
                     Tree_height = 0,
                     Wood_density = 0,
                     Totals = 0),
           epred = TRUE) %>%
  gather_emmeans_draws() %>%
  median_hdi(.width = .9)

# Seed_dry_mass
tree_model_nonthreatened_one %>%
  emtrends(~ Seed_dry_mass,
           var = "Seed_dry_mass",
           at = list(Bark_thickness = 0,
                     Leaf_P_per_mass = 0,
                     Root_depth = 0,
                     Specific_leaf_area = 0,
                     Stem_conduit_diameter = 0,
                     Tree_height = 0,
                     Wood_density = 0,
                     Totals = 0),
           epred = TRUE) %>%
  gather_emmeans_draws() %>%
  median_hdi(.width = .9)

# Stem_conduit_diameter
tree_model_nonthreatened_one %>%
  emtrends(~ Stem_conduit_diameter,
           var = "Stem_conduit_diameter",
           at = list(Bark_thickness = 0,
                     Leaf_P_per_mass = 0,
                     Root_depth = 0,
                     Seed_dry_mass = 0,
                     Specific_leaf_area = 0,
                     Tree_height = 0,
                     Wood_density = 0,
                     Totals = 0),
           epred = TRUE) %>%
  gather_emmeans_draws() %>%
  median_hdi(.width = .9)

# Wood_density
tree_model_nonthreatened_one %>%
  emtrends(~ Wood_density,
           var = "Wood_density",
           at = list(Bark_thickness = 0,
                     Leaf_P_per_mass = 0,
                     Root_depth = 0,
                     Seed_dry_mass = 0,
                     Specific_leaf_area = 0,
                     Stem_conduit_diameter = 0,
                     Tree_height = 0,
                     Totals = 0),
           epred = TRUE) %>%
  gather_emmeans_draws() %>%
  median_hdi(.width = .9)

# Leaf P per mass
newdata_leafppermass_nonthreatened <- expand_grid(Bark_thickness = 0,
                                                  Leaf_P_per_mass = seq(from = min(treetrait_modeldata_nonthreatened$Leaf_P_per_mass),
                                                                        to = max(treetrait_modeldata_nonthreatened$Leaf_P_per_mass),
                                                                        length.out = 100),
                                                  Root_depth = 0,
                                                  Seed_dry_mass = 0,
                                                  Specific_leaf_area = 0,
                                                  Stem_conduit_diameter = 0,
                                                  Tree_height = 0,
                                                  Wood_density = 0,
                                                  Totals = 0)

grand_mean_leafppermass_nonthreatened <- tree_model_nonthreatened_one %>%
  epred_draws(newdata = newdata_leafppermass_nonthreatened,
              re_formula = NA)

save(grand_mean_leafppermass_nonthreatened, file = paste0(path_folder,"05_trait_data/grand_mean_leafppermass_nonthreatened.Rdata"))

load(paste0(path_folder, "05_trait_data/grand_mean_leafppermass_nonthreatened.Rdata"))

tree_model_nonthreatened_one %>%
  emtrends(~ Leaf_P_per_mass,
           var = "Leaf_P_per_mass",
           at = list(Bark_thickness = 0,
                     Root_depth = 0,
                     Seed_dry_mass = 0,
                     Specific_leaf_area = 0,
                     Stem_conduit_diameter = 0,
                     Tree_height = 0,
                     Wood_density = 0,
                     Totals = 0),
           epred = TRUE) %>%
  gather_emmeans_draws() %>%
  median_hdi(.width = .9)

# Specific leaf area
newdata_specificleafarea_nonthreatened <- expand_grid(Bark_thickness = 0,
                                                      Leaf_P_per_mass = 0,
                                                      Root_depth = 0,
                                                      Seed_dry_mass = 0,
                                                      Specific_leaf_area = seq(from = min(treetrait_modeldata_nonthreatened$Specific_leaf_area),
                                                                               to = max(treetrait_modeldata_nonthreatened$Specific_leaf_area),
                                                                               length.out = 100),
                                                      Stem_conduit_diameter = 0,
                                                      Tree_height = 0,
                                                      Wood_density = 0,
                                                      Totals = 0)

grand_mean_specificleafarea_nonthreatened <- tree_model_nonthreatened_one %>%
  epred_draws(newdata = newdata_specificleafarea_nonthreatened,
              re_formula = NA)

save(grand_mean_specificleafarea_nonthreatened, file = paste0(path_folder,"05_trait_data/grand_mean_specificleafarea_nonthreatened.Rdata"))

load(paste0(path_folder, "05_trait_data/grand_mean_specificleafarea_nonthreatened.Rdata"))

tree_model_nonthreatened_one %>%
  emtrends(~ Specific_leaf_area,
           var = "Specific_leaf_area",
           at = list(Bark_thickness = 0,
                     Leaf_P_per_mass = 0,
                     Root_depth = 0,
                     Seed_dry_mass = 0,
                     Stem_conduit_diameter = 0,
                     Tree_height = 0,
                     Wood_density = 0,
                     Totals = 0),
           epred = TRUE) %>%
  gather_emmeans_draws() %>%
  median_hdi(.width = .9)

# Tree height
newdata_treeheight_nonthreatened <- expand_grid(Bark_thickness = 0,
                                                Leaf_P_per_mass = 0,
                                                Root_depth = 0,
                                                Seed_dry_mass = 0,
                                                Specific_leaf_area = 0,
                                                Stem_conduit_diameter = 0,
                                                Tree_height = seq(from = min(treetrait_modeldata_nonthreatened$Tree_height),
                                                                  to = max(treetrait_modeldata_nonthreatened$Tree_height),
                                                                  length.out = 100),
                                                Wood_density = 0,
                                                Totals = 0)

grand_mean_treeheight_nonthreatened <- tree_model_nonthreatened_one %>%
  epred_draws(newdata = newdata_treeheight_nonthreatened,
              re_formula = NA)

save(grand_mean_treeheight_nonthreatened, file = paste0(path_folder,"05_trait_data/grand_mean_treeheight_nonthreatened.Rdata"))

load(paste0(path_folder, "05_trait_data/grand_mean_treeheight_nonthreatened.Rdata"))

tree_model_nonthreatened_one %>%
  emtrends(~ Tree_height,
           var = "Tree_height",
           at = list(Bark_thickness = 0,
                     Leaf_P_per_mass = 0,
                     Root_depth = 0,
                     Seed_dry_mass = 0,
                     Specific_leaf_area = 0,
                     Stem_conduit_diameter = 0,
                     Wood_density = 0,
                     Totals = 0),
           epred = TRUE) %>%
  gather_emmeans_draws() %>%
  median_hdi(.width = .9)

setwd(path_folder)
png(file = paste0(path_folder, "05_trait_data/probability_being_introduced_variables_nonthreatened.png"),
    bg = "transparent", width = 8, height = 8, units = "cm", res = 300)
ggplot() +
  stat_lineribbon(mapping = aes(x = Leaf_P_per_mass, y = .epred), data = grand_mean_leafppermass_nonthreatened,
                  .width = .9, color = "#5D74A5", fill = alpha("#5D74A5", 0.1), linewidth = 0.5) +
  stat_lineribbon(mapping = aes(x = Specific_leaf_area, y = .epred), data = grand_mean_specificleafarea_nonthreatened,
                  .width = .9, color = "#F6D4A9", fill = alpha("#F6D4A9", 0.15), linewidth = 0.5) +
  stat_lineribbon(mapping = aes(x = Tree_height, y = .epred), data = grand_mean_treeheight_nonthreatened,
                  .width = .9, color = "#A8554E", fill = alpha("#A8554E", 0.1), linewidth = 0.5) +
  labs(x = "Standardized variables",
       y = "Probability of being introduced") +
  theme(panel.background = element_blank(),
        axis.line = element_line(),
        axis.text = element_text(family = "Times New Roman", size = 10),
        axis.title = element_text(family = "Times New Roman", size = 12)) +
  scale_y_continuous(labels = c(0, sprintf("%.2f", c(0.25, 0.5, 0.75, 1))),
                     breaks = c(0, 0.25, 0.50, 0.75, 1.00),
                     limits = c(0, 1))
dev.off()

##########################
# plot the economic uses
##########################
# overall tree species
data_Totals <- grand_mean_Totals %>% median_hdi(.width = .9) %>% as.data.frame()
data_Totals <- data_Totals[c(1:4, 6:9, 11, 13: 14), ]

grand_mean_ecouses <- data.frame(plantuse = data_Totals$Totals,
                                 median = data_Totals$.value,
                                 lower_hdi = data_Totals$.lower,
                                 upper_hdi = data_Totals$.upper)

setwd(path_folder)
png(file = paste0(path_folder, "05_trait_data/probability_being_introduced_economic_use.png"),
    bg = "transparent", width = 7, height = 8, units = "cm", res = 300)

ggplot(grand_mean_ecouses) +
  geom_errorbar(aes(x = plantuse, ymin = lower_hdi, ymax = upper_hdi), color = "#C57C7B",
                linewidth = 0.4, width = 0.2) +
  geom_point(aes(x = plantuse, y = median), fill = "#C57C7B", color = "#C57C7B", shape = 21, size = 3) +
  labs(x = "Number of economic uses",
       y = "Probability of being introduced") +
  theme(panel.background = element_blank(),
        axis.line = element_line(),
        axis.text = element_text(size = 10, family = "Times New Roman"),
        axis.title = element_text(size = 12, family = "Times New Roman")) +
  scale_y_continuous(labels = c(0,sprintf("%.2f",c(0.25,0.5,0.75,1))),
                     breaks = c(0,0.25,0.50,0.75,1.00),
                     limits = c(0,1)) +
  theme(legend.position = "none")

dev.off()

###########################
# threatened tree species
data_Totals_threatened <- grand_mean_Totals_threatened %>% median_hdi(.width = .9) %>% as.data.frame()
data_Totals_threatened <- data_Totals_threatened[c(1: 4, 6: 12), ]

grand_mean_ecouses_threatened <- data.frame(plantuse = data_Totals_threatened$Totals,
                                            median = data_Totals_threatened$.value,
                                            lower_hdi = data_Totals_threatened$.lower,
                                            upper_hdi = data_Totals_threatened$.upper)

setwd(path_folder)
png(file = paste0(path_folder, "05_trait_data/probability_being_introduced_economic_use_threatened.png"),
    bg = "transparent", width = 7, height = 8, units = "cm", res = 300)

ggplot(grand_mean_ecouses_threatened) +
  geom_errorbar(aes(x = plantuse, ymin = lower_hdi, ymax = upper_hdi), color = "#C57C7B",
                linewidth = 0.4, width = 0.2) +
  geom_point(aes(x = plantuse, y = median), fill = "#C57C7B", color = "#C57C7B", shape = 21, size = 3) +
  labs(x = "Number of economic uses",
       y = "Probability of being introduced") +
  theme(panel.background = element_blank(),
        axis.line = element_line(),
        axis.text = element_text(size = 10, family = "Times New Roman"),
        axis.title = element_text(size = 12, family = "Times New Roman")) +
  scale_y_continuous(labels = c(0,sprintf("%.2f",c(0.25,0.5,0.75,1))),
                     breaks = c(0,0.25,0.50,0.75,1.00),
                     limits = c(0,1)) +
  theme(legend.position = "none")

dev.off()

###############################
# non-threatened tree species
data_Totals_nonthreatened <- grand_mean_Totals_nonthreatened %>% median_hdi(.width = .9) %>% as.data.frame()
data_Totals_nonthreatened <- data_Totals_nonthreatened[c(1: 9, 11: 12), ]

grand_mean_ecouses_nonthreatened <- data.frame(plantuse = data_Totals_nonthreatened$Totals,
                                               median = data_Totals_nonthreatened$.value,
                                               lower_hdi = data_Totals_nonthreatened$.lower,
                                               upper_hdi = data_Totals_nonthreatened$.upper)

setwd(path_folder)
png(file = paste0(path_folder, "05_trait_data/probability_being_introduced_economic_use_nonthreatened.png"),
    bg = "transparent", width = 7, height = 8, units = "cm", res = 300)

ggplot(grand_mean_ecouses_nonthreatened) +
  geom_errorbar(aes(x = plantuse, ymin = lower_hdi, ymax = upper_hdi), color = "#C57C7B",
                linewidth = 0.4, width = 0.2) +
  geom_point(aes(x = plantuse, y = median), fill = "#C57C7B", color = "#C57C7B", shape = 21, size = 3) +
  labs(x = "Number of economic uses",
       y = "Probability of being introduced") +
  theme(panel.background = element_blank(),
        axis.line = element_line(),
        axis.text = element_text(size = 10, family = "Times New Roman"),
        axis.title = element_text(size = 12, family = "Times New Roman")) +
  scale_y_continuous(labels = c(0,sprintf("%.2f",c(0.25,0.5,0.75,1))),
                     breaks = c(0,0.25,0.50,0.75,1.00),
                     limits = c(0,1)) +
  theme(legend.position = "none")

dev.off()







