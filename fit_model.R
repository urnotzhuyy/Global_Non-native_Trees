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

# path
path_folder <- "D:/non_native_tree/"

# purn the phylogenetic tree for the analysis species
# load the trait file
tree_phylo <- read.tree(paste0(path_folder, "04_phylogenetic_distance/global_tree.tre"))

treetrait_modeldata <- read.table(paste0(path_folder, "05_trait_data/functionaltrait_plantuse_tree_modeldata.txt"),
                                  header = TRUE)
treetrait_modeldata$spacc_phylo <- gsub(" ", "_", treetrait_modeldata$spacc, fixed = TRUE)

# read in phylogenetic tree and estimate variance-covariance matrix (i.e. correlation matrix)
treetrait_phylo <- keep.tip(tree_phylo, tip = treetrait_modeldata$spacc_phylo)
treetrait_vcv <- vcv.phylo(treetrait_phylo, corr = TRUE)

# scale the ecological traits for overall tree species
treetrait_modeldata$Bark_thickness <- scale(treetrait_modeldata$Bark_thickness)
treetrait_modeldata$Leaf_P_per_mass <- scale(treetrait_modeldata$Leaf_P_per_mass)
treetrait_modeldata$Root_depth <- scale(treetrait_modeldata$Root_depth)
treetrait_modeldata$Seed_dry_mass <- scale(treetrait_modeldata$Seed_dry_mass)
treetrait_modeldata$Specific_leaf_area <- scale(treetrait_modeldata$Specific_leaf_area)
treetrait_modeldata$Stem_conduit_diameter <- scale(treetrait_modeldata$Stem_conduit_diameter)
treetrait_modeldata$Tree_height <- scale(treetrait_modeldata$Tree_height)
treetrait_modeldata$Wood_density <- scale(treetrait_modeldata$Wood_density)

# Make sure binary traits are factors 
treetrait_modeldata$AnimalFood <- as.factor(treetrait_modeldata$AnimalFood)
treetrait_modeldata$EnvironmentalUses <- as.factor(treetrait_modeldata$EnvironmentalUses)
treetrait_modeldata$Fuels <- as.factor(treetrait_modeldata$Fuels)
treetrait_modeldata$GeneSources <- as.factor(treetrait_modeldata$GeneSources)
treetrait_modeldata$HumanFood <- as.factor(treetrait_modeldata$HumanFood)
treetrait_modeldata$InvertebrateFood <- as.factor(treetrait_modeldata$InvertebrateFood)
treetrait_modeldata$Materials <- as.factor(treetrait_modeldata$Materials)
treetrait_modeldata$Medicines <- as.factor(treetrait_modeldata$Medicines)
treetrait_modeldata$Poisons <- as.factor(treetrait_modeldata$Poisons)
treetrait_modeldata$SocialUses <- as.factor(treetrait_modeldata$SocialUses)

treetrait_modeldata$Totals <- as.factor(treetrait_modeldata$Totals)

# set the maximum allowed size
options(future.globals.maxSize = 7000*1024^2)
set_cmdstan_path(file.path(Sys.getenv("HOME"), ".cmdstan/cmdstan-2.35.0/"))
treetrait_fit <- brm(
  formula = distribution ~ Bark_thickness + Leaf_P_per_mass + Root_depth + Seed_dry_mass +
    Specific_leaf_area + Stem_conduit_diameter + Tree_height + Wood_density + Totals +
    (1|gr(spacc_phylo, cov = treetrait_vcv)),
  data = treetrait_modeldata,
  family = bernoulli(link = "logit"), sample_prior = TRUE,
  prior = c(
    prior(normal(0,1), "b"),
    prior(normal(0,1), "Intercept"),
    prior(normal(0,.5), "sd")),
  cores = 8, chains = 4, warmup = 2000, iter = 4000,
  backend = "cmdstanr", threads = threading(2),
  data2 = list(treetrait_vcv = treetrait_vcv),
  stan_model_args = list(stanc_options = list("O1")),
  file = paste0(path_folder, "05_trait_data/treetrait_model_ecological_scale.rds"))

# load the fitted model
tree_model <- readRDS(paste0(path_folder, "05_trait_data/treetrait_model_ecological_scale.rds"))

# Probability of Direction, pd (also known as the Maximum Probability of Effect, MPE)
tree_pd <- p_direction(tree_model) %>% as.data.frame()

##################
# Bark thickness
tree_model %>%
  emtrends(~ Bark_thickness,
           var = "Bark_thickness",
           at = list(Leaf_P_per_mass = 0,
                     Root_depth = 0,
                     Seed_dry_mass = 0,
                     Specific_leaf_area = 0,
                     Stem_conduit_diameter = 0,
                     Tree_height = 0,
                     Wood_density = 0,
                     Totals = 1),
           epred = TRUE) %>%
  gather_emmeans_draws() %>%
  median_hdi(.width = .9)

###################
# Leaf P per mass
tree_model %>%
  emtrends(~ Leaf_P_per_mass,
           var = "Leaf_P_per_mass",
           at = list(Bark_thickness = 0,
                     Root_depth = 0,
                     Seed_dry_mass = 0,
                     Specific_leaf_area = 0,
                     Stem_conduit_diameter = 0,
                     Tree_height = 0,
                     Wood_density = 0,
                     Totals = 1),
           epred = TRUE) %>%
  gather_emmeans_draws() %>%
  median_hdi(.width = .9)

##############
# Root_depth
tree_model %>%
  emtrends(~ Root_depth,
           var = "Root_depth",
           at = list(Bark_thickness = 0,
                     Leaf_P_per_mass = 0,
                     Seed_dry_mass = 0,
                     Specific_leaf_area = 0,
                     Stem_conduit_diameter = 0,
                     Tree_height = 0,
                     Wood_density = 0,
                     Totals = 1),
           epred = TRUE) %>%
  gather_emmeans_draws() %>%
  median_hdi(.width = .9)

#########################
# Stem_conduit_diameter
tree_model %>%
  emtrends(~ Stem_conduit_diameter,
           var = "Stem_conduit_diameter",
           at = list(Bark_thickness = 0,
                     Leaf_P_per_mass = 0,
                     Root_depth = 0,
                     Seed_dry_mass = 0,
                     Specific_leaf_area = 0,
                     Tree_height = 0,
                     Wood_density = 0,
                     Totals = 1),
           epred = TRUE) %>%
  gather_emmeans_draws() %>%
  median_hdi(.width = .9)

################
# Wood_density
tree_model %>%
  emtrends(~ Wood_density,
           var = "Wood_density",
           at = list(Bark_thickness = 0,
                     Leaf_P_per_mass = 0,
                     Root_depth = 0,
                     Seed_dry_mass = 0,
                     Specific_leaf_area = 0,
                     Stem_conduit_diameter = 0,
                     Tree_height = 0,
                     Totals = 1),
           epred = TRUE) %>%
  gather_emmeans_draws() %>%
  median_hdi(.width = .9)

#################
# Seed dry mass
newdata_seeddrymass <- expand_grid(Bark_thickness = 0,
                                   Leaf_P_per_mass = 0,
                                   Root_depth = 0,
                                   Seed_dry_mass = seq(from = min(treetrait_modeldata$Seed_dry_mass),
                                                       to = max(treetrait_modeldata$Seed_dry_mass),
                                                       length.out = 100),
                                   Specific_leaf_area = 0,
                                   Stem_conduit_diameter = 0,
                                   Tree_height = 0,
                                   Wood_density = 0,
                                   Totals = 1)
                                   
grand_mean_seeddrymass <- tree_model %>%
  epred_draws(newdata = newdata_seeddrymass,
              re_formula = NA)

save(grand_mean_seeddrymass, file = paste0(path_folder,"05_trait_data/grand_mean_seeddrymass.Rdata"))

load(paste0(path_folder, "05_trait_data/grand_mean_seeddrymass.Rdata"))

tree_model %>%
  emtrends(~ Seed_dry_mass,
           var = "Seed_dry_mass",
           at = list(Bark_thickness = 0,
                     Leaf_P_per_mass = 0,
                     Root_depth = 0,
                     Specific_leaf_area = 0,
                     Stem_conduit_diameter = 0,
                     Tree_height = 0,
                     Wood_density = 0,
                     Totals = 1),
           epred = TRUE) %>%
  gather_emmeans_draws() %>%
  median_hdi(.width = .9)

######################
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
                                        Totals = 1)

grand_mean_specificleafarea <- tree_model %>%
  epred_draws(newdata = newdata_specificleafarea,
              re_formula = NA)

save(grand_mean_specificleafarea, file = paste0(path_folder,"05_trait_data/grand_mean_specificleafarea.Rdata"))

load(paste0(path_folder, "05_trait_data/grand_mean_specificleafarea.Rdata"))

tree_model %>%
  emtrends(~ Specific_leaf_area,
           var = "Specific_leaf_area",
           at = list(Bark_thickness = 0,
                     Leaf_P_per_mass = 0,
                     Root_depth = 0,
                     Seed_dry_mass = 0,
                     Stem_conduit_diameter = 0,
                     Tree_height = 0,
                     Wood_density = 0,
                     Totals = 1),
           epred = TRUE) %>%
  gather_emmeans_draws() %>%
  median_hdi(.width = .9)

###############
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
                                  Totals = 1)

grand_mean_treeheight <- tree_model %>%
  epred_draws(newdata = newdata_treeheight,
              re_formula = NA)

save(grand_mean_treeheight, file = paste0(path_folder,"05_trait_data/grand_mean_treeheight.Rdata"))

load(paste0(path_folder, "05_trait_data/grand_mean_treeheight.Rdata"))

tree_model %>%
  emtrends(~ Tree_height,
           var = "Tree_height",
           at = list(Bark_thickness = 0,
                     Leaf_P_per_mass = 0,
                     Root_depth = 0,
                     Seed_dry_mass = 0,
                     Specific_leaf_area = 0,
                     Stem_conduit_diameter = 0,
                     Wood_density = 0,
                     Totals = 1),
           epred = TRUE) %>%
  gather_emmeans_draws() %>%
  median_hdi(.width = .9)

setwd(path_folder)
png(file=paste0(path_folder, "05_trait_data/probability_being_introduced_variables.png"),
    bg="transparent", width = 8, height = 8, units = "cm", res = 300)
ggplot() +
  stat_lineribbon(mapping = aes(x = Seed_dry_mass, y = .epred), data = grand_mean_seeddrymass,
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
# purn the phylogenetic tree for the analysis species
# load the trait file
tree_phylo <- read.tree(paste0(path_folder, "04_phylogenetic_distance/global_tree.tre"))

treetrait_modeldata_threatened <- read.table(paste0(path_folder, "05_trait_data/functionaltrait_plantuse_tree_modeldata_threatened.txt"),
                                             header = TRUE)
treetrait_modeldata_threatened$spacc_phylo <- gsub(" ", "_", treetrait_modeldata_threatened$spacc, fixed = TRUE)

# read in phylogenetic tree and estimate variance-covariance matrix (i.e. correlation matrix)
treetrait_phylo_threatened <- keep.tip(tree_phylo, tip = treetrait_modeldata_threatened$spacc_phylo)
treetrait_vcv_threatened <- vcv.phylo(treetrait_phylo_threatened, corr = TRUE)

# scale the ecological traits for threatened tree species
treetrait_modeldata_threatened$Bark_thickness <- scale(treetrait_modeldata_threatened$Bark_thickness)
treetrait_modeldata_threatened$Leaf_P_per_mass <- scale(treetrait_modeldata_threatened$Leaf_P_per_mass)
treetrait_modeldata_threatened$Root_depth <- scale(treetrait_modeldata_threatened$Root_depth)
treetrait_modeldata_threatened$Seed_dry_mass <- scale(treetrait_modeldata_threatened$Seed_dry_mass)
treetrait_modeldata_threatened$Specific_leaf_area <- scale(treetrait_modeldata_threatened$Specific_leaf_area)
treetrait_modeldata_threatened$Stem_conduit_diameter <- scale(treetrait_modeldata_threatened$Stem_conduit_diameter)
treetrait_modeldata_threatened$Tree_height <- scale(treetrait_modeldata_threatened$Tree_height)
treetrait_modeldata_threatened$Wood_density <- scale(treetrait_modeldata_threatened$Wood_density)

# make sure binary traits are factors
treetrait_modeldata_threatened$AnimalFood <- as.factor(treetrait_modeldata_threatened$AnimalFood)
treetrait_modeldata_threatened$EnvironmentalUses <- as.factor(treetrait_modeldata_threatened$EnvironmentalUses)
treetrait_modeldata_threatened$Fuels <- as.factor(treetrait_modeldata_threatened$Fuels)
treetrait_modeldata_threatened$GeneSources <- as.factor(treetrait_modeldata_threatened$GeneSources)
treetrait_modeldata_threatened$HumanFood <- as.factor(treetrait_modeldata_threatened$HumanFood)
treetrait_modeldata_threatened$InvertebrateFood <- as.factor(treetrait_modeldata_threatened$InvertebrateFood)
treetrait_modeldata_threatened$Materials <- as.factor(treetrait_modeldata_threatened$Materials)
treetrait_modeldata_threatened$Medicines <- as.factor(treetrait_modeldata_threatened$Medicines)
treetrait_modeldata_threatened$Poisons <- as.factor(treetrait_modeldata_threatened$Poisons)
treetrait_modeldata_threatened$SocialUses <- as.factor(treetrait_modeldata_threatened$SocialUses)

treetrait_modeldata_threatened$Totals <- as.factor(treetrait_modeldata_threatened$Totals)

# just the tree height and economic uses
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
  cores = 8, chains = 4, warmup = 2000, iter = 4000,
  backend = "cmdstanr", threads = threading(2),
  data2 = list(treetrait_vcv_threatened = treetrait_vcv_threatened),
  stan_model_args = list(stanc_options = list("O1")),
  file = paste0(path_folder, "05_trait_data/treetrait_model_ecological_scale_threatened.rds"))

# load the fitted model
tree_model_threatened <- readRDS(paste0(path_folder,"05_trait_data/treetrait_model_ecological_scale_threatened.rds"))

# Probability of Direction, pd (also known as the Maximum Probability of Effect, MPE)
tree_pd_threatened <- p_direction(tree_model_threatened) %>% as.data.frame()

##################
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
                     Totals = 1),
           epred = TRUE) %>%
  gather_emmeans_draws() %>%
  median_hdi(.width = .9)

###################
# Leaf P per mass
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
                     Totals = 1),
           epred = TRUE) %>%
  gather_emmeans_draws() %>%
  median_hdi(.width = .9)

##############
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
                     Totals = 1),
           epred = TRUE) %>%
  gather_emmeans_draws() %>%
  median_hdi(.width = .9)

#########################
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
                     Totals = 1),
           epred = TRUE) %>%
  gather_emmeans_draws() %>%
  median_hdi(.width = .9)

################
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
                     Totals = 1),
           epred = TRUE) %>%
  gather_emmeans_draws() %>%
  median_hdi(.width = .9)

#################
# Seed dry mass
newdata_seeddrymass_threatened <- expand_grid(Bark_thickness = 0,
                                              Leaf_P_per_mass = 0,
                                              Root_depth = 0,
                                              Seed_dry_mass = seq(from = min(treetrait_modeldata_threatened$Seed_dry_mass),
                                                                  to = max(treetrait_modeldata_threatened$Seed_dry_mass),
                                                                  length.out = 100),
                                              Specific_leaf_area = 0,
                                              Stem_conduit_diameter = 0,
                                              Tree_height = 0,
                                              Wood_density = 0,
                                              Totals = 1)

grand_mean_seeddrymass_threatened <- tree_model_threatened %>%
  epred_draws(newdata = newdata_seeddrymass_threatened,
              re_formula = NA)

save(grand_mean_seeddrymass_threatened, file = paste0(path_folder,"05_trait_data/grand_mean_seeddrymass_threatened.Rdata"))

load(paste0(path_folder, "05_trait_data/grand_mean_seeddrymass_threatened.Rdata"))

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
                     Totals = 1),
           epred = TRUE) %>%
  gather_emmeans_draws() %>%
  median_hdi(.width = .9)

######################
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
                                                   Totals = 1)

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
                     Totals = 1),
           epred = TRUE) %>%
  gather_emmeans_draws() %>%
  median_hdi(.width = .9)

###############
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
                                             Totals = 1)

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
                     Totals = 1),
           epred = TRUE) %>%
  gather_emmeans_draws() %>%
  median_hdi(.width = .9)

setwd(path_folder)
png(file=paste0(path_folder, "05_trait_data/probability_being_introduced_variables_threatened.png"),
    bg="transparent", width = 8, height = 8, units = "cm", res = 300)
ggplot() +
  stat_lineribbon(mapping = aes(x = Seed_dry_mass, y = .epred), data = grand_mean_seeddrymass_threatened,
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
# purn the phylogenetic tree for the analysis species
# load the trait file
tree_phylo <- read.tree(paste0(path_folder, "04_phylogenetic_distance/global_tree.tre"))

treetrait_modeldata_nonthreatened <- read.table(paste0(path_folder, "05_trait_data/functionaltrait_plantuse_tree_modeldata_nonthreatened.txt"),
                                                header = TRUE)
treetrait_modeldata_nonthreatened$spacc_phylo <- gsub(" ", "_", treetrait_modeldata_nonthreatened$spacc, fixed = TRUE)

# read in phylogenetic tree and estimate variance-covariance matrix (i.e. correlation matrix)
treetrait_phylo_nonthreatened <- keep.tip(tree_phylo, tip = treetrait_modeldata_nonthreatened$spacc_phylo)
treetrait_vcv_nonthreatened <- vcv.phylo(treetrait_phylo_nonthreatened, corr = TRUE)

# scale the ecological traits for non-threatened tree species
treetrait_modeldata_nonthreatened$Bark_thickness <- scale(treetrait_modeldata_nonthreatened$Bark_thickness)
treetrait_modeldata_nonthreatened$Leaf_P_per_mass <- scale(treetrait_modeldata_nonthreatened$Leaf_P_per_mass)
treetrait_modeldata_nonthreatened$Root_depth <- scale(treetrait_modeldata_nonthreatened$Root_depth)
treetrait_modeldata_nonthreatened$Seed_dry_mass <- scale(treetrait_modeldata_nonthreatened$Seed_dry_mass)
treetrait_modeldata_nonthreatened$Specific_leaf_area <- scale(treetrait_modeldata_nonthreatened$Specific_leaf_area)
treetrait_modeldata_nonthreatened$Stem_conduit_diameter <- scale(treetrait_modeldata_nonthreatened$Stem_conduit_diameter)
treetrait_modeldata_nonthreatened$Tree_height <- scale(treetrait_modeldata_nonthreatened$Tree_height)
treetrait_modeldata_nonthreatened$Wood_density <- scale(treetrait_modeldata_nonthreatened$Wood_density)

# make sure binary traits are factors
treetrait_modeldata_nonthreatened$AnimalFood <- as.factor(treetrait_modeldata_nonthreatened$AnimalFood)
treetrait_modeldata_nonthreatened$EnvironmentalUses <- as.factor(treetrait_modeldata_nonthreatened$EnvironmentalUses)
treetrait_modeldata_nonthreatened$Fuels <- as.factor(treetrait_modeldata_nonthreatened$Fuels)
treetrait_modeldata_nonthreatened$GeneSources <- as.factor(treetrait_modeldata_nonthreatened$GeneSources)
treetrait_modeldata_nonthreatened$HumanFood <- as.factor(treetrait_modeldata_nonthreatened$HumanFood)
treetrait_modeldata_nonthreatened$InvertebrateFood <- as.factor(treetrait_modeldata_nonthreatened$InvertebrateFood)
treetrait_modeldata_nonthreatened$Materials <- as.factor(treetrait_modeldata_nonthreatened$Materials)
treetrait_modeldata_nonthreatened$Medicines <- as.factor(treetrait_modeldata_nonthreatened$Medicines)
treetrait_modeldata_nonthreatened$Poisons <- as.factor(treetrait_modeldata_nonthreatened$Poisons)
treetrait_modeldata_nonthreatened$SocialUses <- as.factor(treetrait_modeldata_nonthreatened$SocialUses)

treetrait_modeldata_nonthreatened$Totals <- as.factor(treetrait_modeldata_nonthreatened$Totals)

# just the ecological traits and the number of economic use type
# set the maximum allowed size
options(future.globals.maxSize = 7000*1024^2)
set_cmdstan_path(file.path(Sys.getenv("HOME"), ".cmdstan/cmdstan-2.35.0/"))
treetrait_fit_nonthreatened <- brm(
  formula = distribution ~ Bark_thickness + Leaf_P_per_mass + Root_depth + Seed_dry_mass +
    Specific_leaf_area + Stem_conduit_diameter + Tree_height + Wood_density + Totals +
    (1|gr(spacc_phylo, cov = treetrait_vcv_nonthreatened)), 
  data = treetrait_modeldata_nonthreatened, 
  family = bernoulli(link = "logit"), sample_prior = TRUE,
  prior = c(
    prior(normal(0,1), "b"),
    prior(normal(0,1), "Intercept"),
    prior(normal(0,.5), "sd")),
  cores = 8, chains = 4, warmup = 2000, iter = 4000,
  backend = "cmdstanr", threads = threading(2),
  data2 = list(treetrait_vcv_nonthreatened = treetrait_vcv_nonthreatened),
  stan_model_args = list(stanc_options = list("O1")),
  file = paste0(path_folder, "05_trait_data/treetrait_model_ecological_scale_nonthreatened.rds"))

# load the fitted model
tree_model_nonthreatened <- readRDS(paste0(path_folder,"05_trait_data/treetrait_model_ecological_scale_nonthreatened.rds"))

# Probability of Direction, pd (also known as the Maximum Probability of Effect, MPE)
tree_pd_nonthreatened <- p_direction(tree_model_nonthreatened) %>% as.data.frame()

##################
# Bark thickness
tree_model_nonthreatened %>%
  emtrends(~ Bark_thickness,
           var = "Bark_thickness",
           at = list(Leaf_P_per_mass = 0,
                     Root_depth = 0,
                     Seed_dry_mass = 0,
                     Specific_leaf_area = 0,
                     Stem_conduit_diameter = 0,
                     Tree_height = 0,
                     Wood_density = 0,
                     Totals = 1),
           epred = TRUE) %>%
  gather_emmeans_draws() %>%
  median_hdi(.width = .9)

###################
# Leaf P per mass
tree_model_nonthreatened %>%
  emtrends(~ Leaf_P_per_mass,
           var = "Leaf_P_per_mass",
           at = list(Bark_thickness = 0,
                     Root_depth = 0,
                     Seed_dry_mass = 0,
                     Specific_leaf_area = 0,
                     Stem_conduit_diameter = 0,
                     Tree_height = 0,
                     Wood_density = 0,
                     Totals = 1),
           epred = TRUE) %>%
  gather_emmeans_draws() %>%
  median_hdi(.width = .9)

##############
# Root_depth
tree_model_nonthreatened %>%
  emtrends(~ Root_depth,
           var = "Root_depth",
           at = list(Bark_thickness = 0,
                     Leaf_P_per_mass = 0,
                     Seed_dry_mass = 0,
                     Specific_leaf_area = 0,
                     Stem_conduit_diameter = 0,
                     Tree_height = 0,
                     Wood_density = 0,
                     Totals = 1),
           epred = TRUE) %>%
  gather_emmeans_draws() %>%
  median_hdi(.width = .9)

#########################
# Stem_conduit_diameter
tree_model_nonthreatened %>%
  emtrends(~ Stem_conduit_diameter,
           var = "Stem_conduit_diameter",
           at = list(Bark_thickness = 0,
                     Leaf_P_per_mass = 0,
                     Root_depth = 0,
                     Seed_dry_mass = 0,
                     Specific_leaf_area = 0,
                     Tree_height = 0,
                     Wood_density = 0,
                     Totals = 1),
           epred = TRUE) %>%
  gather_emmeans_draws() %>%
  median_hdi(.width = .9)

################
# Wood_density
tree_model_nonthreatened %>%
  emtrends(~ Wood_density,
           var = "Wood_density",
           at = list(Bark_thickness = 0,
                     Leaf_P_per_mass = 0,
                     Root_depth = 0,
                     Seed_dry_mass = 0,
                     Specific_leaf_area = 0,
                     Stem_conduit_diameter = 0,
                     Tree_height = 0,
                     Totals = 1),
           epred = TRUE) %>%
  gather_emmeans_draws() %>%
  median_hdi(.width = .9)

#################
# Seed dry mass
newdata_seeddrymass_nonthreatened <- expand_grid(Bark_thickness = 0,
                                                 Leaf_P_per_mass = 0,
                                                 Root_depth = 0,
                                                 Seed_dry_mass = seq(from = min(treetrait_modeldata_nonthreatened$Seed_dry_mass),
                                                                     to = max(treetrait_modeldata_threatened$Seed_dry_mass),
                                                                     length.out = 100),
                                                 Specific_leaf_area = 0,
                                                 Stem_conduit_diameter = 0,
                                                 Tree_height = 0,
                                                 Wood_density = 0,
                                                 Totals = 1)

grand_mean_seeddrymass_nonthreatened <- tree_model_nonthreatened %>%
  epred_draws(newdata = newdata_seeddrymass_nonthreatened,
              re_formula = NA)

save(grand_mean_seeddrymass_nonthreatened, file = paste0(path_folder,"05_trait_data/grand_mean_seeddrymass_nonthreatened.Rdata"))

load(paste0(path_folder, "05_trait_data/grand_mean_seeddrymass_nonthreatened.Rdata"))

tree_model_nonthreatened %>%
  emtrends(~ Seed_dry_mass,
           var = "Seed_dry_mass",
           at = list(Bark_thickness = 0,
                     Leaf_P_per_mass = 0,
                     Root_depth = 0,
                     Specific_leaf_area = 0,
                     Stem_conduit_diameter = 0,
                     Tree_height = 0,
                     Wood_density = 0,
                     Totals = 1),
           epred = TRUE) %>%
  gather_emmeans_draws() %>%
  median_hdi(.width = .9)

######################
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
                                                      Totals = 1)

grand_mean_specificleafarea_nonthreatened <- tree_model_nonthreatened %>%
  epred_draws(newdata = newdata_specificleafarea_nonthreatened,
              re_formula = NA)

save(grand_mean_specificleafarea_nonthreatened, file = paste0(path_folder,"05_trait_data/grand_mean_specificleafarea_nonthreatened.Rdata"))

load(paste0(path_folder, "05_trait_data/grand_mean_specificleafarea_nonthreatened.Rdata"))

tree_model_nonthreatened %>%
  emtrends(~ Specific_leaf_area,
           var = "Specific_leaf_area",
           at = list(Bark_thickness = 0,
                     Leaf_P_per_mass = 0,
                     Root_depth = 0,
                     Seed_dry_mass = 0,
                     Stem_conduit_diameter = 0,
                     Tree_height = 0,
                     Wood_density = 0,
                     Totals = 1),
           epred = TRUE) %>%
  gather_emmeans_draws() %>%
  median_hdi(.width = .9)

###############
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
                                                Totals = 1)

grand_mean_treeheight_nonthreatened <- tree_model_nonthreatened %>%
  epred_draws(newdata = newdata_treeheight_nonthreatened,
              re_formula = NA)

save(grand_mean_treeheight_nonthreatened, file = paste0(path_folder,"05_trait_data/grand_mean_treeheight_nonthreatened.Rdata"))

load(paste0(path_folder, "05_trait_data/grand_mean_treeheight_nonthreatened.Rdata"))

tree_model_nonthreatened %>%
  emtrends(~ Tree_height,
           var = "Tree_height",
           at = list(Bark_thickness = 0,
                     Leaf_P_per_mass = 0,
                     Root_depth = 0,
                     Seed_dry_mass = 0,
                     Specific_leaf_area = 0,
                     Stem_conduit_diameter = 0,
                     Wood_density = 0,
                     Totals = 1),
           epred = TRUE) %>%
  gather_emmeans_draws() %>%
  median_hdi(.width = .9)

setwd(path_folder)
png(file=paste0(path_folder, "05_trait_data/probability_being_introduced_variables_nonthreatened.png"),
    bg="transparent", width = 8, height = 8, units = "cm", res = 300)
ggplot() +
  stat_lineribbon(mapping = aes(x = Seed_dry_mass, y = .epred), data = grand_mean_seeddrymass_nonthreatened,
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







