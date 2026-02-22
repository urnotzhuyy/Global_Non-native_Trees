rm(list=ls())

# loading required libraries
library(circlize)
library(extrafont)
library(ggplot2)
library(cowplot)

# path
path_folder <- "D:/non_native_tree/"

#####################################################
# form the file of native continent of tree species
#####################################################
# load the GADM and TDWG matched file
match_mod <- read.csv(paste0(path_folder, "03_migration_flow/gadm_tdwg_match_bgci_nonnativetree_modify.csv"))

# form the file of native continent of tree species
for (id in 1: length(native_all)) {
  
  # species name
  sp <- native_all[id]
  
  # extract native continent of species
  if(sp %in% native_list){
    
    # POWO and IUCN distribution information
    # native data (origin continent)
    sp_file_ori <- read.table(paste0(path_folder, "aggregate_data/tree_data_native_final_update/", sp, "_native.txt"), header = TRUE)
    # origin
    sp_ori <- unique(sp_file_ori$Level1_cod)
    
  }else{
    
    # BGCI distribution information
    # native data (origin continent)
    sp_file_ori <- read.table(paste0(path_folder, "aggregate_data/tree_data_native_final_update/bgci_part/", sp, "_native.txt"), header = TRUE)
    # origin
    sp_ori <- unique(as.character(sp_file_ori$GID_0))
    sp_ori <- unique(match_mod$TDWG_level1_cod[which(match_mod$GID_0 %in% sp_ori)])
    
  }
  
  sp_res <- cbind(sp, sp_ori)
  write.table(sp_res, file = paste0(path_folder, "03_migration_flow/tree_native_continent.txt"),
              col.names = FALSE, row.names = FALSE, append = TRUE)
  
  print(id)
  
}

# modify deal with São Tomé and Príncipe, Fiji, Mayotte
res_file <- read.table(paste0(path_folder, "03_migration_flow/tree_native_continent.txt"))
colnames(res_file) <- c("spacc", "origin")
write.table(res_file, file = paste0(path_folder, "03_migration_flow/tree_native_continent.txt"),
            col.names = TRUE, row.names = FALSE)

res_file <- read.table(paste0(path_folder, "03_migration_flow/tree_native_continent.txt"), header = TRUE)

###################
# random sampling
###################
# overall non-native tree
# load the migration flow file
tree_flow <- read.table(paste0(path_folder, "03_migration_flow/tree_flow_all.txt"), header = TRUE)
# load the observed flow file
flow_data <- read.table(paste0(path_folder, "03_migration_flow/migration_flow_among_continent.txt"), header = TRUE)

# calculate flows for each source-recipient combination
observed <- flow_data$flow

# the total number of non-native tree species in each continent
europe_dest <- length(unique(tree_flow$spacc[which(tree_flow$destination == 1)])) 
africa_dest <- length(unique(tree_flow$spacc[which(tree_flow$destination == 2)])) 
temperateasia_dest <- length(unique(tree_flow$spacc[which(tree_flow$destination == 3)]))
tropicalasia_dest <- length(unique(tree_flow$spacc[which(tree_flow$destination == 4)])) 
australasia_dest <- length(unique(tree_flow$spacc[which(tree_flow$destination == 5)])) 
pacific_dest <- length(unique(tree_flow$spacc[which(tree_flow$destination == 6)])) 
northamerica_dest <- length(unique(tree_flow$spacc[which(tree_flow$destination == 7)])) 
southamerica_dest <- length(unique(tree_flow$spacc[which(tree_flow$destination == 8)])) 

# number of random draws
n_random <- 9999

# create vectors for each source-receptor combination
europe_europe <- rep(0, n_random + 1)
europe_africa <- rep(0, n_random + 1)
europe_temperateasia <- rep(0, n_random + 1)
europe_tropicalasia <- rep(0, n_random + 1)
europe_australasia <- rep(0, n_random + 1)
europe_pacific <- rep(0, n_random + 1)
europe_northamerica <- rep(0, n_random + 1)
europe_southamerica <- rep(0, n_random + 1)

africa_europe <- rep(0, n_random + 1)
africa_africa <- rep(0, n_random + 1)
africa_temperateasia <- rep(0, n_random + 1)
africa_tropicalasia <- rep(0, n_random + 1)
africa_australasia <- rep(0, n_random + 1)
africa_pacific <- rep(0, n_random + 1)
africa_northamerica <- rep(0, n_random + 1)
africa_southamerica <- rep(0, n_random + 1)

temperateasia_europe <- rep(0, n_random + 1)
temperateasia_africa <- rep(0, n_random + 1)
temperateasia_temperateasia <- rep(0, n_random + 1)
temperateasia_tropicalasia <- rep(0, n_random + 1)
temperateasia_australasia <- rep(0, n_random + 1)
temperateasia_pacific <- rep(0, n_random + 1)
temperateasia_northamerica <- rep(0, n_random + 1)
temperateasia_southamerica <- rep(0, n_random + 1)

tropicalasia_europe <- rep(0, n_random + 1)
tropicalasia_africa <- rep(0, n_random + 1)
tropicalasia_temperateasia <- rep(0, n_random + 1)
tropicalasia_tropicalasia <- rep(0, n_random + 1)
tropicalasia_australasia <- rep(0, n_random + 1)
tropicalasia_pacific <- rep(0, n_random + 1)
tropicalasia_northamerica <- rep(0, n_random + 1)
tropicalasia_southamerica <- rep(0, n_random + 1)

australasia_europe <- rep(0, n_random + 1)
australasia_africa <- rep(0, n_random + 1)
australasia_temperateasia <- rep(0, n_random + 1)
australasia_tropicalasia <- rep(0, n_random + 1)
australasia_australasia <- rep(0, n_random + 1)
australasia_pacific <- rep(0, n_random + 1)
australasia_northamerica <- rep(0, n_random + 1)
australasia_southamerica <- rep(0, n_random + 1)

pacific_europe <- rep(0, n_random + 1)
pacific_africa <- rep(0, n_random + 1)
pacific_temperateasia <- rep(0, n_random + 1)
pacific_tropicalasia <- rep(0, n_random + 1)
pacific_australasia <- rep(0, n_random + 1)
pacific_pacific <- rep(0, n_random + 1)
pacific_northamerica <- rep(0, n_random + 1)
pacific_southamerica <- rep(0, n_random + 1)

northamerica_europe <- rep(0, n_random + 1)
northamerica_africa <- rep(0, n_random + 1)
northamerica_temperateasia <- rep(0, n_random + 1)
northamerica_tropicalasia <- rep(0, n_random + 1)
northamerica_australasia <- rep(0, n_random + 1)
northamerica_pacific <- rep(0, n_random + 1)
northamerica_northamerica <- rep(0, n_random + 1)
northamerica_southamerica <- rep(0, n_random + 1)

southamerica_europe <- rep(0, n_random + 1)
southamerica_africa <- rep(0, n_random + 1)
southamerica_temperateasia <- rep(0, n_random + 1)
southamerica_tropicalasia <- rep(0, n_random + 1)
southamerica_australasia <- rep(0, n_random + 1)
southamerica_pacific <- rep(0, n_random + 1)
southamerica_northamerica <- rep(0, n_random + 1)
southamerica_southamerica <- rep(0, n_random + 1)

# create an empty data frame from the vectors above
count_sourcereceptor <- as.data.frame(cbind(europe_europe,
                                            europe_africa,
                                            europe_temperateasia,
                                            europe_tropicalasia,
                                            europe_australasia,
                                            europe_pacific,
                                            europe_northamerica,
                                            europe_southamerica,
                                            africa_europe,
                                            africa_africa,
                                            africa_temperateasia,
                                            africa_tropicalasia,
                                            africa_australasia,
                                            africa_pacific,
                                            africa_northamerica,
                                            africa_southamerica,
                                            temperateasia_europe,
                                            temperateasia_africa,
                                            temperateasia_temperateasia,
                                            temperateasia_tropicalasia,
                                            temperateasia_australasia,
                                            temperateasia_pacific,
                                            temperateasia_northamerica,
                                            temperateasia_southamerica,
                                            tropicalasia_europe,
                                            tropicalasia_africa,
                                            tropicalasia_temperateasia,
                                            tropicalasia_tropicalasia,
                                            tropicalasia_australasia,
                                            tropicalasia_pacific,
                                            tropicalasia_northamerica,
                                            tropicalasia_southamerica,
                                            australasia_europe,
                                            australasia_africa,
                                            australasia_temperateasia,
                                            australasia_tropicalasia,
                                            australasia_australasia,
                                            australasia_pacific,
                                            australasia_northamerica,
                                            australasia_southamerica,
                                            pacific_europe,
                                            pacific_africa,
                                            pacific_temperateasia,
                                            pacific_tropicalasia,
                                            pacific_australasia,
                                            pacific_pacific,
                                            pacific_northamerica,
                                            pacific_southamerica,
                                            northamerica_europe,
                                            northamerica_africa,
                                            northamerica_temperateasia,
                                            northamerica_tropicalasia,
                                            northamerica_australasia,
                                            northamerica_pacific,
                                            northamerica_northamerica,
                                            northamerica_southamerica,
                                            southamerica_europe,
                                            southamerica_africa,
                                            southamerica_temperateasia,
                                            southamerica_tropicalasia,
                                            southamerica_australasia,
                                            southamerica_pacific,
                                            southamerica_northamerica,
                                            southamerica_southamerica))

# put observed values in first row of data frame
count_sourcereceptor[1, ] <- observed

for (id in 2: (n_random + 1)) {
  
  print(id)
  
  # draw for each continent separately a number of taxa corresponding to the observed number in that continent
  # make sure the same order as destination
  europe_draw <- sample(native_all, europe_dest)
  africa_draw <- sample(native_all, africa_dest)
  temperateasia_draw <- sample(native_all, temperateasia_dest)
  tropicalasia_draw <- sample(native_all, tropicalasia_dest)
  australasia_draw <- sample(native_all, australasia_dest)
  pacific_draw <- sample(native_all, pacific_dest)
  northamerica_draw <- sample(native_all, northamerica_dest)
  southamerica_draw <- sample(native_all, southamerica_dest)
  
  # extract the migration flow according to the random species
  europe_draw_data <- res_file[which(res_file$spacc %in% europe_draw), ]
  africa_draw_data <- res_file[which(res_file$spacc %in% africa_draw), ]
  temperateasia_draw_data <- res_file[which(res_file$spacc %in% temperateasia_draw), ]
  tropicalasia_draw_data <- res_file[which(res_file$spacc %in% tropicalasia_draw), ]
  australasia_draw_data <- res_file[which(res_file$spacc %in% australasia_draw), ]
  pacific_draw_data <- res_file[which(res_file$spacc %in% pacific_draw), ]
  northamerica_draw_data <- res_file[which(res_file$spacc %in% northamerica_draw), ]
  southamerica_draw_data <- res_file[which(res_file$spacc %in% southamerica_draw), ]
  
  # add the sums for each flow to the columns in the count_sourcereceptor dataframe
  count_sourcereceptor[id, 1] <- length(unique(europe_draw_data$spacc[which(europe_draw_data$origin == 1)]))
  count_sourcereceptor[id, 2] <- length(unique(africa_draw_data$spacc[which(africa_draw_data$origin == 1)]))
  count_sourcereceptor[id, 3] <- length(unique(temperateasia_draw_data$spacc[which(temperateasia_draw_data$origin == 1)]))
  count_sourcereceptor[id, 4] <- length(unique(tropicalasia_draw_data$spacc[which(tropicalasia_draw_data$origin == 1)]))
  count_sourcereceptor[id, 5] <- length(unique(australasia_draw_data$spacc[which(australasia_draw_data$origin == 1)]))
  count_sourcereceptor[id, 6] <- length(unique(pacific_draw_data$spacc[which(pacific_draw_data$origin == 1)]))
  count_sourcereceptor[id, 7] <- length(unique(northamerica_draw_data$spacc[which(northamerica_draw_data$origin == 1)]))
  count_sourcereceptor[id, 8] <- length(unique(southamerica_draw_data$spacc[which(southamerica_draw_data$origin == 1)]))
  
  count_sourcereceptor[id, 9] <- length(unique(europe_draw_data$spacc[which(europe_draw_data$origin == 2)]))
  count_sourcereceptor[id, 10] <- length(unique(africa_draw_data$spacc[which(africa_draw_data$origin == 2)]))
  count_sourcereceptor[id, 11] <- length(unique(temperateasia_draw_data$spacc[which(temperateasia_draw_data$origin == 2)]))
  count_sourcereceptor[id, 12] <- length(unique(tropicalasia_draw_data$spacc[which(tropicalasia_draw_data$origin == 2)]))
  count_sourcereceptor[id, 13] <- length(unique(australasia_draw_data$spacc[which(australasia_draw_data$origin == 2)]))
  count_sourcereceptor[id, 14] <- length(unique(pacific_draw_data$spacc[which(pacific_draw_data$origin == 2)]))
  count_sourcereceptor[id, 15] <- length(unique(northamerica_draw_data$spacc[which(northamerica_draw_data$origin == 2)]))
  count_sourcereceptor[id, 16] <- length(unique(southamerica_draw_data$spacc[which(southamerica_draw_data$origin == 2)]))
  
  count_sourcereceptor[id, 17] <- length(unique(europe_draw_data$spacc[which(europe_draw_data$origin == 3)]))
  count_sourcereceptor[id, 18] <- length(unique(africa_draw_data$spacc[which(africa_draw_data$origin == 3)]))
  count_sourcereceptor[id, 19] <- length(unique(temperateasia_draw_data$spacc[which(temperateasia_draw_data$origin == 3)]))
  count_sourcereceptor[id, 20] <- length(unique(tropicalasia_draw_data$spacc[which(tropicalasia_draw_data$origin == 3)]))
  count_sourcereceptor[id, 21] <- length(unique(australasia_draw_data$spacc[which(australasia_draw_data$origin == 3)]))
  count_sourcereceptor[id, 22] <- length(unique(pacific_draw_data$spacc[which(pacific_draw_data$origin == 3)]))
  count_sourcereceptor[id, 23] <- length(unique(northamerica_draw_data$spacc[which(northamerica_draw_data$origin == 3)]))
  count_sourcereceptor[id, 24] <- length(unique(southamerica_draw_data$spacc[which(southamerica_draw_data$origin == 3)]))
  
  count_sourcereceptor[id, 25] <- length(unique(europe_draw_data$spacc[which(europe_draw_data$origin == 4)]))
  count_sourcereceptor[id, 26] <- length(unique(africa_draw_data$spacc[which(africa_draw_data$origin == 4)]))
  count_sourcereceptor[id, 27] <- length(unique(temperateasia_draw_data$spacc[which(temperateasia_draw_data$origin == 4)]))
  count_sourcereceptor[id, 28] <- length(unique(tropicalasia_draw_data$spacc[which(tropicalasia_draw_data$origin == 4)]))
  count_sourcereceptor[id, 29] <- length(unique(australasia_draw_data$spacc[which(australasia_draw_data$origin == 4)]))
  count_sourcereceptor[id, 30] <- length(unique(pacific_draw_data$spacc[which(pacific_draw_data$origin == 4)]))
  count_sourcereceptor[id, 31] <- length(unique(northamerica_draw_data$spacc[which(northamerica_draw_data$origin == 4)]))
  count_sourcereceptor[id, 32] <- length(unique(southamerica_draw_data$spacc[which(southamerica_draw_data$origin == 4)]))
  
  count_sourcereceptor[id, 33] <- length(unique(europe_draw_data$spacc[which(europe_draw_data$origin == 5)]))
  count_sourcereceptor[id, 34] <- length(unique(africa_draw_data$spacc[which(africa_draw_data$origin == 5)]))
  count_sourcereceptor[id, 35] <- length(unique(temperateasia_draw_data$spacc[which(temperateasia_draw_data$origin == 5)]))
  count_sourcereceptor[id, 36] <- length(unique(tropicalasia_draw_data$spacc[which(tropicalasia_draw_data$origin == 5)]))
  count_sourcereceptor[id, 37] <- length(unique(australasia_draw_data$spacc[which(australasia_draw_data$origin == 5)]))
  count_sourcereceptor[id, 38] <- length(unique(pacific_draw_data$spacc[which(pacific_draw_data$origin == 5)]))
  count_sourcereceptor[id, 39] <- length(unique(northamerica_draw_data$spacc[which(northamerica_draw_data$origin == 5)]))
  count_sourcereceptor[id, 40] <- length(unique(southamerica_draw_data$spacc[which(southamerica_draw_data$origin == 5)]))
  
  count_sourcereceptor[id, 41] <- length(unique(europe_draw_data$spacc[which(europe_draw_data$origin == 6)]))
  count_sourcereceptor[id, 42] <- length(unique(africa_draw_data$spacc[which(africa_draw_data$origin == 6)]))
  count_sourcereceptor[id, 43] <- length(unique(temperateasia_draw_data$spacc[which(temperateasia_draw_data$origin == 6)]))
  count_sourcereceptor[id, 44] <- length(unique(tropicalasia_draw_data$spacc[which(tropicalasia_draw_data$origin == 6)]))
  count_sourcereceptor[id, 45] <- length(unique(australasia_draw_data$spacc[which(australasia_draw_data$origin == 6)]))
  count_sourcereceptor[id, 46] <- length(unique(pacific_draw_data$spacc[which(pacific_draw_data$origin == 6)]))
  count_sourcereceptor[id, 47] <- length(unique(northamerica_draw_data$spacc[which(northamerica_draw_data$origin == 6)]))
  count_sourcereceptor[id, 48] <- length(unique(southamerica_draw_data$spacc[which(southamerica_draw_data$origin == 6)]))
  
  count_sourcereceptor[id, 49] <- length(unique(europe_draw_data$spacc[which(europe_draw_data$origin == 7)]))
  count_sourcereceptor[id, 50] <- length(unique(africa_draw_data$spacc[which(africa_draw_data$origin == 7)]))
  count_sourcereceptor[id, 51] <- length(unique(temperateasia_draw_data$spacc[which(temperateasia_draw_data$origin == 7)]))
  count_sourcereceptor[id, 52] <- length(unique(tropicalasia_draw_data$spacc[which(tropicalasia_draw_data$origin == 7)]))
  count_sourcereceptor[id, 53] <- length(unique(australasia_draw_data$spacc[which(australasia_draw_data$origin == 7)]))
  count_sourcereceptor[id, 54] <- length(unique(pacific_draw_data$spacc[which(pacific_draw_data$origin == 7)]))
  count_sourcereceptor[id, 55] <- length(unique(northamerica_draw_data$spacc[which(northamerica_draw_data$origin == 7)]))
  count_sourcereceptor[id, 56] <- length(unique(southamerica_draw_data$spacc[which(southamerica_draw_data$origin == 7)]))
  
  count_sourcereceptor[id, 57] <- length(unique(europe_draw_data$spacc[which(europe_draw_data$origin == 8)]))
  count_sourcereceptor[id, 58] <- length(unique(africa_draw_data$spacc[which(africa_draw_data$origin == 8)]))
  count_sourcereceptor[id, 59] <- length(unique(temperateasia_draw_data$spacc[which(temperateasia_draw_data$origin == 8)]))
  count_sourcereceptor[id, 60] <- length(unique(tropicalasia_draw_data$spacc[which(tropicalasia_draw_data$origin == 8)]))
  count_sourcereceptor[id, 61] <- length(unique(australasia_draw_data$spacc[which(australasia_draw_data$origin == 8)]))
  count_sourcereceptor[id, 62] <- length(unique(pacific_draw_data$spacc[which(pacific_draw_data$origin == 8)]))
  count_sourcereceptor[id, 63] <- length(unique(northamerica_draw_data$spacc[which(northamerica_draw_data$origin == 8)]))
  count_sourcereceptor[id, 64] <- length(unique(southamerica_draw_data$spacc[which(southamerica_draw_data$origin == 8)]))
  
}

write.table(count_sourcereceptor, file = paste0(path_folder, "03_migration_flow/random_sample_count_sourcereceptor.txt"),
            col.names = TRUE, row.names = FALSE)

count_sourcereceptor <- read.table(paste0(path_folder, "03_migration_flow/random_sample_count_sourcereceptor.txt"),
                                        header = TRUE)

# create dataframe with all combinations, observed values and statistics from the randomization test
source_receptor <- c("europe_europe",
                     "europe_africa",
                     "europe_temperateasia",
                     "europe_tropicalasia",
                     "europe_australasia",
                     "europe_pacific",
                     "europe_northamerica",
                     "europe_southamerica",
                     "africa_europe",
                     "africa_africa",
                     "africa_temperateasia",
                     "africa_tropicalasia",
                     "africa_australasia",
                     "africa_pacific",
                     "africa_northamerica",
                     "africa_southamerica",
                     "temperateasia_europe",
                     "temperateasia_africa",
                     "temperateasia_temperateasia",
                     "temperateasia_tropicalasia",
                     "temperateasia_australasia",
                     "temperateasia_pacific",
                     "temperateasia_northamerica",
                     "temperateasia_southamerica",
                     "tropicalasia_europe",
                     "tropicalasia_africa",
                     "tropicalasia_temperateasia",
                     "tropicalasia_tropicalasia",
                     "tropicalasia_australasia",
                     "tropicalasia_pacific",
                     "tropicalasia_northamerica",
                     "tropicalasia_southamerica",
                     "australasia_europe",
                     "australasia_africa",
                     "australasia_temperateasia",
                     "australasia_tropicalasia",
                     "australasia_australasia",
                     "australasia_pacific",
                     "australasia_northamerica",
                     "australasia_southamerica",
                     "pacific_europe",
                     "pacific_africa",
                     "pacific_temperateasia",
                     "pacific_tropicalasia",
                     "pacific_australasia",
                     "pacific_pacific",
                     "pacific_northamerica",
                     "pacific_southamerica",
                     "northamerica_europe",
                     "northamerica_africa",
                     "northamerica_temperateasia",
                     "northamerica_tropicalasia",
                     "northamerica_australasia",
                     "northamerica_pacific",
                     "northamerica_northamerica",
                     "northamerica_southamerica",
                     "southamerica_europe",
                     "southamerica_africa",
                     "southamerica_temperateasia",
                     "southamerica_tropicalasia",
                     "southamerica_australasia",
                     "southamerica_pacific",
                     "southamerica_northamerica",
                     "southamerica_southamerica")

random_res <- as.data.frame(cbind(source_receptor, observed))
# create empty column in random_res for the proportion of random values that is smaller than the observed value
random_res$proportion <- rep(NA, length(random_res$observed))
# create empty column in random_res for mean value of the random draws
random_res$mean <- rep(NA, length(random_res$observed))
# create empty column in random_res for median of the random draws
random_res$median <- rep(NA, length(random_res$observed))

# calculate whether observed flows are higher (proportion close to 1) or lower (proportion close to 0) than expected
for (id in 1: dim(random_res)[1]) {
  
  random_data <- count_sourcereceptor[-1, id]
  observed_data <- count_sourcereceptor[1, id]
  
  random_res$proportion[id] <- length(random_data[random_data < observed_data])/9999
  random_res$mean[id] <- mean(count_sourcereceptor[-1, id])
  random_res$median[id] <- median(count_sourcereceptor[-1, id])
  
  print(id)
  
}

# define line colors for observed-value lines in histograms
random_res$line_color <- rep("black", dim(random_res)[1]) 
random_res$line_color <- replace(random_res$line_color, random_res$proportion > 0.975, "red")
random_res$line_color <- replace(random_res$line_color, random_res$proportion < 0.025, "blue")

###########################
# plot the expected flows
# form the dataframe
for (id in 1: dim(random_res)[1]) {
  
  sour_rece_name <- random_res[id, 1]
  
  origin_name <- strsplit(sour_rece_name, "_", fixed = TRUE)[[1]][1]
  destination_name <- strsplit(sour_rece_name, "_", fixed = TRUE)[[1]][2]
  
  sour_rece_flow <- random_res[id, 5]
  sour_rece_res <- cbind(origin_name, destination_name, sour_rece_flow)
  write.table(sour_rece_res, file = paste0(path_folder, "03_migration_flow/random_sample_migration_flow_median.txt"),
              col.names = FALSE, row.names = FALSE, append = TRUE)
  
  print(id)
  
}

flow_data_random <- read.table(paste0(path_folder, "03_migration_flow/random_sample_migration_flow_median.txt"))
colnames(flow_data_random) <- c("origin", "destination", "flow")
write.table(flow_data_random, file = paste0(path_folder, "03_migration_flow/random_sample_migration_flow_median.txt"),
            col.names = TRUE, row.names = FALSE)

# manually modify the name of continents
flow_data_random <- read.table(paste0(path_folder, "03_migration_flow/random_sample_migration_flow_median.txt"),
                               header = TRUE)

# select the top 50% of migration flows
flow_data_random <- flow_data_random[order(flow_data_random$flow), ]
flow_data_random <- flow_data_random[33: 64, ]

# plot
setwd(path_folder)
png(file = paste0(path_folder, "03_migration_flow/migrationflow_nonnativetree_randomsamplemedian.png"),
    bg = "transparent", width = 8, height = 8, units = "cm", res = 300)

circos.clear()

circos.par(start.degree = 90,
           gap.after = 5,
           track.margin = c(0.01, 0.01))

chordDiagram(x = flow_data_random, 
             order = c("South America",
                       "North America",
                       "Pacific",
                       "Tropical Asia",
                       "Australasia",
                       "Temperate Asia",
                       "Africa",
                       "Europe"),
             grid.col = c("South America" = "#5D74A5",
                          "North America" = "#8CA5CA",
                          "Pacific" = "#BBD1E2",
                          "Tropical Asia" = "#E7EAD0",
                          "Australasia" = "#F8DEB2",
                          "Temperate Asia" = "#EDAC88",
                          "Africa" = "#CE7F69",
                          "Europe" = "#A8554E"),
             directional = 1, 
             annotationTrack = c("grid"),
             transparency = 0,
             annotationTrackHeight = c(0.06, 0.07),
             preAllocateTracks = list(track.height = mm_h(0.1)),
             direction.type = "diffHeight",
             diffHeight  = mm_h(4),
             link.target.prop = FALSE)

# add labels and axis
circos.track(bg.border = NA, ylim = c(0,1),
             panel.fun = function(x, y) {
               xcenter = get.cell.meta.data("xcenter")
               sector.index = get.cell.meta.data("sector.index")
               circos.text(x = xcenter, y = 2.1, labels = sector.index, 
                           facing = "bending", cex = 0.7,font = 1, family = "Times New Roman")
               circos.axis(h = "top", labels = FALSE, track.index = 2,
                           minor.ticks = 0,major.at = seq(0, 35000, 1000),
                           lwd = 1,major.tick.length = mm_y(0.5))
             })

dev.off()

#############################################
# threatened non-native tree migration flow
#############################################

##############################
# threatened non-native tree
# load the migration flow file
tree_flow <- read.table(paste0(path_folder, "03_migration_flow/tree_flow_all.txt"), header = TRUE)
threatenedtree_flow <- tree_flow[which(tree_flow$spacc %in% nonnative_threatened_sppre), ]
# load the observed flow file
threatenedflow_data <- read.table(paste0(path_folder, "03_migration_flow/migration_flow_among_continent_threatened.txt"), header = TRUE)

# calculate flows for each source-recipient combination
observed_threatened <- threatenedflow_data$flow

# the total number of threatened non-native tree species in each continent
europe_dest_threatened <- length(unique(threatenedtree_flow$spacc[which(threatenedtree_flow$destination == 1)])) 
africa_dest_threatened <- length(unique(threatenedtree_flow$spacc[which(threatenedtree_flow$destination == 2)])) 
temperateasia_dest_threatened <- length(unique(threatenedtree_flow$spacc[which(threatenedtree_flow$destination == 3)])) 
tropicalasia_dest_threatened <- length(unique(threatenedtree_flow$spacc[which(threatenedtree_flow$destination == 4)])) 
australasia_dest_threatened <- length(unique(threatenedtree_flow$spacc[which(threatenedtree_flow$destination == 5)])) 
pacific_dest_threatened <- length(unique(threatenedtree_flow$spacc[which(threatenedtree_flow$destination == 6)])) 
northamerica_dest_threatened <- length(unique(threatenedtree_flow$spacc[which(threatenedtree_flow$destination == 7)])) 
southamerica_dest_threatened <- length(unique(threatenedtree_flow$spacc[which(threatenedtree_flow$destination == 8)])) 

# number of random draws
n_random <- 9999

# create vectors for each source-receptor combination
europe_europe <- rep(0, n_random + 1)
europe_africa <- rep(0, n_random + 1)
europe_temperateasia <- rep(0, n_random + 1)
europe_tropicalasia <- rep(0, n_random + 1)
europe_australasia <- rep(0, n_random + 1)
europe_pacific <- rep(0, n_random + 1)
europe_northamerica <- rep(0, n_random + 1)
europe_southamerica <- rep(0, n_random + 1)

africa_europe <- rep(0, n_random + 1)
africa_africa <- rep(0, n_random + 1)
africa_temperateasia <- rep(0, n_random + 1)
africa_tropicalasia <- rep(0, n_random + 1)
africa_australasia <- rep(0, n_random + 1)
africa_pacific <- rep(0, n_random + 1)
africa_northamerica <- rep(0, n_random + 1)
africa_southamerica <- rep(0, n_random + 1)

temperateasia_europe <- rep(0, n_random + 1)
temperateasia_africa <- rep(0, n_random + 1)
temperateasia_temperateasia <- rep(0, n_random + 1)
temperateasia_tropicalasia <- rep(0, n_random + 1)
temperateasia_australasia <- rep(0, n_random + 1)
temperateasia_pacific <- rep(0, n_random + 1)
temperateasia_northamerica <- rep(0, n_random + 1)
temperateasia_southamerica <- rep(0, n_random + 1)

tropicalasia_europe <- rep(0, n_random + 1)
tropicalasia_africa <- rep(0, n_random + 1)
tropicalasia_temperateasia <- rep(0, n_random + 1)
tropicalasia_tropicalasia <- rep(0, n_random + 1)
tropicalasia_australasia <- rep(0, n_random + 1)
tropicalasia_pacific <- rep(0, n_random + 1)
tropicalasia_northamerica <- rep(0, n_random + 1)
tropicalasia_southamerica <- rep(0, n_random + 1)

australasia_europe <- rep(0, n_random + 1)
australasia_africa <- rep(0, n_random + 1)
australasia_temperateasia <- rep(0, n_random + 1)
australasia_tropicalasia <- rep(0, n_random + 1)
australasia_australasia <- rep(0, n_random + 1)
australasia_pacific <- rep(0, n_random + 1)
australasia_northamerica <- rep(0, n_random + 1)
australasia_southamerica <- rep(0, n_random + 1)

pacific_europe <- rep(0, n_random + 1)
pacific_africa <- rep(0, n_random + 1)
pacific_temperateasia <- rep(0, n_random + 1)
pacific_tropicalasia <- rep(0, n_random + 1)
pacific_australasia <- rep(0, n_random + 1)
pacific_pacific <- rep(0, n_random + 1)
pacific_northamerica <- rep(0, n_random + 1)
pacific_southamerica <- rep(0, n_random + 1)

northamerica_europe <- rep(0, n_random + 1)
northamerica_africa <- rep(0, n_random + 1)
northamerica_temperateasia <- rep(0, n_random + 1)
northamerica_tropicalasia <- rep(0, n_random + 1)
northamerica_australasia <- rep(0, n_random + 1)
northamerica_pacific <- rep(0, n_random + 1)
northamerica_northamerica <- rep(0, n_random + 1)
northamerica_southamerica <- rep(0, n_random + 1)

southamerica_europe <- rep(0, n_random + 1)
southamerica_africa <- rep(0, n_random + 1)
southamerica_temperateasia <- rep(0, n_random + 1)
southamerica_tropicalasia <- rep(0, n_random + 1)
southamerica_australasia <- rep(0, n_random + 1)
southamerica_pacific <- rep(0, n_random + 1)
southamerica_northamerica <- rep(0, n_random + 1)
southamerica_southamerica <- rep(0, n_random + 1)

# create an empty data frame from the vectors above
count_sourcereceptor_threatened <- as.data.frame(cbind(europe_europe,
                                                       europe_africa,
                                                       europe_temperateasia,
                                                       europe_tropicalasia,
                                                       europe_australasia,
                                                       europe_pacific,
                                                       europe_northamerica,
                                                       europe_southamerica,
                                                       africa_europe,
                                                       africa_africa,
                                                       africa_temperateasia,
                                                       africa_tropicalasia,
                                                       africa_australasia,
                                                       africa_pacific,
                                                       africa_northamerica,
                                                       africa_southamerica,
                                                       temperateasia_europe,
                                                       temperateasia_africa,
                                                       temperateasia_temperateasia,
                                                       temperateasia_tropicalasia,
                                                       temperateasia_australasia,
                                                       temperateasia_pacific,
                                                       temperateasia_northamerica,
                                                       temperateasia_southamerica,
                                                       tropicalasia_europe,
                                                       tropicalasia_africa,
                                                       tropicalasia_temperateasia,
                                                       tropicalasia_tropicalasia,
                                                       tropicalasia_australasia,
                                                       tropicalasia_pacific,
                                                       tropicalasia_northamerica,
                                                       tropicalasia_southamerica,
                                                       australasia_europe,
                                                       australasia_africa,
                                                       australasia_temperateasia,
                                                       australasia_tropicalasia,
                                                       australasia_australasia,
                                                       australasia_pacific,
                                                       australasia_northamerica,
                                                       australasia_southamerica,
                                                       pacific_europe,
                                                       pacific_africa,
                                                       pacific_temperateasia,
                                                       pacific_tropicalasia,
                                                       pacific_australasia,
                                                       pacific_pacific,
                                                       pacific_northamerica,
                                                       pacific_southamerica,
                                                       northamerica_europe,
                                                       northamerica_africa,
                                                       northamerica_temperateasia,
                                                       northamerica_tropicalasia,
                                                       northamerica_australasia,
                                                       northamerica_pacific,
                                                       northamerica_northamerica,
                                                       northamerica_southamerica,
                                                       southamerica_europe,
                                                       southamerica_africa,
                                                       southamerica_temperateasia,
                                                       southamerica_tropicalasia,
                                                       southamerica_australasia,
                                                       southamerica_pacific,
                                                       southamerica_northamerica,
                                                       southamerica_southamerica))

# put observed values in first row of data frame
count_sourcereceptor_threatened[1, ] <- observed_threatened

for (id in 2: (n_random + 1)) {
  
  print(id)
  
  # draw for each continent separately a number of threatened taxa corresponding to the observed threatened number in that continent
  # make sure the same order as destination
  europe_draw <- sample(native_threatened_sppre, europe_dest_threatened)
  africa_draw <- sample(native_threatened_sppre, africa_dest_threatened)
  temperateasia_draw <- sample(native_threatened_sppre, temperateasia_dest_threatened)
  tropicalasia_draw <- sample(native_threatened_sppre, tropicalasia_dest_threatened)
  australasia_draw <- sample(native_threatened_sppre, australasia_dest_threatened)
  pacific_draw <- sample(native_threatened_sppre, pacific_dest_threatened)
  northamerica_draw <- sample(native_threatened_sppre, northamerica_dest_threatened)
  southamerica_draw <- sample(native_threatened_sppre, southamerica_dest_threatened)
  
  # extract the migration flow according to the random species
  europe_draw_data <- res_file[which(res_file$spacc %in% europe_draw), ]
  africa_draw_data <- res_file[which(res_file$spacc %in% africa_draw), ]
  temperateasia_draw_data <- res_file[which(res_file$spacc %in% temperateasia_draw), ]
  tropicalasia_draw_data <- res_file[which(res_file$spacc %in% tropicalasia_draw), ]
  australasia_draw_data <- res_file[which(res_file$spacc %in% australasia_draw), ]
  pacific_draw_data <- res_file[which(res_file$spacc %in% pacific_draw), ]
  northamerica_draw_data <- res_file[which(res_file$spacc %in% northamerica_draw), ]
  southamerica_draw_data <- res_file[which(res_file$spacc %in% southamerica_draw), ]
  
  # add the sums for each flow to the columns in the count_sourcereceptor dataframe
  count_sourcereceptor_threatened[id, 1] <- length(unique(europe_draw_data$spacc[which(europe_draw_data$origin == 1)]))
  count_sourcereceptor_threatened[id, 2] <- length(unique(africa_draw_data$spacc[which(africa_draw_data$origin == 1)]))
  count_sourcereceptor_threatened[id, 3] <- length(unique(temperateasia_draw_data$spacc[which(temperateasia_draw_data$origin == 1)]))
  count_sourcereceptor_threatened[id, 4] <- length(unique(tropicalasia_draw_data$spacc[which(tropicalasia_draw_data$origin == 1)]))
  count_sourcereceptor_threatened[id, 5] <- length(unique(australasia_draw_data$spacc[which(australasia_draw_data$origin == 1)]))
  count_sourcereceptor_threatened[id, 6] <- length(unique(pacific_draw_data$spacc[which(pacific_draw_data$origin == 1)]))
  count_sourcereceptor_threatened[id, 7] <- length(unique(northamerica_draw_data$spacc[which(northamerica_draw_data$origin == 1)]))
  count_sourcereceptor_threatened[id, 8] <- length(unique(southamerica_draw_data$spacc[which(southamerica_draw_data$origin == 1)]))
  
  count_sourcereceptor_threatened[id, 9] <- length(unique(europe_draw_data$spacc[which(europe_draw_data$origin == 2)]))
  count_sourcereceptor_threatened[id, 10] <- length(unique(africa_draw_data$spacc[which(africa_draw_data$origin == 2)]))
  count_sourcereceptor_threatened[id, 11] <- length(unique(temperateasia_draw_data$spacc[which(temperateasia_draw_data$origin == 2)]))
  count_sourcereceptor_threatened[id, 12] <- length(unique(tropicalasia_draw_data$spacc[which(tropicalasia_draw_data$origin == 2)]))
  count_sourcereceptor_threatened[id, 13] <- length(unique(australasia_draw_data$spacc[which(australasia_draw_data$origin == 2)]))
  count_sourcereceptor_threatened[id, 14] <- length(unique(pacific_draw_data$spacc[which(pacific_draw_data$origin == 2)]))
  count_sourcereceptor_threatened[id, 15] <- length(unique(northamerica_draw_data$spacc[which(northamerica_draw_data$origin == 2)]))
  count_sourcereceptor_threatened[id, 16] <- length(unique(southamerica_draw_data$spacc[which(southamerica_draw_data$origin == 2)]))
  
  count_sourcereceptor_threatened[id, 17] <- length(unique(europe_draw_data$spacc[which(europe_draw_data$origin == 3)]))
  count_sourcereceptor_threatened[id, 18] <- length(unique(africa_draw_data$spacc[which(africa_draw_data$origin == 3)]))
  count_sourcereceptor_threatened[id, 19] <- length(unique(temperateasia_draw_data$spacc[which(temperateasia_draw_data$origin == 3)]))
  count_sourcereceptor_threatened[id, 20] <- length(unique(tropicalasia_draw_data$spacc[which(tropicalasia_draw_data$origin == 3)]))
  count_sourcereceptor_threatened[id, 21] <- length(unique(australasia_draw_data$spacc[which(australasia_draw_data$origin == 3)]))
  count_sourcereceptor_threatened[id, 22] <- length(unique(pacific_draw_data$spacc[which(pacific_draw_data$origin == 3)]))
  count_sourcereceptor_threatened[id, 23] <- length(unique(northamerica_draw_data$spacc[which(northamerica_draw_data$origin == 3)]))
  count_sourcereceptor_threatened[id, 24] <- length(unique(southamerica_draw_data$spacc[which(southamerica_draw_data$origin == 3)]))
  
  count_sourcereceptor_threatened[id, 25] <- length(unique(europe_draw_data$spacc[which(europe_draw_data$origin == 4)]))
  count_sourcereceptor_threatened[id, 26] <- length(unique(africa_draw_data$spacc[which(africa_draw_data$origin == 4)]))
  count_sourcereceptor_threatened[id, 27] <- length(unique(temperateasia_draw_data$spacc[which(temperateasia_draw_data$origin == 4)]))
  count_sourcereceptor_threatened[id, 28] <- length(unique(tropicalasia_draw_data$spacc[which(tropicalasia_draw_data$origin == 4)]))
  count_sourcereceptor_threatened[id, 29] <- length(unique(australasia_draw_data$spacc[which(australasia_draw_data$origin == 4)]))
  count_sourcereceptor_threatened[id, 30] <- length(unique(pacific_draw_data$spacc[which(pacific_draw_data$origin == 4)]))
  count_sourcereceptor_threatened[id, 31] <- length(unique(northamerica_draw_data$spacc[which(northamerica_draw_data$origin == 4)]))
  count_sourcereceptor_threatened[id, 32] <- length(unique(southamerica_draw_data$spacc[which(southamerica_draw_data$origin == 4)]))
  
  count_sourcereceptor_threatened[id, 33] <- length(unique(europe_draw_data$spacc[which(europe_draw_data$origin == 5)]))
  count_sourcereceptor_threatened[id, 34] <- length(unique(africa_draw_data$spacc[which(africa_draw_data$origin == 5)]))
  count_sourcereceptor_threatened[id, 35] <- length(unique(temperateasia_draw_data$spacc[which(temperateasia_draw_data$origin == 5)]))
  count_sourcereceptor_threatened[id, 36] <- length(unique(tropicalasia_draw_data$spacc[which(tropicalasia_draw_data$origin == 5)]))
  count_sourcereceptor_threatened[id, 37] <- length(unique(australasia_draw_data$spacc[which(australasia_draw_data$origin == 5)]))
  count_sourcereceptor_threatened[id, 38] <- length(unique(pacific_draw_data$spacc[which(pacific_draw_data$origin == 5)]))
  count_sourcereceptor_threatened[id, 39] <- length(unique(northamerica_draw_data$spacc[which(northamerica_draw_data$origin == 5)]))
  count_sourcereceptor_threatened[id, 40] <- length(unique(southamerica_draw_data$spacc[which(southamerica_draw_data$origin == 5)]))
  
  count_sourcereceptor_threatened[id, 41] <- length(unique(europe_draw_data$spacc[which(europe_draw_data$origin == 6)]))
  count_sourcereceptor_threatened[id, 42] <- length(unique(africa_draw_data$spacc[which(africa_draw_data$origin == 6)]))
  count_sourcereceptor_threatened[id, 43] <- length(unique(temperateasia_draw_data$spacc[which(temperateasia_draw_data$origin == 6)]))
  count_sourcereceptor_threatened[id, 44] <- length(unique(tropicalasia_draw_data$spacc[which(tropicalasia_draw_data$origin == 6)]))
  count_sourcereceptor_threatened[id, 45] <- length(unique(australasia_draw_data$spacc[which(australasia_draw_data$origin == 6)]))
  count_sourcereceptor_threatened[id, 46] <- length(unique(pacific_draw_data$spacc[which(pacific_draw_data$origin == 6)]))
  count_sourcereceptor_threatened[id, 47] <- length(unique(northamerica_draw_data$spacc[which(northamerica_draw_data$origin == 6)]))
  count_sourcereceptor_threatened[id, 48] <- length(unique(southamerica_draw_data$spacc[which(southamerica_draw_data$origin == 6)]))
  
  count_sourcereceptor_threatened[id, 49] <- length(unique(europe_draw_data$spacc[which(europe_draw_data$origin == 7)]))
  count_sourcereceptor_threatened[id, 50] <- length(unique(africa_draw_data$spacc[which(africa_draw_data$origin == 7)]))
  count_sourcereceptor_threatened[id, 51] <- length(unique(temperateasia_draw_data$spacc[which(temperateasia_draw_data$origin == 7)]))
  count_sourcereceptor_threatened[id, 52] <- length(unique(tropicalasia_draw_data$spacc[which(tropicalasia_draw_data$origin == 7)]))
  count_sourcereceptor_threatened[id, 53] <- length(unique(australasia_draw_data$spacc[which(australasia_draw_data$origin == 7)]))
  count_sourcereceptor_threatened[id, 54] <- length(unique(pacific_draw_data$spacc[which(pacific_draw_data$origin == 7)]))
  count_sourcereceptor_threatened[id, 55] <- length(unique(northamerica_draw_data$spacc[which(northamerica_draw_data$origin == 7)]))
  count_sourcereceptor_threatened[id, 56] <- length(unique(southamerica_draw_data$spacc[which(southamerica_draw_data$origin == 7)]))
  
  count_sourcereceptor_threatened[id, 57] <- length(unique(europe_draw_data$spacc[which(europe_draw_data$origin == 8)]))
  count_sourcereceptor_threatened[id, 58] <- length(unique(africa_draw_data$spacc[which(africa_draw_data$origin == 8)]))
  count_sourcereceptor_threatened[id, 59] <- length(unique(temperateasia_draw_data$spacc[which(temperateasia_draw_data$origin == 8)]))
  count_sourcereceptor_threatened[id, 60] <- length(unique(tropicalasia_draw_data$spacc[which(tropicalasia_draw_data$origin == 8)]))
  count_sourcereceptor_threatened[id, 61] <- length(unique(australasia_draw_data$spacc[which(australasia_draw_data$origin == 8)]))
  count_sourcereceptor_threatened[id, 62] <- length(unique(pacific_draw_data$spacc[which(pacific_draw_data$origin == 8)]))
  count_sourcereceptor_threatened[id, 63] <- length(unique(northamerica_draw_data$spacc[which(northamerica_draw_data$origin == 8)]))
  count_sourcereceptor_threatened[id, 64] <- length(unique(southamerica_draw_data$spacc[which(southamerica_draw_data$origin == 8)]))
  
}

write.table(count_sourcereceptor_threatened, file = paste0(path_folder, "03_migration_flow/random_sample_count_sourcereceptor_threatened.txt"),
            col.names = TRUE, row.names = FALSE)

count_sourcereceptor_threatened <- read.table(paste0(path_folder, "03_migration_flow/random_sample_count_sourcereceptor_threatened.txt"),
                                              header = TRUE)

# create dataframe with all combinations, observed values and statistics from the randomization test
source_receptor <- c("europe_europe",
                     "europe_africa",
                     "europe_temperateasia",
                     "europe_tropicalasia",
                     "europe_australasia",
                     "europe_pacific",
                     "europe_northamerica",
                     "europe_southamerica",
                     "africa_europe",
                     "africa_africa",
                     "africa_temperateasia",
                     "africa_tropicalasia",
                     "africa_australasia",
                     "africa_pacific",
                     "africa_northamerica",
                     "africa_southamerica",
                     "temperateasia_europe",
                     "temperateasia_africa",
                     "temperateasia_temperateasia",
                     "temperateasia_tropicalasia",
                     "temperateasia_australasia",
                     "temperateasia_pacific",
                     "temperateasia_northamerica",
                     "temperateasia_southamerica",
                     "tropicalasia_europe",
                     "tropicalasia_africa",
                     "tropicalasia_temperateasia",
                     "tropicalasia_tropicalasia",
                     "tropicalasia_australasia",
                     "tropicalasia_pacific",
                     "tropicalasia_northamerica",
                     "tropicalasia_southamerica",
                     "australasia_europe",
                     "australasia_africa",
                     "australasia_temperateasia",
                     "australasia_tropicalasia",
                     "australasia_australasia",
                     "australasia_pacific",
                     "australasia_northamerica",
                     "australasia_southamerica",
                     "pacific_europe",
                     "pacific_africa",
                     "pacific_temperateasia",
                     "pacific_tropicalasia",
                     "pacific_australasia",
                     "pacific_pacific",
                     "pacific_northamerica",
                     "pacific_southamerica",
                     "northamerica_europe",
                     "northamerica_africa",
                     "northamerica_temperateasia",
                     "northamerica_tropicalasia",
                     "northamerica_australasia",
                     "northamerica_pacific",
                     "northamerica_northamerica",
                     "northamerica_southamerica",
                     "southamerica_europe",
                     "southamerica_africa",
                     "southamerica_temperateasia",
                     "southamerica_tropicalasia",
                     "southamerica_australasia",
                     "southamerica_pacific",
                     "southamerica_northamerica",
                     "southamerica_southamerica")

random_res_threatened <- as.data.frame(cbind(source_receptor, observed_threatened))
# create empty column in random_res for the proportion of random values that is smaller than the observed value
random_res_threatened$proportion <- rep(NA, length(random_res_threatened$observed_threatened))
# create empty column in random_res for mean value of the random draws
random_res_threatened$mean <- rep(NA, length(random_res_threatened$observed_threatened))
# create empty column in random_res for median of the random draws
random_res_threatened$median <- rep(NA, length(random_res_threatened$observed_threatened))

# calculate whether observed flows are higher (proportion close to 1) or lower (proportion close to 0) than expected
for (id in 1: dim(random_res_threatened)[1]) {
  
  random_data_threatened <- count_sourcereceptor_threatened[-1, id]
  observed_data_threatened <- count_sourcereceptor_threatened[1, id]
  
  random_res_threatened$proportion[id] <- length(random_data_threatened[random_data_threatened < observed_data_threatened])/9999
  random_res_threatened$mean[id] <- mean(count_sourcereceptor_threatened[-1, id])
  random_res_threatened$median[id] <- median(count_sourcereceptor_threatened[-1, id])
  
  print(id)
  
}

# define line colors for observed-value lines in histograms
random_res_threatened$line_color <- rep("black", dim(random_res_threatened)[1]) 
random_res_threatened$line_color <- replace(random_res_threatened$line_color, random_res_threatened$proportion > 0.975, "red")
random_res_threatened$line_color <- replace(random_res_threatened$line_color, random_res_threatened$proportion < 0.025, "blue")

###########################
# plot the expected flows
# form the dataframe
for (id in 1: dim(random_res_threatened)[1]) {
  
  sour_rece_name <- random_res_threatened[id, 1]
  
  origin_name <- strsplit(sour_rece_name, "_", fixed = TRUE)[[1]][1]
  destination_name <- strsplit(sour_rece_name, "_", fixed = TRUE)[[1]][2]
  
  sour_rece_flow <- random_res_threatened[id, 5]
  sour_rece_res <- cbind(origin_name, destination_name, sour_rece_flow)
  write.table(sour_rece_res, file = paste0(path_folder, "03_migration_flow/random_sample_migration_flow_median_threatened.txt"),
              col.names = FALSE, row.names = FALSE, append = TRUE)
  
  print(id)
  
}

flow_data_random_threatened <- read.table(paste0(path_folder, "03_migration_flow/random_sample_migration_flow_median_threatened.txt"))
colnames(flow_data_random_threatened) <- c("origin", "destination", "flow")
write.table(flow_data_random_threatened, file = paste0(path_folder, "03_migration_flow/random_sample_migration_flow_median_threatened.txt"),
            col.names = TRUE, row.names = FALSE)

# manually modify the name of continents
flow_data_random_threatened <- read.table(paste0(path_folder, "03_migration_flow/random_sample_migration_flow_median_threatened.txt"),
                                          header = TRUE)

# select the top 50% of migration flows
flow_data_random_threatened <- flow_data_random_threatened[order(flow_data_random_threatened$flow), ]
flow_data_random_threatened <- flow_data_random_threatened[33: 64, ]

# plot
setwd(path_folder)
png(file = paste0(path_folder, "03_migration_flow/migrationflow_nonnativetree_randomsamplemedian_threatened.png"),
    bg = "transparent", width = 8, height = 8, units = "cm", res = 300)

circos.clear()

circos.par(start.degree = 90,
           gap.after = 5,
           track.margin = c(0.01, 0.01))

chordDiagram(x = flow_data_random_threatened, 
             order = c("South America",
                       "North America",
                       "Pacific",
                       "Tropical Asia",
                       "Australasia",
                       "Temperate Asia",
                       "Africa",
                       "Europe"),
             grid.col = c("South America" = "#5D74A5",
                          "North America" = "#8CA5CA",
                          "Pacific" = "#BBD1E2",
                          "Tropical Asia" = "#E7EAD0",
                          "Australasia" = "#F8DEB2",
                          "Temperate Asia" = "#EDAC88",
                          "Africa" = "#CE7F69",
                          "Europe" = "#A8554E"),
             directional = 1, 
             annotationTrack = c("grid"),
             transparency = 0,
             annotationTrackHeight = c(0.06, 0.07),
             preAllocateTracks = list(track.height = mm_h(0.1)),
             direction.type = "diffHeight",
             diffHeight  = mm_h(4),
             link.target.prop = FALSE)

# add labels and axis
circos.track(bg.border = NA, ylim = c(0,1),
             panel.fun = function(x, y) {
               xcenter = get.cell.meta.data("xcenter")
               sector.index = get.cell.meta.data("sector.index")
               circos.text(x = xcenter, y = 2.1, labels = sector.index, 
                           facing = "bending", cex = 0.7,font = 1, family = "Times New Roman")
               circos.axis(h = "top", labels = FALSE, track.index = 2,
                           minor.ticks = 0,major.at = seq(0, 10000, 100),
                           lwd = 1,major.tick.length = mm_y(0.5))
             })

dev.off()

#################################################
# non-threatened non-native tree migration flow
#################################################

##################################
# non-threatened non-native tree
# load the migration flow file
tree_flow <- read.table(paste0(path_folder, "03_migration_flow/tree_flow_all.txt"), header = TRUE)
nonthreatenedtree_flow <- tree_flow[which(tree_flow$spacc %in% nonnative_nonthreatened_sppre), ]
# load the observed flow file
nonthreatenedflow_data <- read.table(paste0(path_folder, "03_migration_flow/migration_flow_among_continent_nonthreatened.txt"), header = TRUE)

# calculate flows for each source-recipient combination
observed_nonthreatened <- nonthreatenedflow_data$flow

# the total number of non-threatened non-native tree species in each continent
europe_dest_nonthreatened <- length(unique(nonthreatenedtree_flow$spacc[which(nonthreatenedtree_flow$destination == 1)])) 
africa_dest_nonthreatened <- length(unique(nonthreatenedtree_flow$spacc[which(nonthreatenedtree_flow$destination == 2)])) 
temperateasia_dest_nonthreatened <- length(unique(nonthreatenedtree_flow$spacc[which(nonthreatenedtree_flow$destination == 3)])) 
tropicalasia_dest_nonthreatened <- length(unique(nonthreatenedtree_flow$spacc[which(nonthreatenedtree_flow$destination == 4)])) 
australasia_dest_nonthreatened <- length(unique(nonthreatenedtree_flow$spacc[which(nonthreatenedtree_flow$destination == 5)])) 
pacific_dest_nonthreatened <- length(unique(nonthreatenedtree_flow$spacc[which(nonthreatenedtree_flow$destination == 6)])) 
northamerica_dest_nonthreatened <- length(unique(nonthreatenedtree_flow$spacc[which(nonthreatenedtree_flow$destination == 7)])) 
southamerica_dest_nonthreatened <- length(unique(nonthreatenedtree_flow$spacc[which(nonthreatenedtree_flow$destination == 8)])) 

# number of random draws
n_random <- 9999

# create vectors for each source-receptor combination
europe_europe <- rep(0, n_random + 1)
europe_africa <- rep(0, n_random + 1)
europe_temperateasia <- rep(0, n_random + 1)
europe_tropicalasia <- rep(0, n_random + 1)
europe_australasia <- rep(0, n_random + 1)
europe_pacific <- rep(0, n_random + 1)
europe_northamerica <- rep(0, n_random + 1)
europe_southamerica <- rep(0, n_random + 1)

africa_europe <- rep(0, n_random + 1)
africa_africa <- rep(0, n_random + 1)
africa_temperateasia <- rep(0, n_random + 1)
africa_tropicalasia <- rep(0, n_random + 1)
africa_australasia <- rep(0, n_random + 1)
africa_pacific <- rep(0, n_random + 1)
africa_northamerica <- rep(0, n_random + 1)
africa_southamerica <- rep(0, n_random + 1)

temperateasia_europe <- rep(0, n_random + 1)
temperateasia_africa <- rep(0, n_random + 1)
temperateasia_temperateasia <- rep(0, n_random + 1)
temperateasia_tropicalasia <- rep(0, n_random + 1)
temperateasia_australasia <- rep(0, n_random + 1)
temperateasia_pacific <- rep(0, n_random + 1)
temperateasia_northamerica <- rep(0, n_random + 1)
temperateasia_southamerica <- rep(0, n_random + 1)

tropicalasia_europe <- rep(0, n_random + 1)
tropicalasia_africa <- rep(0, n_random + 1)
tropicalasia_temperateasia <- rep(0, n_random + 1)
tropicalasia_tropicalasia <- rep(0, n_random + 1)
tropicalasia_australasia <- rep(0, n_random + 1)
tropicalasia_pacific <- rep(0, n_random + 1)
tropicalasia_northamerica <- rep(0, n_random + 1)
tropicalasia_southamerica <- rep(0, n_random + 1)

australasia_europe <- rep(0, n_random + 1)
australasia_africa <- rep(0, n_random + 1)
australasia_temperateasia <- rep(0, n_random + 1)
australasia_tropicalasia <- rep(0, n_random + 1)
australasia_australasia <- rep(0, n_random + 1)
australasia_pacific <- rep(0, n_random + 1)
australasia_northamerica <- rep(0, n_random + 1)
australasia_southamerica <- rep(0, n_random + 1)

pacific_europe <- rep(0, n_random + 1)
pacific_africa <- rep(0, n_random + 1)
pacific_temperateasia <- rep(0, n_random + 1)
pacific_tropicalasia <- rep(0, n_random + 1)
pacific_australasia <- rep(0, n_random + 1)
pacific_pacific <- rep(0, n_random + 1)
pacific_northamerica <- rep(0, n_random + 1)
pacific_southamerica <- rep(0, n_random + 1)

northamerica_europe <- rep(0, n_random + 1)
northamerica_africa <- rep(0, n_random + 1)
northamerica_temperateasia <- rep(0, n_random + 1)
northamerica_tropicalasia <- rep(0, n_random + 1)
northamerica_australasia <- rep(0, n_random + 1)
northamerica_pacific <- rep(0, n_random + 1)
northamerica_northamerica <- rep(0, n_random + 1)
northamerica_southamerica <- rep(0, n_random + 1)

southamerica_europe <- rep(0, n_random + 1)
southamerica_africa <- rep(0, n_random + 1)
southamerica_temperateasia <- rep(0, n_random + 1)
southamerica_tropicalasia <- rep(0, n_random + 1)
southamerica_australasia <- rep(0, n_random + 1)
southamerica_pacific <- rep(0, n_random + 1)
southamerica_northamerica <- rep(0, n_random + 1)
southamerica_southamerica <- rep(0, n_random + 1)

# create an empty data frame from the vectors above
count_sourcereceptor_nonthreatened <- as.data.frame(cbind(europe_europe,
                                                          europe_africa,
                                                          europe_temperateasia,
                                                          europe_tropicalasia,
                                                          europe_australasia,
                                                          europe_pacific,
                                                          europe_northamerica,
                                                          europe_southamerica,
                                                          africa_europe,
                                                          africa_africa,
                                                          africa_temperateasia,
                                                          africa_tropicalasia,
                                                          africa_australasia,
                                                          africa_pacific,
                                                          africa_northamerica,
                                                          africa_southamerica,
                                                          temperateasia_europe,
                                                          temperateasia_africa,
                                                          temperateasia_temperateasia,
                                                          temperateasia_tropicalasia,
                                                          temperateasia_australasia,
                                                          temperateasia_pacific,
                                                          temperateasia_northamerica,
                                                          temperateasia_southamerica,
                                                          tropicalasia_europe,
                                                          tropicalasia_africa,
                                                          tropicalasia_temperateasia,
                                                          tropicalasia_tropicalasia,
                                                          tropicalasia_australasia,
                                                          tropicalasia_pacific,
                                                          tropicalasia_northamerica,
                                                          tropicalasia_southamerica,
                                                          australasia_europe,
                                                          australasia_africa,
                                                          australasia_temperateasia,
                                                          australasia_tropicalasia,
                                                          australasia_australasia,
                                                          australasia_pacific,
                                                          australasia_northamerica,
                                                          australasia_southamerica,
                                                          pacific_europe,
                                                          pacific_africa,
                                                          pacific_temperateasia,
                                                          pacific_tropicalasia,
                                                          pacific_australasia,
                                                          pacific_pacific,
                                                          pacific_northamerica,
                                                          pacific_southamerica,
                                                          northamerica_europe,
                                                          northamerica_africa,
                                                          northamerica_temperateasia,
                                                          northamerica_tropicalasia,
                                                          northamerica_australasia,
                                                          northamerica_pacific,
                                                          northamerica_northamerica,
                                                          northamerica_southamerica,
                                                          southamerica_europe,
                                                          southamerica_africa,
                                                          southamerica_temperateasia,
                                                          southamerica_tropicalasia,
                                                          southamerica_australasia,
                                                          southamerica_pacific,
                                                          southamerica_northamerica,
                                                          southamerica_southamerica))

# put observed values in first row of data frame
count_sourcereceptor_nonthreatened[1, ] <- observed_nonthreatened

for (id in 2: (n_random + 1)) {
  
  print(id)
  
  # draw for each continent separately a number of threatened taxa corresponding to the observed threatened number in that continent
  # make sure the same order as destination
  europe_draw <- sample(native_nonthreatened_sppre, europe_dest_nonthreatened)
  africa_draw <- sample(native_nonthreatened_sppre, africa_dest_nonthreatened)
  temperateasia_draw <- sample(native_nonthreatened_sppre, temperateasia_dest_nonthreatened)
  tropicalasia_draw <- sample(native_nonthreatened_sppre, tropicalasia_dest_nonthreatened)
  australasia_draw <- sample(native_nonthreatened_sppre, australasia_dest_nonthreatened)
  pacific_draw <- sample(native_nonthreatened_sppre, pacific_dest_nonthreatened)
  northamerica_draw <- sample(native_nonthreatened_sppre, northamerica_dest_nonthreatened)
  southamerica_draw <- sample(native_nonthreatened_sppre, southamerica_dest_nonthreatened)
  
  # extract the migration flow according to the random species
  europe_draw_data <- res_file[which(res_file$spacc %in% europe_draw), ]
  africa_draw_data <- res_file[which(res_file$spacc %in% africa_draw), ]
  temperateasia_draw_data <- res_file[which(res_file$spacc %in% temperateasia_draw), ]
  tropicalasia_draw_data <- res_file[which(res_file$spacc %in% tropicalasia_draw), ]
  australasia_draw_data <- res_file[which(res_file$spacc %in% australasia_draw), ]
  pacific_draw_data <- res_file[which(res_file$spacc %in% pacific_draw), ]
  northamerica_draw_data <- res_file[which(res_file$spacc %in% northamerica_draw), ]
  southamerica_draw_data <- res_file[which(res_file$spacc %in% southamerica_draw), ]
  
  # add the sums for each flow to the columns in the count_sourcereceptor dataframe
  count_sourcereceptor_nonthreatened[id, 1] <- length(unique(europe_draw_data$spacc[which(europe_draw_data$origin == 1)]))
  count_sourcereceptor_nonthreatened[id, 2] <- length(unique(africa_draw_data$spacc[which(africa_draw_data$origin == 1)]))
  count_sourcereceptor_nonthreatened[id, 3] <- length(unique(temperateasia_draw_data$spacc[which(temperateasia_draw_data$origin == 1)]))
  count_sourcereceptor_nonthreatened[id, 4] <- length(unique(tropicalasia_draw_data$spacc[which(tropicalasia_draw_data$origin == 1)]))
  count_sourcereceptor_nonthreatened[id, 5] <- length(unique(australasia_draw_data$spacc[which(australasia_draw_data$origin == 1)]))
  count_sourcereceptor_nonthreatened[id, 6] <- length(unique(pacific_draw_data$spacc[which(pacific_draw_data$origin == 1)]))
  count_sourcereceptor_nonthreatened[id, 7] <- length(unique(northamerica_draw_data$spacc[which(northamerica_draw_data$origin == 1)]))
  count_sourcereceptor_nonthreatened[id, 8] <- length(unique(southamerica_draw_data$spacc[which(southamerica_draw_data$origin == 1)]))
  
  count_sourcereceptor_nonthreatened[id, 9] <- length(unique(europe_draw_data$spacc[which(europe_draw_data$origin == 2)]))
  count_sourcereceptor_nonthreatened[id, 10] <- length(unique(africa_draw_data$spacc[which(africa_draw_data$origin == 2)]))
  count_sourcereceptor_nonthreatened[id, 11] <- length(unique(temperateasia_draw_data$spacc[which(temperateasia_draw_data$origin == 2)]))
  count_sourcereceptor_nonthreatened[id, 12] <- length(unique(tropicalasia_draw_data$spacc[which(tropicalasia_draw_data$origin == 2)]))
  count_sourcereceptor_nonthreatened[id, 13] <- length(unique(australasia_draw_data$spacc[which(australasia_draw_data$origin == 2)]))
  count_sourcereceptor_nonthreatened[id, 14] <- length(unique(pacific_draw_data$spacc[which(pacific_draw_data$origin == 2)]))
  count_sourcereceptor_nonthreatened[id, 15] <- length(unique(northamerica_draw_data$spacc[which(northamerica_draw_data$origin == 2)]))
  count_sourcereceptor_nonthreatened[id, 16] <- length(unique(southamerica_draw_data$spacc[which(southamerica_draw_data$origin == 2)]))
  
  count_sourcereceptor_nonthreatened[id, 17] <- length(unique(europe_draw_data$spacc[which(europe_draw_data$origin == 3)]))
  count_sourcereceptor_nonthreatened[id, 18] <- length(unique(africa_draw_data$spacc[which(africa_draw_data$origin == 3)]))
  count_sourcereceptor_nonthreatened[id, 19] <- length(unique(temperateasia_draw_data$spacc[which(temperateasia_draw_data$origin == 3)]))
  count_sourcereceptor_nonthreatened[id, 20] <- length(unique(tropicalasia_draw_data$spacc[which(tropicalasia_draw_data$origin == 3)]))
  count_sourcereceptor_nonthreatened[id, 21] <- length(unique(australasia_draw_data$spacc[which(australasia_draw_data$origin == 3)]))
  count_sourcereceptor_nonthreatened[id, 22] <- length(unique(pacific_draw_data$spacc[which(pacific_draw_data$origin == 3)]))
  count_sourcereceptor_nonthreatened[id, 23] <- length(unique(northamerica_draw_data$spacc[which(northamerica_draw_data$origin == 3)]))
  count_sourcereceptor_nonthreatened[id, 24] <- length(unique(southamerica_draw_data$spacc[which(southamerica_draw_data$origin == 3)]))
  
  count_sourcereceptor_nonthreatened[id, 25] <- length(unique(europe_draw_data$spacc[which(europe_draw_data$origin == 4)]))
  count_sourcereceptor_nonthreatened[id, 26] <- length(unique(africa_draw_data$spacc[which(africa_draw_data$origin == 4)]))
  count_sourcereceptor_nonthreatened[id, 27] <- length(unique(temperateasia_draw_data$spacc[which(temperateasia_draw_data$origin == 4)]))
  count_sourcereceptor_nonthreatened[id, 28] <- length(unique(tropicalasia_draw_data$spacc[which(tropicalasia_draw_data$origin == 4)]))
  count_sourcereceptor_nonthreatened[id, 29] <- length(unique(australasia_draw_data$spacc[which(australasia_draw_data$origin == 4)]))
  count_sourcereceptor_nonthreatened[id, 30] <- length(unique(pacific_draw_data$spacc[which(pacific_draw_data$origin == 4)]))
  count_sourcereceptor_nonthreatened[id, 31] <- length(unique(northamerica_draw_data$spacc[which(northamerica_draw_data$origin == 4)]))
  count_sourcereceptor_nonthreatened[id, 32] <- length(unique(southamerica_draw_data$spacc[which(southamerica_draw_data$origin == 4)]))
  
  count_sourcereceptor_nonthreatened[id, 33] <- length(unique(europe_draw_data$spacc[which(europe_draw_data$origin == 5)]))
  count_sourcereceptor_nonthreatened[id, 34] <- length(unique(africa_draw_data$spacc[which(africa_draw_data$origin == 5)]))
  count_sourcereceptor_nonthreatened[id, 35] <- length(unique(temperateasia_draw_data$spacc[which(temperateasia_draw_data$origin == 5)]))
  count_sourcereceptor_nonthreatened[id, 36] <- length(unique(tropicalasia_draw_data$spacc[which(tropicalasia_draw_data$origin == 5)]))
  count_sourcereceptor_nonthreatened[id, 37] <- length(unique(australasia_draw_data$spacc[which(australasia_draw_data$origin == 5)]))
  count_sourcereceptor_nonthreatened[id, 38] <- length(unique(pacific_draw_data$spacc[which(pacific_draw_data$origin == 5)]))
  count_sourcereceptor_nonthreatened[id, 39] <- length(unique(northamerica_draw_data$spacc[which(northamerica_draw_data$origin == 5)]))
  count_sourcereceptor_nonthreatened[id, 40] <- length(unique(southamerica_draw_data$spacc[which(southamerica_draw_data$origin == 5)]))
  
  count_sourcereceptor_nonthreatened[id, 41] <- length(unique(europe_draw_data$spacc[which(europe_draw_data$origin == 6)]))
  count_sourcereceptor_nonthreatened[id, 42] <- length(unique(africa_draw_data$spacc[which(africa_draw_data$origin == 6)]))
  count_sourcereceptor_nonthreatened[id, 43] <- length(unique(temperateasia_draw_data$spacc[which(temperateasia_draw_data$origin == 6)]))
  count_sourcereceptor_nonthreatened[id, 44] <- length(unique(tropicalasia_draw_data$spacc[which(tropicalasia_draw_data$origin == 6)]))
  count_sourcereceptor_nonthreatened[id, 45] <- length(unique(australasia_draw_data$spacc[which(australasia_draw_data$origin == 6)]))
  count_sourcereceptor_nonthreatened[id, 46] <- length(unique(pacific_draw_data$spacc[which(pacific_draw_data$origin == 6)]))
  count_sourcereceptor_nonthreatened[id, 47] <- length(unique(northamerica_draw_data$spacc[which(northamerica_draw_data$origin == 6)]))
  count_sourcereceptor_nonthreatened[id, 48] <- length(unique(southamerica_draw_data$spacc[which(southamerica_draw_data$origin == 6)]))
  
  count_sourcereceptor_nonthreatened[id, 49] <- length(unique(europe_draw_data$spacc[which(europe_draw_data$origin == 7)]))
  count_sourcereceptor_nonthreatened[id, 50] <- length(unique(africa_draw_data$spacc[which(africa_draw_data$origin == 7)]))
  count_sourcereceptor_nonthreatened[id, 51] <- length(unique(temperateasia_draw_data$spacc[which(temperateasia_draw_data$origin == 7)]))
  count_sourcereceptor_nonthreatened[id, 52] <- length(unique(tropicalasia_draw_data$spacc[which(tropicalasia_draw_data$origin == 7)]))
  count_sourcereceptor_nonthreatened[id, 53] <- length(unique(australasia_draw_data$spacc[which(australasia_draw_data$origin == 7)]))
  count_sourcereceptor_nonthreatened[id, 54] <- length(unique(pacific_draw_data$spacc[which(pacific_draw_data$origin == 7)]))
  count_sourcereceptor_nonthreatened[id, 55] <- length(unique(northamerica_draw_data$spacc[which(northamerica_draw_data$origin == 7)]))
  count_sourcereceptor_nonthreatened[id, 56] <- length(unique(southamerica_draw_data$spacc[which(southamerica_draw_data$origin == 7)]))
  
  count_sourcereceptor_nonthreatened[id, 57] <- length(unique(europe_draw_data$spacc[which(europe_draw_data$origin == 8)]))
  count_sourcereceptor_nonthreatened[id, 58] <- length(unique(africa_draw_data$spacc[which(africa_draw_data$origin == 8)]))
  count_sourcereceptor_nonthreatened[id, 59] <- length(unique(temperateasia_draw_data$spacc[which(temperateasia_draw_data$origin == 8)]))
  count_sourcereceptor_nonthreatened[id, 60] <- length(unique(tropicalasia_draw_data$spacc[which(tropicalasia_draw_data$origin == 8)]))
  count_sourcereceptor_nonthreatened[id, 61] <- length(unique(australasia_draw_data$spacc[which(australasia_draw_data$origin == 8)]))
  count_sourcereceptor_nonthreatened[id, 62] <- length(unique(pacific_draw_data$spacc[which(pacific_draw_data$origin == 8)]))
  count_sourcereceptor_nonthreatened[id, 63] <- length(unique(northamerica_draw_data$spacc[which(northamerica_draw_data$origin == 8)]))
  count_sourcereceptor_nonthreatened[id, 64] <- length(unique(southamerica_draw_data$spacc[which(southamerica_draw_data$origin == 8)]))
  
}

write.table(count_sourcereceptor_nonthreatened, file = paste0(path_folder, "03_migration_flow/random_sample_count_sourcereceptor_nonthreatened.txt"),
            col.names = TRUE, row.names = FALSE)

count_sourcereceptor_nonthreatened <- read.table(paste0(path_folder, "03_migration_flow/random_sample_count_sourcereceptor_nonthreatened.txt"),
                                                 header = TRUE)

# create dataframe with all combinations, observed values and statistics from the randomization test
source_receptor <- c("europe_europe",
                     "europe_africa",
                     "europe_temperateasia",
                     "europe_tropicalasia",
                     "europe_australasia",
                     "europe_pacific",
                     "europe_northamerica",
                     "europe_southamerica",
                     "africa_europe",
                     "africa_africa",
                     "africa_temperateasia",
                     "africa_tropicalasia",
                     "africa_australasia",
                     "africa_pacific",
                     "africa_northamerica",
                     "africa_southamerica",
                     "temperateasia_europe",
                     "temperateasia_africa",
                     "temperateasia_temperateasia",
                     "temperateasia_tropicalasia",
                     "temperateasia_australasia",
                     "temperateasia_pacific",
                     "temperateasia_northamerica",
                     "temperateasia_southamerica",
                     "tropicalasia_europe",
                     "tropicalasia_africa",
                     "tropicalasia_temperateasia",
                     "tropicalasia_tropicalasia",
                     "tropicalasia_australasia",
                     "tropicalasia_pacific",
                     "tropicalasia_northamerica",
                     "tropicalasia_southamerica",
                     "australasia_europe",
                     "australasia_africa",
                     "australasia_temperateasia",
                     "australasia_tropicalasia",
                     "australasia_australasia",
                     "australasia_pacific",
                     "australasia_northamerica",
                     "australasia_southamerica",
                     "pacific_europe",
                     "pacific_africa",
                     "pacific_temperateasia",
                     "pacific_tropicalasia",
                     "pacific_australasia",
                     "pacific_pacific",
                     "pacific_northamerica",
                     "pacific_southamerica",
                     "northamerica_europe",
                     "northamerica_africa",
                     "northamerica_temperateasia",
                     "northamerica_tropicalasia",
                     "northamerica_australasia",
                     "northamerica_pacific",
                     "northamerica_northamerica",
                     "northamerica_southamerica",
                     "southamerica_europe",
                     "southamerica_africa",
                     "southamerica_temperateasia",
                     "southamerica_tropicalasia",
                     "southamerica_australasia",
                     "southamerica_pacific",
                     "southamerica_northamerica",
                     "southamerica_southamerica")

random_res_nonthreatened <- as.data.frame(cbind(source_receptor, observed_nonthreatened))
# create empty column in random_res for the proportion of random values that is smaller than the observed value
random_res_nonthreatened$proportion <- rep(NA, length(random_res_nonthreatened$observed_nonthreatened))
# create empty column in random_res for mean value of the random draws
random_res_nonthreatened$mean <- rep(NA, length(random_res_nonthreatened$observed_nonthreatened))
# create empty column in random_res for median of the random draws
random_res_nonthreatened$median <- rep(NA, length(random_res_nonthreatened$observed_nonthreatened))

# calculate whether observed flows are higher (proportion close to 1) or lower (proportion close to 0) than expected
for (id in 1: dim(random_res_nonthreatened)[1]) {
  
  random_data_nonthreatened <- count_sourcereceptor_nonthreatened[-1, id]
  observed_data_nonthreatened <- count_sourcereceptor_nonthreatened[1, id]
  
  random_res_nonthreatened$proportion[id] <- length(random_data_nonthreatened[random_data_nonthreatened < observed_data_nonthreatened])/9999
  random_res_nonthreatened$mean[id] <- mean(count_sourcereceptor_nonthreatened[-1, id])
  random_res_nonthreatened$median[id] <- median(count_sourcereceptor_nonthreatened[-1, id])
  
  print(id)
  
}

# define line colors for observed-value lines in histograms
random_res_nonthreatened$line_color <- rep("black", dim(random_res_nonthreatened)[1]) 
random_res_nonthreatened$line_color <- replace(random_res_nonthreatened$line_color, random_res_nonthreatened$proportion > 0.975, "red")
random_res_nonthreatened$line_color <- replace(random_res_nonthreatened$line_color, random_res_nonthreatened$proportion < 0.025, "blue")

###########################
# plot the expected flows
# form the dataframe
for (id in 1: dim(random_res_nonthreatened)[1]) {
  
  sour_rece_name <- random_res_nonthreatened[id, 1]
  
  origin_name <- strsplit(sour_rece_name, "_", fixed = TRUE)[[1]][1]
  destination_name <- strsplit(sour_rece_name, "_", fixed = TRUE)[[1]][2]
  
  sour_rece_flow <- random_res_nonthreatened[id, 5]
  sour_rece_res <- cbind(origin_name, destination_name, sour_rece_flow)
  write.table(sour_rece_res, file = paste0(path_folder, "03_migration_flow/random_sample_migration_flow_median_nonthreatened.txt"),
              col.names = FALSE, row.names = FALSE, append = TRUE)
  
  print(id)
  
}

flow_data_random_nonthreatened <- read.table(paste0(path_folder, "03_migration_flow/random_sample_migration_flow_median_nonthreatened.txt"))
colnames(flow_data_random_nonthreatened) <- c("origin", "destination", "flow")
write.table(flow_data_random_nonthreatened, file = paste0(path_folder, "03_migration_flow/random_sample_migration_flow_median_nonthreatened.txt"),
            col.names = TRUE, row.names = FALSE)

# manually modify the name of continents
flow_data_random_nonthreatened <- read.table(paste0(path_folder, "03_migration_flow/random_sample_migration_flow_median_nonthreatened.txt"),
                                             header = TRUE)

# select the top 50% of migration flows
flow_data_random_nonthreatened <- flow_data_random_nonthreatened[order(flow_data_random_nonthreatened$flow), ]
flow_data_random_nonthreatened <- flow_data_random_nonthreatened[33: 64, ]

# plot
setwd(path_folder)
png(file = paste0(path_folder, "03_migration_flow/migrationflow_nonnativetree_randomsamplemedian_nonthreatened.png"),
    bg = "transparent", width = 8, height = 8, units = "cm", res = 300)

circos.clear()

circos.par(start.degree = 90,
           gap.after = 5,
           track.margin = c(0.01, 0.01))

chordDiagram(x = flow_data_random_nonthreatened, 
             order = c("South America",
                       "North America",
                       "Pacific",
                       "Tropical Asia",
                       "Australasia",
                       "Temperate Asia",
                       "Africa",
                       "Europe"),
             grid.col = c("South America" = "#5D74A5",
                          "North America" = "#8CA5CA",
                          "Pacific" = "#BBD1E2",
                          "Tropical Asia" = "#E7EAD0",
                          "Australasia" = "#F8DEB2",
                          "Temperate Asia" = "#EDAC88",
                          "Africa" = "#CE7F69",
                          "Europe" = "#A8554E"),
             directional = 1, 
             annotationTrack = c("grid"),
             transparency = 0,
             annotationTrackHeight = c(0.06, 0.07),
             preAllocateTracks = list(track.height = mm_h(0.1)),
             direction.type = "diffHeight",
             diffHeight  = mm_h(4),
             link.target.prop = FALSE)

# add labels and axis
circos.track(bg.border = NA, ylim = c(0,1),
             panel.fun = function(x, y) {
               xcenter = get.cell.meta.data("xcenter")
               sector.index = get.cell.meta.data("sector.index")
               circos.text(x = xcenter, y = 2.1, labels = sector.index, 
                           facing = "bending", cex = 0.7,font = 1, family = "Times New Roman")
               circos.axis(h = "top", labels = FALSE, track.index = 2,
                           minor.ticks = 0,major.at = seq(0, 25000, 1000),
                           lwd = 1,major.tick.length = mm_y(0.5))
             })

dev.off()







