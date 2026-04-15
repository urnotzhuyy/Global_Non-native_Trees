rm(list=ls())

# loading required libraries
library(raster)
library(ggplot2)
library(MetBrewer)
library(scales)

# path
path_folder <- "D:/non_native_tree/"

###############################
# generate the template layer
world_shp <- shapefile(paste0(path_folder, "01_distribution_pattern/level4_remove_Antarctica/level4_remove_Antarctica.shp"))

raster_temp <- raster(resolution = c(1, 1),
                      crs = "EPSG:4326",
                      vals = 0)

world_raster <- rasterize(world_shp, raster_temp, getCover = TRUE)

# remove the land area less than 30% of each cell
world_raster_land <- world_raster
values(world_raster_land)[which(values(world_raster_land) < 0.3)] <- NA
values(world_raster_land)[which(!is.na(values(world_raster_land)))] <- 0
writeRaster(world_raster_land, filename = paste0(path_folder, "01_distribution_pattern/global_land_onedegree.tif"))

########################
# distribution pattern
########################

# load the null map
null_map <- raster(paste0(path_folder, "01_distribution_pattern/global_land_onedegree.tif"))

cell_data <- matrix(data = 0, nrow = 64800, ncol = 2)
for (cellid in 1: 64800) {
  
  cell_data[cellid, 1] <- cellid
  print(cellid)
  
}

###############
# native tree
native_celldata <- cell_data

for (id in 1: length(native_all_sppre)) {
  
  sp <- native_all_sppre[id]
  
  if(sp %in% native_list){
    # powo and iucn distribution information
    sp_file <- read.table(paste0(path_folder, "aggregate_data/tree_data_native_final_update/", sp, "_native.txt"), 
                          header = TRUE)
  }else{
    # bgci distribution information
    sp_file <- read.table(paste0(path_folder, "aggregate_data/tree_data_native_final_update/bgci_part/", sp, "_native.txt"),
                          header = TRUE)
  }
  
  sp_coor <- sp_file[, c(2, 3)]
  
  cell_num <- cellFromXY(null_map, sp_coor)
  cell_uni <- unique(cell_num)
  
  for (cell in 1: length(cell_uni)) {
    
    cell_sp <- cell_uni[cell]
    
    cell_oldval <- native_celldata[cell_sp, 2]
    cell_newval <- cell_oldval + 1
    
    native_celldata[cell_sp, 2] <- cell_newval
    
  }
  
  print(id)
  
}

# form the distribution pattern of native trees
native_map <- null_map

for (id in 1: 64800) {
  
  cellid <- id
  
  if(!is.na(native_map[cellid])){
    
    native_val <- native_celldata[cellid, 2]
    native_map[cellid] <- native_val
    
  }
  
  print(id)
  
}

writeRaster(native_map, filename = paste0(path_folder, "01_distribution_pattern/distribution_pattern_nativetree.tif"))

###################
# non-native tree
nonnative_celldata <- cell_data

for (id in 1: length(nonnative_all_sppre)) {
  
  sp <- nonnative_all_sppre[id]
  
  if(sp %in% nonnative_list){
    # powo and iucn distribution information
    sp_file <- read.table(paste0(path_folder, "aggregate_data/tree_data_nonnative_final_update/", sp, "_nonnative.txt"),
                          header = TRUE)
  }else{
    # bgci distribution information
    sp_file <- read.table(paste0(path_folder, "aggregate_data/tree_data_nonnative_final_update/bgci_part/", sp, "_nonnative.txt"),
                          header = TRUE)
  }
  
  sp_coor <- sp_file[, c(2, 3)]
  
  cell_num <- cellFromXY(null_map, sp_coor)
  cell_uni <- unique(cell_num)
  
  for (cell in 1: length(cell_uni)) {
    
    cell_sp <- cell_uni[cell]
    
    cell_oldval <- nonnative_celldata[cell_sp, 2]
    cell_newval <- cell_oldval + 1
    
    nonnative_celldata[cell_sp, 2] <- cell_newval
    
  }
  
  print(id)
  
}

# form the distribution pattern of non-native trees
nonnative_map <- null_map

for (id in 1: 64800) {
  
  cellid <- id
  
  if(!is.na(nonnative_map[cellid])){
    
    nonnative_val <- nonnative_celldata[cellid, 2]
    nonnative_map[cellid] <- nonnative_val
    
  }
  
  print(id)
  
}

writeRaster(nonnative_map, filename = paste0(path_folder, "01_distribution_pattern/distribution_pattern_nonnativetree.tif"))

###################
# threatened tree
###################

# distribution patterns of threatened trees
# threatened native tree
threatenednative_celldata <- cell_data

for (id in 1: length(native_threatened_sppre)) {
  
  sp <- native_threatened_sppre[id]
  
  if(sp %in% native_list){
    # powo and iucn distribution information
    sp_file <- read.table(paste0(path_folder, "aggregate_data/tree_data_native_final_update/", sp, "_native.txt"),
                          header = TRUE)
  }else{
    # bgci distribution information
    sp_file <- read.table(paste0(path_folder, "aggregate_data/tree_data_native_final_update/bgci_part/", sp, "_native.txt"),
                          header = TRUE)
  }
  
  sp_coor <- sp_file[, c(2, 3)]
  
  cell_num <- cellFromXY(null_map, sp_coor)
  cell_uni <- unique(cell_num)
  
  for (cell in 1: length(cell_uni)) {
    
    cell_sp <- cell_uni[cell]
    
    cell_oldval <- threatenednative_celldata[cell_sp, 2]
    cell_newval <- cell_oldval + 1
    
    threatenednative_celldata[cell_sp, 2] <- cell_newval
    
  }
  
  print(id)
  
}

# form the distribution pattern of threatened native trees
threatenednative_map <- null_map

for (id in 1: 64800) {
  
  cellid <- id
  
  if(!is.na(threatenednative_map[cellid])){
    
    threatenednative_val <- threatenednative_celldata[cellid, 2]
    threatenednative_map[cellid] <- threatenednative_val
    
  }
  
  print(id)
  
}

writeRaster(threatenednative_map, filename = paste0(path_folder, "01_distribution_pattern/distribution_pattern_nativetree_threatened.tif"))

##############################
# threatened non-native tree
threatenednonnative_celldata <- cell_data

for (id in 1: length(nonnative_threatened_sppre)) {
  
  sp <- nonnative_threatened_sppre[id]
  
  if(sp %in% nonnative_list){
    # powo and iucn distribution information
    sp_file <- read.table(paste0(path_folder, "aggregate_data/tree_data_nonnative_final_update/", sp, "_nonnative.txt"),
                          header = TRUE)
  }else{
    # bgci distribution information
    sp_file <- read.table(paste0(path_folder, "aggregate_data/tree_data_nonnative_final_update/bgci_part/", sp, "_nonnative.txt"),
                          header = TRUE)
  }
  
  sp_coor <- sp_file[, c(2, 3)]
  
  cell_num <- cellFromXY(null_map, sp_coor)
  cell_uni <- unique(cell_num)
  
  for (cell in 1: length(cell_uni)) {
    
    cell_sp <- cell_uni[cell]
    
    cell_oldval <- threatenednonnative_celldata[cell_sp, 2]
    cell_newval <- cell_oldval + 1
    
    threatenednonnative_celldata[cell_sp, 2] <- cell_newval
    
  }
  
  print(id)
  
}

# form the distribution pattern of threatened non-native trees
threatenednonnative_map <- null_map

for (id in 1: 64800) {
  
  cellid <- id
  
  if(!is.na(threatenednonnative_map[cellid])){
    
    threatenednonnative_val <- threatenednonnative_celldata[cellid, 2]
    threatenednonnative_map[cellid] <- threatenednonnative_val
    
  }
  
  print(id)
  
}

writeRaster(threatenednonnative_map, filename = paste0(path_folder, "01_distribution_pattern/distribution_pattern_nonnativetree_threatened.tif"))

#######################
# non-threatened tree
#######################

# distribution patterns of non-threatened trees
# non-threatened native tree
nonthreatenednative_celldata <- cell_data

for (id in 1: length(native_nonthreatened_sppre)) {
  
  sp <- native_nonthreatened_sppre[id]
  
  if(sp %in% native_list){
    # powo and iucn distribution information
    sp_file <- read.table(paste0(path_folder, "aggregate_data/tree_data_native_final_update/", sp, "_native.txt"),
                          header = TRUE)
  }else{
    # bgci distribution information
    sp_file <- read.table(paste0(path_folder, "aggregate_data/tree_data_native_final_update/bgci_part/", sp, "_native.txt"),
                          header = TRUE)
  }
  
  sp_coor <- sp_file[, c(2, 3)]
  
  cell_num <- cellFromXY(null_map, sp_coor)
  cell_uni <- unique(cell_num)
  
  for (cell in 1: length(cell_uni)) {
    
    cell_sp <- cell_uni[cell]
    
    cell_oldval <- nonthreatenednative_celldata[cell_sp, 2]
    cell_newval <- cell_oldval + 1
    
    nonthreatenednative_celldata[cell_sp, 2] <- cell_newval
    
  }
  
  print(id)
  
}

# form the distribution pattern of non-threatened native trees
nonthreatenednative_map <- null_map

for (id in 1: 64800) {
  
  cellid <- id
  
  if(!is.na(nonthreatenednative_map[cellid])){
    
    nonthreatenednative_val <- nonthreatenednative_celldata[cellid, 2]
    nonthreatenednative_map[cellid] <- nonthreatenednative_val
    
  }
  
  print(id)
  
}

writeRaster(nonthreatenednative_map, filename = paste0(path_folder, "01_distribution_pattern/distribution_pattern_nativetree_nonthreatened.tif"))

##################################
# non-threatened non-native tree
nonthreatenednonnative_celldata <- cell_data

for (id in 1: length(nonnative_nonthreatened_sppre)) {
  
  sp <- nonnative_nonthreatened_sppre[id]
  
  if(sp %in% nonnative_list){
    # powo and iucn distribution information
    sp_file <- read.table(paste0(path_folder, "aggregate_data/tree_data_nonnative_final_update/", sp, "_nonnative.txt"),
                          header = TRUE)
  }else{
    # bgci distribution information
    sp_file <- read.table(paste0(path_folder, "aggregate_data/tree_data_nonnative_final_update/bgci_part/", sp, "_nonnative.txt"),
                          header = TRUE)
  }
  
  sp_coor <- sp_file[, c(2, 3)]
  
  cell_num <- cellFromXY(null_map, sp_coor)
  cell_uni <- unique(cell_num)
  
  for (cell in 1: length(cell_uni)) {
    
    cell_sp <- cell_uni[cell]
    
    cell_oldval <- nonthreatenednonnative_celldata[cell_sp, 2]
    cell_newval <- cell_oldval + 1
    
    nonthreatenednonnative_celldata[cell_sp, 2] <- cell_newval
    
  }
  
  print(id)
  
}

# form the distribution pattern of non-threatened non-native trees
nonthreatenednonnative_map <- null_map

for (id in 1: 64800) {
  
  cellid <- id
  
  if(!is.na(nonthreatenednonnative_map[cellid])){
    
    nonthreatenednonnative_val <- nonthreatenednonnative_celldata[cellid, 2]
    nonthreatenednonnative_map[cellid] <- nonthreatenednonnative_val
    
  }
  
  print(id)
  
}

writeRaster(nonthreatenednonnative_map, filename = paste0(path_folder, "01_distribution_pattern/distribution_pattern_nonnativetree_nonthreatened.tif"))

###########################################################
# relationship between native and non-native tree species
###########################################################

native_map <- raster(paste0(path_folder, "01_distribution_pattern/distribution_pattern_nativetree.tif"))
nonnative_map <- raster(paste0(path_folder, "01_distribution_pattern/distribution_pattern_nonnativetree.tif"))

threatenednative_map <- raster(paste0(path_folder, "01_distribution_pattern/distribution_pattern_nativetree_threatened.tif"))
threatenednonnative_map <- raster(paste0(path_folder, "01_distribution_pattern/distribution_pattern_nonnativetree_threatened.tif"))

nonthreatenednative_map <- raster(paste0(path_folder, "01_distribution_pattern/distribution_pattern_nativetree_nonthreatened.tif"))
nonthreatenednonnative_map <- raster(paste0(path_folder, "01_distribution_pattern/distribution_pattern_nonnativetree_nonthreatened.tif"))

for (cellid in 1: 64800) {
  
  # grid id
  gridid <- cellid
  
  if(!is.na(null_map[cellid])){
    
    # cell xy
    cellx <- xFromCell(null_map, cellid)
    celly <- yFromCell(null_map, cellid)
    
    # richness
    native_ri <- native_map[cellid]
    nonnative_ri <- nonnative_map[cellid]
    
    threatenednative_ri <- threatenednative_map[cellid]
    threatenednonnative_ri <- threatenednonnative_map[cellid]
    
    nonthreatenednative_ri <- nonthreatenednative_map[cellid]
    nonthreatenednonnative_ri <- nonthreatenednonnative_map[cellid]
    
    cell_res <- cbind(gridid, cellx, celly, native_ri, nonnative_ri,
                      threatenednative_ri, threatenednonnative_ri,
                      nonthreatenednative_ri, nonthreatenednonnative_ri)
    write.table(cell_res, file = paste0(path_folder, "01_distribution_pattern/richness_coordinate.txt"),
                col.names = FALSE, row.names = FALSE, append = TRUE)
    
  }
  
  print(cellid)
  
}

cell_res <- read.table(paste0(path_folder, "01_distribution_pattern/richness_coordinate.txt"))
colnames(cell_res) <- c("gridid", "cellx", "celly", "native_ri", "nonnative_ri",
                        "threatenednative_ri", "threatenednonnative_ri",
                        "nonthreatenednative_ri", "nonthreatenednonnative_ri")

cell_remove <- cell_res$gridid[which(cell_res$native_ri == 0 & cell_res$nonnative_ri == 0 &
                                       cell_res$threatenednative_ri == 0 & cell_res$threatenednonnative_ri == 0 &
                                       cell_res$nonthreatenednative_ri == 0 & cell_res$nonthreatenednonnative_ri == 0)]

cell_keep <- cell_res[which(!cell_res$gridid %in% cell_remove), ]
write.table(cell_keep, file = paste0(path_folder, "01_distribution_pattern/richness_coordinate.txt"),
            col.names = TRUE, row.names = FALSE)

# TDWG continents
richness_tdwg <- read.csv(paste0(path_folder, "01_distribution_pattern/richness_coordinate_continent.csv"))
richness_tdwg <- richness_tdwg[which(richness_tdwg$LEVEL1_COD != 9), ] # remove one grid cell that occurred in Antarctic islands

richness_tdwg_plot <- richness_tdwg[, c(4: 12, 16)]
colnames(richness_tdwg_plot) <- c("gridid", "cellx", "celly",
                                   "native_ri", "nonnative_ri",
                                   "threatenednative_ri", "threatenednonnative_ri",
                                   "nonthreatenednative_ri", "nonthreatenednonnative_ri",
                                   "tdwg_continent")

richness_tdwg_plot$tdwg_continent <- as.factor(richness_tdwg_plot$tdwg_continent)

# non-native trees
setwd(path_folder)
tiff(filename = paste0(path_folder, "01_distribution_pattern/relationship_native_nonnative.tif"),
     bg = "transparent", width = 7.5, height = 7.5, units = "cm", res = 300)

ggplot(richness_tdwg_plot, aes(x = native_ri, y = nonnative_ri)) +
  geom_point(aes(color = tdwg_continent), size = 0.5) +
  geom_smooth(method = "loess", linewidth = 0.1, color = "black") +
  scale_color_manual(values = c("1" = "#A8554E",
                                "2" = "#CE7F69",
                                "3" = "#EDAC88",
                                "4" = "#E7EAD0",
                                "5" = "#F8DEB2",
                                "6" = "#BBD1E2",
                                "7" = "#8CA5CA",
                                "8" = "#5D74A5")) +
  theme(panel.background = element_blank(),
        axis.line = element_line(),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 11)) +
  xlab("Overall native richness") +
  ylab("Overall non-native richness") +
  theme(legend.position = "none")
  
dev.off()

# threatened non-native trees
setwd(path_folder)
tiff(filename = paste0(path_folder, "01_distribution_pattern/relationship_threatenednative_threatenednonnative.tif"),
     bg = "transparent", width = 7.5, height = 7.5, units = "cm", res = 300)

ggplot(richness_tdwg_plot, aes(x = threatenednative_ri, y = threatenednonnative_ri)) +
  geom_point(aes(color = tdwg_continent), size = 0.5) +
  geom_smooth(method = "loess", linewidth = 0.1, color = "black") +
  scale_color_manual(values = c("1" = "#A8554E",
                                "2" = "#CE7F69",
                                "3" = "#EDAC88",
                                "4" = "#E7EAD0",
                                "5" = "#F8DEB2",
                                "6" = "#BBD1E2",
                                "7" = "#8CA5CA",
                                "8" = "#5D74A5")) +
  scale_x_continuous(breaks = c(0, 100, 200, 300, 400),
                     labels = c(0, 100, 200, 300, 400)) +
  theme(panel.background = element_blank(),
        axis.line = element_line(),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 11)) +
  xlab("Threatened native richness") +
  ylab("Threatened non-native richness") +
  theme(legend.position = "none")

dev.off()

# non-threatened non-native trees
setwd(path_folder)
tiff(filename = paste0(path_folder, "01_distribution_pattern/relationship_nonthreatenednative_nonthreatenednonnative.tif"),
     bg = "transparent", width = 7.5, height = 7.5, units = "cm", res = 300)

ggplot(richness_tdwg_plot, aes(x = nonthreatenednative_ri, y = nonthreatenednonnative_ri)) +
  geom_point(aes(color = tdwg_continent), size = 0.5) +
  geom_smooth(method = "loess", linewidth = 0.1, color = "black") +
  scale_color_manual(values = c("1" = "#A8554E",
                                "2" = "#CE7F69",
                                "3" = "#EDAC88",
                                "4" = "#E7EAD0",
                                "5" = "#F8DEB2",
                                "6" = "#BBD1E2",
                                "7" = "#8CA5CA",
                                "8" = "#5D74A5")) +
  theme(panel.background = element_blank(),
        axis.line = element_line(),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 11)) +
  xlab("Non-threatened native richness") +
  ylab("Non-threatened non-native richness") +
  theme(legend.position = "none")

dev.off()





