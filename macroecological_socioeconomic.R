rm(list=ls())

# loading required libraries
library(raster)
library(nlme)
library(car)
library(MuMIn)
library(ggplot2)
library(ncf)
library(PNWColors)
library(scales)

# path
path_folder <- "D:/non_native_tree/"

#################################
# extract the climate variables
#################################
# load the list of climate layers

clim_1 <- raster(paste0(path_folder, "02_macroecological_socioeconomic/raw_data/clim/wc2.1_30s_bio_1.tif"))
clim_7 <- raster(paste0(path_folder, "02_macroecological_socioeconomic/raw_data/clim/wc2.1_30s_bio_7.tif"))
clim_12 <- raster(paste0(path_folder, "02_macroecological_socioeconomic/raw_data/clim/wc2.1_30s_bio_12.tif"))
clim_15 <- raster(paste0(path_folder, "02_macroecological_socioeconomic/raw_data/clim/wc2.1_30s_bio_15.tif"))

clim_stack <- stack(clim_1, clim_7, clim_12, clim_15)

# load the grid cell
cell_data <- read.table(paste0(path_folder, "01_distribution_pattern/richness_coordinate.txt"), header = TRUE)

# load the null map
null_map <- raster(paste0(path_folder, "01_distribution_pattern/global_land_onedegree.tif"))

for (id in 1: dim(cell_data)[1]) {
  
  # grid id
  gridid <- cell_data[id, 1]
  
  # the raster of grid cell
  cell <- rasterFromCells(null_map, gridid, values = FALSE)
  
  # the value of grid cell is zero
  cell_poly <- rasterToPolygons(cell, na.rm = TRUE)
  
  # extract
  cell_value <- extract(clim_stack, cell_poly, fun = mean, na.rm = TRUE)
  
  cell_all <- cbind(gridid, cell_value)
  write.table(cell_all, file = paste0(path_folder, "02_macroecological_socioeconomic/richness_coordinate_clim.txt"),
              col.names = FALSE, row.names = FALSE, append = TRUE)
  
  print(id)
  
}

cell_clim <- read.table(paste0(path_folder, "02_macroecological_socioeconomic/richness_coordinate_clim.txt"))
colnames(cell_clim) <- c("gridid", "clim_1", "clim_7", "clim_12", "clim_15")
write.table(cell_clim, file = paste0(path_folder, "02_macroecological_socioeconomic/richness_coordinate_clim.txt"),
            col.names = TRUE, row.names = FALSE)

####################################################
# match the socioeconomic and insularity variables
####################################################

# extract the socioeconomic variables
gdp <- raster(paste0(path_folder, "02_macroecological_socioeconomic/raw_data/GDP/2019GDP-WGS84/2019GDP-WGS84.tif"))
ghm <- raster(paste0(path_folder, "02_macroecological_socioeconomic/raw_data/GHM/GHM-WGS84/GHM-WGS84.tif"))

# load the grid cell
cell_data <- read.table(paste0(path_folder, "01_distribution_pattern/richness_coordinate.txt"), header = TRUE)

# load the null map
null_map <- raster(paste0(path_folder, "01_distribution_pattern/global_land_onedegree.tif"))

for (id in 1: dim(cell_data)[1]) {
  
  # grid id
  gridid <- cell_data[id, 1]
  
  # the raster of grid cell
  cell <- rasterFromCells(null_map, gridid, values = FALSE)
  
  # the value of grid cell is zero
  cell_poly <- rasterToPolygons(cell, na.rm = TRUE)
  
  # extract
  cell_value_gdp <- extract(gdp, cell_poly, fun = mean, na.rm = TRUE)
  cell_value_ghm <- extract(ghm, cell_poly, fun = mean, na.rm = TRUE)
  
  cell_all <- cbind(gridid, cell_value_gdp, cell_value_ghm)
  write.table(cell_all, file = paste0(path_folder, "02_macroecological_socioeconomic/richness_coordinate_gdpghm.txt"),
              col.names = FALSE, row.names = FALSE, append = TRUE)
  
  print(id)
  
}

cell_gdpghm <- read.table(paste0(path_folder, "02_macroecological_socioeconomic/richness_coordinate_gdpghm.txt"))
colnames(cell_gdpghm) <- c("gridid", "gdp", "ghm")
write.table(cell_gdpghm, file = paste0(path_folder, "02_macroecological_socioeconomic/richness_coordinate_gdpghm.txt"),
            col.names = TRUE, row.names = FALSE)

####################################
# load the socioeconomic variables
gpi <- read.csv(paste0(path_folder, "02_macroecological_socioeconomic/raw_data/GPI/GADM_Country_GPI.csv"))

# load the GADM file of grid cell (need to manually deal with the gadm code of XCA to IRN)
cell_gadm <- read.csv(paste0(path_folder, "02_macroecological_socioeconomic/richness_coordinate_gadm.csv"))

# match the richness coordinates
cell_total <- merge(cell_gadm, gpi, by = "GID_gadm")
cell_total <- cell_total[, c(1: 10, 13: 14)]
cell_total <- cell_total[order(cell_total$gridid), ]
colnames(cell_total) <- c("GID_gadm", "gridid", "cellx", "celly",
                          "native_ri", "nonnative_ri", "threatenednative_ri", "threatenednonnative_ri",
                          "nonthreatenednative_ri", "nonthreatenednonnative_ri", "Country", "GPI")
write.table(cell_total, file = paste0(path_folder, "02_macroecological_socioeconomic/richness_coordinate_gpi.txt"),
            col.names = TRUE, row.names = FALSE)

#######################################
# TDWG level-4 region (random effect)
#######################################
# load the TDWG file of grid cell (need to manually deal with the part of coding problem)
cell_tdwg <- read.csv(paste0(path_folder, "02_macroecological_socioeconomic/richness_coordinate_tdwg.csv"))
cell_tdwg <- cell_tdwg[, c(1: 9, 12, 14: 16)]

write.table(cell_tdwg, file = paste0(path_folder, "02_macroecological_socioeconomic/richness_coordinate_tdwglevelall.txt"),
            col.names = TRUE, row.names = FALSE)

#############################
# linear mixed-effect model
#############################
# load the files
clim_data <- read.table(paste0(path_folder, "02_macroecological_socioeconomic/richness_coordinate_clim.txt"), header = TRUE)
gdpghm_data <- read.table(paste0(path_folder, "02_macroecological_socioeconomic/richness_coordinate_gdpghm.txt"), header = TRUE)
gpi_data <- read.table(paste0(path_folder, "02_macroecological_socioeconomic/richness_coordinate_gpi.txt"), header = TRUE)
insu_data <- read.csv(paste0(path_folder, "02_macroecological_socioeconomic/richness_coordinate_insularity.csv"), header = TRUE)
tdwg_data <- read.table(paste0(path_folder, "02_macroecological_socioeconomic/richness_coordinate_tdwglevelall.txt"), header = TRUE)

lmm_data <- cbind(tdwg_data, clim_data[, c(2: 5)], gdpghm_data[, c(2: 3)], gpi_data[, 12], insu_data[, 11])
colnames(lmm_data)[18] <- c("GDP")
colnames(lmm_data)[19] <- c("GHM")
colnames(lmm_data)[20] <- c("GPI")
colnames(lmm_data)[21] <- c("insularity")
write.table(lmm_data, file = paste0(path_folder, "02_macroecological_socioeconomic/linearmixedeffectmodel_data.txt"),
            col.names = TRUE, row.names = FALSE)

lmm_data <- read.table(paste0(path_folder, "02_macroecological_socioeconomic/linearmixedeffectmodel_data.txt"), header = TRUE)

lmm_data$Level4_cod <- as.factor(lmm_data$Level4_cod)
lmm_data$Level3_cod <- as.factor(lmm_data$Level3_cod)
lmm_data$Level2_cod <- as.factor(lmm_data$Level2_cod)
lmm_data$Level1_cod <- as.factor(lmm_data$Level1_cod)
lmm_data$GDP <- as.numeric(lmm_data$GDP)
lmm_data$insularity <- factor(lmm_data$insularity, levels = c("1", "0")) # island is the reference level

##################
# nonnative tree
##################
lmm_nonnative <- lmm_data[, c(2: 3, 5, 10: 21)]
lmm_nonnative <- lmm_nonnative[which(lmm_nonnative$nonnative_ri != 0), ]

# ln transform for richness
lmm_nonnative$nonnative_ri_ln <- log(lmm_nonnative$nonnative_ri + 1)

# z-score standardized
lmm_nonnative$clim_1_sta <- scale(lmm_nonnative$clim_1)
lmm_nonnative$clim_7_sta <- scale(lmm_nonnative$clim_7)
lmm_nonnative$clim_12_sta <- scale(lmm_nonnative$clim_12)
lmm_nonnative$clim_15_sta <- scale(lmm_nonnative$clim_15)
lmm_nonnative$GDP_sta <- scale(lmm_nonnative$GDP)
lmm_nonnative$GHM_sta <- scale(lmm_nonnative$GHM)
lmm_nonnative$GPI_sta <- scale(lmm_nonnative$GPI)

# fit model
fit_nonnative <- lme(nonnative_ri_ln ~ clim_1_sta + clim_7_sta + clim_12_sta + clim_15_sta + GDP_sta + GHM_sta + GPI_sta + insularity,
                     random = ~1| Level1_cod/ Level2_cod/ Level4_cod, data = lmm_nonnative)

###################
# model selection
options(na.action = "na.fail")

select_nonnative <- dredge(fit_nonnative)

# model average models with delta AICc < 2
average_nonnative <- model.avg(select_nonnative, subset = delta < 2)
best_nonnative <- summary(average_nonnative)

##############################
# threatened non-native tree
##############################
lmm_threatenednonnative <- lmm_data[, c(2: 3, 7, 10: 21)]
lmm_threatenednonnative <- lmm_threatenednonnative[which(lmm_threatenednonnative$threatenednonnative_ri != 0), ]

# ln transform for richness
lmm_threatenednonnative$threatenednonnative_ri_ln <- log(lmm_threatenednonnative$threatenednonnative_ri + 1)

# z-score standardized
lmm_threatenednonnative$clim_1_sta <- scale(lmm_threatenednonnative$clim_1)
lmm_threatenednonnative$clim_7_sta <- scale(lmm_threatenednonnative$clim_7)
lmm_threatenednonnative$clim_12_sta <- scale(lmm_threatenednonnative$clim_12)
lmm_threatenednonnative$clim_15_sta <- scale(lmm_threatenednonnative$clim_15)
lmm_threatenednonnative$GDP_sta <- scale(lmm_threatenednonnative$GDP)
lmm_threatenednonnative$GHM_sta <- scale(lmm_threatenednonnative$GHM)
lmm_threatenednonnative$GPI_sta <- scale(lmm_threatenednonnative$GPI)

# fit model
fit_threatenednonnative <- lme(threatenednonnative_ri_ln ~ clim_1_sta + clim_7_sta + clim_12_sta + clim_15_sta + GDP_sta + GHM_sta + GPI_sta + insularity,
                               random = ~1| Level1_cod/ Level2_cod/ Level4_cod, data = lmm_threatenednonnative)

###################
# model selection
options(na.action = "na.fail")

select_threatenednonnative <- dredge(fit_threatenednonnative)

# model average models with delta AICc < 2
# average_threatenednonnative <- model.avg(select_threatenednonnative, subset = delta < 2)
# best_threatenednonnative <- summary(average_threatenednonnative)

summary(get.models(select_threatenednonnative, 1)[[1]])

##################################
# non-threatened non-native tree
##################################
lmm_nonthreatenednonnative <- lmm_data[, c(2: 3, 9, 10: 21)]
lmm_nonthreatenednonnative <- lmm_nonthreatenednonnative[which(lmm_nonthreatenednonnative$nonthreatenednonnative_ri != 0), ]

# ln transform for richness
lmm_nonthreatenednonnative$nonthreatenednonnative_ri_ln <- log(lmm_nonthreatenednonnative$nonthreatenednonnative_ri + 1)

# z-score standardized
lmm_nonthreatenednonnative$clim_1_sta <- scale(lmm_nonthreatenednonnative$clim_1)
lmm_nonthreatenednonnative$clim_7_sta <- scale(lmm_nonthreatenednonnative$clim_7)
lmm_nonthreatenednonnative$clim_12_sta <- scale(lmm_nonthreatenednonnative$clim_12)
lmm_nonthreatenednonnative$clim_15_sta <- scale(lmm_nonthreatenednonnative$clim_15)
lmm_nonthreatenednonnative$GDP_sta <- scale(lmm_nonthreatenednonnative$GDP)
lmm_nonthreatenednonnative$GHM_sta <- scale(lmm_nonthreatenednonnative$GHM)
lmm_nonthreatenednonnative$GPI_sta <- scale(lmm_nonthreatenednonnative$GPI)

# fit model
fit_nonthreatenednonnative <- lme(nonthreatenednonnative_ri_ln ~ clim_1_sta + clim_7_sta + clim_12_sta + clim_15_sta + GDP_sta + GHM_sta + GPI_sta + insularity,
                                  random = ~1| Level1_cod/ Level2_cod/ Level4_cod, data = lmm_nonthreatenednonnative)

###################
# model selection
options(na.action = "na.fail")

select_nonthreatenednonnative <- dredge(fit_nonthreatenednonnative)

# model average models with delta AICc < 2
average_nonthreatenednonnative <- model.avg(select_nonthreatenednonnative, subset = delta < 2)
best_nonthreatenednonnative <- summary(average_nonthreatenednonnative)

##############################################
# spatial autocorrelation (of the residuals)
##############################################

# non-native tree species
# fit the best model
fit_nonnative <- lme(nonnative_ri_ln ~ clim_1_sta + clim_7_sta + clim_12_sta + clim_15_sta + GDP_sta + GHM_sta + GPI_sta + insularity,
                     random = ~1| Level1_cod/ Level2_cod/ Level4_cod, data = lmm_nonnative)

resi_nonnative <- residuals(fit_nonnative)

# Residual correlograms
resi_nonnative <- data.frame(longitude = lmm_nonnative$cellx,
                             latitude = lmm_nonnative$celly,
                             model_resi = residuals(fit_nonnative))

# plot the correlograms of Moran's I
cor_nonnative <- spline.correlog(x = resi_nonnative$longitude,
                                 y = resi_nonnative$latitude,
                                 z = resi_nonnative$model_resi,
                                 resamp = 200,
                                 latlon = TRUE)

save(cor_nonnative,
     file = paste0(path_folder, "02_macroecological_socioeconomic/correlograms_nonnativetree.Rdata"))

load(paste0(path_folder, "02_macroecological_socioeconomic/correlograms_nonnativetree.Rdata"))

setwd(path_folder)
png(file = paste0(path_folder, "02_macroecological_socioeconomic/sac_nonnativetree.png"),
    bg = "transparent", width = 8, height = 6, units = "cm", res = 300)

real_nonnative <- as.data.frame(cbind(as.numeric(cor_nonnative$real$predicted$x),
                                      as.numeric(cor_nonnative$real$predicted$y)))
se_nonnative <- as.data.frame(cbind(as.numeric(cor_nonnative$boot$boot.summary$predicted$x),
                                    as.numeric(cor_nonnative$boot$boot.summary$predicted$y[2, ]),
                                    as.numeric(cor_nonnative$boot$boot.summary$predicted$y[10, ])))

ggplot() +
  geom_line(data = real_nonnative, aes(V1, V2),linewidth = 0.1) +
  geom_line(data = se_nonnative, aes(V1, V2),linewidth = 0.1) +
  geom_line(data = se_nonnative, aes(V1, V3),linewidth = 0.1) +
  geom_hline(yintercept = 0, linewidth = 0.1) +
  scale_y_continuous(labels = c(sprintf("%.1f", c(-1, -0.5)), 0, sprintf("%.1f", c(0.5, 1))),
                     limits = c(-1, 1)) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA)) +
  ylab("Correlation") +
  xlab("Distance")

dev.off()

# significant level
# width of the correlogram increment [km]
increment = 1
# how many times to resample for significance testing
resample = 200

cor_nonnative_significant <- ncf::correlog(x = resi_nonnative$longitude,
                                           y = resi_nonnative$latitude,
                                           z = resi_nonnative$model_resi,
                                           latlon = TRUE,
                                           resamp = resample,
                                           increment = increment)

save(cor_nonnative_significant,
     file = paste0(path_folder, "02_macroecological_socioeconomic/correlograms_significant_nonnativetree.Rdata"))

######################################
# threatened non-native tree species
# fit the best model
fit_threatenednonnative <- lme(threatenednonnative_ri_ln ~ clim_1_sta + clim_7_sta + clim_12_sta + clim_15_sta + GDP_sta + GHM_sta,
                               random = ~1| Level1_cod/ Level2_cod/ Level4_cod, data = lmm_threatenednonnative)

resi_threatenednonnative <- residuals(fit_threatenednonnative)

# Residual correlograms
resi_threatenednonnative <- data.frame(longitude = lmm_threatenednonnative$cellx,
                                       latitude = lmm_threatenednonnative$celly,
                                       model_resi = residuals(fit_threatenednonnative))

# plot the correlograms of Moran's I
cor_threatenednonnative <- spline.correlog(x = resi_threatenednonnative$longitude,
                                           y = resi_threatenednonnative$latitude,
                                           z = resi_threatenednonnative$model_resi,
                                           resamp = 200,
                                           latlon = TRUE)

save(cor_threatenednonnative,
     file = paste0(path_folder, "02_macroecological_socioeconomic/correlograms_nonnativetree_threatened.Rdata"))

load(paste0(path_folder, "02_macroecological_socioeconomic/correlograms_nonnativetree_threatened.Rdata"))

setwd(path_folder)
png(file = paste0(path_folder, "02_macroecological_socioeconomic/sac_nonnativetree_threatened.png"),
    bg = "transparent", width = 8, height = 6, units = "cm", res = 300)

real_threatenednonnative <- as.data.frame(cbind(as.numeric(cor_threatenednonnative$real$predicted$x),
                                                as.numeric(cor_threatenednonnative$real$predicted$y)))
se_threatenednonnative <- as.data.frame(cbind(as.numeric(cor_threatenednonnative$boot$boot.summary$predicted$x),
                                              as.numeric(cor_threatenednonnative$boot$boot.summary$predicted$y[2, ]),
                                              as.numeric(cor_threatenednonnative$boot$boot.summary$predicted$y[10, ])))

ggplot() +
  geom_line(data = real_threatenednonnative, aes(V1, V2),linewidth = 0.1) +
  geom_line(data = se_threatenednonnative, aes(V1, V2),linewidth = 0.1) +
  geom_line(data = se_threatenednonnative, aes(V1, V3),linewidth = 0.1) +
  geom_hline(yintercept = 0, linewidth = 0.1) +
  scale_y_continuous(labels = c(sprintf("%.1f", c(-1, -0.5)), 0, sprintf("%.1f", c(0.5, 1))),
                     limits = c(-1, 1)) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA)) +
  ylab("Correlation") +
  xlab("Distance")

dev.off()

# significant level
# width of the correlogram increment [km]
increment = 1
# how many times to resample for significance testing
resample = 200

cor_threatenednonnative_significant <- ncf::correlog(x = resi_threatenednonnative$longitude,
                                                     y = resi_threatenednonnative$latitude,
                                                     z = resi_threatenednonnative$model_resi,
                                                     latlon = TRUE,
                                                     resamp = resample,
                                                     increment = increment)

save(cor_threatenednonnative_significant,
     file = paste0(path_folder, "02_macroecological_socioeconomic/correlograms_significant_nonnativetree_threatened.Rdata"))

##########################################
# non-threatened non-native tree species
# fit the best model
fit_nonthreatenednonnative <- lme(nonthreatenednonnative_ri_ln ~ clim_1_sta + clim_7_sta + clim_12_sta + clim_15_sta + GDP_sta + GHM_sta + GPI_sta + insularity,
                                  random = ~1| Level1_cod/ Level2_cod/ Level4_cod, data = lmm_nonthreatenednonnative)

resi_nonthreatenednonnative <- residuals(fit_nonthreatenednonnative)

# Residual correlograms
resi_nonthreatenednonnative <- data.frame(longitude = lmm_nonthreatenednonnative$cellx,
                                          latitude = lmm_nonthreatenednonnative$celly,
                                          model_resi = residuals(fit_nonthreatenednonnative))

# plot the correlograms of Moran's I
cor_nonthreatenednonnative <- spline.correlog(x = resi_nonthreatenednonnative$longitude,
                                              y = resi_nonthreatenednonnative$latitude,
                                              z = resi_nonthreatenednonnative$model_resi,
                                              resamp = 200,
                                              latlon = TRUE)

save(cor_nonthreatenednonnative,
     file = paste0(path_folder, "02_macroecological_socioeconomic/correlograms_nonnativetree_nonthreatened.Rdata"))

load(paste0(path_folder, "02_macroecological_socioeconomic/correlograms_nonnativetree_nonthreatened.Rdata"))

setwd(path_folder)
png(file = paste0(path_folder, "02_macroecological_socioeconomic/sac_nonnativetree_nonthreatened.png"),
    bg = "transparent", width = 8, height = 6, units = "cm", res = 300)

real_nonthreatenednonnative <- as.data.frame(cbind(as.numeric(cor_nonthreatenednonnative$real$predicted$x),
                                                   as.numeric(cor_nonthreatenednonnative$real$predicted$y)))
se_nonthreatenednonnative <- as.data.frame(cbind(as.numeric(cor_nonthreatenednonnative$boot$boot.summary$predicted$x),
                                                 as.numeric(cor_nonthreatenednonnative$boot$boot.summary$predicted$y[2, ]),
                                                 as.numeric(cor_nonthreatenednonnative$boot$boot.summary$predicted$y[10, ])))

ggplot() +
  geom_line(data = real_nonthreatenednonnative, aes(V1, V2),linewidth = 0.1) +
  geom_line(data = se_nonthreatenednonnative, aes(V1, V2),linewidth = 0.1) +
  geom_line(data = se_nonthreatenednonnative, aes(V1, V3),linewidth = 0.1) +
  geom_hline(yintercept = 0, linewidth = 0.1) +
  scale_y_continuous(labels = c(sprintf("%.1f", c(-1, -0.5)), 0, sprintf("%.1f", c(0.5, 1))),
                     limits = c(-1, 1)) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA)) +
  ylab("Correlation") +
  xlab("Distance")

dev.off()

# significant level
# width of the correlogram increment [km]
increment = 1
# how many times to resample for significance testing
resample = 200

cor_nonthreatenednonnative_significant <- ncf::correlog(x = resi_nonthreatenednonnative$longitude,
                                                        y = resi_nonthreatenednonnative$latitude,
                                                        z = resi_nonthreatenednonnative$model_resi,
                                                        latlon = TRUE,
                                                        resamp = resample,
                                                        increment = increment)

save(cor_nonthreatenednonnative_significant,
     file = paste0(path_folder, "02_macroecological_socioeconomic/correlograms_significant_nonnativetree_nonthreatened.Rdata"))







