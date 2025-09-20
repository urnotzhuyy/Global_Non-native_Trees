rm(list=ls())

# loading required libraries
library(ggplot2)
library(rcompanion)
library(extrafont)

# path
path_folder <- "D:/non_native_tree/"

# load the phylogenetic distance file
pd_list<-list.files(paste0(path_folder,"04_phylogenetic_distance/pairwise_phylogenetic_distance/"),full.names = TRUE)

for (id in 1:length(pd_list)) {
  
  pd_file<-read.table(pd_list[id],header = TRUE)
  
  write.table(pd_file,file = paste0(path_folder,"04_phylogenetic_distance/phylogenetic_distance_global.txt"),
              col.names = FALSE,row.names = FALSE,append = TRUE)
  
  print(id)
  
}

pd_global<-read.table(paste0(path_folder,"04_phylogenetic_distance/phylogenetic_distance_global.txt"))
colnames(pd_global)<-c("species","probability","pdmean","pdmin","region")
write.table(pd_global,file = paste0(path_folder,"04_phylogenetic_distance/phylogenetic_distance_global.txt"),
            col.names = TRUE,row.names = FALSE)

pd_global<-read.table(paste0(path_folder,"04_phylogenetic_distance/phylogenetic_distance_global.txt"),header = TRUE)

# z-score standardized
pd_global$pdmean_sta<-scale(pd_global$pdmean)
pd_global$pdmin_sta<-scale(pd_global$pdmin)

# binomial glm (clog-log link function), including the quadratic term for phylogenetic distance
# PD mean
pdmean_glm<-glm(probability ~ pdmean_sta + I(pdmean_sta^2),
            family = binomial(link = "cloglog"),data = pd_global)

summary(pdmean_glm)

nagelkerke(pdmean_glm)

# plot
setwd(path_folder)
png(file = paste0(path_folder, "04_phylogenetic_distance/phylogenetic_distance_mean.png"),
    bg = "transparent", width = 7, height = 7, units = "cm", res = 300)

ggplot(pd_global, aes(x = pdmean_sta,y = probability)) +
  geom_smooth(method = "glm",
              formula = y ~ x + I(x^2),
              method.args = list(family = binomial(link = "cloglog")),
              se = TRUE,
              colour = "black",
              linewidth = 0.5) +
  scale_y_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5),
                     labels = c(0, sprintf("%.3f", c(0.1, 0.2, 0.3, 0.4, 0.5)))) +
  scale_x_continuous(breaks = seq(-4, 6, 2)) +
  theme(panel.background = element_blank(),
        axis.line = element_line(),
        axis.text = element_text(family = "Times New Roman", size = 10),
        axis.title = element_text(family = "Times New Roman", size = 12)) +
  xlab("MDNS") +
  ylab("Probability of being introduced")

dev.off()

##########
# PD min
pdmin_glm<-glm(probability ~ pdmin_sta + I(pdmin_sta^2),
            family = binomial(link = "cloglog"),data = pd_global)

summary(pdmin_glm)

nagelkerke(pdmin_glm)

# plot
setwd(path_folder)
png(file = paste0(path_folder, "04_phylogenetic_distance/phylogenetic_distance_min.png"),
    bg = "transparent", width = 7, height = 7, units = "cm", res = 300)

ggplot(pd_global, aes(x = pdmin_sta,y = probability)) +
  geom_smooth(method = "glm",
              formula = y ~ x + I(x^2),
              method.args = list(family = binomial(link = "cloglog")),
              se = TRUE,
              colour = "black",
              linewidth = 0.5) +
  scale_y_continuous(breaks = c(0, 0.005, 0.010, 0.015, 0.020, 0.025),
                     labels = c(0, sprintf("%.3f", c(0.005, 0.010, 0.015, 0.020, 0.025)))) +
  theme(panel.background = element_blank(),
        axis.line = element_line(),
        axis.text = element_text(family = "Times New Roman", size = 10),
        axis.title = element_text(family = "Times New Roman", size = 12)) +
  xlab("DNNS") +
  ylab("Probability of being introduced")

dev.off()

###################
# threatened tree
###################

# load the phylogenetic distance file for threatened tree
threatenedpd_list<-list.files(paste0(path_folder,"04_phylogenetic_distance/pairwise_phylogenetic_distance_threatened/"),full.names = TRUE)

for (id in 1:length(threatenedpd_list)) {
  
  threatenedpd_file<-read.table(threatenedpd_list[id],header = TRUE)
  
  write.table(threatenedpd_file,file = paste0(path_folder,"04_phylogenetic_distance/phylogenetic_distance_global_threatened.txt"),
              col.names = FALSE,row.names = FALSE,append = TRUE)
  
  print(id)
  
}

threatenedpd_global<-read.table(paste0(path_folder,"04_phylogenetic_distance/phylogenetic_distance_global_threatened.txt"))
colnames(threatenedpd_global)<-c("species","probability","pdmean","pdmin","region")
write.table(threatenedpd_global,file = paste0(path_folder,"04_phylogenetic_distance/phylogenetic_distance_global_threatened.txt"),
            col.names = TRUE,row.names = FALSE)

threatenedpd_global<-read.table(paste0(path_folder,"04_phylogenetic_distance/phylogenetic_distance_global_threatened.txt"),header = TRUE)

# z-score standardized
threatenedpd_global$pdmean_sta<-scale(threatenedpd_global$pdmean)
threatenedpd_global$pdmin_sta<-scale(threatenedpd_global$pdmin)

# binomial glm (clog-log link function), including the quadratic term for phylogenetic distance
# PD mean
threatenedpdmean_glm<-glm(probability ~ pdmean_sta + I(pdmean_sta^2),
                family = binomial(link = "cloglog"),data = threatenedpd_global)

summary(threatenedpdmean_glm)

nagelkerke(threatenedpdmean_glm)

# plot
setwd(path_folder)
png(file = paste0(path_folder, "04_phylogenetic_distance/phylogenetic_distance_mean_threatened.png"),
    bg = "transparent", width = 7, height = 7, units = "cm", res = 300)

ggplot(threatenedpd_global, aes(x = pdmean_sta,y = probability)) +
  geom_smooth(method = "glm",
              formula = y ~ x + I(x^2),
              method.args = list(family = binomial(link = "cloglog")),
              se = TRUE,
              colour = "black",
              linewidth = 0.5) +
  scale_y_continuous(breaks = c(0, 0.01, 0.02, 0.03, 0.04),
                     labels = c(0, sprintf("%.3f", c(0.01, 0.02, 0.03, 0.04)))) +
  scale_x_continuous(breaks = c(-2, 0, 2, 4, 6)) +
  theme(panel.background = element_blank(),
        axis.line = element_line(),
        axis.text = element_text(family = "Times New Roman", size = 10),
        axis.title = element_text(family = "Times New Roman", size = 12)) +
  xlab("MDNS") +
  ylab("Probability of being introduced")

dev.off()

##########
# PD min

threatenedpdmin_glm<-glm(probability ~ pdmin_sta + I(pdmin_sta^2),
               family = binomial(link = "cloglog"),data = threatenedpd_global)

summary(threatenedpdmin_glm)

nagelkerke(threatenedpdmin_glm)

# plot
setwd(path_folder)
png(file = paste0(path_folder, "04_phylogenetic_distance/phylogenetic_distance_min_threatened.png"),
    bg = "transparent", width = 7, height = 7, units = "cm", res = 300)

ggplot(threatenedpd_global, aes(x = pdmin_sta,y = probability)) +
  geom_smooth(method = "glm",
              formula = y ~ x + I(x^2),
              method.args = list(family = binomial(link = "cloglog")),
              se = TRUE,
              colour = "black",
              linewidth = 0.5) +
  scale_y_continuous(breaks = c(0, 0.005, 0.01, 0.015),
                     labels = c(0, sprintf("%.3f", c(0.005, 0.01, 0.015)))) +
  theme(panel.background = element_blank(),
        axis.line = element_line(),
        axis.text = element_text(family = "Times New Roman", size = 10),
        axis.title = element_text(family = "Times New Roman", size = 12)) +
  xlab("DNNS") +
  ylab("Probability of being introduced")

dev.off()

#######################
# non-threatened tree
#######################

# load the phylogenetic distance file for non-threatened tree
nonthreatenedpd_list<-list.files(paste0(path_folder,"04_phylogenetic_distance/pairwise_phylogenetic_distance_nonthreatened/"),full.names = TRUE)

for (id in 1:length(nonthreatenedpd_list)) {
  
  nonthreatenedpd_file<-read.table(nonthreatenedpd_list[id],header = TRUE)
  
  write.table(nonthreatenedpd_file,file = paste0(path_folder,"04_phylogenetic_distance/phylogenetic_distance_global_nonthreatened.txt"),
              col.names = FALSE,row.names = FALSE,append = TRUE)
  
  print(id)
  
}

nonthreatenedpd_global<-read.table(paste0(path_folder,"04_phylogenetic_distance/phylogenetic_distance_global_nonthreatened.txt"))
colnames(nonthreatenedpd_global)<-c("species","probability","pdmean","pdmin","region")
write.table(nonthreatenedpd_global,file = paste0(path_folder,"04_phylogenetic_distance/phylogenetic_distance_global_nonthreatened.txt"),
            col.names = TRUE,row.names = FALSE)

nonthreatenedpd_global<-read.table(paste0(path_folder,"04_phylogenetic_distance/phylogenetic_distance_global_nonthreatened.txt"),header = TRUE)

# z-score standardized
nonthreatenedpd_global$pdmean_sta<-scale(nonthreatenedpd_global$pdmean)
nonthreatenedpd_global$pdmin_sta<-scale(nonthreatenedpd_global$pdmin)

# binomial glm (clog-log link function), including the quadratic term for phylogenetic distance
# PD mean
nonthreatenedpdmean_glm<-glm(probability ~ pdmean_sta + I(pdmean_sta^2),
                          family = binomial(link = "cloglog"),data = nonthreatenedpd_global)

summary(nonthreatenedpdmean_glm)

nagelkerke(nonthreatenedpdmean_glm)

# plot
setwd(path_folder)
png(file = paste0(path_folder, "04_phylogenetic_distance/phylogenetic_distance_mean_nonthreatened.png"),
    bg = "transparent", width = 7, height = 7, units = "cm", res = 300)

ggplot(nonthreatenedpd_global, aes(x = pdmean_sta,y = probability)) +
  geom_smooth(method = "glm",
              formula = y ~ x + I(x^2),
              method.args = list(family = binomial(link = "cloglog")),
              se = TRUE,
              colour = "black",
              linewidth = 0.5) +
  scale_y_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5),
                     labels = c(0, sprintf("%.3f", c(0.1, 0.2, 0.3, 0.4, 0.5)))) +
  scale_x_continuous(breaks = seq(-4, 6, 2)) +
  theme(panel.background = element_blank(),
        axis.line = element_line(),
        axis.text = element_text(family = "Times New Roman", size = 10),
        axis.title = element_text(family = "Times New Roman", size = 12)) +
  xlab("MDNS") +
  ylab("Probability of being introduced")

dev.off()

##########
# PD min

nonthreatenedpdmin_glm<-glm(probability ~ pdmin_sta + I(pdmin_sta^2),
                         family = binomial(link = "cloglog"),data = nonthreatenedpd_global)

summary(nonthreatenedpdmin_glm)

nagelkerke(nonthreatenedpdmin_glm)

# plot
setwd(path_folder)
png(file = paste0(path_folder, "04_phylogenetic_distance/phylogenetic_distance_min_nonthreatened.png"),
    bg = "transparent", width = 7, height = 7, units = "cm", res = 300)

ggplot(nonthreatenedpd_global, aes(x = pdmin_sta,y = probability)) +
  geom_smooth(method = "glm",
              formula = y ~ x + I(x^2),
              method.args = list(family = binomial(link = "cloglog")),
              se = TRUE,
              colour = "black",
              linewidth = 0.5) +
  scale_y_continuous(breaks = c(0, 0.01, 0.02, 0.03),
                     labels = c(0, sprintf("%.3f", c(0.01, 0.02, 0.03)))) +
  theme(panel.background = element_blank(),
        axis.line = element_line(),
        axis.text = element_text(family = "Times New Roman", size = 10),
        axis.title = element_text(family = "Times New Roman", size = 12)) +
  xlab("DNNS") +
  ylab("Probability of being introduced")

dev.off()




