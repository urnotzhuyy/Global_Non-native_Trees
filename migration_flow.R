rm(list=ls())

# loading required libraries
library(circlize)
library(extrafont)

# path
path_folder <- "D:/non_native_tree/"

####################
# data preparation
####################

# native tree
native_list<-list.files(paste0(path_folder,"aggregate_data/tree_data_native_final/"),pattern = ".txt")
native_list<-sub("_native.txt","",native_list,fixed = TRUE)

native_bgci_list<-list.files(paste0(path_folder,"aggregate_data/tree_data_native_final/bgci_part/"))
native_bgci_list<-sub("_native.txt","",native_bgci_list,fixed = TRUE)

native_all<-c(native_list,native_bgci_list) # 51218

# non-native tree
nonnative_list<-list.files(paste0(path_folder,"aggregate_data/tree_data_nonnative_final/"),pattern = ".txt")
nonnative_list<-sub("_nonnative.txt","",nonnative_list,fixed = TRUE)

nonnative_bgci_list<-list.files(paste0(path_folder,"aggregate_data/tree_data_nonnative_final/bgci_part/"))
nonnative_bgci_list<-sub("_nonnative.txt","",nonnative_bgci_list,fixed = TRUE)

nonnative_all<-c(nonnative_list,nonnative_bgci_list) # 27417

# 921 non-native tree donnot have native records, subsequent analysis only use the 26496 non-native tree
nonnative_all<-nonnative_all[nonnative_all%in%native_all] # 26496

##################################################################
# match the GADM and TDWG region of non-native tree in BGCI list
##################################################################

for (id in 1:length(nonnative_bgci_list)) {
  
  sp<-nonnative_bgci_list[id]
  
  if(sp%in%nonnative_all){
    
    sp_file_nonnative<-read.table(paste0(path_folder,"aggregate_data/tree_data_nonnative_final/bgci_part/",sp,"_nonnative.txt"),header = TRUE)
    sp_file_native<-read.table(paste0(path_folder,"aggregate_data/tree_data_native_final/bgci_part/",sp,"_native.txt"),header = TRUE)
    
    sp_file<-rbind(sp_file_nonnative,sp_file_native)
    
    sp_region<-sp_file[,c(4,5)]
    sp_region<-sp_region[!duplicated(sp_region),]
    
    write.table(sp_region,file = paste0(path_folder,"03_migration_flow/gadm_tdwg_match_bgci_nonnativetree.txt"),
                col.names = FALSE,row.names = FALSE,append = TRUE)
    
  }
  
  print(id)
  
}

match_file<-read.table(paste0(path_folder,"03_migration_flow/gadm_tdwg_match_bgci_nonnativetree.txt"))
colnames(match_file)<-c("GID_0","NAME_0")
match_file<-match_file[!duplicated(match_file),]

write.table(match_file,file = paste0(path_folder,"03_migration_flow/gadm_tdwg_match_bgci_nonnativetree.txt"),
            col.names = TRUE,row.names = FALSE)

##################################
# non-native tree migration flow
##################################
# load the GADM and TDWG matched file
match_mod<-read.csv(paste0(path_folder,"03_migration_flow/gadm_tdwg_match_bgci_nonnativetree_modify.csv"))

# form the migration flow of each non-native tree
for (id in 1:length(nonnative_all)) {
  
  # species name
  sp<-nonnative_all[id]
  
  # extract non-native continent (destination) of species
  if(sp%in%nonnative_list){
    # POWO and IUCN distribution information
    # native data (origin continent)
    sp_file_ori<-read.table(paste0(path_folder,"aggregate_data/tree_data_native_final/",sp,"_native.txt"),header = TRUE)
    # origin
    sp_ori<-unique(sp_file_ori$Level1_cod)
  
    # non-native data (destination continent)
    sp_file_des<-read.table(paste0(path_folder,"aggregate_data/tree_data_nonnative_final/",sp,"_nonnative.txt"),header = TRUE)
    # destination
    sp_des<-unique(sp_file_des$Level1_cod)
    
    # form the migration flow information
    for (oriid in 1:length(sp_ori)) {
      
      ori_name<-sp_ori[oriid]
      
      for (desid in 1:length(sp_des)) {
        
        des_name<-sp_des[desid]
        
        sp_flow<-cbind(sp,ori_name,des_name)
        
        write.table(sp_flow,file = paste0(path_folder,"03_migration_flow/tree_flow/",sp,"_nonnativetree_flow.txt"),
                    col.names = FALSE,row.names = FALSE,append = TRUE)
        
      }
      
    }
    
  }else{
    # BGCI distribution information
    # native data (origin continent)
    sp_file_ori<-read.table(paste0(path_folder,"aggregate_data/tree_data_native_final/bgci_part/",sp,"_native.txt"),header = TRUE)
    # origin
    sp_ori<-unique(as.character(sp_file_ori$GID_0))
    sp_ori<-unique(match_mod$TDWG_level1_cod[which(match_mod$GID_0%in%sp_ori)])
    
    # non-native data (destination continent)
    sp_file_des<-read.table(paste0(path_folder,"aggregate_data/tree_data_nonnative_final/bgci_part/",sp,"_nonnative.txt"),header = TRUE)
    # destination
    sp_des<-unique(as.character(sp_file_des$GID_0))
    sp_des<-unique(match_mod$TDWG_level1_cod[which(match_mod$GID_0%in%sp_des)])
    
    # form the migration flow information
    for (oriid in 1:length(sp_ori)) {
      
      ori_name<-sp_ori[oriid]
      
      for (desid in 1:length(sp_des)) {
        
        des_name<-sp_des[desid]
        
        sp_flow<-cbind(sp,ori_name,des_name)
        
        write.table(sp_flow,file = paste0(path_folder,"03_migration_flow/tree_flow/",sp,"_nonnativetree_flow.txt"),
                    col.names = FALSE,row.names = FALSE,append = TRUE)
        
      }
      
    }
    
  }
  
  print(id)
  
}

# add the column names of species flow files
for (id in 1:length(nonnative_all)) {
  
  sp<-nonnative_all[id]
  
  sp_flow<-read.table(paste0(path_folder,"03_migration_flow/tree_flow/",sp,"_nonnativetree_flow.txt"))
  colnames(sp_flow)<-c("spacc","origin","destination")
  
  write.table(sp_flow,file = paste0(path_folder,"03_migration_flow/tree_flow/",sp,"_nonnativetree_flow.txt"),
              col.names = TRUE,row.names = FALSE)
  
  print(id)
  
}

##################################
# form the continental flow file
for (id in 1:length(nonnative_all)) {
  
  sp<-nonnative_all[id]
  
  sp_flow<-read.table(paste0(path_folder,"03_migration_flow/tree_flow/",sp,"_nonnativetree_flow.txt"),header = TRUE)
  
  write.table(sp_flow,file = paste0(path_folder,"03_migration_flow/tree_flow_all.txt"),
              col.names = FALSE,row.names = FALSE,append = TRUE)
  
  print(id)
  
}

tree_flow<-read.table(paste0(path_folder,"03_migration_flow/tree_flow_all.txt"))
colnames(tree_flow)<-c("spacc","origin","destination")
write.table(tree_flow,file = paste0(path_folder,"03_migration_flow/tree_flow_all.txt"),
            col.names = TRUE,row.names = FALSE)

tree_flow<-read.table(paste0(path_folder,"03_migration_flow/tree_flow_all.txt"),header = TRUE)

# TDWG Level 1 code
region_list<-c("Europe","Africa","Temperate Asia","Tropical Asia","Australasia","Pacific","North America","South America")

for (oriid in 1:length(region_list)) {
  
  ori_code<-oriid
  
  ori_reg<-region_list[oriid]
  
  ori_file<-tree_flow[which(tree_flow$origin==ori_code),]
  
  for (desid in 1:length(region_list)) {
    
    des_code<-desid
    
    des_reg<-region_list[desid]
    
    des_file<-ori_file[which(ori_file$destination==des_code),]
    
    if(dim(des_file)[1]!=0){
      
      flow_num<-dim(des_file)[1]
      
    }else{
      
      flow_num<-0
      
    }
    
    flow_data<-cbind(ori_reg,des_reg,flow_num)
    
    write.table(flow_data,file = paste0(path_folder,"03_migration_flow/migration_flow_among_continent.txt"),
                col.names = FALSE,row.names = FALSE,append = TRUE)
    
  }
  
  print(oriid)
  
}

flow_data<-read.table(paste0(path_folder,"03_migration_flow/migration_flow_among_continent.txt"))
colnames(flow_data)<-c("origin","destination","flow")
write.table(flow_data,file = paste0(path_folder,"03_migration_flow/migration_flow_among_continent.txt"),
            col.names = TRUE,row.names = FALSE)

flow_data<-read.table(paste0(path_folder,"03_migration_flow/migration_flow_among_continent.txt"),header = TRUE)

# select the top 50% of migration flows
flow_data<-flow_data[order(flow_data$flow),]
flow_data<-flow_data[33:64,]

# plot
setwd(path_folder)
png(file = paste0(path_folder,"03_migration_flow/migrationflow_nonnativetree.png"),
    bg="transparent",width = 8,height = 8,units = "cm",res = 300)

circos.clear()

circos.par(start.degree=90,
           gap.after= 5,
           track.margin=c(0.01,0.01))

chordDiagram(x = flow_data, 
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
             directional = 1, #1 means the direction is from the first column to the second column 
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
                           minor.ticks = 0,major.at = seq(0,33000,1000),
                           lwd = 1,major.tick.length = mm_y(0.5))
             })

dev.off()

#############################################
# threatened non-native tree migration flow
#############################################

# load the threatened list file (Critically Endangered, Endangered and Vulnerable)
threatened_list<-read.table(paste0(path_folder,"threatened_data/threatened_species_list.txt"),header = TRUE)
# 23950 threatened trees
threatened_all<-unique(threatened_list$species_accept[which(threatened_list$iucn_category%in%c("Critically Endangered","Endangered","Vulnerable"))])

native_threatened<-threatened_all[threatened_all%in%native_all] # 11385
nonnative_threatened<-threatened_all[threatened_all%in%nonnative_all] # 3626

#########################################################
# extract the threatened non-native tree migration flow

tree_flow<-read.table(paste0(path_folder,"03_migration_flow/tree_flow_all.txt"),header = TRUE)

threatenedtree_flow<-tree_flow[which(tree_flow$spacc%in%nonnative_threatened),]

# TDWG Level 1 code
region_list<-c("Europe","Africa","Temperate Asia","Tropical Asia","Australasia","Pacific","North America","South America")

for (oriid in 1:length(region_list)) {
  
  ori_code<-oriid
  
  ori_reg<-region_list[oriid]
  
  ori_file<-threatenedtree_flow[which(threatenedtree_flow$origin==ori_code),]
  
  for (desid in 1:length(region_list)) {
    
    des_code<-desid
    
    des_reg<-region_list[desid]
    
    des_file<-ori_file[which(ori_file$destination==des_code),]
    
    if(dim(des_file)[1]!=0){
      
      flow_num<-dim(des_file)[1]
      
    }else{
      
      flow_num<-0
      
    }
    
    flow_data<-cbind(ori_reg,des_reg,flow_num)
    
    write.table(flow_data,file = paste0(path_folder,"03_migration_flow/migration_flow_among_continent_threatened.txt"),
                col.names = FALSE,row.names = FALSE,append = TRUE)
    
  }
  
  print(oriid)
  
}

threatenedflow_data<-read.table(paste0(path_folder,"03_migration_flow/migration_flow_among_continent_threatened.txt"))
colnames(threatenedflow_data)<-c("origin","destination","flow")
write.table(threatenedflow_data,file = paste0(path_folder,"03_migration_flow/migration_flow_among_continent_threatened.txt"),
            col.names = TRUE,row.names = FALSE)

threatenedflow_data<-read.table(paste0(path_folder,"03_migration_flow/migration_flow_among_continent_threatened.txt"),header = TRUE)

# select the top 50% of migration flows
threatenedflow_data<-threatenedflow_data[order(threatenedflow_data$flow),]
threatenedflow_data<-threatenedflow_data[33:64,]

# plot
setwd(path_folder)
png(file=paste0(path_folder,"03_migration_flow/migrationflow_nonnativetree_threatened.png"),
    bg="transparent",width = 8,height = 8,units = "cm",res = 300)

circos.clear()

circos.par(start.degree=90,
           gap.after= 5,
           track.margin=c(0.01,0.01))

chordDiagram(x = threatenedflow_data,
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
                           minor.ticks = 0,major.at = seq(0,4000,100),
                           lwd = 1,major.tick.length = mm_y(0.5))
             })

dev.off()

#################################################
# non-threatened non-native tree migration flow
#################################################

# load the non-threatened trees
native_nonthreatened<-native_all[!native_all%in%native_threatened] # 39833
nonnative_nonthreatened<-nonnative_all[!nonnative_all%in%nonnative_threatened] # 22870

#############################################################
# extract the non-threatened non-native tree migration flow

tree_flow<-read.table(paste0(path_folder,"03_migration_flow/tree_flow_all.txt"),header = TRUE)

nonthreatenedtree_flow<-tree_flow[which(tree_flow$spacc%in%nonnative_nonthreatened),]

# TDWG Level 1 code
region_list<-c("Europe","Africa","Temperate Asia","Tropical Asia","Australasia","Pacific","North America","South America")

for (oriid in 1:length(region_list)) {
  
  ori_code<-oriid
  
  ori_reg<-region_list[oriid]
  
  ori_file<-nonthreatenedtree_flow[which(nonthreatenedtree_flow$origin==ori_code),]
  
  for (desid in 1:length(region_list)) {
    
    des_code<-desid
    
    des_reg<-region_list[desid]
    
    des_file<-ori_file[which(ori_file$destination==des_code),]
    
    if(dim(des_file)[1]!=0){
      
      flow_num<-dim(des_file)[1]
      
    }else{
      
      flow_num<-0
      
    }
    
    flow_data<-cbind(ori_reg,des_reg,flow_num)
    
    write.table(flow_data,file = paste0(path_folder,"03_migration_flow/migration_flow_among_continent_nonthreatened.txt"),
                col.names = FALSE,row.names = FALSE,append = TRUE)
    
  }
  
  print(oriid)
  
}

nonthreatenedflow_data<-read.table(paste0(path_folder,"03_migration_flow/migration_flow_among_continent_nonthreatened.txt"))
colnames(nonthreatenedflow_data)<-c("origin","destination","flow")
write.table(nonthreatenedflow_data,file = paste0(path_folder,"03_migration_flow/migration_flow_among_continent_nonthreatened.txt"),
            col.names = TRUE,row.names = FALSE)

nonthreatenedflow_data<-read.table(paste0(path_folder,"03_migration_flow/migration_flow_among_continent_nonthreatened.txt"),header = TRUE)

# select the top 50% of migration flows
nonthreatenedflow_data<-nonthreatenedflow_data[order(nonthreatenedflow_data$flow),]
nonthreatenedflow_data<-nonthreatenedflow_data[33:64,]

# plot
setwd(path_folder)
png(file=paste0(path_folder,"03_migration_flow/migrationflow_nonnativetree_nonthreatened.png"),
    bg="transparent",width = 8,height = 8,units = "cm",res = 300)

circos.clear()

circos.par(start.degree=90,
           gap.after= 5,
           track.margin=c(0.01,0.01))

chordDiagram(x = nonthreatenedflow_data,
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
             annotationTrackHeight = c(0.06, 0.07),#外面一圈的宽度
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
                           minor.ticks = 0,major.at = seq(0,33000,1000),
                           lwd = 1,major.tick.length = mm_y(0.5))
             })

dev.off()










