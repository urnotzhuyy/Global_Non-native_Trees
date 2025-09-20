rm(list=ls())

# loading required libraries
library(ape)
library(ggplot2) 
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(tidytree)
library(PNWColors)
library(scales)
library(caper)

# path
path_folder <- "D:/non_native_tree/"

# load the phylogeny and information file
tree_phylo<-read.tree(paste0(path_folder,"04_phylogenetic_distance/global_tree.tre"))
tree_info<-read.table(paste0(path_folder,"04_phylogenetic_distance/global_tree.txt"))

# grouping the species by family
tree_fam<-unique(tree_info$family) # 255 families

tree_cls<-list()
length(tree_cls)<-length(tree_fam)
names(tree_cls)<-tree_fam

for (id in 1:length(tree_fam)) {
  
  fam<-tree_fam[id]
  
  fam_sp<-tree_info$species[which(tree_info$family==fam)]
  fam_sp<-gsub(" ","_",fam_sp,fixed = TRUE)
  
  tree_cls[[id]]<-fam_sp
  
  print(id)
  
}

tree_phylo_info<-groupOTU(tree_phylo,tree_cls)
tree_phylo_info<-as.phylo(tree_phylo_info)

# calculate the percent of non-native tree species per family
# load the list of native and non-native tree species
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

for (id in 1:length(tree_fam)) {
  
  fam<-tree_fam[id]
  
  fam_sp<-tree_info$species[which(tree_info$family==fam)]
  
  # native trees for this family
  fam_sp_native<-length(fam_sp[fam_sp%in%native_all])
  
  # non-native trees for this family
  fam_sp_nonnative<-length(fam_sp[fam_sp%in%nonnative_all])
  
  # percent of clade (for family) that is non-native
  fam_propor<-fam_sp_nonnative/fam_sp_native*100
  
  fam_res<-cbind(fam,fam_sp_native,fam_sp_nonnative,fam_propor)
  write.table(fam_res,file = paste0(path_folder,"04_phylogenetic_distance/percent_nonnative_per_family.txt"),
              col.names = FALSE,row.names = FALSE,append = TRUE)
  
  print(id)
  
}

fam_res<-read.table(paste0(path_folder,"04_phylogenetic_distance/percent_nonnative_per_family.txt"))
colnames(fam_res)<-c("family","native_richness","nonnative_richness","porportion")
write.table(fam_res,file = paste0(path_folder,"04_phylogenetic_distance/percent_nonnative_per_family.txt"),
            col.names = TRUE,row.names = FALSE)

fam_res<-read.table(paste0(path_folder,"04_phylogenetic_distance/percent_nonnative_per_family.txt"),
                    header = TRUE)

###########################
# match the tree and data
tree_tibble<-as_tibble(tree_phylo_info)
tree_tibble_df<-as.data.frame(tree_tibble)

for (id in 1:dim(tree_tibble_df)[1]) {
  
  node_id<-tree_tibble_df[id,2]
  
  node_fam<-as.character(tree_tibble_df[id,5])  
  
  if(node_fam==0){
    
    percent_nonnative<-0
    
  }else{
    
    percent_nonnative<-fam_res$porportion[which(fam_res$family==node_fam)]
    
  }
  
  node_res<-cbind(node_id,node_fam,percent_nonnative)
  write.table(node_res,file = paste0(path_folder,"04_phylogenetic_distance/phylogeny_percent_nonnative_per_family.txt"),
              col.names = FALSE,row.names = FALSE,append = TRUE)
  
  print(id)
  
}

node_info<-read.table(paste0(path_folder,"04_phylogenetic_distance/phylogeny_percent_nonnative_per_family.txt"))
colnames(node_info)<-c("node","family","percent_nonnative")
write.table(node_info,file = paste0(path_folder,"04_phylogenetic_distance/phylogeny_percent_nonnative_per_family.txt"),
            col.names = TRUE,row.names = FALSE)

node_info<-read.table(paste0(path_folder,"04_phylogenetic_distance/phylogeny_percent_nonnative_per_family.txt"),
                      header = TRUE)

node_data<-data.frame(node = node_info$node,trait = node_info$percent_nonnative)
tree_phylo_nonnative<-full_join(tree_phylo,node_data,by = "node")

##############################################
# match the data frame used for outer-circle
# load the threatened list file (Critically Endangered, Endangered and Vulnerable)
threatened_list<-read.table(paste0(path_folder,"threatened_data/threatened_species_list.txt"),header = TRUE)
# 23950 threatened trees
threatened_all<-unique(threatened_list$species_accept[which(threatened_list$iucn_category%in%c("Critically Endangered","Endangered","Vulnerable"))])

native_threatened<-threatened_all[threatened_all%in%native_all] # 11385
nonnative_threatened<-threatened_all[threatened_all%in%nonnative_all] # 3626

# extract the non-threatened trees
native_nonthreatened<-native_all[!native_all%in%native_threatened] # 39833
nonnative_nonthreatened<-nonnative_all[!nonnative_all%in%nonnative_threatened] # 22870

for (id in 1:dim(tree_info)[1]) {
  
  sp<-tree_info[id,1]
  sp_tip<-gsub(" ","_",sp,fixed = TRUE)
  
  # a species is native or nonnative
  if(sp%in%nonnative_all){
    
    sp_dis<-"nonnative"
    
  }else{
    
    sp_dis<-"native"
    
  }
  
  # a species is threatened or not
  if(sp%in%nonnative_all){
    
    # non-native trees
    
    if(sp%in%nonnative_threatened){
      
      sp_status<-"threatened_nonnative"
      
    }else{
      
      sp_status<-"nonthreatened_nonnative"
      
    }
    
  }else{
    
    # native trees
    
    if(sp%in%native_threatened){
      
      sp_status<-"threatened_native"
      
    }else{
      
      sp_status<-"nonthreatened_native"
      
    }
    
  }
  
  sp_fam<-tree_info[id,3]
  
  sp_res<-cbind(sp_tip,sp_dis,sp_status,sp_fam)
  write.table(sp_res,file = paste0(path_folder,"04_phylogenetic_distance/phylogeny_outer_circle_dataframe.txt"),
              col.names = FALSE,row.names = FALSE,append = TRUE)
  
  print(id)
  
}

sp_circle<-read.table(paste0(path_folder,"04_phylogenetic_distance/phylogeny_outer_circle_dataframe.txt"))
colnames(sp_circle)<-c("species","distribution","status","family")
write.table(sp_circle,file = paste0(path_folder,"04_phylogenetic_distance/phylogeny_outer_circle_dataframe.txt"),
            col.names = TRUE,row.names = FALSE)

sp_circle<-read.table(paste0(path_folder,"04_phylogenetic_distance/phylogeny_outer_circle_dataframe.txt"),
                      header = TRUE)

##################
# plot phylogeny
setwd(path_folder)
png(filename = paste0(path_folder,"04_phylogenetic_distance/phylogeny_percent_nonnative_per_family.png"), 
     bg="transparent",width = 20,height = 20,units = "cm",res = 300)

p1<-ggtree(tree_phylo_nonnative,aes(color = trait),layout = "circular",size = 0.1) +
  scale_color_gradientn(colours = c("#5d74a5","#b0cbe7","#fef7c7","#eba07e","#a8554e"))

p2<-p1+new_scale_fill()+
  geom_fruit(data = sp_circle,
             geom = geom_tile,
             mapping = aes(y = species, fill = distribution),
             width = 12,
             offset = 0.03)+
  scale_fill_manual(values = c("native"="#dec5da","nonnative"="#574571"))

p3<-p2+new_scale_fill()+
  geom_fruit(data = sp_circle,
             geom = geom_tile,
             mapping = aes(y = species, fill = status),
             width = 12,
             offset = 0.05)+
  scale_fill_manual(values = c("threatened_native" = "#133e7e", 
                               "threatened_nonnative" = "#133e7e",
                               "nonthreatened_native" = "#bad6f9", 
                               "nonthreatened_nonnative" = "#bad6f9"))

p4<-p3+new_scale_fill()+
  geom_fruit(data = sp_circle,
             geom = geom_tile,
             mapping = aes(y = species, fill = family),
             width = 5,
             offset = 0.04)+
  scale_fill_manual(values = c("Fabaceae" = "black","Myrtaceae" = "black",
                               "Rubiaceae" = "grey","Lauraceae" = "black",
                               "Euphorbiaceae" = "black","Malvaceae" = "grey",
                               "Melastomataceae" = "grey","Annonaceae" = "grey",
                               "Arecaceae" = "black","Sapindaceae" = "black",
                               "Fagaceae" = "black","Moraceae" = "black",
                               "Sapotaceae" = "grey","Rosaceae" = "black",
                               "Salicaceae" = "black","Rutaceae" = "grey",
                               "Phyllanthaceae" = "grey","Apocynaceae" = "black",
                               "Araliaceae" = "grey","Clusiaceae" = "black"),
                    na.value = "white")+
  theme(legend.position = "none")

p4

dev.off()

#######################
# phylogenetic signal
#######################

# overall non-native tree species
tree_phylo<-read.tree(paste0(path_folder,"04_phylogenetic_distance/global_tree.tre"))
sp_circle<-read.table(paste0(path_folder,"04_phylogenetic_distance/phylogeny_outer_circle_dataframe.txt"),
                      header = TRUE)

# transfer to binary variable
sp_circle$bivar<-0
sp_circle$bivar[which(sp_circle$distribution=="nonnative")]<-1

tree_phylo$node.label<-NULL
sp_circle_d<-data.frame(species = sp_circle$species,distribution = as.integer(sp_circle$bivar))

# calculate the phylogenetic D statistic
treephylo_data<-comparative.data(tree_phylo, sp_circle_d, species)
treephylo_res<-phylo.d(treephylo_data, binvar = distribution)

######################################
# threatened non-native tree species
# transfer to binary variable
sp_circle$bivar_threatened<-0
sp_circle$bivar_threatened[which(sp_circle$status=="threatened_nonnative")]<-1

tree_phylo$node.label<-NULL
sp_circle_d_threatened<-data.frame(species = sp_circle$species,status = as.integer(sp_circle$bivar_threatened))

# calculate the phylogenetic D statistic
treephylo_data_threatened<-comparative.data(tree_phylo, sp_circle_d_threatened, species)
treephylo_res_threatened<-phylo.d(treephylo_data_threatened, binvar = status)









