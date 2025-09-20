rm(list=ls())

# loading required libraries
library(plantlist)
library(taxize)
library(ape)
library(V.PhyloMaker)

# path
path_folder <- "D:/non_native_tree/"

# load the all tree list (i.e. the native tree due to the containing of non-native tree)
# native tree
native_list<-list.files(paste0(path_folder,"aggregate_data/tree_data_native_final/"),pattern = ".txt")
native_list<-sub("_native.txt","",native_list,fixed = TRUE)

native_bgci_list<-list.files(paste0(path_folder,"aggregate_data/tree_data_native_final/bgci_part/"))
native_bgci_list<-sub("_native.txt","",native_bgci_list,fixed = TRUE)

native_all<-c(native_list,native_bgci_list) # 51218

# matching the genus and family of native tree
for (id in 1:length(native_all)) {
  
  spacc<-native_all[id]

  # confirm genus and family via TPL function
  sp_tpl<-TPL(spacc)
  
  if(dim(sp_tpl)[1]==1){
    
    sp_gen<-sp_tpl[1,2]
    sp_fam<-sp_tpl[1,3]
    
    # result
    sp_inf<-cbind(spacc,sp_gen,sp_fam)
    
    write.table(sp_inf,file = paste0(path_folder,"04_phylogenetic_distance/tree_confirm_genus_family_native_tpl.txt"),
                row.names = FALSE,col.names = FALSE,append = TRUE)
    
  }else{
    
    print(spacc)
    
  }
  
}

native_confirm<-read.table(paste0(path_folder,"04_phylogenetic_distance/tree_confirm_genus_family_native_tpl.txt"))
colnames(native_confirm)<-c("speciesaccepted","genus","family")
write.table(native_confirm,file = paste0(path_folder,"04_phylogenetic_distance/tree_confirm_genus_family_native_tpl.txt"),
            col.names = TRUE,row.names = FALSE)

native_confirm<-read.table(paste0(path_folder,"04_phylogenetic_distance/tree_confirm_genus_family_native_tpl.txt"),header = TRUE)

native_final<-read.table(paste0(path_folder,"04_phylogenetic_distance/tree_confirm_genus_family_native.txt"),header = TRUE)

# finally forming the phylogenetic tree 
sp_all<-native_final$speciesaccepted
gen_all<-native_final$genus
fam_all<-native_final$family
data_sp<-data.frame(species = sp_all, genus = gen_all, family = fam_all)

tree_base<-phylo.maker(sp.list = data_sp, scenarios="S1",output.tree = TRUE)

tree_base_onlysp<-tree_base$scenario.1

write.tree(tree_base_onlysp,file=paste0(path_folder,"04_phylogenetic_distance/global_tree.tre"))
write.table(tree_base$species.list,file = paste0(path_folder,"04_phylogenetic_distance/global_tree.txt"))

tree_phylo<-read.tree(paste0(path_folder,"04_phylogenetic_distance/global_tree.tre"))

################################################
# calculate the pairwise phylogenetic distance
################################################

# data preparation
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

#######################################################
# match the GADM and TDWG region of tree in BGCI list
for (id in 1:length(native_bgci_list)) {
  
  sp<-native_bgci_list[id]
  
  # native tree file
  sp_file_native<-read.table(paste0(path_folder,"aggregate_data/tree_data_native_final/bgci_part/",sp,"_native.txt"),header = TRUE)
  
  sp_region<-sp_file_native[,c(4,5)]
  sp_region<-sp_region[!duplicated(sp_region),]
  
  write.table(sp_region,file = paste0(path_folder,"04_phylogenetic_distance/gadm_tdwg_match_bgci_tree.txt"),
              col.names = FALSE,row.names = FALSE,append = TRUE)
  
  if(sp%in%nonnative_bgci_list){
    
    # non-native tree file
    sp_file_nonnative<-read.table(paste0(path_folder,"aggregate_data/tree_data_nonnative_final/bgci_part/",sp,"_nonnative.txt"),header = TRUE)
    
    sp_region<-sp_file_nonnative[,c(4,5)]
    sp_region<-sp_region[!duplicated(sp_region),]
    
    write.table(sp_region,file = paste0(path_folder,"04_phylogenetic_distance/gadm_tdwg_match_bgci_tree.txt"),
                col.names = FALSE,row.names = FALSE,append = TRUE)
    
  }
  
  print(id)
  
}

match_file<-read.table(paste0(path_folder,"04_phylogenetic_distance/gadm_tdwg_match_bgci_tree.txt"))
colnames(match_file)<-c("GID_0","NAME_0")
match_file<-match_file[!duplicated(match_file),]

write.table(match_file,file = paste0(path_folder,"04_phylogenetic_distance/gadm_tdwg_match_bgci_tree.txt"),
            col.names = TRUE,row.names = FALSE)

#########################
# phylogenetic distance
native_tdwg<-read.table(paste0(path_folder,"04_phylogenetic_distance/nativetree_tdwglevel3.txt"),header = TRUE)
nonnative_tdwg<-read.table(paste0(path_folder,"04_phylogenetic_distance/nonnativetree_tdwglevel3.txt"),header = TRUE)

# matching the phylogenetic distance among each TDWG Level-3 region
native_region<-unique(native_tdwg$Level3_cod) # 336 regions
nonnative_region<-unique(nonnative_tdwg$Level3_cod) # 338 regions

tree_region<-c(native_region,nonnative_region)
tree_region<-tree_region[duplicated(tree_region)] # 333 regions

# now total tree pool is native_all file
tree_pool<-native_all

for (id in 1:length(tree_region)) {
  
  reg_cod<-tree_region[id]
  
  # native tree
  reg_native<-native_tdwg$spacc[which(native_tdwg$Level3_cod==reg_cod)]
  
  # non-native tree
  reg_nonnative<-nonnative_tdwg$spacc[which(nonnative_tdwg$Level3_cod==reg_cod)]
  # due to the matching of BGCI region with TDWG Level-3 region, some new TDWG Level-3 cod are overlapped or some from TDWG Level-4
  reg_nonnative<-reg_nonnative[!reg_nonnative%in%reg_native]
  
  # failed invasion tree
  reg_failedinvasion<-tree_pool[!tree_pool%in%c(reg_native,reg_nonnative)]
  
  # pariwise phylogenetic distance
  sp_nati<-reg_native
  sp_nati<-gsub(" ","_",sp_nati,fixed = TRUE)
  
  # invasion species (non-native tree)
  divi_nonn<-length(reg_nonnative)%/%1000
  mod_nonn<-length(reg_nonnative)%%1000
  
  if(divi_nonn==0){
    
    # reg_nonnative is from 1 to 999
    sp_nonn<-reg_nonnative
    sp_nonn<-gsub(" ","_",sp_nonn,fixed = TRUE)
    
    sp_group<-c(sp_nonn,sp_nati)
    
    sp_phylo<-keep.tip(tree_phylo,tip = sp_group)
    sp_pd<-cophenetic.phylo(sp_phylo)
    
    # form the matrix that rownames are nonnative trees and colnames are native trees
    sp_pd_nonn<-sp_pd[rownames(sp_pd)%in%sp_nonn,colnames(sp_pd)%in%sp_nati]
    sp_pd_mean<-apply(sp_pd_nonn, 1, mean)
    sp_pd_min<-apply(sp_pd_nonn, 1, min)
    sp_pd_pro<-rep(1,length(sp_nonn))
    
    sp_pd_rowname<-gsub("_"," ",rownames(sp_pd_nonn),fixed = TRUE)
    sp_pd_reg<-rep(reg_cod,length(sp_nonn))
    sp_pd_final<-cbind(sp_pd_rowname,sp_pd_pro,sp_pd_mean,sp_pd_min,sp_pd_reg)
    
    write.table(sp_pd_final,file = paste0(path_folder,"04_phylogenetic_distance/pairwise_phylogenetic_distance/phylogenetic_distance_",reg_cod,".txt"),
                col.names = FALSE,row.names = FALSE,append = TRUE)
    
  }else{
    
    if(mod_nonn!=0){
      
      # reg_nonnative is larger than 1000 (including 1000) and its mod is unequal to 0
      for (subid in 1:divi_nonn) {
        
        sp_nonn<-reg_nonnative[((subid-1)*1000+1):(subid*1000)]
        sp_nonn<-gsub(" ","_",sp_nonn,fixed = TRUE)
        
        sp_group<-c(sp_nonn,sp_nati)
        
        sp_phylo<-keep.tip(tree_phylo,tip = sp_group)
        sp_pd<-cophenetic.phylo(sp_phylo)
        
        # form the matrix that rownames are nonnative trees and colnames are native trees
        sp_pd_nonn<-sp_pd[rownames(sp_pd)%in%sp_nonn,colnames(sp_pd)%in%sp_nati]
        sp_pd_mean<-apply(sp_pd_nonn, 1, mean)
        sp_pd_min<-apply(sp_pd_nonn, 1, min)
        sp_pd_pro<-rep(1,length(sp_nonn))
        
        sp_pd_rowname<-gsub("_"," ",rownames(sp_pd_nonn),fixed = TRUE)
        sp_pd_reg<-rep(reg_cod,length(sp_nonn))
        sp_pd_final<-cbind(sp_pd_rowname,sp_pd_pro,sp_pd_mean,sp_pd_min,sp_pd_reg)
        
        write.table(sp_pd_final,file = paste0(path_folder,"04_phylogenetic_distance/pairwise_phylogenetic_distance/phylogenetic_distance_",reg_cod,".txt"),
                    col.names = FALSE,row.names = FALSE,append = TRUE)
        
        print(subid)
        
      }
      
      # extract the mod part
      sp_nonn<-reg_nonnative[(divi_nonn*1000+1):(divi_nonn*1000+mod_nonn)]
      sp_nonn<-gsub(" ","_",sp_nonn,fixed = TRUE)
      
      sp_group<-c(sp_nonn,sp_nati)
      
      sp_phylo<-keep.tip(tree_phylo,tip = sp_group)
      sp_pd<-cophenetic.phylo(sp_phylo)
      
      # form the matrix that rownames are nonnative trees and colnames are native trees
      sp_pd_nonn<-sp_pd[rownames(sp_pd)%in%sp_nonn,colnames(sp_pd)%in%sp_nati]
      sp_pd_mean<-apply(sp_pd_nonn, 1, mean)
      sp_pd_min<-apply(sp_pd_nonn, 1, min)
      sp_pd_pro<-rep(1,length(sp_nonn))
      
      sp_pd_rowname<-gsub("_"," ",rownames(sp_pd_nonn),fixed = TRUE)
      sp_pd_reg<-rep(reg_cod,length(sp_nonn))
      sp_pd_final<-cbind(sp_pd_rowname,sp_pd_pro,sp_pd_mean,sp_pd_min,sp_pd_reg)
      
      write.table(sp_pd_final,file = paste0(path_folder,"04_phylogenetic_distance/pairwise_phylogenetic_distance/phylogenetic_distance_",reg_cod,".txt"),
                  col.names = FALSE,row.names = FALSE,append = TRUE)
      
    }else{
      
      # reg_nonnative is larger than 1000 (including 1000) and its mod is equal to 0
      for (subid in 1:divi_nonn) {
        
        sp_nonn<-reg_nonnative[((subid-1)*1000+1):(subid*1000)]
        sp_nonn<-gsub(" ","_",sp_nonn,fixed = TRUE)
        
        sp_group<-c(sp_nonn,sp_nati)
        
        sp_phylo<-keep.tip(tree_phylo,tip = sp_group)
        sp_pd<-cophenetic.phylo(sp_phylo)
        
        # form the matrix that rownames are nonnative trees and colnames are native trees
        sp_pd_nonn<-sp_pd[rownames(sp_pd)%in%sp_nonn,colnames(sp_pd)%in%sp_nati]
        sp_pd_mean<-apply(sp_pd_nonn, 1, mean)
        sp_pd_min<-apply(sp_pd_nonn, 1, min)
        sp_pd_pro<-rep(1,length(sp_nonn))
        
        sp_pd_rowname<-gsub("_"," ",rownames(sp_pd_nonn),fixed = TRUE)
        sp_pd_reg<-rep(reg_cod,length(sp_nonn))
        sp_pd_final<-cbind(sp_pd_rowname,sp_pd_pro,sp_pd_mean,sp_pd_min,sp_pd_reg)
        
        write.table(sp_pd_final,file = paste0(path_folder,"04_phylogenetic_distance/pairwise_phylogenetic_distance/phylogenetic_distance_",reg_cod,".txt"),
                    col.names = FALSE,row.names = FALSE,append = TRUE)
        
        print(subid)
        
      }
      
    }
    
  }
  
  # failed invasion
  divi_fail<-length(reg_failedinvasion)%/%1000
  mod_fail<-length(reg_failedinvasion)%%1000
  
  if(divi_fail==0){
    
    # reg_failedinvasion is from 1 to 999
    sp_fail<-reg_failedinvasion
    sp_fail<-gsub(" ","_",sp_fail,fixed = TRUE)
    
    sp_group<-c(sp_fail,sp_nati)
    
    sp_phylo<-keep.tip(tree_phylo,tip = sp_group)
    sp_pd<-cophenetic.phylo(sp_phylo)
    
    # form the matrix that rownames are failed invasion trees and colnames are native trees
    sp_pd_fail<-sp_pd[rownames(sp_pd)%in%sp_fail,colnames(sp_pd)%in%sp_nati]
    sp_pd_mean<-apply(sp_pd_fail, 1, mean)
    sp_pd_min<-apply(sp_pd_fail, 1, min)
    sp_pd_pro<-rep(0,length(sp_fail))
    
    sp_pd_rowname<-gsub("_"," ",rownames(sp_pd_fail),fixed = TRUE)
    sp_pd_reg<-rep(reg_cod,length(sp_fail))
    sp_pd_final<-cbind(sp_pd_rowname,sp_pd_pro,sp_pd_mean,sp_pd_min,sp_pd_reg)
    
    write.table(sp_pd_final,file = paste0(path_folder,"04_phylogenetic_distance/pairwise_phylogenetic_distance/phylogenetic_distance_",reg_cod,".txt"),
                col.names = FALSE,row.names = FALSE,append = TRUE)
    
  }else{
    
    if(mod_fail!=0){
      
      # reg_failedinvasion is larger than 1000 (including 1000) and its mod is unequal to 0
      for (subid in 1:divi_fail) {
        
        sp_fail<-reg_failedinvasion[((subid-1)*1000+1):(subid*1000)]
        sp_fail<-gsub(" ","_",sp_fail,fixed = TRUE)
        
        sp_group<-c(sp_fail,sp_nati)
        
        sp_phylo<-keep.tip(tree_phylo,tip = sp_group)
        sp_pd<-cophenetic.phylo(sp_phylo)
        
        # form the matrix that rownames are failed invasion trees and colnames are native trees
        sp_pd_fail<-sp_pd[rownames(sp_pd)%in%sp_fail,colnames(sp_pd)%in%sp_nati]
        sp_pd_mean<-apply(sp_pd_fail, 1, mean)
        sp_pd_min<-apply(sp_pd_fail, 1, min)
        sp_pd_pro<-rep(0,length(sp_fail))
        
        sp_pd_rowname<-gsub("_"," ",rownames(sp_pd_fail),fixed = TRUE)
        sp_pd_reg<-rep(reg_cod,length(sp_fail))
        sp_pd_final<-cbind(sp_pd_rowname,sp_pd_pro,sp_pd_mean,sp_pd_min,sp_pd_reg)
        
        write.table(sp_pd_final,file = paste0(path_folder,"04_phylogenetic_distance/pairwise_phylogenetic_distance/phylogenetic_distance_",reg_cod,".txt"),
                    col.names = FALSE,row.names = FALSE,append = TRUE)
        
        print(subid)
        
      }
      
      # extract the mod part
      sp_fail<-reg_failedinvasion[(divi_fail*1000+1):(divi_fail*1000+mod_fail)]
      sp_fail<-gsub(" ","_",sp_fail,fixed = TRUE)
      
      sp_group<-c(sp_fail,sp_nati)
      
      sp_phylo<-keep.tip(tree_phylo,tip = sp_group)
      sp_pd<-cophenetic.phylo(sp_phylo)
      
      # form the matrix that rownames are failed invasion trees and colnames are native trees
      sp_pd_fail<-sp_pd[rownames(sp_pd)%in%sp_fail,colnames(sp_pd)%in%sp_nati]
      sp_pd_mean<-apply(sp_pd_fail, 1, mean)
      sp_pd_min<-apply(sp_pd_fail, 1, min)
      sp_pd_pro<-rep(0,length(sp_fail))
      
      sp_pd_rowname<-gsub("_"," ",rownames(sp_pd_fail),fixed = TRUE)
      sp_pd_reg<-rep(reg_cod,length(sp_fail))
      sp_pd_final<-cbind(sp_pd_rowname,sp_pd_pro,sp_pd_mean,sp_pd_min,sp_pd_reg)
      
      write.table(sp_pd_final,file = paste0(path_folder,"04_phylogenetic_distance/pairwise_phylogenetic_distance/phylogenetic_distance_",reg_cod,".txt"),
                  col.names = FALSE,row.names = FALSE,append = TRUE)
      
    }else{
      
      # reg_failedinvasion is larger than 1000 (including 1000) and its mod is equal to 0
      for (subid in 1:divi_fail) {
        
        sp_fail<-reg_failedinvasion[((subid-1)*1000+1):(subid*1000)]
        sp_fail<-gsub(" ","_",sp_fail,fixed = TRUE)
        
        sp_group<-c(sp_fail,sp_nati)
        
        sp_phylo<-keep.tip(tree_phylo,tip = sp_group)
        sp_pd<-cophenetic.phylo(sp_phylo)
        
        # form the matrix that rownames are failed invasion trees and colnames are native trees
        sp_pd_fail<-sp_pd[rownames(sp_pd)%in%sp_fail,colnames(sp_pd)%in%sp_nati]
        sp_pd_mean<-apply(sp_pd_fail, 1, mean)
        sp_pd_min<-apply(sp_pd_fail, 1, min)
        sp_pd_pro<-rep(0,length(sp_fail))
        
        sp_pd_rowname<-gsub("_"," ",rownames(sp_pd_fail),fixed = TRUE)
        sp_pd_reg<-rep(reg_cod,length(sp_fail))
        sp_pd_final<-cbind(sp_pd_rowname,sp_pd_pro,sp_pd_mean,sp_pd_min,sp_pd_reg)
        
        write.table(sp_pd_final,file = paste0(path_folder,"04_phylogenetic_distance/pairwise_phylogenetic_distance/phylogenetic_distance_",reg_cod,".txt"),
                    col.names = FALSE,row.names = FALSE,append = TRUE)
        
        print(subid)
        
      }
      
    }
    
  }
  
  print(id)
  
}

# "FOR", "MRS", "NUE", "CPI", "MDV", "MPE" need manually deal with

# adding the column names
for (id in 1:length(tree_region)) {
  
  reg_name<-tree_region[id]
  
  reg_file<-read.table(paste0(path_folder,"04_phylogenetic_distance/pairwise_phylogenetic_distance/phylogenetic_distance_",reg_name,".txt"))
  colnames(reg_file)<-c("species","probability","pdmean","pdmin","region")
  
  write.table(reg_file,file = paste0(path_folder,"04_phylogenetic_distance/pairwise_phylogenetic_distance/phylogenetic_distance_",reg_name,".txt"),
              col.names = TRUE,row.names = FALSE)
  
  print(id)
  
}

###################
# threatened tree
###################

# load the threatened list file (Critically Endangered, Endangered and Vulnerable)
threatened_list<-read.table(paste0(path_folder,"threatened_data/threatened_species_list.txt"),header = TRUE)
# 23950 threatened trees
threatened_all<-unique(threatened_list$species_accept[which(threatened_list$iucn_category%in%c("Critically Endangered","Endangered","Vulnerable"))])

native_threatened<-threatened_all[threatened_all%in%native_all] # 11385
nonnative_threatened<-threatened_all[threatened_all%in%nonnative_all] # 3626

###########################################
# phylogeny (the same as all trees above)
tree_phylo<-read.tree(paste0(path_folder,"04_phylogenetic_distance/global_tree.tre"))

# load the data files (the same as all trees above)
native_tdwg<-read.table(paste0(path_folder,"04_phylogenetic_distance/nativetree_tdwglevel3.txt"),header = TRUE)
nonnative_tdwg<-read.table(paste0(path_folder,"04_phylogenetic_distance/nonnativetree_tdwglevel3.txt"),header = TRUE)

# extract the region having threatened non-native tree (keeping the same native tree list)
nonnative_tdwg<-nonnative_tdwg[which(nonnative_tdwg$spacc%in%nonnative_threatened),]

# matching the phylogenetic distance among each TDWG Level-3 region
native_region<-unique(native_tdwg$Level3_cod) # 336 regions
nonnative_region<-unique(nonnative_tdwg$Level3_cod) # 301 regions
# select the region having both threatened non-native tree and native species (i.e. threatened or not)
tree_region<-nonnative_region[nonnative_region%in%native_region] # 299 regions

# now potential introduction (invasion) tree pool is native_threatened file
tree_pool<-native_threatened

for (id in 1:length(tree_region)) {
  
  reg_cod<-tree_region[id]
  
  # native tree
  reg_native<-native_tdwg$spacc[which(native_tdwg$Level3_cod==reg_cod)]
  
  # non-native tree
  reg_nonnative<-nonnative_tdwg$spacc[which(nonnative_tdwg$Level3_cod==reg_cod)]
  # due to the matching of BGCI region with TDWG Level-3 region, some new TDWG Level-3 cod are overlapped or some from TDWG Level-4
  reg_nonnative<-reg_nonnative[!reg_nonnative%in%reg_native]
  
  # failed invasion tree
  reg_failedinvasion<-tree_pool[!tree_pool%in%c(reg_native,reg_nonnative)]
  
  # pariwise phylogenetic distance
  sp_nati<-reg_native
  sp_nati<-gsub(" ","_",sp_nati,fixed = TRUE)
  
  # invasion species (non-native tree)
  divi_nonn<-length(reg_nonnative)%/%1000
  mod_nonn<-length(reg_nonnative)%%1000
  
  if(divi_nonn==0){
    
    # reg_nonnative is from 1 to 999
    sp_nonn<-reg_nonnative
    sp_nonn<-gsub(" ","_",sp_nonn,fixed = TRUE)
    
    sp_group<-c(sp_nonn,sp_nati)
    
    sp_phylo<-keep.tip(tree_phylo,tip = sp_group)
    sp_pd<-cophenetic.phylo(sp_phylo)
    
    # form the matrix that rownames are nonnative trees and colnames are native trees
    sp_pd_nonn<-sp_pd[rownames(sp_pd)%in%sp_nonn,colnames(sp_pd)%in%sp_nati]
    sp_pd_mean<-apply(sp_pd_nonn, 1, mean)
    sp_pd_min<-apply(sp_pd_nonn, 1, min)
    sp_pd_pro<-rep(1,length(sp_nonn))
    
    sp_pd_rowname<-gsub("_"," ",rownames(sp_pd_nonn),fixed = TRUE)
    sp_pd_reg<-rep(reg_cod,length(sp_nonn))
    sp_pd_final<-cbind(sp_pd_rowname,sp_pd_pro,sp_pd_mean,sp_pd_min,sp_pd_reg)
    
    write.table(sp_pd_final,file = paste0(path_folder,"04_phylogenetic_distance/pairwise_phylogenetic_distance_threatened/phylogenetic_distance_",reg_cod,"_threatened.txt"),
                col.names = FALSE,row.names = FALSE,append = TRUE)
    
  }else{
    
    if(mod_nonn!=0){
      
      # reg_nonnative is larger than 1000 (including 1000) and its mod is unequal to 0
      for (subid in 1:divi_nonn) {
        
        sp_nonn<-reg_nonnative[((subid-1)*1000+1):(subid*1000)]
        sp_nonn<-gsub(" ","_",sp_nonn,fixed = TRUE)
        
        sp_group<-c(sp_nonn,sp_nati)
        
        sp_phylo<-keep.tip(tree_phylo,tip = sp_group)
        sp_pd<-cophenetic.phylo(sp_phylo)
        
        # form the matrix that rownames are nonnative trees and colnames are native trees
        sp_pd_nonn<-sp_pd[rownames(sp_pd)%in%sp_nonn,colnames(sp_pd)%in%sp_nati]
        sp_pd_mean<-apply(sp_pd_nonn, 1, mean)
        sp_pd_min<-apply(sp_pd_nonn, 1, min)
        sp_pd_pro<-rep(1,length(sp_nonn))
        
        sp_pd_rowname<-gsub("_"," ",rownames(sp_pd_nonn),fixed = TRUE)
        sp_pd_reg<-rep(reg_cod,length(sp_nonn))
        sp_pd_final<-cbind(sp_pd_rowname,sp_pd_pro,sp_pd_mean,sp_pd_min,sp_pd_reg)
        
        write.table(sp_pd_final,file = paste0(path_folder,"04_phylogenetic_distance/pairwise_phylogenetic_distance_threatened/phylogenetic_distance_",reg_cod,"_threatened.txt"),
                    col.names = FALSE,row.names = FALSE,append = TRUE)
        
        print(subid)
        
      }
      
      # extract the mod part
      sp_nonn<-reg_nonnative[(divi_nonn*1000+1):(divi_nonn*1000+mod_nonn)]
      sp_nonn<-gsub(" ","_",sp_nonn,fixed = TRUE)
      
      sp_group<-c(sp_nonn,sp_nati)
      
      sp_phylo<-keep.tip(tree_phylo,tip = sp_group)
      sp_pd<-cophenetic.phylo(sp_phylo)
      
      # form the matrix that rownames are nonnative trees and colnames are native trees
      sp_pd_nonn<-sp_pd[rownames(sp_pd)%in%sp_nonn,colnames(sp_pd)%in%sp_nati]
      sp_pd_mean<-apply(sp_pd_nonn, 1, mean)
      sp_pd_min<-apply(sp_pd_nonn, 1, min)
      sp_pd_pro<-rep(1,length(sp_nonn))
      
      sp_pd_rowname<-gsub("_"," ",rownames(sp_pd_nonn),fixed = TRUE)
      sp_pd_reg<-rep(reg_cod,length(sp_nonn))
      sp_pd_final<-cbind(sp_pd_rowname,sp_pd_pro,sp_pd_mean,sp_pd_min,sp_pd_reg)
      
      write.table(sp_pd_final,file = paste0(path_folder,"04_phylogenetic_distance/pairwise_phylogenetic_distance_threatened/phylogenetic_distance_",reg_cod,"_threatened.txt"),
                  col.names = FALSE,row.names = FALSE,append = TRUE)
      
    }else{
      
      # reg_nonnative is larger than 1000 (including 1000) and its mod is equal to 0
      for (subid in 1:divi_nonn) {
        
        sp_nonn<-reg_nonnative[((subid-1)*1000+1):(subid*1000)]
        sp_nonn<-gsub(" ","_",sp_nonn,fixed = TRUE)
        
        sp_group<-c(sp_nonn,sp_nati)
        
        sp_phylo<-keep.tip(tree_phylo,tip = sp_group)
        sp_pd<-cophenetic.phylo(sp_phylo)
        
        # form the matrix that rownames are nonnative trees and colnames are native trees
        sp_pd_nonn<-sp_pd[rownames(sp_pd)%in%sp_nonn,colnames(sp_pd)%in%sp_nati]
        sp_pd_mean<-apply(sp_pd_nonn, 1, mean)
        sp_pd_min<-apply(sp_pd_nonn, 1, min)
        sp_pd_pro<-rep(1,length(sp_nonn))
        
        sp_pd_rowname<-gsub("_"," ",rownames(sp_pd_nonn),fixed = TRUE)
        sp_pd_reg<-rep(reg_cod,length(sp_nonn))
        sp_pd_final<-cbind(sp_pd_rowname,sp_pd_pro,sp_pd_mean,sp_pd_min,sp_pd_reg)
        
        write.table(sp_pd_final,file = paste0(path_folder,"04_phylogenetic_distance/pairwise_phylogenetic_distance_threatened/phylogenetic_distance_",reg_cod,"_threatened.txt"),
                    col.names = FALSE,row.names = FALSE,append = TRUE)
        
        print(subid)
        
      }
      
    }
    
  }
  
  # failed invasion
  divi_fail<-length(reg_failedinvasion)%/%1000
  mod_fail<-length(reg_failedinvasion)%%1000
  
  if(divi_fail==0){
    
    # reg_failedinvasion is from 1 to 999
    sp_fail<-reg_failedinvasion
    sp_fail<-gsub(" ","_",sp_fail,fixed = TRUE)
    
    sp_group<-c(sp_fail,sp_nati)
    
    sp_phylo<-keep.tip(tree_phylo,tip = sp_group)
    sp_pd<-cophenetic.phylo(sp_phylo)
    
    # form the matrix that rownames are failed invasion trees and colnames are native trees
    sp_pd_fail<-sp_pd[rownames(sp_pd)%in%sp_fail,colnames(sp_pd)%in%sp_nati]
    sp_pd_mean<-apply(sp_pd_fail, 1, mean)
    sp_pd_min<-apply(sp_pd_fail, 1, min)
    sp_pd_pro<-rep(0,length(sp_fail))
    
    sp_pd_rowname<-gsub("_"," ",rownames(sp_pd_fail),fixed = TRUE)
    sp_pd_reg<-rep(reg_cod,length(sp_fail))
    sp_pd_final<-cbind(sp_pd_rowname,sp_pd_pro,sp_pd_mean,sp_pd_min,sp_pd_reg)
    
    write.table(sp_pd_final,file = paste0(path_folder,"04_phylogenetic_distance/pairwise_phylogenetic_distance_threatened/phylogenetic_distance_",reg_cod,"_threatened.txt"),
                col.names = FALSE,row.names = FALSE,append = TRUE)
    
  }else{
    
    if(mod_fail!=0){
      
      # reg_failedinvasion is larger than 1000 (including 1000) and its mod is unequal to 0
      for (subid in 1:divi_fail) {
        
        sp_fail<-reg_failedinvasion[((subid-1)*1000+1):(subid*1000)]
        sp_fail<-gsub(" ","_",sp_fail,fixed = TRUE)
        
        sp_group<-c(sp_fail,sp_nati)
        
        sp_phylo<-keep.tip(tree_phylo,tip = sp_group)
        sp_pd<-cophenetic.phylo(sp_phylo)
        
        # form the matrix that rownames are failed invasion trees and colnames are native trees
        sp_pd_fail<-sp_pd[rownames(sp_pd)%in%sp_fail,colnames(sp_pd)%in%sp_nati]
        sp_pd_mean<-apply(sp_pd_fail, 1, mean)
        sp_pd_min<-apply(sp_pd_fail, 1, min)
        sp_pd_pro<-rep(0,length(sp_fail))
        
        sp_pd_rowname<-gsub("_"," ",rownames(sp_pd_fail),fixed = TRUE)
        sp_pd_reg<-rep(reg_cod,length(sp_fail))
        sp_pd_final<-cbind(sp_pd_rowname,sp_pd_pro,sp_pd_mean,sp_pd_min,sp_pd_reg)
        
        write.table(sp_pd_final,file = paste0(path_folder,"04_phylogenetic_distance/pairwise_phylogenetic_distance_threatened/phylogenetic_distance_",reg_cod,"_threatened.txt"),
                    col.names = FALSE,row.names = FALSE,append = TRUE)
        
        print(subid)
        
      }
      
      # extract the mod part
      sp_fail<-reg_failedinvasion[(divi_fail*1000+1):(divi_fail*1000+mod_fail)]
      sp_fail<-gsub(" ","_",sp_fail,fixed = TRUE)
      
      sp_group<-c(sp_fail,sp_nati)
      
      sp_phylo<-keep.tip(tree_phylo,tip = sp_group)
      sp_pd<-cophenetic.phylo(sp_phylo)
      
      # form the matrix that rownames are failed invasion trees and colnames are native trees
      sp_pd_fail<-sp_pd[rownames(sp_pd)%in%sp_fail,colnames(sp_pd)%in%sp_nati]
      sp_pd_mean<-apply(sp_pd_fail, 1, mean)
      sp_pd_min<-apply(sp_pd_fail, 1, min)
      sp_pd_pro<-rep(0,length(sp_fail))
      
      sp_pd_rowname<-gsub("_"," ",rownames(sp_pd_fail),fixed = TRUE)
      sp_pd_reg<-rep(reg_cod,length(sp_fail))
      sp_pd_final<-cbind(sp_pd_rowname,sp_pd_pro,sp_pd_mean,sp_pd_min,sp_pd_reg)
      
      write.table(sp_pd_final,file = paste0(path_folder,"04_phylogenetic_distance/pairwise_phylogenetic_distance_threatened/phylogenetic_distance_",reg_cod,"_threatened.txt"),
                  col.names = FALSE,row.names = FALSE,append = TRUE)
      
    }else{
      
      # reg_failedinvasion is larger than 1000 (including 1000) and its mod is equal to 0
      for (subid in 1:divi_fail) {
        
        sp_fail<-reg_failedinvasion[((subid-1)*1000+1):(subid*1000)]
        sp_fail<-gsub(" ","_",sp_fail,fixed = TRUE)
        
        sp_group<-c(sp_fail,sp_nati)
        
        sp_phylo<-keep.tip(tree_phylo,tip = sp_group)
        sp_pd<-cophenetic.phylo(sp_phylo)
        
        # form the matrix that rownames are failed invasion trees and colnames are native trees
        sp_pd_fail<-sp_pd[rownames(sp_pd)%in%sp_fail,colnames(sp_pd)%in%sp_nati]
        sp_pd_mean<-apply(sp_pd_fail, 1, mean)
        sp_pd_min<-apply(sp_pd_fail, 1, min)
        sp_pd_pro<-rep(0,length(sp_fail))
        
        sp_pd_rowname<-gsub("_"," ",rownames(sp_pd_fail),fixed = TRUE)
        sp_pd_reg<-rep(reg_cod,length(sp_fail))
        sp_pd_final<-cbind(sp_pd_rowname,sp_pd_pro,sp_pd_mean,sp_pd_min,sp_pd_reg)
        
        write.table(sp_pd_final,file = paste0(path_folder,"04_phylogenetic_distance/pairwise_phylogenetic_distance_threatened/phylogenetic_distance_",reg_cod,"_threatened.txt"),
                    col.names = FALSE,row.names = FALSE,append = TRUE)
        
        print(subid)
        
      }
      
    }
    
  }
  
  print(id)
  
}

# "CRL", "CVI", "GNB", "AND", "YAK", "KUR", "YUK", "NDA", "TUA", "CAY", "SAM" need manually deal with

# adding the column names
for (id in 1:length(tree_region)) {
  
  reg_name<-tree_region[id]
  
  reg_file<-read.table(paste0(path_folder,"04_phylogenetic_distance/pairwise_phylogenetic_distance_threatened/phylogenetic_distance_",reg_name,"_threatened.txt"))
  colnames(reg_file)<-c("species","probability","pdmean","pdmin","region")
  
  write.table(reg_file,file = paste0(path_folder,"04_phylogenetic_distance/pairwise_phylogenetic_distance_threatened/phylogenetic_distance_",reg_name,"_threatened.txt"),
              col.names = TRUE,row.names = FALSE)
  
  print(id)
  
}

#######################
# non-threatened tree
#######################

# extract the non-threatened trees
native_nonthreatened<-native_all[!native_all%in%native_threatened] # 39833
nonnative_nonthreatened<-nonnative_all[!nonnative_all%in%nonnative_threatened] # 22870

###########################################
# phylogeny (the same as all trees above)
tree_phylo<-read.tree(paste0(path_folder,"04_phylogenetic_distance/global_tree.tre"))

# load the data files (the same as all trees above)
native_tdwg<-read.table(paste0(path_folder,"04_phylogenetic_distance/nativetree_tdwglevel3.txt"),header = TRUE)
nonnative_tdwg<-read.table(paste0(path_folder,"04_phylogenetic_distance/nonnativetree_tdwglevel3.txt"),header = TRUE)

# extract the region having non-threatened non-native tree (keeping the same native tree list)
nonnative_tdwg<-nonnative_tdwg[which(nonnative_tdwg$spacc%in%nonnative_nonthreatened),]

# matching the phylogenetic distance among each TDWG Level-3 region
native_region<-unique(native_tdwg$Level3_cod) # 336 regions
nonnative_region<-unique(nonnative_tdwg$Level3_cod) # 337 regions
# select the region having both non-threatened non-native tree and native species (i.e. threatened or not)
tree_region<-nonnative_region[nonnative_region%in%native_region] # 333 regions

# now potential introduction (invasion) tree pool is native_threatened file
tree_pool<-native_nonthreatened

for (id in 1:length(tree_region)) {
  
  reg_cod<-tree_region[id]
  
  # native tree
  reg_native<-native_tdwg$spacc[which(native_tdwg$Level3_cod==reg_cod)]
  
  # non-native tree
  reg_nonnative<-nonnative_tdwg$spacc[which(nonnative_tdwg$Level3_cod==reg_cod)]
  # due to the matching of BGCI region with TDWG Level-3 region, some new TDWG Level-3 cod are overlapped or some from TDWG Level-4
  reg_nonnative<-reg_nonnative[!reg_nonnative%in%reg_native]
  
  # failed invasion tree
  reg_failedinvasion<-tree_pool[!tree_pool%in%c(reg_native,reg_nonnative)]
  
  # pariwise phylogenetic distance
  sp_nati<-reg_native
  sp_nati<-gsub(" ","_",sp_nati,fixed = TRUE)
  
  # invasion species (non-native tree)
  divi_nonn<-length(reg_nonnative)%/%1000
  mod_nonn<-length(reg_nonnative)%%1000
  
  if(divi_nonn==0){
    
    # reg_nonnative is from 1 to 999
    sp_nonn<-reg_nonnative
    sp_nonn<-gsub(" ","_",sp_nonn,fixed = TRUE)
    
    sp_group<-c(sp_nonn,sp_nati)
    
    sp_phylo<-keep.tip(tree_phylo,tip = sp_group)
    sp_pd<-cophenetic.phylo(sp_phylo)
    
    # form the matrix that rownames are nonnative trees and colnames are native trees
    sp_pd_nonn<-sp_pd[rownames(sp_pd)%in%sp_nonn,colnames(sp_pd)%in%sp_nati]
    sp_pd_mean<-apply(sp_pd_nonn, 1, mean)
    sp_pd_min<-apply(sp_pd_nonn, 1, min)
    sp_pd_pro<-rep(1,length(sp_nonn))
    
    sp_pd_rowname<-gsub("_"," ",rownames(sp_pd_nonn),fixed = TRUE)
    sp_pd_reg<-rep(reg_cod,length(sp_nonn))
    sp_pd_final<-cbind(sp_pd_rowname,sp_pd_pro,sp_pd_mean,sp_pd_min,sp_pd_reg)
    
    write.table(sp_pd_final,file = paste0(path_folder,"04_phylogenetic_distance/pairwise_phylogenetic_distance_nonthreatened/phylogenetic_distance_",reg_cod,"_nonthreatened.txt"),
                col.names = FALSE,row.names = FALSE,append = TRUE)
    
  }else{
    
    if(mod_nonn!=0){
      
      # reg_nonnative is larger than 1000 (including 1000) and its mod is unequal to 0
      for (subid in 1:divi_nonn) {
        
        sp_nonn<-reg_nonnative[((subid-1)*1000+1):(subid*1000)]
        sp_nonn<-gsub(" ","_",sp_nonn,fixed = TRUE)
        
        sp_group<-c(sp_nonn,sp_nati)
        
        sp_phylo<-keep.tip(tree_phylo,tip = sp_group)
        sp_pd<-cophenetic.phylo(sp_phylo)
        
        # form the matrix that rownames are nonnative trees and colnames are native trees
        sp_pd_nonn<-sp_pd[rownames(sp_pd)%in%sp_nonn,colnames(sp_pd)%in%sp_nati]
        sp_pd_mean<-apply(sp_pd_nonn, 1, mean)
        sp_pd_min<-apply(sp_pd_nonn, 1, min)
        sp_pd_pro<-rep(1,length(sp_nonn))
        
        sp_pd_rowname<-gsub("_"," ",rownames(sp_pd_nonn),fixed = TRUE)
        sp_pd_reg<-rep(reg_cod,length(sp_nonn))
        sp_pd_final<-cbind(sp_pd_rowname,sp_pd_pro,sp_pd_mean,sp_pd_min,sp_pd_reg)
        
        write.table(sp_pd_final,file = paste0(path_folder,"04_phylogenetic_distance/pairwise_phylogenetic_distance_nonthreatened/phylogenetic_distance_",reg_cod,"_nonthreatened.txt"),
                    col.names = FALSE,row.names = FALSE,append = TRUE)
        
        print(subid)
        
      }
      
      # extract the mod part
      sp_nonn<-reg_nonnative[(divi_nonn*1000+1):(divi_nonn*1000+mod_nonn)]
      sp_nonn<-gsub(" ","_",sp_nonn,fixed = TRUE)
      
      sp_group<-c(sp_nonn,sp_nati)
      
      sp_phylo<-keep.tip(tree_phylo,tip = sp_group)
      sp_pd<-cophenetic.phylo(sp_phylo)
      
      # form the matrix that rownames are nonnative trees and colnames are native trees
      sp_pd_nonn<-sp_pd[rownames(sp_pd)%in%sp_nonn,colnames(sp_pd)%in%sp_nati]
      sp_pd_mean<-apply(sp_pd_nonn, 1, mean)
      sp_pd_min<-apply(sp_pd_nonn, 1, min)
      sp_pd_pro<-rep(1,length(sp_nonn))
      
      sp_pd_rowname<-gsub("_"," ",rownames(sp_pd_nonn),fixed = TRUE)
      sp_pd_reg<-rep(reg_cod,length(sp_nonn))
      sp_pd_final<-cbind(sp_pd_rowname,sp_pd_pro,sp_pd_mean,sp_pd_min,sp_pd_reg)
      
      write.table(sp_pd_final,file = paste0(path_folder,"04_phylogenetic_distance/pairwise_phylogenetic_distance_nonthreatened/phylogenetic_distance_",reg_cod,"_nonthreatened.txt"),
                  col.names = FALSE,row.names = FALSE,append = TRUE)
      
    }else{
      
      # reg_nonnative is larger than 1000 (including 1000) and its mod is equal to 0
      for (subid in 1:divi_nonn) {
        
        sp_nonn<-reg_nonnative[((subid-1)*1000+1):(subid*1000)]
        sp_nonn<-gsub(" ","_",sp_nonn,fixed = TRUE)
        
        sp_group<-c(sp_nonn,sp_nati)
        
        sp_phylo<-keep.tip(tree_phylo,tip = sp_group)
        sp_pd<-cophenetic.phylo(sp_phylo)
        
        # form the matrix that rownames are nonnative trees and colnames are native trees
        sp_pd_nonn<-sp_pd[rownames(sp_pd)%in%sp_nonn,colnames(sp_pd)%in%sp_nati]
        sp_pd_mean<-apply(sp_pd_nonn, 1, mean)
        sp_pd_min<-apply(sp_pd_nonn, 1, min)
        sp_pd_pro<-rep(1,length(sp_nonn))
        
        sp_pd_rowname<-gsub("_"," ",rownames(sp_pd_nonn),fixed = TRUE)
        sp_pd_reg<-rep(reg_cod,length(sp_nonn))
        sp_pd_final<-cbind(sp_pd_rowname,sp_pd_pro,sp_pd_mean,sp_pd_min,sp_pd_reg)
        
        write.table(sp_pd_final,file = paste0(path_folder,"04_phylogenetic_distance/pairwise_phylogenetic_distance_nonthreatened/phylogenetic_distance_",reg_cod,"_nonthreatened.txt"),
                    col.names = FALSE,row.names = FALSE,append = TRUE)
        
        print(subid)
        
      }
      
    }
    
  }
  
  # failed invasion
  divi_fail<-length(reg_failedinvasion)%/%1000
  mod_fail<-length(reg_failedinvasion)%%1000
  
  if(divi_fail==0){
    
    # reg_failedinvasion is from 1 to 999
    sp_fail<-reg_failedinvasion
    sp_fail<-gsub(" ","_",sp_fail,fixed = TRUE)
    
    sp_group<-c(sp_fail,sp_nati)
    
    sp_phylo<-keep.tip(tree_phylo,tip = sp_group)
    sp_pd<-cophenetic.phylo(sp_phylo)
    
    # form the matrix that rownames are failed invasion trees and colnames are native trees
    sp_pd_fail<-sp_pd[rownames(sp_pd)%in%sp_fail,colnames(sp_pd)%in%sp_nati]
    sp_pd_mean<-apply(sp_pd_fail, 1, mean)
    sp_pd_min<-apply(sp_pd_fail, 1, min)
    sp_pd_pro<-rep(0,length(sp_fail))
    
    sp_pd_rowname<-gsub("_"," ",rownames(sp_pd_fail),fixed = TRUE)
    sp_pd_reg<-rep(reg_cod,length(sp_fail))
    sp_pd_final<-cbind(sp_pd_rowname,sp_pd_pro,sp_pd_mean,sp_pd_min,sp_pd_reg)
    
    write.table(sp_pd_final,file = paste0(path_folder,"04_phylogenetic_distance/pairwise_phylogenetic_distance_nonthreatened/phylogenetic_distance_",reg_cod,"_nonthreatened.txt"),
                col.names = FALSE,row.names = FALSE,append = TRUE)
    
  }else{
    
    if(mod_fail!=0){
      
      # reg_failedinvasion is larger than 1000 (including 1000) and its mod is unequal to 0
      for (subid in 1:divi_fail) {
        
        sp_fail<-reg_failedinvasion[((subid-1)*1000+1):(subid*1000)]
        sp_fail<-gsub(" ","_",sp_fail,fixed = TRUE)
        
        sp_group<-c(sp_fail,sp_nati)
        
        sp_phylo<-keep.tip(tree_phylo,tip = sp_group)
        sp_pd<-cophenetic.phylo(sp_phylo)
        
        # form the matrix that rownames are failed invasion trees and colnames are native trees
        sp_pd_fail<-sp_pd[rownames(sp_pd)%in%sp_fail,colnames(sp_pd)%in%sp_nati]
        sp_pd_mean<-apply(sp_pd_fail, 1, mean)
        sp_pd_min<-apply(sp_pd_fail, 1, min)
        sp_pd_pro<-rep(0,length(sp_fail))
        
        sp_pd_rowname<-gsub("_"," ",rownames(sp_pd_fail),fixed = TRUE)
        sp_pd_reg<-rep(reg_cod,length(sp_fail))
        sp_pd_final<-cbind(sp_pd_rowname,sp_pd_pro,sp_pd_mean,sp_pd_min,sp_pd_reg)
        
        write.table(sp_pd_final,file = paste0(path_folder,"04_phylogenetic_distance/pairwise_phylogenetic_distance_nonthreatened/phylogenetic_distance_",reg_cod,"_nonthreatened.txt"),
                    col.names = FALSE,row.names = FALSE,append = TRUE)
        
        print(subid)
        
      }
      
      # extract the mod part
      sp_fail<-reg_failedinvasion[(divi_fail*1000+1):(divi_fail*1000+mod_fail)]
      sp_fail<-gsub(" ","_",sp_fail,fixed = TRUE)
      
      sp_group<-c(sp_fail,sp_nati)
      
      sp_phylo<-keep.tip(tree_phylo,tip = sp_group)
      sp_pd<-cophenetic.phylo(sp_phylo)
      
      # form the matrix that rownames are failed invasion trees and colnames are native trees
      sp_pd_fail<-sp_pd[rownames(sp_pd)%in%sp_fail,colnames(sp_pd)%in%sp_nati]
      sp_pd_mean<-apply(sp_pd_fail, 1, mean)
      sp_pd_min<-apply(sp_pd_fail, 1, min)
      sp_pd_pro<-rep(0,length(sp_fail))
      
      sp_pd_rowname<-gsub("_"," ",rownames(sp_pd_fail),fixed = TRUE)
      sp_pd_reg<-rep(reg_cod,length(sp_fail))
      sp_pd_final<-cbind(sp_pd_rowname,sp_pd_pro,sp_pd_mean,sp_pd_min,sp_pd_reg)
      
      write.table(sp_pd_final,file = paste0(path_folder,"04_phylogenetic_distance/pairwise_phylogenetic_distance_nonthreatened/phylogenetic_distance_",reg_cod,"_nonthreatened.txt"),
                  col.names = FALSE,row.names = FALSE,append = TRUE)
      
    }else{
      
      # reg_failedinvasion is larger than 1000 (including 1000) and its mod is equal to 0
      for (subid in 1:divi_fail) {
        
        sp_fail<-reg_failedinvasion[((subid-1)*1000+1):(subid*1000)]
        sp_fail<-gsub(" ","_",sp_fail,fixed = TRUE)
        
        sp_group<-c(sp_fail,sp_nati)
        
        sp_phylo<-keep.tip(tree_phylo,tip = sp_group)
        sp_pd<-cophenetic.phylo(sp_phylo)
        
        # form the matrix that rownames are failed invasion trees and colnames are native trees
        sp_pd_fail<-sp_pd[rownames(sp_pd)%in%sp_fail,colnames(sp_pd)%in%sp_nati]
        sp_pd_mean<-apply(sp_pd_fail, 1, mean)
        sp_pd_min<-apply(sp_pd_fail, 1, min)
        sp_pd_pro<-rep(0,length(sp_fail))
        
        sp_pd_rowname<-gsub("_"," ",rownames(sp_pd_fail),fixed = TRUE)
        sp_pd_reg<-rep(reg_cod,length(sp_fail))
        sp_pd_final<-cbind(sp_pd_rowname,sp_pd_pro,sp_pd_mean,sp_pd_min,sp_pd_reg)
        
        write.table(sp_pd_final,file = paste0(path_folder,"04_phylogenetic_distance/pairwise_phylogenetic_distance_nonthreatened/phylogenetic_distance_",reg_cod,"_nonthreatened.txt"),
                    col.names = FALSE,row.names = FALSE,append = TRUE)
        
        print(subid)
        
      }
      
    }
    
  }
  
  print(id)
  
}

# "FOR", "MRS", "NUE", "CPI", "MDV", "MPE" need manually deal with

# adding the column names
for (id in 1:length(tree_region)) {
  
  reg_name<-tree_region[id]
  
  reg_file<-read.table(paste0(path_folder,"04_phylogenetic_distance/pairwise_phylogenetic_distance_nonthreatened/phylogenetic_distance_",reg_name,"_nonthreatened.txt"))
  colnames(reg_file)<-c("species","probability","pdmean","pdmin","region")
  
  write.table(reg_file,file = paste0(path_folder,"04_phylogenetic_distance/pairwise_phylogenetic_distance_nonthreatened/phylogenetic_distance_",reg_name,"_nonthreatened.txt"),
              col.names = TRUE,row.names = FALSE)
  
  print(id)
  
}




