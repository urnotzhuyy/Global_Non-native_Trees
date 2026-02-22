rm(list=ls())

# loading required library
library(CoordinateCleaner)

# path
path_folder <- "D:/non_native_tree/"

#######################################################
# the Global Biodiversity Information Facility (GBIF)
#######################################################

# load the gbif rawdata
gbif.file<-list.files(paste0(path_folder,"gbif/rawdata/"))

for (id in 1:length(gbif.file)) {
  
  spname<-sub(".txt","",gbif.file[id],fixed = TRUE)
  
  sp.file<-read.table(paste0(path_folder,"gbif/rawdata/",spname,".txt"),header = TRUE)
  
  # remove locations without validity
  sp.file.val<-cc_val(sp.file,lon = "decimalLongitude",lat = "decimalLatitude")
  
  if(dim(sp.file.val)[1]!=0){
    
    # remove the zero coordinates
    sp.file.zero<-cc_zero(sp.file.val,lon = "decimalLongitude",lat = "decimalLatitude")
    
    if(dim(sp.file.zero)[1]!=0){
      
      # remove duplicated locations,because of validating of scientific name, species="kingdom"
      sp.file.dupl<-cc_dupl(sp.file.zero,lon = "decimalLongitude",lat = "decimalLatitude",species = "kingdom")
      
      if(dim(sp.file.dupl)[1]!=0){
        
        # remove locations located in sea
        sp.file.sea<-cc_sea(sp.file.dupl,lon = "decimalLongitude",lat = "decimalLatitude")
        
        if(dim(sp.file.sea)[1]!=0){
          
          # remove locations located in the centroids of country and province
          sp.file.centroid<-cc_cen(sp.file.sea,lon = "decimalLongitude",lat = "decimalLatitude",test = "both")
          
          if(dim(sp.file.centroid)[1]!=0){
            
            # remove locations located in the capitals of country
            sp.file.capital<-cc_cap(sp.file.centroid,lon = "decimalLongitude",lat = "decimalLatitude")
            
            if(dim(sp.file.capital)[1]!=0){
              
              # extract the information
              sp.file.basic<-sp.file.capital[,c("species","decimalLongitude","decimalLatitude","basisOfRecord")]
              write.table(sp.file.basic,file = paste0(path_folder,"gbif/cleaning/",spname,"_cleaning.txt"),
                          col.names = TRUE,row.names = FALSE)
              
            }
            
          }
          
        }
        
      }
      
    }
    
  }
  
  print(id)
  
}

# check the categories of basisofrecord
gbif.cleaning<-list.files(paste0(path_folder,"gbif/cleaning/"))

bor.all<-c()

for (id in 1:length(gbif.cleaning)) {
  
  spname<-sub("_cleaning.txt","",gbif.cleaning[id],fixed = TRUE)
  
  sp.file<-read.table(paste0(path_folder,"gbif/cleaning/",spname,"_cleaning.txt"),header = TRUE)
  
  sp.bor<-unique(as.character(sp.file$basisOfRecord))
  
  bor.all<-c(bor.all,sp.bor)
  
  print(id)
  
}

bor.uni<-unique(bor.all)

# remove the fossil records
gbif.cleaning<-list.files(paste0(path_folder,"gbif/cleaning/"))

for (id in 1:length(gbif.cleaning)) {
  
  spname<-sub("_cleaning.txt","",gbif.cleaning[id],fixed = TRUE)
  
  sp.file<-read.table(paste0(path_folder,"gbif/cleaning/",spname,"_cleaning.txt"),header = TRUE)
  
  sp.file.ext<-sp.file[which(sp.file$basisOfRecord!="FOSSIL_SPECIMEN"),]
  
  if(dim(sp.file.ext)[1]!=0){
    write.table(sp.file.ext,file = paste0(path_folder,"gbif/cleaning/",spname,"_cleaning.txt"),
                col.names = TRUE,row.names = FALSE)
  }else{
    file.remove(paste0(path_folder,"gbif/cleaning/",spname,"_cleaning.txt"))
  }
  
  print(id)
  
}

# arranging the cleaning files
gbif.cleaning<-list.files(paste0(path_folder,"gbif/cleaning/"))

for (id in 1:length(gbif.cleaning)) {
  
  spname<-sub("_cleaning.txt","",gbif.cleaning[id],fixed = TRUE)
  
  sp.file<-read.table(paste0(path_folder,"gbif/cleaning/",spname,"_cleaning.txt"),header = TRUE)
  sp.file$confirmed_name<-spname
  
  sp.file.arrange<-sp.file[,c("confirmed_name","decimalLongitude","decimalLatitude")]
  colnames(sp.file.arrange)<-c("acceptedname","longitude","latitude")
  
  write.table(sp.file.arrange,file = paste0(path_folder,"gbif/cleaning/",spname,"_cleaning.txt"),
              col.names = TRUE,row.names = FALSE)
  
  print(id)
  
}

# removing the location flagged in the GBIF headquarters
gbif.cleaning<-list.files(paste0(path_folder,"gbif/cleaning/"))

for (id in 1:length(gbif.cleaning)) {
  
  spname<-gbif.cleaning[id]
  spname<-sub("_cleaning.txt","",spname,fixed = TRUE)
  
  sp.file<-read.table(paste0(path_folder,"gbif/cleaning/",spname,"_cleaning.txt"),header = TRUE)
  
  sp.file.gbif<-cc_gbif(sp.file,lon = "longitude",lat = "latitude")
  
  if(dim(sp.file.gbif)[1]!=0){
    write.table(sp.file.gbif,file = paste0(path_folder,"gbif/cleaning/",spname,"_cleaning.txt"),
                col.names = TRUE,row.names = FALSE)
  }else{
    file.remove(paste0(path_folder,"gbif/cleaning/",spname,"_cleaning.txt"))
  }
  
  print(id)
  
}

#######################################
# the Atlas of Living Australia (ALA)
#######################################

# load the ala rawdata
ala.file<-list.files(paste0(path_folder,"ala/rawdata/"))

for (id in 1:length(ala.file)) {
  
  spname<-sub(".txt","",ala.file[id],fixed = TRUE)
  
  sp.file<-read.table(paste0(path_folder,"ala/rawdata/",spname,".txt"),header = TRUE)
  
  # adding column of confirmed name
  sp.file$confirmed_name<-spname
  
  # remove locations without validity
  sp.file.val<-cc_val(sp.file,lon = "decimalLongitude",lat = "decimalLatitude")
  
  if(dim(sp.file.val)[1]!=0){
    
    # remove the zero coordinates
    sp.file.zero<-cc_zero(sp.file.val,lon = "decimalLongitude",lat = "decimalLatitude")
    
    if(dim(sp.file.zero)[1]!=0){
      
      # remove duplicated locations
      sp.file.dupl<-cc_dupl(sp.file.zero,lon = "decimalLongitude",lat = "decimalLatitude",species = "confirmed_name")
      
      if(dim(sp.file.dupl)[1]!=0){
        
        # remove locations located in sea
        sp.file.sea<-cc_sea(sp.file.dupl,lon = "decimalLongitude",lat = "decimalLatitude")
        
        if(dim(sp.file.sea)[1]!=0){
          
          # remove locations located in the centroids of country and province
          sp.file.centroid<-cc_cen(sp.file.sea,lon = "decimalLongitude",lat = "decimalLatitude",test = "both")
          
          if(dim(sp.file.centroid)[1]!=0){
            
            # remove locations located in the capitals of country
            sp.file.capital<-cc_cap(sp.file.centroid,lon = "decimalLongitude",lat = "decimalLatitude")
            
            if(dim(sp.file.capital)[1]!=0){
              
              # extract the information
              sp.file.basic<-sp.file.capital[,c("scientificName","decimalLongitude","decimalLatitude","basisOfRecord","confirmed_name")]
              write.table(sp.file.basic,file = paste0(path_folder,"ala/cleaning/",spname,"_cleaning.txt"),
                          col.names = TRUE,row.names = FALSE)
              
            }
            
          }
          
        }
        
      }
      
    }
    
  }
  
  print(id)
  
}

# check the categories of basisofrecord
ala.cleaning<-list.files(paste0(path_folder,"ala/cleaning/"))

bor.all<-c()

for (id in 1:length(ala.cleaning)) {
  
  spname<-sub("_cleaning.txt","",ala.cleaning[id],fixed = TRUE)
  
  sp.file<-read.table(paste0(path_folder,"ala/cleaning/",spname,"_cleaning.txt"),header = TRUE)
  
  sp.bor<-unique(as.character(sp.file$basisOfRecord))
  
  bor.all<-c(bor.all,sp.bor)
  
  print(id)
  
}

bor.uni<-unique(bor.all)

# arranging the cleaning files
ala.cleaning<-list.files(paste0(path_folder,"ala/cleaning/"))

for (id in 1:length(ala.cleaning)) {
  
  spname<-sub("_cleaning.txt","",ala.cleaning[id],fixed = TRUE)
  
  sp.file<-read.table(paste0(path_folder,"ala/cleaning/",spname,"_cleaning.txt"),header = TRUE)
  
  sp.file.arrange<-sp.file[,c("confirmed_name","decimalLongitude","decimalLatitude")]
  colnames(sp.file.arrange)<-c("acceptedname","longitude","latitude")
  
  write.table(sp.file.arrange,file = paste0(path_folder,"ala/cleaning/",spname,"_cleaning.txt"),
              col.names = TRUE,row.names = FALSE)
  
  print(id)
  
}

#############################################################################
# the public domain Botanical Information and Ecological Network v.4 (BIEN)
#############################################################################

# load the bien rawdata
bien.file<-list.files(paste0(path_folder,"bien/rawdata/"))

for (id in 1:length(bien.file)) {
  
  spname<-sub(".txt","",bien.file[id],fixed = TRUE)
  sp.file<-read.table(paste0(path_folder,"bien/rawdata/",spname,".txt"),header = TRUE)
  
  # adding column of confirmed name
  sp.file$confirmed_name<-spname
  
  # remove locations without validity
  sp.file.val<-cc_val(sp.file,lon = "longitude",lat = "latitude")
  
  if(dim(sp.file.val)[1]!=0){
    
    # remove the zero coordinates
    sp.file.zero<-cc_zero(sp.file.val,lon = "longitude",lat = "latitude")
    
    if(dim(sp.file.zero)[1]!=0){
      
      # remove duplicated locations
      sp.file.dupl<-cc_dupl(sp.file.zero,lon = "longitude",lat = "latitude",species = "confirmed_name")
      
      if(dim(sp.file.dupl)[1]!=0){
        
        # remove locations located in sea
        sp.file.sea<-cc_sea(sp.file.dupl,lon = "longitude",lat = "latitude")
        
        if(dim(sp.file.sea)[1]!=0){
          
          # remove locations located in the centroids of country and province
          sp.file.centroid<-cc_cen(sp.file.sea,lon = "longitude",lat = "latitude",test = "both")
          
          if(dim(sp.file.centroid)[1]!=0){
            
            # remove locations located in the capitals of country
            sp.file.capital<-cc_cap(sp.file.centroid,lon = "longitude",lat = "latitude")
            
            if(dim(sp.file.capital)[1]!=0){
              
              # extract the information
              sp.file.basic<-sp.file.capital[,c("scrubbed_species_binomial","longitude","latitude","confirmed_name")]
              write.table(sp.file.basic,file = paste0(path_folder,"bien/cleaning/",spname,"_cleaning.txt"),
                          col.names = TRUE,row.names = FALSE)
              
            }
            
          }
          
        }
        
      }
      
    }
    
  }
  
  print(id)
  
}

# arranging the cleaning files
bien.cleaning<-list.files(paste0(path_folder,"bien/cleaning/"))

for (id in 1:length(bien.cleaning)) {
  
  spname<-sub("_cleaning.txt","",bien.cleaning[id],fixed = TRUE)
  
  sp.file<-read.table(paste0(path_folder,"bien/cleaning/",spname,"_cleaning.txt"),header = TRUE)
  
  sp.file.arrange<-sp.file[,c("confirmed_name","longitude","latitude")]
  colnames(sp.file.arrange)<-c("acceptedname","longitude","latitude")
  
  write.table(sp.file.arrange,file = paste0(path_folder,"bien/cleaning/",spname,"_cleaning.txt"),
              col.names = TRUE,row.names = FALSE)
  
  print(id)
  
}

###########################################################
# the Biodiversity Information Serving Our Nation (BISON)
###########################################################

# load the bison rawdata
bison.file<-list.files(paste0(path_folder,"bison/rawdata/"))

for (id in 1:length(bison.file)) {
  
  spname<-sub(".txt","",bison.file[id],fixed = TRUE)
  
  sp.file<-read.table(paste0(path_folder,"bison/rawdata/",spname,".txt"),header = TRUE)
  
  # adding column of confirmed name
  sp.file$confirmed_name<-spname
  
  # remove locations without validity
  sp.file.val<-cc_val(sp.file,lon = "decimalLongitude",lat = "decimalLatitude")
  
  if(dim(sp.file.val)[1]!=0){
    
    # remove the zero coordinates
    sp.file.zero<-cc_zero(sp.file.val,lon = "decimalLongitude",lat = "decimalLatitude")
    
    if(dim(sp.file.zero)[1]!=0){
      
      # remove duplicated locations
      sp.file.dupl<-cc_dupl(sp.file.zero,lon = "decimalLongitude",lat = "decimalLatitude",species = "confirmed_name")
      
      if(dim(sp.file.dupl)[1]!=0){
        
        # remove locations located in sea
        sp.file.sea<-cc_sea(sp.file.dupl,lon = "decimalLongitude",lat = "decimalLatitude")
        
        if(dim(sp.file.sea)[1]!=0){
          
          # remove locations located in the centroids of country and province
          sp.file.centroid<-cc_cen(sp.file.sea,lon = "decimalLongitude",lat = "decimalLatitude",test = "both")
          
          if(dim(sp.file.centroid)[1]!=0){
            
            # remove locations located in the capitals of country
            sp.file.capital<-cc_cap(sp.file.centroid,lon = "decimalLongitude",lat = "decimalLatitude")
            
            if(dim(sp.file.capital)[1]!=0){
              
              # extract the information
              sp.file.basic<-sp.file.capital[,c("scientificName","decimalLongitude","decimalLatitude","basisOfRecord","confirmed_name")]
              write.table(sp.file.basic,file = paste0(path_folder,"bison/cleaning/",spname,"_cleaning.txt"),
                          col.names = TRUE,row.names = FALSE)
              
            }
            
          }
          
        }
        
      }
      
    }
    
  }
  
  print(id)
  
}

# check the categories of basisofrecord
bison.cleaning<-list.files(paste0(path_folder,"bison/cleaning/"))

bor.all<-c()

for (id in 1:length(bison.cleaning)) {
  
  spname<-sub("_cleaning.txt","",bison.cleaning[id],fixed = TRUE)
  
  sp.file<-read.table(paste0(path_folder,"bison/cleaning/",spname,"_cleaning.txt"),header = TRUE)
  
  sp.bor<-unique(as.character(sp.file$basisOfRecord))
  
  bor.all<-c(bor.all,sp.bor)
  
  print(id)
  
}

bor.uni<-unique(bor.all)

# remove the fossil records
bison.cleaning<-list.files(paste0(path_folder,"bison/cleaning/"))

for (id in 1:length(bison.cleaning)) {
  
  spname<-sub("_cleaning.txt","",bison.cleaning[id],fixed = TRUE)
  
  sp.file<-read.table(paste0(path_folder,"bison/cleaning/",spname,"_cleaning.txt"),header = TRUE)
  
  sp.file.ext<-sp.file[which(sp.file$basisOfRecord!="FOSSIL_SPECIMEN"),]
  
  if(dim(sp.file.ext)[1]!=0){
    write.table(sp.file.ext,file = paste0(path_folder,"bison/cleaning/",spname,"_cleaning.txt"),
                col.names = TRUE,row.names = FALSE)
  }else{
    file.remove(paste0(path_folder,"bison/cleaning/",spname,"_cleaning.txt"))
  }
  
  print(id)
  
}

# arranging the cleaning files
bison.cleaning<-list.files(paste0(path_folder,"bison/cleaning/"))

for (id in 1:length(bison.cleaning)) {
  
  spname<-sub("_cleaning.txt","",bison.cleaning[id],fixed = TRUE)
  
  sp.file<-read.table(paste0(path_folder,"bison/cleaning/",spname,"_cleaning.txt"),header = TRUE)
  
  sp.file.arrange<-sp.file[,c("confirmed_name","decimalLongitude","decimalLatitude")]
  colnames(sp.file.arrange)<-c("acceptedname","longitude","latitude")
  
  write.table(sp.file.arrange,file = paste0(path_folder,"bison/cleaning/",spname,"_cleaning.txt"),
              col.names = TRUE,row.names = FALSE)
  
  print(id)
  
}

#############################################################################
# Latin American Seasonally Dry Tropical Forest Floristic Network (DRYFLOR)
#############################################################################

# load the dryflor rawdata
dryflor.file<-list.files(paste0(path_folder,"dryflor/rawdata/"))

for (id in 1:length(dryflor.file)) {
  
  spname<-sub(".txt","",dryflor.file[id],fixed = TRUE)
  
  sp.file<-read.table(paste0(path_folder,"dryflor/rawdata/",spname,".txt"),header = TRUE)
  
  # adding column of confirmed name
  sp.file$confirmed_name<-spname
  
  # remove locations without validity
  sp.file.val<-cc_val(sp.file,lon = "Long",lat = "Lat")
  
  if(dim(sp.file.val)[1]!=0){
    
    # remove the zero coordinates
    sp.file.zero<-cc_zero(sp.file.val,lon = "Long",lat = "Lat")
    
    if(dim(sp.file.zero)[1]!=0){
      
      # remove duplicated locations
      sp.file.dupl<-cc_dupl(sp.file.zero,lon = "Long",lat = "Lat",species = "confirmed_name")
      
      if(dim(sp.file.dupl)[1]!=0){
        
        # remove locations located in sea
        sp.file.sea<-cc_sea(sp.file.dupl,lon = "Long",lat = "Lat")
        
        if(dim(sp.file.sea)[1]!=0){
          
          # remove locations located in the centroids of country and province
          sp.file.centroid<-cc_cen(sp.file.sea,lon = "Long",lat = "Lat",test = "both")
          
          if(dim(sp.file.centroid)[1]!=0){
            
            # remove locations located in the capitals of country
            sp.file.capital<-cc_cap(sp.file.centroid,lon = "Long",lat = "Lat")
            
            if(dim(sp.file.capital)[1]!=0){
              
              # extract the information
              sp.file.basic<-sp.file.capital[,c("sp_","Long","Lat","confirmed_name")]
              write.table(sp.file.basic,file = paste0(path_folder,"dryflor/cleaning/",spname,"_cleaning.txt"),
                          col.names = TRUE,row.names = FALSE)
              
            }
            
          }
          
        }
        
      }
      
    }
    
  }
  
  print(id)
  
}

# arranging the cleaning files
dryflor.cleaning<-list.files(paste0(path_folder,"dryflor/cleaning/"))

for (id in 1:length(dryflor.cleaning)) {
  
  spname<-sub("_cleaning.txt","",dryflor.cleaning[id],fixed = TRUE)
  
  sp.file<-read.table(paste0(path_folder,"dryflor/cleaning/",spname,"_cleaning.txt"),header = TRUE)
  
  sp.file.arrange<-sp.file[,c("confirmed_name","Long","Lat")]
  colnames(sp.file.arrange)<-c("acceptedname","longitude","latitude")
  
  write.table(sp.file.arrange,file = paste0(path_folder,"dryflor/cleaning/",spname,"_cleaning.txt"),
              col.names = TRUE,row.names = FALSE)
  
  print(id)
  
}

#####################################################
# the Integrated Digitized Biocollections (iDigBio)
#####################################################

# load the idigbio rawdata
idigbio.file<-list.files(paste0(path_folder,"idigbio/rawdata/"))

for (id in 1:length(idigbio.file)) {
  
  spname<-sub(".txt","",idigbio.file[id],fixed = TRUE)
  
  sp.file<-read.table(paste0(path_folder,"idigbio/rawdata/",spname,".txt"),header = TRUE)
  
  # adding column of confirmed name
  sp.file$confirmed_name<-spname
  
  # remove locations without validity
  sp.file.val<-cc_val(sp.file,lon = "dwc.decimalLongitude",lat = "dwc.decimalLatitude")
  
  if(dim(sp.file.val)[1]!=0){
    
    # remove the zero coordinates
    sp.file.zero<-cc_zero(sp.file.val,lon = "dwc.decimalLongitude",lat = "dwc.decimalLatitude")
    
    if(dim(sp.file.zero)[1]!=0){
      
      # remove duplicated locations
      sp.file.dupl<-cc_dupl(sp.file.zero,lon = "dwc.decimalLongitude",lat = "dwc.decimalLatitude",species = "confirmed_name")
      
      if(dim(sp.file.dupl)[1]!=0){
        
        # remove locations located in sea
        sp.file.sea<-cc_sea(sp.file.dupl,lon = "dwc.decimalLongitude",lat = "dwc.decimalLatitude")
        
        if(dim(sp.file.sea)[1]!=0){
          
          # remove locations located in the centroids of country and province
          sp.file.centroid<-cc_cen(sp.file.sea,lon = "dwc.decimalLongitude",lat = "dwc.decimalLatitude",test = "both")
          
          if(dim(sp.file.centroid)[1]!=0){
            
            # remove locations located in the capitals of country
            sp.file.capital<-cc_cap(sp.file.centroid,lon = "dwc.decimalLongitude",lat = "dwc.decimalLatitude")
            
            if(dim(sp.file.capital)[1]!=0){
              
              # extract the information
              sp.file.basic<-sp.file.capital[,c("dwc.scientificName","dwc.decimalLongitude","dwc.decimalLatitude","dwc.basisOfRecord","confirmed_name")]
              write.table(sp.file.basic,file = paste0(path_folder,"idigbio/cleaning/",spname,"_cleaning.txt"),
                          col.names = TRUE,row.names = FALSE)
              
            }
            
          }
          
        }
        
      }
      
    }
    
  }
  
  print(id)
  
}

# check the categories of basisofrecord
idigbio.cleaning<-list.files(paste0(path_folder,"idigbio/cleaning/"))

bor.all<-c()

for (id in 1:length(idigbio.cleaning)) {
  
  spname<-sub("_cleaning.txt","",idigbio.cleaning[id],fixed = TRUE)
  
  sp.file<-read.table(paste0(path_folder,"idigbio/cleaning/",spname,"_cleaning.txt"),header = TRUE)
  
  sp.bor<-unique(as.character(sp.file$dwc.basisOfRecord))
  
  bor.all<-c(bor.all,sp.bor)
  
  print(id)
  
}

bor.uni<-unique(bor.all)

# remove the fossil records
idigbio.cleaning<-list.files(paste0(path_folder,"idigbio/cleaning/"))

for (id in 1:length(idigbio.cleaning)) {
  
  spname<-sub("_cleaning.txt","",idigbio.cleaning[id],fixed = TRUE)
  
  sp.file<-read.table(paste0(path_folder,"idigbio/cleaning/",spname,"_cleaning.txt"),header = TRUE)
  
  sp.file.ext<-sp.file[which(sp.file$dwc.basisOfRecord!="FossilSpecimen"),]
  
  if(dim(sp.file.ext)[1]!=0){
    write.table(sp.file.ext,file = paste0(path_folder,"idigbio/cleaning/",spname,"_cleaning.txt"),
                col.names = TRUE,row.names = FALSE)
  }else{
    file.remove(paste0(path_folder,"idigbio/cleaning/",spname,"_cleaning.txt"))
  }
  
  print(id)
  
}

# arranging the cleaning files
idigbio.cleaning<-list.files(paste0(path_folder,"idigbio/cleaning/"))

for (id in 1:length(idigbio.cleaning)) {
  
  spname<-sub("_cleaning.txt","",idigbio.cleaning[id],fixed = TRUE)
  
  sp.file<-read.table(paste0(path_folder,"idigbio/cleaning/",spname,"_cleaning.txt"),header = TRUE)
  
  sp.file.arrange<-sp.file[,c("confirmed_name","dwc.decimalLongitude","dwc.decimalLatitude")]
  colnames(sp.file.arrange)<-c("acceptedname","longitude","latitude")
  
  write.table(sp.file.arrange,file = paste0(path_folder,"idigbio/cleaning/",spname,"_cleaning.txt"),
              col.names = TRUE,row.names = FALSE)
  
  print(id)
  
}

#########################################################
# International Union for Conservation of Nature (IUCN)
#########################################################

# load the iucn rawdata
iucn.file<-list.files(paste0(path_folder,"iucn/rawdata/"))

for (id in 1:length(iucn.file)) {
  
  spname<-sub(".txt","",iucn.file[id],fixed = TRUE)
  
  sp.file<-read.table(paste0(path_folder,"iucn/rawdata/",spname,".txt"),header = TRUE)
  
  # adding column of confirmed name
  sp.file$confirmed_name<-spname
  
  # remove locations without validity
  sp.file.val<-cc_val(sp.file,lon = "longitude",lat = "latitude")
  
  if(dim(sp.file.val)[1]!=0){
    
    # remove the zero coordinates
    sp.file.zero<-cc_zero(sp.file.val,lon = "longitude",lat = "latitude")
    
    if(dim(sp.file.zero)[1]!=0){
      
      # remove duplicated locations
      sp.file.dupl<-cc_dupl(sp.file.zero,lon = "longitude",lat = "latitude",species = "confirmed_name")
      
      if(dim(sp.file.dupl)[1]!=0){
        
        # remove locations located in sea
        sp.file.sea<-cc_sea(sp.file.dupl,lon = "longitude",lat = "latitude")
        
        if(dim(sp.file.sea)[1]!=0){
          
          # remove locations located in the centroids of country and province
          sp.file.centroid<-cc_cen(sp.file.sea,lon = "longitude",lat = "latitude",test = "both")
          
          if(dim(sp.file.centroid)[1]!=0){
            
            # remove locations located in the capitals of country
            sp.file.capital<-cc_cap(sp.file.centroid,lon = "longitude",lat = "latitude")
            
            if(dim(sp.file.capital)[1]!=0){
              
              # extract the information
              sp.file.basic<-sp.file.capital[,c("binomial","longitude","latitude","basisofrec","presence","confirmed_name")]
              write.table(sp.file.basic,file = paste0(path_folder,"iucn/cleaning/",spname,"_cleaning.txt"),
                          col.names = TRUE,row.names = FALSE)
              
            }
            
          }
          
        }
        
      }
      
    }
    
  }
  
  print(id)
  
}

# check the categories of basisofrec and presence
iucn.cleaning<-list.files(paste0(path_folder,"iucn/cleaning/"))

bor.all<-c()
pre.all<-c()

for (id in 1:length(iucn.cleaning)) {
  
  spname<-sub("_cleaning.txt","",iucn.cleaning[id],fixed = TRUE)
  
  sp.file<-read.table(paste0(path_folder,"iucn/cleaning/",spname,"_cleaning.txt"),header = TRUE)
  
  sp.bor<-unique(as.character(sp.file$basisofrec))
  sp.pre<-unique(as.character(sp.file$presence))
  
  bor.all<-c(bor.all,sp.bor)
  pre.all<-c(pre.all,sp.pre)
  
  print(id)
  
}

bor.uni<-unique(bor.all)
pre.uni<-unique(pre.all)

##remove the fossil records, only keep the presence is one
iucn.cleaning<-list.files(paste0(path_folder,"iucn/cleaning/"))

for (id in 1:length(iucn.cleaning)) {
  
  spname<-sub("_cleaning.txt","",iucn.cleaning[id],fixed = TRUE)
  
  sp.file<-read.table(paste0(path_folder,"iucn/cleaning/",spname,"_cleaning.txt"),header = TRUE)
  
  sp.file.ext<-sp.file[which(sp.file$basisofrec!="FOSSIL_SPECIMEN"),]
  sp.file.ext<-sp.file.ext[which(sp.file.ext$basisofrec!="FossilSpecimen"),]
  sp.file.ext<-sp.file.ext[which(sp.file.ext$presence==1),]
  
  if(dim(sp.file.ext)[1]!=0){
    write.table(sp.file.ext,file = paste0(path_folder,"iucn/cleaning/",spname,"_cleaning.txt"),
                col.names = TRUE,row.names = FALSE)
  }else{
    file.remove(paste0(path_folder,"iucn/cleaning/",spname,"_cleaning.txt"))
  }
  
  print(id)
  
}

# arranging the cleaning files
iucn.cleaning<-list.files(paste0(path_folder,"iucn/cleaning/"))

for (id in 1:length(iucn.cleaning)) {
  
  spname<-sub("_cleaning.txt","",iucn.cleaning[id],fixed = TRUE)
  
  sp.file<-read.table(paste0(path_folder,"iucn/cleaning/",spname,"_cleaning.txt"),header = TRUE)
  
  sp.file.arrange<-sp.file[,c("confirmed_name","longitude","latitude")]
  colnames(sp.file.arrange)<-c("acceptedname","longitude","latitude")
  
  write.table(sp.file.arrange,file = paste0(path_folder,"iucn/cleaning/",spname,"_cleaning.txt"),
              col.names = TRUE,row.names = FALSE)
  
  print(id)
  
}

#######################################################
# National Specimen Information Infrastructure (NSII)
#######################################################

# load the nsii rawdata
nsii.file<-list.files(paste0(path_folder,"nsii/rawdata_longitude_latitude/"))

for (id in 1:length(nsii.file)) {
  
  spname<-sub(".txt","",nsii.file[id],fixed = TRUE)
  
  sp.file<-read.table(paste0(path_folder,"nsii/rawdata_longitude_latitude/",spname,".txt"),header = TRUE)
  
  # adding column of confirmed name
  sp.file$confirmed_name<-spname
  
  # remove locations without validity
  sp.file.val<-cc_val(sp.file,lon = "longitude",lat = "latitude")
  
  if(dim(sp.file.val)[1]!=0){
    
    # remove the zero coordinates
    sp.file.zero<-cc_zero(sp.file.val,lon = "longitude",lat = "latitude")
    
    if(dim(sp.file.zero)[1]!=0){
      
      # remove duplicated locations
      sp.file.dupl<-cc_dupl(sp.file.zero,lon = "longitude",lat = "latitude",species = "confirmed_name")
      
      if(dim(sp.file.dupl)[1]!=0){
        
        # remove locations located in sea
        sp.file.sea<-cc_sea(sp.file.dupl,lon = "longitude",lat = "latitude")
        
        if(dim(sp.file.sea)[1]!=0){
          
          # remove locations located in the centroids of country and province
          sp.file.centroid<-cc_cen(sp.file.sea,lon = "longitude",lat = "latitude",test = "both")
          
          if(dim(sp.file.centroid)[1]!=0){
            
            # remove locations located in the capitals of country
            sp.file.capital<-cc_cap(sp.file.centroid,lon = "longitude",lat = "latitude")
            
            if(dim(sp.file.capital)[1]!=0){
              
              # extract the information
              sp.file.basic<-sp.file.capital[,c("scientificname","longitude","latitude","confirmed_name")]
              write.table(sp.file.basic,file = paste0(path_folder,"nsii/cleaning/",spname,"_cleaning.txt"),
                          col.names = TRUE,row.names = FALSE)
              
            }
            
          }
          
        }
        
      }
      
    }
    
  }
  
  print(id)
  
}

# arranging the cleaning files
nsii.cleaning<-list.files(paste0(path_folder,"nsii/cleaning/"))

for (id in 1:length(nsii.cleaning)) {
  
  spname<-sub("_cleaning.txt","",nsii.cleaning[id],fixed = TRUE)
  
  sp.file<-read.table(paste0(path_folder,"nsii/cleaning/",spname,"_cleaning.txt"),header = TRUE)
  
  sp.file.arrange<-sp.file[,c("confirmed_name","longitude","latitude")]
  colnames(sp.file.arrange)<-c("acceptedname","longitude","latitude")
  
  write.table(sp.file.arrange,file = paste0(path_folder,"nsii/cleaning/",spname,"_cleaning.txt"),
              col.names = TRUE,row.names = FALSE)
  
  print(id)
  
}

####################################
# Plant Photo Bank of China (PPBC)
####################################

# load the pbcc rawdata
pbcc.file<-list.files(paste0(path_folder,"pbcc/rawdata_longitude_latitude/"))

for (id in 1:length(pbcc.file)) {
  
  spname<-sub(".txt","",pbcc.file[id],fixed = TRUE)
  
  sp.file<-read.table(paste0(path_folder,"pbcc/rawdata_longitude_latitude/",spname,".txt"),header = TRUE)
  
  # adding column of confirmed name
  sp.file$confirmed_name<-spname
  
  # remove locations without validity
  sp.file.val<-cc_val(sp.file,lon = "longitude",lat = "latitude")
  
  if(dim(sp.file.val)[1]!=0){
    
    # remove the zero coordinates
    sp.file.zero<-cc_zero(sp.file.val,lon = "longitude",lat = "latitude")
    
    if(dim(sp.file.zero)[1]!=0){
      
      # remove duplicated locations
      sp.file.dupl<-cc_dupl(sp.file.zero,lon = "longitude",lat = "latitude",species = "confirmed_name")
      
      if(dim(sp.file.dupl)[1]!=0){
        
        # remove locations located in sea
        sp.file.sea<-cc_sea(sp.file.dupl,lon = "longitude",lat = "latitude")
        
        if(dim(sp.file.sea)[1]!=0){
          
          # remove locations located in the centroids of country and province
          sp.file.centroid<-cc_cen(sp.file.sea,lon = "longitude",lat = "latitude",test = "both")
          
          if(dim(sp.file.centroid)[1]!=0){
            
            # remove locations located in the capitals of country
            sp.file.capital<-cc_cap(sp.file.centroid,lon = "longitude",lat = "latitude")
            
            if(dim(sp.file.capital)[1]!=0){
              
              # extract the information
              sp.file.basic<-sp.file.capital[,c("scientificname","longitude","latitude","confirmed_name")]
              write.table(sp.file.basic,file = paste0(path_folder,"pbcc/cleaning/",spname,"_cleaning.txt"),
                          col.names = TRUE,row.names = FALSE)
              
            }
            
          }
          
        }
        
      }
      
    }
    
  }
  
  print(id)
  
}

# arranging the cleaning files
pbcc.cleaning<-list.files(paste0(path_folder,"pbcc/cleaning/"))

for (id in 1:length(pbcc.cleaning)) {
  
  spname<-sub("_cleaning.txt","",pbcc.cleaning[id],fixed = TRUE)
  
  sp.file<-read.table(paste0(path_folder,"pbcc/cleaning/",spname,"_cleaning.txt"),header = TRUE)
  
  sp.file.arrange<-sp.file[,c("confirmed_name","longitude","latitude")]
  colnames(sp.file.arrange)<-c("acceptedname","longitude","latitude")
  
  write.table(sp.file.arrange,file = paste0(path_folder,"pbcc/cleaning/",spname,"_cleaning.txt"),
              col.names = TRUE,row.names = FALSE)
  
  print(id)
  
}

###########
# RAINBIO
###########

# load the rainbio rawdata
rainbio.file<-list.files(paste0(path_folder,"rainbio/rawdata/"))

for (id in 1:length(rainbio.file)) {
  
  spname<-sub(".txt","",rainbio.file[id],fixed = TRUE)
  
  sp.file<-read.table(paste0(path_folder,"rainbio/rawdata/",spname,".txt"),header = TRUE)
  
  # adding column of confirmed name
  sp.file$confirmed_name<-spname
  
  # remove locations without validity
  sp.file.val<-cc_val(sp.file,lon = "decimalLongitude",lat = "decimalLatitude")
  
  if(dim(sp.file.val)[1]!=0){
    
    # remove the zero coordinates
    sp.file.zero<-cc_zero(sp.file.val,lon = "decimalLongitude",lat = "decimalLatitude")
    
    if(dim(sp.file.zero)[1]!=0){
      
      # remove duplicated locations
      sp.file.dupl<-cc_dupl(sp.file.zero,lon = "decimalLongitude",lat = "decimalLatitude",species = "confirmed_name")
      
      if(dim(sp.file.dupl)[1]!=0){
        
        # remove locations located in sea
        sp.file.sea<-cc_sea(sp.file.dupl,lon = "decimalLongitude",lat = "decimalLatitude")
        
        if(dim(sp.file.sea)[1]!=0){
          
          # remove locations located in the centroids of country and province
          sp.file.centroid<-cc_cen(sp.file.sea,lon = "decimalLongitude",lat = "decimalLatitude",test = "both")
          
          if(dim(sp.file.centroid)[1]!=0){
            
            # remove locations located in the capitals of country
            sp.file.capital<-cc_cap(sp.file.centroid,lon = "decimalLongitude",lat = "decimalLatitude")
            
            if(dim(sp.file.capital)[1]!=0){
              
              # extract the information
              sp.file.basic<-sp.file.capital[,c("species","decimalLongitude","decimalLatitude","basisOfRecord","confirmed_name")]
              write.table(sp.file.basic,file = paste0(path_folder,"rainbio/cleaning/",spname,"_cleaning.txt"),
                          col.names = TRUE,row.names = FALSE)
              
            }
            
          }
          
        }
        
      }
      
    }
    
  }
  
  print(id)
  
}

# check the categories of basisofrecord
rainbio.cleaning<-list.files(paste0(path_folder,"rainbio/cleaning/"))

bor.all<-c()

for (id in 1:length(rainbio.cleaning)) {
  
  spname<-sub("_cleaning.txt","",rainbio.cleaning[id],fixed = TRUE)
  
  sp.file<-read.table(paste0(path_folder,"rainbio/cleaning/",spname,"_cleaning.txt"),header = TRUE)
  
  sp.bor<-unique(as.character(sp.file$basisOfRecord))
  
  bor.all<-c(bor.all,sp.bor)
  
  print(id)
  
}

bor.uni<-unique(bor.all)

# arranging the cleaning files
rainbio.cleaning<-list.files(paste0(path_folder,"rainbio/cleaning/"))

for (id in 1:length(rainbio.cleaning)) {
  
  spname<-sub("_cleaning.txt","",rainbio.cleaning[id],fixed = TRUE)
  
  sp.file<-read.table(paste0(path_folder,"rainbio/cleaning/",spname,"_cleaning.txt"),header = TRUE)
  
  sp.file.arrange<-sp.file[,c("confirmed_name","decimalLongitude","decimalLatitude")]
  colnames(sp.file.arrange)<-c("acceptedname","longitude","latitude")
  
  write.table(sp.file.arrange,file = paste0(path_folder,"rainbio/cleaning/",spname,"_cleaning.txt"),
              col.names = TRUE,row.names = FALSE)
  
  print(id)
  
}








