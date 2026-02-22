rm(list=ls())

# loading required library
library(rWCVP)

# path
path_folder <- "D:/non_native_tree/"

# load the list of global trees
sp_list <- read.csv(paste0(path_folder, "speciesname_confirm/global_tree.csv"))

# validating the names via the WCVP
sp_list_match <- wcvp_match_names(names_df = sp_list,
                                  name_col = "TaxonName",
                                  author_col = "Author",
                                  fuzzy = TRUE)

save(sp_list_match, file = paste0(path_folder, "speciesname_confirm/sp_list_match.Rdata"))
load(paste0(path_folder, "speciesname_confirm/sp_list_match.Rdata"))

# exact matching without multiple matches
exact_match <- sp_list_match[which(sp_list_match$match_type %in% c("Exact (with author)", "Exact (without author)")), ]
exact_match_one <- exact_match[which(exact_match$multiple_matches == FALSE), ] 

# exact matching with multiple matches
exact_match_multiple <- exact_match[which(exact_match$multiple_matches == TRUE), ]
write.table(exact_match_multiple, file = paste0(path_folder, "speciesname_confirm/wcvp_exact_match_multiple_manual.txt"),
            row.names = FALSE, col.names = TRUE)

exact_match_multiple_manual <- read.csv(paste0(path_folder, "speciesname_confirm/wcvp_exact_match_multiple_manual.csv"))

# others need manually deal with
other_match <- sp_list_match[which(!sp_list_match$TaxonName %in% exact_match$TaxonName), ] 
write.table(other_match, file = paste0(path_folder, "speciesname_confirm/wcvp_other_match_manual.txt"),
            row.names = FALSE, col.names = TRUE)

other_match_manual <- read.csv(paste0(path_folder, "speciesname_confirm/wcvp_other_match_manual.csv"))

# aggregating the validating result
valid_data <- rbind(exact_match_one, exact_match_multiple_manual, other_match_manual)
write.table(valid_data, file = paste0(path_folder, "speciesname_confirm/tree_list_speciesname_wcvp.txt"),
            row.names = FALSE, col.names = TRUE)







