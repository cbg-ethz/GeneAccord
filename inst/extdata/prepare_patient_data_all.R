
args <- commandArgs(TRUE)

# defines the path to the external data directory
ext_data_dir <- args[1]

library(magrittr)

# source the functions necessary
packageDir <- paste0(ext_data_dir, "/../../")
for (file in list.files(paste0(packageDir, "/R"), full.names = TRUE)) {source(file)}

# define the output RData file
dataToSave_clone_tbl_all_pats_all_trees <- paste0(packageDir, "/data/clone_tbl_all_pats_all_trees.RData")
dataToSave_avg_rates_m <- paste0(packageDir, "/data/avg_rates_m.RData")
dataToSave_list_of_num_trees_all_pats <- paste0(packageDir, "/data/list_of_num_trees_all_pats.RData")
dataToSave_list_of_clon_excl_all_pats <- paste0(packageDir, "/data/list_of_clon_excl_all_pats.RData")
dataToSave_ecdf_list <- paste0(packageDir, "/data/ecdf_list.RData")

# set the patient ids
all_pat_ids <- c()
for(i in seq(1,89)){
  # these patients did not have more than one gene assigned to a clone, 
  # i.e. no GeneAccord analysis possible
  if(i %in% c(3, 34, 51, 53, 58, 65, 87)){
    next
  } 
  all_pat_ids <- c(all_pat_ids, ifelse(i<10, paste0("0", i), as.character(i)))
}

# this will be the vector with the average rates of clonal exclusivity of each patient
avg_rates_m <- c()

# these will be the two lists with an entry for each patient with
# a histogram of how often pairs occur in the collection of trees
list_of_num_trees_all_pats <- list()
# a histogram of how ofthen pairs were clonally exclusive
list_of_clon_excl_all_pats <- list()
# this is the list of tibbles that will contain the gene-to-clone assignments from all patients from all trees
clone_tbl_all_pats_all_trees_list <- list()


# loop over all patients and prepare the input data
cnt <- 0
for (this_pat in all_pat_ids){
  cnt <- cnt + 1

  # find all csv-files from current patient, which are in the clonal_genotypes subdirectory, with seeds 5, 10, 15, ..., 100
  input_files_this_pat <- paste0(ext_data_dir, "/clonal_genotypes/cloe_seed", seq(5, 100, by = 5), "/", this_pat, ".csv")

  # there should be files from 20 different seeds
  stopifnot(length(input_files_this_pat) == 20)
  
  clone_tbl_this_pat_all_trees <- create_tbl_tree_collection(input_files_this_pat)

  rates_clon_excl_this_pat <- compute_rates_clon_excl(clone_tbl_this_pat_all_trees)

  # average rate for current patient across trees
  m_this_pat <- mean(rates_clon_excl_this_pat)

  # compute the histogram of how often gene pairs occur in the trees, and how often they are clonally exclusive
  hist_of_num_trees_clon_excl_this_pat <- get_hist_clon_excl(clone_tbl_this_pat_all_trees)

  avg_rates_m <- c(avg_rates_m, m_this_pat)
  list_of_num_trees_all_pats[[cnt]] <- hist_of_num_trees_clon_excl_this_pat[[1]]
  list_of_clon_excl_all_pats[[cnt]] <- hist_of_num_trees_clon_excl_this_pat[[2]]
  names(avg_rates_m)[cnt] <- names(list_of_num_trees_all_pats)[cnt] <- names(list_of_clon_excl_all_pats)[cnt] <- this_pat
  clone_tbl_all_pats_all_trees_list[[cnt]] <- clone_tbl_this_pat_all_trees 
}

# Now we aggregate the list of tibbles into one big tibble
clone_tbl_all_pats_all_trees <- do.call("rbind", clone_tbl_all_pats_all_trees_list)


# save the data 
save(clone_tbl_all_pats_all_trees, file = dataToSave_clone_tbl_all_pats_all_trees)
save(avg_rates_m, file = dataToSave_avg_rates_m)
save(list_of_num_trees_all_pats, file = dataToSave_list_of_num_trees_all_pats)
save(list_of_clon_excl_all_pats , file = dataToSave_list_of_clon_excl_all_pats)


# generate the ecdf
set.seed(1234)
load(dataToSave_avg_rates_m)
load(dataToSave_list_of_num_trees_all_pats)
load(dataToSave_list_of_clon_excl_all_pats)
num_pat_pair_max <- 3
ecdf_list <- generate_ecdf_test_stat(avg_rates_m, list_of_num_trees_all_pats, list_of_clon_excl_all_pats, num_pat_pair_max, num_pairs_sim = 100000)
save(ecdf_list, file = dataToSave_ecdf_list)

