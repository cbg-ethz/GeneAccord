## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  echo = TRUE,
  collapse = TRUE,
  comment = "#>"
)

## ------------------------------------------------------------------------
library(GeneAccord)

## ------------------------------------------------------------------------
# defines the path to the external data directory
ext_data_dir <- system.file('extdata', package = 'GeneAccord')

# set the patient id
pat_id <- "01"

# find all csv-files from patient 01, which are in the clonal_genotypes subdirectory, with seeds 5, 10, 15, ..., 100
input_files_01 <- paste(ext_data_dir, "/clonal_genotypes/cloe_seed", seq(5, 100, by = 5), "/", pat_id, ".csv", sep = "")
# 'input_files_01' is now a vector with all 20 csv-files for this patient

## ------------------------------------------------------------------------
head(read.csv(input_files_01[1]))

## ------------------------------------------------------------------------
clone_tbl_01_all_trees <- create_tbl_tree_collection(input_files_01)


## ------------------------------------------------------------------------
clone_tbl_01_all_trees

## ------------------------------------------------------------------------
# print all 20 tree id's
unique(as.character(clone_tbl_01_all_trees$tree_id))

## ------------------------------------------------------------------------
rates_clon_excl_01 <- compute_rates_clon_excl(clone_tbl_01_all_trees)

## ------------------------------------------------------------------------
head(rates_clon_excl_01)

# average rate for patient 01 across trees
m_01 <- mean(rates_clon_excl_01)

## ------------------------------------------------------------------------
# compute the histogram of how often gene pairs occur in the trees, and how often they are clonally exclusive
hist_of_num_trees_clon_excl_01 <- get_hist_clon_excl(clone_tbl_01_all_trees)

## ------------------------------------------------------------------------
head(hist_of_num_trees_clon_excl_01[[1]])
head(hist_of_num_trees_clon_excl_01[[2]])

## ---- eval = F-----------------------------------------------------------
#  # set the patient identifiers
#  all_pat_ids <- c()
#  for(i in seq(1, 89)){
#    if(i %in% c(3, 34, 51, 53, 58, 65, 87)) {
#      next
#    }
#    all_pat_ids <- c(all_pat_ids, ifelse(i<10, paste0("0", i), as.character(i)))
#  }
#  
#  # this will be the vector with the average rates of clonal exclusivity of each patient
#  avg_rates_m <- c()
#  
#  # these will be the two lists with an entry for each patient with
#  # a histogram of how often pairs occur in the collection of trees
#  list_of_num_trees_all_pats <- list()
#  # a histogram of how ofthen pairs were clonally exclusive
#  list_of_clon_excl_all_pats <- list()
#  # this is the list of tibbles that will contain the gene-to-clone assignments from all patients from all trees
#  clone_tbl_all_pats_all_trees_list <- list()
#  
#  
#  # loop over all patients and prepare the input data
#  cnt <- 0
#  for (this_pat in all_pat_ids){
#    cnt <- cnt + 1
#  
#    # find all csv-files from current patient, which are in the clonal_genotypes subdirectory, with seeds 5, 10, 15, ..., 100
#    input_files_this_pat <- paste0(ext_data_dir, "/clonal_genotypes/cloe_seed", seq(5, 100, by = 5), "/", this_pat, ".csv")
#  
#    # here, there should be files from 20 different seeds
#    stopifnot(length(input_files_this_pat) == 20)
#  
#    # generate the tibble containing the gene-to-clone-assignments
#    clone_tbl_this_pat_all_trees <- create_tbl_tree_collection(input_files_this_pat)
#  
#    # compute the rates of clonal exclusivity
#    rates_clon_excl_this_pat <- compute_rates_clon_excl(clone_tbl_this_pat_all_trees)
#  
#    # average rate for current patient across trees
#    m_this_pat <- mean(rates_clon_excl_this_pat)
#  
#    # compute the histogram of how often gene pairs occur in the trees, and how often they are clonally exclusive
#    hist_of_num_trees_clon_excl_this_pat <- get_hist_clon_excl(clone_tbl_this_pat_all_trees)
#  
#    # save this information into named vectors/lists
#    avg_rates_m <- c(avg_rates_m, m_this_pat)
#    list_of_num_trees_all_pats[[cnt]] <- hist_of_num_trees_clon_excl_this_pat[[1]]
#    list_of_clon_excl_all_pats[[cnt]] <- hist_of_num_trees_clon_excl_this_pat[[2]]
#    names(avg_rates_m)[cnt] <- names(list_of_num_trees_all_pats)[cnt] <- names(list_of_clon_excl_all_pats)[cnt] <- this_pat
#    clone_tbl_all_pats_all_trees_list[[cnt]] <- clone_tbl_this_pat_all_trees
#  }
#  
#  # Now we aggregate the list of tibbles into one big tibble
#  clone_tbl_all_pats_all_trees <- do.call("rbind", clone_tbl_all_pats_all_trees_list)

## ------------------------------------------------------------------------
# load the data as generated in step 1 from all 82 patients
data("clone_tbl_all_pats_all_trees")
data("avg_rates_m")
data("list_of_num_trees_all_pats")
data("list_of_clon_excl_all_pats")

## ------------------------------------------------------------------------
# the tibble containing the gene-to-clone assignments from all patients and the whole collection of trees
print(clone_tbl_all_pats_all_trees)

## ------------------------------------------------------------------------
# the average rates of clonal exclusivity from each patient
head(avg_rates_m)

## ------------------------------------------------------------------------
# the histograms of how often pairs were mutated across the collection of trees - here printed for patient 02
list_of_num_trees_all_pats[[2]]

## ------------------------------------------------------------------------
# the histograms of how often pairs were clonally exclusive across the collection of trees - here printed for patient 02
list_of_clon_excl_all_pats[[2]]

## ------------------------------------------------------------------------
pairs_in_patients_hist(clone_tbl_all_pats_all_trees)

## ------------------------------------------------------------------------
num_pat_pair_max <- 3

## ------------------------------------------------------------------------
set.seed(1234)

ecdf_list <- generate_ecdf_test_stat(avg_rates_m, list_of_num_trees_all_pats, list_of_clon_excl_all_pats, num_pat_pair_max, num_pairs_sim = 1000)

## ------------------------------------------------------------------------
data("ecdf_list")

## ------------------------------------------------------------------------
# perform the clonal exclusivity test
res_pairs <- GeneAccord(clone_tbl_all_pats_all_trees, avg_rates_m, ecdf_list)

## ------------------------------------------------------------------------
suppressMessages(library(magrittr))
suppressMessages(library(dplyr))

# extract all pairs with an adjusted p-value less than 0.05
sig_pairs <- res_pairs %>% filter(qval < 0.05)

# print the top pairs
sig_pairs[order(sig_pairs$pval),]

## ---- eval = F-----------------------------------------------------------
#  # this will create a tibble with ensembl gene id's and corresponding hgnc gene symbols of the human genes
#  all_genes_tbl <- create_ensembl_gene_tbl_hg()

## ------------------------------------------------------------------------
data("all_genes_tbl")

## ------------------------------------------------------------------------
# get the hgnc symbols of the significant pairs
sig_pairs <- map_pairs_to_hgnc_symbols(sig_pairs, all_genes_tbl)

# print the top pairs - here just showing the hgnc gene symbols and not the ensembl gene id's 
sig_pairs[order(sig_pairs$pval),] %>% select(-entity_A, -entity_B, -num_patients)

## ------------------------------------------------------------------------
# find out in which patients they occur and are clonally exclusive
sig_pairs <- take_pairs_and_get_patients(clone_tbl_all_pats_all_trees, sig_pairs)

# print the top pairs
sig_pairs[order(sig_pairs$pval),] %>% select(-entity_A, -entity_B, -mle_delta, -test_statistic, -num_patients)

## ------------------------------------------------------------------------
# define output file
sig_pairs_tsv <- "GeneAccord_sig_pairs.tsv"

# save the results into the tsv-file
sig_pairs_table <- write_res_pairs_to_disk(sig_pairs, avg_rates_m, sig_pairs_tsv)

## ---- echo=FALSE, results='asis'-----------------------------------------
rownames(sig_pairs_table) <- NULL
knitr::kable(sig_pairs_table)

## ---- echo=F, message=F, results='hide'----------------------------------
if(file.exists(sig_pairs_tsv)){
  suppressMessages(file.remove(sig_pairs_tsv))
}

## ---- fig.height=5, fig.width=12-----------------------------------------
# plot the ECDF of the test statistic under the null
plot_ecdf_test_stat(ecdf_list)


## ------------------------------------------------------------------------
res_pairs %>% select(num_patients) %>% group_by(num_patients) %>% tally()

## ---- fig.height=10, fig.width=5-----------------------------------------
plot_rates_clon_excl(avg_rates_m, clone_tbl_all_pats_all_trees)

## ---- fig.width=8, fig.height=5------------------------------------------
# extract the patient id's in which the pairs are mutated
pairs_of_interest <- sig_pairs %>% filter(qval < 0.02)
pat_ids_of_interest <- unique(unlist(strsplit(as.character(pairs_of_interest$mutated_in), ";")))

# from the gene-to-clone tibble extract the data for one tree
this_tree_id <- 1
clone_tbl_all_pats_tree1 <- clone_tbl_all_pats_all_trees %>% 
                   filter(tree_id == this_tree_id, patient_id %in% pat_ids_of_interest) %>% 
                   select(-tree_id)

# plot the heatmap of clones and genes for the patients and genes of interest
heatmap_clones_gene_pat(pairs_of_interest, clone_tbl_all_pats_tree1, all_genes_tbl, first_clone_is_N = TRUE)


## ------------------------------------------------------------------------
data("ensg_reactome_path_map")

## ------------------------------------------------------------------------
clone_tbl_pat01_tree1 <- clone_tbl_all_pats_all_trees %>%
                            filter(patient_id == "01") %>%
                            filter(tree_id == 1) %>%
                            select(-tree_id)

clone_tbl_pat01_tree1_pw <- convert_ensembl_to_reactome_pw_tbl(clone_tbl_pat01_tree1, ensg_reactome_path_map)


clone_tbl_pat01_tree1_pw

## ------------------------------------------------------------------------
# this tells GeneAccord that both cases should be tested: 
# pairs that tend to co-occur together in the same clone more often than expected,
# and pairs that tend to be clonally exclusive.
alternative <- "two.sided"

# to limit the number of tests done, the genes that should be
# tested can be pre-specified
# here, the genes VHL, PXDN
genes_of_interest <- c("ENSG00000134086", "ENSG00000130508")

# this tells GeneAccord that we only want to test pairs where both gene A AND gene B
# are from the 'genes of interest', so we're only testing the pair {VHL, PXDN}
AND_OR <- "AND"

# run GeneAccord with these options
res_pairs_ts <- suppressMessages(GeneAccord(clone_tbl_all_pats_all_trees, avg_rates_m, ecdf_list, alternative, genes_of_interest, AND_OR))

res_pairs_ts

## ------------------------------------------------------------------------
AND_OR <- "OR"

# Then, to limit the number of tests performed, it is recommended to only specify a small
# number of genes
# here, only SART3
genes_of_interest <- c("ENSG00000075856")

# run GeneAccord testing all pairs that include the gene SART3
res_pairs_ts <- suppressMessages(GeneAccord(clone_tbl_all_pats_all_trees, avg_rates_m, ecdf_list, alternative, genes_of_interest, AND_OR))

res_pairs_ts

## ------------------------------------------------------------------------
sessionInfo()

