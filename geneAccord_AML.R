# load libraries
library(dplyr)
library(purrr)
library(tidyr)

# test if GeneAccord is installed
GeneAccordFlag <- FALSE
try(GeneAccordFlag <- packageVersion("GeneAccord") > '1.0.0', silent = TRUE)

if(GeneAccordFlag) {
  library(GeneAccord)
} else {
  source("./R/geneAccord_fns.R")
}


# read in data
clone_df <- read.table("./data/AMLsubclones_trees.txt", header = TRUE)
# process clone data frame into clonal exclusivity data frame
cx_df <- clone_df_to_cx_df(clone_df)
# select only gene pairs in more than 4 patients
cx_selected <- filter_cx_df(cx_df, 4)

# run the standard geneAccord testing
cx_results <- geneAccord(cx_selected)
cx_results

# run testing of patient selection
pat_results <- geneAccord(cx_selected, test_type = "patient")
pat_results

# run combined testing of patient selection and 
# location of gene pairs in those patients
comb_results <- geneAccord(cx_selected, test_type = "combined")
comb_results
