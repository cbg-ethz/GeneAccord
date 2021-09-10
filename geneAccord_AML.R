# load libraries
library(dplyr)
library(purrr)
library(tidyr)

# test if GeneAccord is installed
GeneAccordFlag <- FALSE
#try(GeneAccordFlag <- packageVersion("GeneAccord") > '1.0.0', silent = TRUE)

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

# run the placement geneAccord testing
(place_results <- geneAccord(cx_selected))

# run occurrence geneAccord testing
(occur_results <- geneAccord(cx_selected, test_type = "occurrence"))

# run combined testing of occurrence in patients and 
# location of gene pairs in those patients
(comb_results <- geneAccord(cx_selected, test_type = "combined"))
