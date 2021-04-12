
# test if GeneAccord is installed
GeneAccordFlag <- FALSE
try(GeneAccordFlag <- packageVersion("GeneAccord") > 0, silent = TRUE)

if(GeneAccordFlag) {
  library(GeneAccord)
} else {
  library(dplyr)
  library(purrr)
  library(tidyr)
  source("./R/geneAccord_fns.R")
}

# read in data
clone_df <- read.table("./data/AMLsubclones_trees.txt", header = TRUE)
# process clone data frame into clonal exclusivity data frame
cx_df <- clone_df_to_cx_df(clone_df)
# select only gene pairs in more than 4 patients
cx_selected <- filter_cx_df(cx_df, 4)
# run the geneAccord testing
cx_results <- geneAccord(cx_selected)
# compare to the chi-squared asymptotics
#cx_results$p_chisq <- pchisq(cx_results$LR, df = 1, lower.tail = FALSE)
#cx_results$q_chisq <- p.adjust(cx_results$p_chisq, method = "BH")

cx_results


