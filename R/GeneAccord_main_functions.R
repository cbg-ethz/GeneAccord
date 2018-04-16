#' Method to detect clonally exclusive gene or pathway pairs in a cohort of cancer patients
#' 
#' After running a tool such as \code{Cloe} that identifies clones in a tumor and infers the phylogenetic history,
#' the user has for each tumor a list of alterations and their clone assignments. Since the tree inference includes
#' uncertainty, it may be run several times. Given a tibble containing the information of which genes/pathways are mutated
#' in which patient and clone and from which tree, this function systematically tests the data for significant
#' clonal exclusivities. That is, it checks for each gene/pathway pair whether the number of clonal exclusivities
#' is significantly different from what would be expected by chance. Such a tibble can be generated with
#' \code{\link{create_tbl_tree_collection}}, and then adding the additional column 'tree_id' to indicate which tree of the tree
#' inference was used. For instance, if the tree inference tool was run several times using different seeds, the column
#' 'tree_id' may contain the seed of the respective tree. Hence, the tibble is expected to have the columns 'file_name',
#' 'patient_id', 'altered_entity', 'clone1', 'clone2', ... up to the maximal number of clones (Default: until 'clone7'), and
#' 'tree_id'. Note that the labelling of the clones does not matter and only needs to stay fixed within each patient and
#' tree inference. There is also the option to test two-sided, meaning that also pairs will be tested that tend to occur more
#' often together in the same clones. Hence it also allows to detect significant clonal co-occurrence. An additional option is
#' to test only a specific subset of genes.
#'
#' @title Detection of clonally exclusive gene or pathway pairs in a cohort of cancer patients
#' @param clone_tbl The tibble containing the information of which gene/pathway is mutated in which
#' clone from which patient and in which tree from the collection of trees. Can be generated with \code{\link{create_tbl_tree_collection}}
#' for each patient separately and then appended.
#' @param avg_rates_m The average rates of clonal exclusivity for each patient as computed with
#' \code{\link{compute_rates_clon_excl}}. The name of each rate is the respective patient id. The rates are assumed to be the average
#' over all tree inferences from a patient.
#' @param ecdf_list The list of ECDF's of the test statistic under the null distribution. Can be generated with
#' \code{\link{generate_ecdf_test_stat}}.
#' @param alternative The character indicating whether pairs should only be tested if delta > 0 or if all pairs should be tested.
#' Can be one of "greater" or "two.sided". Default: "greater".
#' @param genes_of_interest A character vector of genes to test for clonal exclusivity. The genes have to be in the same identifier as 
#' the one in the tibble. Per default, all genes are tested. Default: "ALL".
#' @param AND_OR If \code{genes_of_interest} is specified, this indicator tells whether to test only pairs within the \code{genes_of_interest} (AND),
#' or whether all pairs involving at least one of these genes should be tested (OR). I.e. can be one of "AND", "OR". Default: "OR". If 
#' \code{genes_of_interest} is "ALL", then all gene pairs will be tested and this parameter is ignored.
#' @return A tibble containing the test result for each pair of mutated genes/pathways that was tested. More precisely, it contains the columns
#' 'entity_A', 'entity_B', 'num_patients', 'pval', 'mle_delta', 'test_statistic', and 'qval'.
#' Each row is then a gene or pathway pair which is specified with 'entity_A', and 'entity_B'. Note that the test is symmetric,
#' hence switching the labels A and B does not change the results. The column 'num_patients' contains the information in how
#' many patients both of the genes/pathways were mutated and hence how many patients' rates were used for the test. The 'pval' is 
#' the p-value of the clonal exclusivity test. The 'mlde_delta' is the maximum likelihood estimate of the delta for the elevated clonal exclusivity 
#' rate in the alternative model. The column 'test_statistic' is the likelihood ratio test statistic. The 'qval' is the adjusted p-value after multiple testing
#' correction with Benjamini-Hochberg.
#' @author Ariane L. Moore, \email{ariane.moore@@bsse.ethz.ch} 
#' @import 
#' caTools
#' dplyr
#' tibble
#' @examples
#' clone_tbl <- tibble::as_tibble(cbind("file_name" =
#'                                                rep(c(rep(c("fn1", "fn2"), each = 3)), 2),
#'                                    "patient_id" = rep(c(rep(c("pat1", "pat2"), each = 3)), 2),
#'                                    "altered_entity" = c(rep(c("geneA", "geneB", "geneC"), 4)),
#'                                    "clone1" = c(0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0),
#'                                    "clone2" = c(1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1),
#'                                    "tree_id" = c(rep(5, 6), rep(10, 6))))
#' clone_tbl_pat1 <- dplyr::filter(clone_tbl, patient_id == "pat1")
#' clone_tbl_pat2 <- dplyr::filter(clone_tbl, patient_id == "pat2")
#' rates_exmpl_1 <- compute_rates_clon_excl(clone_tbl_pat1)
#' rates_exmpl_2 <- compute_rates_clon_excl(clone_tbl_pat2)
#' avg_rates_m <- apply(cbind(rates_exmpl_1, rates_exmpl_2), 2, mean)
#' names(avg_rates_m) <- c(names(rates_exmpl_1)[1], names(rates_exmpl_2)[1])
#' values_clon_excl_num_trees_pat1 <- get_hist_clon_excl(clone_tbl_pat1)
#' values_clon_excl_num_trees_pat2 <- get_hist_clon_excl(clone_tbl_pat2)
#' list_of_num_trees_all_pats <- list(pat1=values_clon_excl_num_trees_pat1[[1]], 
#'                                    pat2=values_clon_excl_num_trees_pat2[[1]])
#' list_of_clon_excl_all_pats <- list(pat1=values_clon_excl_num_trees_pat1[[2]],
#'                                    pat2=values_clon_excl_num_trees_pat2[[2]])
#' num_pat_pair_max <- 2
#' num_pairs_sim <- 10
#' ecdf_list <- generate_ecdf_test_stat(avg_rates_m, list_of_num_trees_all_pats, 
#'                              list_of_clon_excl_all_pats, num_pat_pair_max, num_pairs_sim)
#' alternative <- "greater"
#' GeneAccord(clone_tbl, avg_rates_m, ecdf_list, alternative)
#' alternative <- "two.sided"
#' GeneAccord(clone_tbl, avg_rates_m, ecdf_list, alternative)
#' genes_of_interest <- c("geneB", "geneC")
#' GeneAccord(clone_tbl, avg_rates_m, ecdf_list, alternative, genes_of_interest)
#' AND_OR = "AND"
#' GeneAccord(clone_tbl, avg_rates_m, ecdf_list, alternative, genes_of_interest, AND_OR)
GeneAccord <- function(clone_tbl, avg_rates_m, ecdf_list, 
                       alternative = "greater", genes_of_interest = "ALL", AND_OR = "OR"){
  patient_id <- altered_entity <- file_name <- tree_id <- qval <- pval <- mle_delta <- num_patients <- pairs <- n <- NULL
  stopifnot(dplyr::is.tbl(clone_tbl))
  stopifnot("file_name" %in% colnames(clone_tbl))
  stopifnot("patient_id" %in% colnames(clone_tbl))
  stopifnot("altered_entity" %in% colnames(clone_tbl))
  stopifnot("clone1" %in% colnames(clone_tbl))
  stopifnot("clone2" %in% colnames(clone_tbl))
  stopifnot("tree_id" %in% colnames(clone_tbl))
  stopifnot(is.numeric(avg_rates_m))
  stopifnot(is.list(ecdf_list))
  stopifnot(is.character(alternative))
  stopifnot(alternative == "greater" || alternative == "two.sided")
  stopifnot(is.character(genes_of_interest))
  stopifnot(is.character(AND_OR))
  stopifnot(AND_OR == "AND" || AND_OR == "OR") 
  
  # message to user
  num_entries <- dim(clone_tbl)[1]
  all_tree_ids <- unique(as.character(clone_tbl$tree_id))
  num_trees <- length(all_tree_ids)
  message(paste0("There were ", num_trees, " different tree inferences per patient.\n",
                 num_entries, " alterations in total found in the provided tibble."))
  num_genes_of_interest <- length(genes_of_interest)
  spec_genes <- FALSE # should only a specific set of genes be tested? 
  if(num_genes_of_interest == 1 && genes_of_interest != "ALL"){
    message(paste("There is one gene of interest provided, namely ", genes_of_interest, ".", sep = ""))
    if(AND_OR == "AND"){
      message(paste0("Since there is only one gene of interest provided, the AND_OR parameter has to be 'OR',",
                     " otherwise no gene pairs are left"))
      AND_OR <- "OR"
    }
    spec_genes <- TRUE
  } else if (num_genes_of_interest > 1){
    message(paste("There are ", num_genes_of_interest, " genes of interest provided.", sep = ""))
    spec_genes <- TRUE
  }
  
  # get all genes/pathways, here referred to as entity (ent)
  all_ents <- unique(as.character(clone_tbl$altered_entity))
  num_ents <- length(all_ents)
  
  # get all patients
  all_pats <- unique(as.character(clone_tbl$patient_id))
  num_pats <- length(all_pats)
  stopifnot(length(avg_rates_m) == num_pats)
  stopifnot(length(setdiff(names(avg_rates_m), all_pats)) == 0)
  
  # get the pairs from each patient separately, and then check which pairs re-occur across at least two patients
  list_all_pairs_all_pats <- lapply(all_pats, function(x){
    clone_tbl_this_pat <- clone_tbl %>% filter(patient_id == x)
    all_ents_this_pat <- unique(as.character(clone_tbl_this_pat$altered_entity))
    num_ents_this_pat <- length(all_ents_this_pat)
    num_pairs_this_pat <- num_ents_this_pat*(num_ents_this_pat-1)/2
    # get all pairs of entities of that patient efficiently
    all_pairs_this_pat <- caTools::combs(all_ents_this_pat, 2)
    # to make sure the gene pairs are all in the same order such that
    # later, when we check which pairs occur in multiple patients,
    # that we do not miss pairs just because they are in a different
    # order
    all_pairs_this_pat <- t(apply(all_pairs_this_pat, 1, function(pair){sort(pair)}))
    stopifnot(num_pairs_this_pat == dim(all_pairs_this_pat)[1])
    return(all_pairs_this_pat)
  })
  all_pairs_all_pats <- do.call("rbind", list_all_pairs_all_pats)
  unique_all_pairs_all_pats <- unique(all_pairs_all_pats)
  num_pairs <- dim(unique_all_pairs_all_pats)[1]
  all_pairs <- unique_all_pairs_all_pats
  
  
  # already try to filter such that only those pairs remain which are mutated in two or more
  min_num_pat <- 2
  all_pairs_all_pats_one_string <- apply(all_pairs_all_pats, 1, function(x){paste0(x, collapse = "__")})
  all_pairs_all_pats_one_string_minPat <- tibble::tibble(pairs = all_pairs_all_pats_one_string) %>% 
                                               dplyr::group_by(pairs) %>% 
                                               dplyr::tally() %>% 
                                               dplyr::filter(n>=min_num_pat)
  all_pairs <- t(vapply(all_pairs_all_pats_one_string_minPat$pairs, function(x){unlist(strsplit(x, "__"))}, character(2)))
  rownames(all_pairs) <- NULL
  num_pairs <- dim(all_pairs)[1]
  
  message(paste0("Found in total ", num_ents, " different mutated genes/pathways in the provided tibble.\n",
                 "Found ", num_pats, " different patient id's in the provided tibble.\n",
                 "A total of", num_pairs, " pairs of genes/pathways is mutated in at least one patient.\n",
                 "And ", num_pairs, " pairs of genes/pathways are mutated in at least ", min_num_pat," patients."))
  
  # as an additional safety mechanism, we check what is the maximum n that pairs are mutated in how many patients
  max_num_pat_pairs <- max(as.numeric(as.character(all_pairs_all_pats_one_string_minPat$n)))
  if(length(ecdf_list) < max_num_pat_pairs){
    message(paste0("The maximum number of patients that a pair is mutated in is ", max_num_pat_pairs,
                   ", which is higher than the number of ECDF's provided. This may be problematic.\n",
                   "If you have more time, please re-run the function 'generate_ecdf_test_stat' with ",
                   "num_pat_pair_max = ", max_num_pat_pairs, ".\nOtherwise, you can also wait and see if ",
                   "one of the pairs that will be tested is actually mutated in so many patients.\n",
                   "Maybe all pairs tested will be mutated in less patients, and you're fine.\n"))
  }
  
  ## get all pairs of entities efficiently
  #num_pairs <- num_ents*(num_ents-1)/2
  #all_pairs <- caTools::combs(all_ents, 2)
  #stopifnot(num_pairs == dim(all_pairs)[1])
  
  # distinguish between the case where all genes should be tested or where only a specific subset should be tested
  pairs_to_test <- c()
  if(spec_genes){ # if TRUE, only test genes of interest
    # sanity check: make sure that the genes of interest are among the genes in the tibble
    for (this_gene in genes_of_interest){
      if(! this_gene %in% all_ents){
        warning(paste0("The gene " , this_gene, " is among the genes of interest, but is not among the genes in the clone tibble.",
                       " Are you sure that there is no typo in this gene?"))
      }
    }
    # now filter the pairs according to the genes of interest and the AND_OR
    if(AND_OR == "AND"){
      message(paste("Only pairs will be tested, where both genes are in the set of genes of interest.", sep = ""))
      pairs_to_test <- apply(all_pairs, 1, function(x){
        this_ent_A <- x[1]
        this_ent_B <- x[2]
        if(this_ent_A %in% genes_of_interest && this_ent_B %in% genes_of_interest){
          return(x)
        }
      })
    } else if(AND_OR == "OR"){
      message(paste("Only pairs will be tested, where at least one of the genes is in the set of genes of interest.", sep = ""))
      pairs_to_test <- apply(all_pairs, 1, function(x){
        this_ent_A <- x[1]
        this_ent_B <- x[2]
        if(this_ent_A %in% genes_of_interest || this_ent_B %in% genes_of_interest){
          return(x)
        }
      })
    }
  } else {
    pairs_to_test <- all_pairs
  }
  
  # remove the NULL entries
  if(is.null(pairs_to_test)){
    message("There are no pairs that fulfill the restriction with genes_of_interest and the AND_OR.")
    message("No test will be done.")
    res_all_pairs <- tibble::tibble(entity_A = character(), 
                                    entity_B = character(),
                                    num_patients = numeric(),
                                    pval = numeric(),
                                    mle_delta = numeric(),
                                    test_statistic = numeric(),
                                    qval = numeric())
    
    return(res_all_pairs)
  } else {
    # num pairs to test
    if(spec_genes){ # if only a subset of all pairs is being tested
      pairs_to_test <- matrix(unlist(pairs_to_test[unlist(lapply(pairs_to_test, function(x){!is.null(x)}))]), ncol = 2, byrow = TRUE)
    }
    num_pairs_to_test <- dim(pairs_to_test)[1]
    message(paste("", num_pairs_to_test, " gene/pathway pairs will be processed.", sep = ""))
    
    # then use apply on the pairs with a function that takes both entities, the clone tbl, and
    # the rates of average clonal exclusivity
    # the function invokes the clonal exclusivity test for the pair
    # returns a list of lists: list(p_val = p_val, num_patients = num_shared_pats, mle_delta = mle_delta,
    # test_statistic = test_stat) for each pair
    all_tests <- list()
    
    all_tests <- apply(pairs_to_test, 1, function(x){
      suppressMessages(ecdf_lr_test_clon_excl_avg_rate(x[1], x[2], clone_tbl, avg_rates_m, ecdf_list, alternative))})
    these_deltas <- unlist(lapply(all_tests, function(x){x[[3]]}))
    these_test_stats <- unlist(lapply(all_tests, function(x){x[[4]]}))
    stopifnot(is.list(all_tests))
    
    # extract the p-values
    all_pairs_pval <- unlist(lapply(all_tests, function(x){x[[1]]}))
    stopifnot(length(all_pairs_pval) == num_pairs_to_test)
    
    # extract num_patients, deltas, and test statistic
    all_pairs_num_pat <- unlist(lapply(all_tests, function(x){x[[2]]}))
    stopifnot(length(all_pairs_num_pat) == num_pairs_to_test)
    stopifnot(length(these_deltas) == num_pairs_to_test)
    stopifnot(length(these_test_stats) == num_pairs_to_test)

    # create final tibble to return
    # contains the columns 'entity_A', 'entity_B', 'num_patients', 'pval', 'mle_delta', 'test_statistic', 'qval'
    res_all_pairs <- dplyr::as.tbl(data.frame(entity_A = pairs_to_test[,1], 
                                              entity_B = pairs_to_test[,2],
                                              num_patients = all_pairs_num_pat,
                                              pval = all_pairs_pval,
                                              mle_delta = these_deltas,
                                              test_statistic = these_test_stats))
    
    # now filter out the pairs where no test was performed
    num_pairs_in_zero_pats <- dim(res_all_pairs %>% dplyr::filter(num_patients == 0))[1]
    num_pairs_in_one_pat <- dim(res_all_pairs %>% dplyr::filter(num_patients == 1))[1]
    #message(paste("There were ", num_pairs_in_zero_pats, " pairs that were mutated together in zero patients.", sep = ""))
    #message(paste("And ", num_pairs_in_one_pat, " pairs that were mutated together in one patient.", sep = ""))
    res_all_pairs_minNumPat2 <- res_all_pairs %>% dplyr::filter(num_patients >= 2)
    num_pairs_in_two_or_more_pats <- dim(res_all_pairs_minNumPat2)[1]
    # convert the delta to numeric such that it can be filtered according to delta
    res_all_pairs_minNumPat2$mle_delta <- as.numeric(as.character(res_all_pairs_minNumPat2$mle_delta))
    
    if (alternative == "greater"){
      # filter out the pairs where delta <= 0, because they were not tested
      final_res_pairs <- res_all_pairs_minNumPat2 %>% dplyr::filter(mle_delta > 0)
      num_pairs_deltaGT0 <- dim(final_res_pairs)[1]
      message(paste0(num_pairs_in_two_or_more_pats ," pairs were mutated in more than one patient and for them the delta was computed.\n",
                     "There were ", num_pairs_deltaGT0, " pairs where delta is > 0 and for which the p-value was computed."))
    } else {
      message(paste(num_pairs_in_two_or_more_pats ,
                    " pairs were mutated in more than one patient and for them the delta and the p-value was computed.", sep = ""))
      final_res_pairs <- res_all_pairs_minNumPat2
    }
    
    # perform correction for multiple testing
    all_pvals <- as.numeric(as.vector(as.data.frame(final_res_pairs %>% dplyr::select(pval))[,1]))
    qvals_BH <- p.adjust(all_pvals, method = "BH")
    stopifnot(length(qvals_BH) == dim(final_res_pairs)[1])
    final_res_pairs_pval_corrected <- final_res_pairs %>% dplyr::mutate(qval = qvals_BH)
    
    return(final_res_pairs_pval_corrected)
  }
}


