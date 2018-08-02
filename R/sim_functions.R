#' Generate the ECDF of the test statistic under the null distribution - 
#' taking the average rates of clonal exclusivity, as well as
#' sampling from the real data for each patient, in how many trees a pair 
#' occurs and is clonally excl.
#'
#' This function takes the computed average rates of clonal exclusivity from 
#' the data (m1, ... mN), which are specific to each
#' patient and averaged over several trees from the collection of tree 
#' inferences. It also takes the histogram for each patient, of the values of 
#' how often a pair was clonally exclusive over the number of trees it was 
#' mutated in. It then simulates the test statistic under the
#' null for each number of patients a pair is be mutated in from 2, 3, ... 
#' 'num_pat_pair_max'. Afterwards, it generates the empirical cumulative
#' distribution function (ECDF) using the \code{ecdf} function of the stats 
#' package and returns the list with the ECDF's for the
#' number of patients n=2, 3, ..., N. This step is necessary for each new 
#' data set before the clonal exclusivity test can be
#' done. In the clonal exclusivity test, the observed test statistics are 
#' compared to the ECDF.
#'
#' @title Generate the ECDF of the test statistic under the null distribution
#'  - taking the average rates of clonal exclusivity
#' @param avg_rates_m The average rates of clonal exclusivity from all the 
#' patients in the cohort, and averaged over several trees 
#' from the collection of tree inferences.
#' @param list_of_num_trees_all_pats A named list that contains an entry for 
#' each patient which is the vector with the values of the 
#' information from each pair in a patient of how often it was mutated across 
#' trees. The patient odering in the list has to be the 
#' same as in \code{avg_rates_m}.
#' @param list_of_clon_excl_all_pats A named list with an entry for each 
#' patient that is a vector with the values of in how many trees 
#' a pair was clonally exclusive. The patient ordering in the list has to 
#' be the same as in \code{avg_rates_m}.
#' @param num_pat_pair_max The maximum number of patients a pair is 
#' mutated in.
#' @param num_pairs_sim The number of simulated gene/pathway pairs to be 
#' generated, i.e. the number of times the test statistic
#' is computed. Recommended to choose a big number, e.g. 100000.
#' @param beta_distortion The value \code{M=alpha + beta} for the beta 
#' distribution, with which the average rates will be distorted. The bigger 
#' the M the higher the distribution is peaked around the actual rate. 
#' Therefore, the lesser the M, the more distorted the rates will be. 
#' Default: 1000.
#' @author Ariane L. Moore
#' @return The return value is a list with ECDF's. The first list entry 
#' is just set to NULL for technical reasons.
#' @export
#' @importFrom magrittr "%>%"
#' @importFrom dplyr tibble filter is.tbl pull
#' @importFrom stats ecdf
#' @examples
#' clone_tbl <- dplyr::tibble("file_name" =
#'    rep(c(rep(c("fn1", "fn2"), each=3)), 2),
#'    "patient_id"=rep(c(rep(c("pat1", "pat2"), each=3)), 2),
#'    "altered_entity"=c(rep(c("geneA", "geneB", "geneC"), 4)),
#'    "clone1"=c(0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0),
#'    "clone2"=c(1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1),
#'    "tree_id"=c(rep(5, 6), rep(10, 6)))
#' clone_tbl_pat1 <- dplyr::filter(clone_tbl, patient_id == "pat1")
#' clone_tbl_pat2 <- dplyr::filter(clone_tbl, patient_id == "pat2")
#' rates_exmpl_1 <- compute_rates_clon_excl(clone_tbl_pat1)
#' rates_exmpl_2 <- compute_rates_clon_excl(clone_tbl_pat2)
#' avg_rates_m <- apply(cbind(rates_exmpl_1, rates_exmpl_2), 2, mean)
#' names(avg_rates_m) <- c(names(rates_exmpl_1)[1], names(rates_exmpl_2)[1])
#' values_clon_excl_num_trees_pat1 <- get_hist_clon_excl(clone_tbl_pat1)
#' values_clon_excl_num_trees_pat2 <- get_hist_clon_excl(clone_tbl_pat2)
#' list_of_num_trees_all_pats <-
#'     list(pat1=values_clon_excl_num_trees_pat1[[1]], 
#'     pat2=values_clon_excl_num_trees_pat2[[1]])
#' list_of_clon_excl_all_pats <-
#'     list(pat1=values_clon_excl_num_trees_pat1[[2]],
#'     pat2=values_clon_excl_num_trees_pat2[[2]])
#' num_pat_pair_max <- 2
#' num_pairs_sim <- 10
#' ecdf_list <- generate_ecdf_test_stat(avg_rates_m, 
#'                 list_of_num_trees_all_pats, list_of_clon_excl_all_pats,
#'                 num_pat_pair_max, num_pairs_sim)
#' plot(ecdf_list[[2]])
generate_ecdf_test_stat <- function(avg_rates_m, list_of_num_trees_all_pats, 
                                    list_of_clon_excl_all_pats, 
                                    num_pat_pair_max, num_pairs_sim, 
                                    beta_distortion=1000){
    test_statistic <- NULL
    stopifnot(is.numeric(avg_rates_m))
    stopifnot(is.list(list_of_num_trees_all_pats))
    stopifnot(is.list(list_of_clon_excl_all_pats))
    stopifnot(is.numeric(num_pat_pair_max))
    stopifnot(is.numeric(num_pairs_sim))
    for(this_rate_m in avg_rates_m){
        stopifnot(this_rate_m <= 1 && this_rate_m >= 0)
    }
    num_pats_total=length(avg_rates_m)
    list_of_clon_excl_frac_trees_all_pats <- 
        list(list_of_num_trees_all_pats, list_of_clon_excl_all_pats)
    stopifnot(length(list_of_clon_excl_frac_trees_all_pats) == 2)
    stopifnot(length(list_of_clon_excl_frac_trees_all_pats[[1]]) == 
        num_pats_total)
    stopifnot(length(list_of_clon_excl_frac_trees_all_pats[[2]]) == 
        num_pats_total)
    stopifnot(is.numeric(beta_distortion))
    if(beta_distortion < 50){
        warning("A beta_distortion value of ", beta_distortion, 
        " is not recommended as it will distort the rates a lot.",
        " Please consider choosing a higher value, e.g. 100,",
        " or 1000.")
    }
    num_pat_pair_max <- min(c(num_pat_pair_max, num_pats_total)) 
    ## the number of patients that the simulated pairs are mutated in 
    ## are 1, 2, ..., num_pat_pair_max
    
    ## also make sure that the rates and the lists with histogram are in 
    ## the same patient order:
    for (i in seq_len(num_pats_total) ){
        if((names(avg_rates_m)[i] != 
            names(list_of_clon_excl_frac_trees_all_pats[[1]])[i]) ||
            (names(avg_rates_m)[i] != 
            names(list_of_clon_excl_frac_trees_all_pats[[2]])[i]) ){
            stop("[Function: generate_ecdf_test_stat]: The rates seem to",
            " be\nin a different patient order than the ",
            "lists with patient's histograms\nlist_of_num_trees_all_pats,",
            " list_of_clon_excl_all_pats. See for example\nthe entry",
            " number ",
            i , " of the rates.\nMake sure that patients rates and the",
            " lists of histogram are ",
            "in the same patient order!")
        }
    }
    
    ## message to user
    message("There were ", num_pats_total, " rates of ",
    "average clonal exclusivity provided.\n",
    "The average rates will be distorted with a beta",
    " distribution and M=alpha+beta=", beta_distortion,
    ".\n", 
    "Additionally, from the ", num_pats_total, 
    " patients, the values from the data\nfor all pairs of",
    " # clon. excl./ # occurrence across all trees are",
    " provided.\n",
    "From those, it will be sampled ", num_pairs_sim, 
    " times to obtain the test statistic under the null.\n",
    "This procedure is done for num_patients=2,...",
    num_pat_pair_max , ".")
    ecdf_lrtest_stat <- list()
    ecdf_lrtest_stat[[1]] <- NULL
    
    ## loop over the number of patients 2, 3, ..., num_pat_pair_max
    for (i in seq(2, num_pat_pair_max)){
        message("Generate ECDF for the case where a gene pair",
        " is mutated together in ", i, " patient(s).")
      
        ## generate samples of the test statistic under null for i patients
        this_test_stat_hist <- 
            suppressMessages(build_null_test_statistic(avg_rates_m, 
            list_of_clon_excl_frac_trees_all_pats, 
            num_pat_pair=i,
            num_pairs_sim=num_pairs_sim,
            beta_distortion=beta_distortion))
        stopifnot(dim(this_test_stat_hist)[2] == (2+i*4))
      
        stopifnot(is.tbl(this_test_stat_hist))
        stopifnot("test_statistic" %in% colnames(this_test_stat_hist))
        stopifnot(dim(this_test_stat_hist)[1] == num_pairs_sim)
        
        ## build the ecdf
        ecdf_lrtest_stat[[i]] <- 
            ecdf(this_test_stat_hist %>% pull(test_statistic))
    }
    
    return(ecdf_lrtest_stat)
}



#' Generate the values of the test statistic under the null, and also 
#' p-values of the clonal exclusivity test under the null.
#' Taking the average rates of clonal exclusivity, as well as sampling 
#' from the real data for each patient, in how many trees a 
#' pair occurs and is clonally exclusive.
#'
#' This function takes the computed average rates of clonal exclusivity 
#' from the data (m1, ... mN), which are specific to each
#' patient and averaged over several trees from the collection of tree 
#' inferences. It also takes the histogram for each patient, 
#' of the values of how often a pair was clonally exclusive over the number 
#' of trees it was mutated in. It also takes the empirical 
#' cumulative distribution function (ECDF) which was generated with 
#' \code{\link{generate_ecdf_test_stat}}. It then computes the p-value of
#' the simulated pairs under the null.
#'
#' @title Generate the test statistic and p-values under the null 
#' distribution
#' @param avg_rates_m The average rates of clonal exclusivity from all 
#' the patients in the cohort, and averaged over several trees 
#' from the collection of tree inferences.
#' @param list_of_num_trees_all_pats A named list that contains an entry 
#' for each patient which is the vector with the values of the 
#' information from each pair in a patient of how often it was mutated 
#' across trees. The patient odering in the list has to be the 
#' same as in \code{avg_rates_m}.
#' @param list_of_clon_excl_all_pats A named list with an entry for each 
#' patient that is a vector with the values of in how many trees 
#' a pair was clonally exclusive. The patient ordering in the list has to 
#' be the same as in \code{avg_rates_m}.
#' @param ecdf_list The list with ECDF's as generated with 
#' \code{\link{generate_ecdf_test_stat}}. 
#' @param num_pat_pair_max The maximum number of patients a pair is 
#' mutated in.
#' @param num_pairs_sim The number of simulated gene/pathway pairs to be 
#' generated, i.e. the number of times the test statistic
#' is computed.
#' @param beta_distortion The value \code{M=alpha + beta} for the beta 
#' distribution, with which the average rates will be distorted. The bigger 
#' the M
#' the higher the distribution is peaked around the actual rate. Therefore, 
#' the lesser the M, the more distorted the rates will be. Default: 1000.
#' @author Ariane L. Moore
#' @return The return value is a list of tibbles with a tibble for each 
#' number of patients, a pair can be mutated in. Each tibble contains the 
#' columns 'test_statistic', 'mle_delta', and then num_pat_pair columns 
#' of the rates of each patient 'pat1', 'pat2', ...; as well as 
#' num_pat_pair columns with the information about each patient, in how 
#' many trees the pair was occurring and in how many trees the pair was 
#' clonally exclusive. The tibble also contains a column 'pval' with the
#' p-value of the simulated pair. The list of tibbles is of length 
#' min{num_pat_pair_max, length(avg_rates_m)}.
#' @export
#' @importFrom magrittr "%>%"
#' @importFrom dplyr tibble filter is.tbl pull mutate
#' @examples
#' clone_tbl <- dplyr::tibble("file_name" =
#'        rep(c(rep(c("fn1", "fn2"), each=3)), 2),
#'        "patient_id"=rep(c(rep(c("pat1", "pat2"), each=3)), 2),
#'        "altered_entity"=c(rep(c("geneA", "geneB", "geneC"), 4)),
#'        "clone1"=c(0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0),
#'        "clone2"=c(1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1),
#'        "tree_id"=c(rep(5, 6), rep(10, 6)))
#' clone_tbl_pat1 <- dplyr::filter(clone_tbl, patient_id == "pat1")
#' clone_tbl_pat2 <- dplyr::filter(clone_tbl, patient_id == "pat2")
#' rates_exmpl_1 <- compute_rates_clon_excl(clone_tbl_pat1)
#' rates_exmpl_2 <- compute_rates_clon_excl(clone_tbl_pat2)
#' avg_rates_m <- apply(cbind(rates_exmpl_1, rates_exmpl_2), 2, mean)
#' names(avg_rates_m) <- c(names(rates_exmpl_1)[1], names(rates_exmpl_2)[1])
#' values_clon_excl_num_trees_pat1 <- get_hist_clon_excl(clone_tbl_pat1)
#' values_clon_excl_num_trees_pat2 <- get_hist_clon_excl(clone_tbl_pat2)
#' list_of_num_trees_all_pats <-
#'  list(pat1=values_clon_excl_num_trees_pat1[[1]], 
#'       pat2=values_clon_excl_num_trees_pat2[[1]])
#' list_of_clon_excl_all_pats <- 
#'   list(pat1=values_clon_excl_num_trees_pat1[[2]],
#'        pat2=values_clon_excl_num_trees_pat2[[2]])
#' num_pat_pair_max <- 2
#' num_pairs_sim <- 10
#' ecdf_list <- generate_ecdf_test_stat(avg_rates_m, 
#'                            list_of_num_trees_all_pats, 
#'                            list_of_clon_excl_all_pats,
#'                            num_pat_pair_max, num_pairs_sim)
#' sim_res <- generate_test_stat_hist(avg_rates_m, 
#'                                   list_of_num_trees_all_pats, 
#'                                   list_of_clon_excl_all_pats, 
#'                                   ecdf_list, 
#'                                   num_pat_pair_max, 
#'                                   num_pairs_sim)
generate_test_stat_hist <- function(avg_rates_m, list_of_num_trees_all_pats, 
    list_of_clon_excl_all_pats, 
    ecdf_list, num_pat_pair_max, 
    num_pairs_sim, beta_distortion=1000){
    test_statistic <- pval <- num_pat_pair <- NULL
    stopifnot(is.numeric(avg_rates_m))
    stopifnot(is.list(list_of_num_trees_all_pats))
    stopifnot(is.list(list_of_clon_excl_all_pats))
    stopifnot(is.list(ecdf_list))
    stopifnot(is.numeric(num_pat_pair_max))
    stopifnot(is.numeric(num_pairs_sim))
    for(this_rate_m in avg_rates_m){
        stopifnot(this_rate_m <= 1 && this_rate_m >= 0)
    }
    num_rates=length(avg_rates_m)
    list_of_clon_excl_frac_trees_all_pats <- list(list_of_num_trees_all_pats, 
        list_of_clon_excl_all_pats)
    stopifnot(length(list_of_clon_excl_frac_trees_all_pats) == 2)
    stopifnot(length(list_of_clon_excl_frac_trees_all_pats[[1]]) == num_rates)
    stopifnot(length(list_of_clon_excl_frac_trees_all_pats[[2]]) == num_rates)
    stopifnot(is.numeric(beta_distortion))
    if(beta_distortion < 50){
        warning("A beta_distortion value of ", beta_distortion, 
        " is not recommended as it will distort the rates a lot.",
        " Please consider choosing a higher value, e.g. 100,",
        " or 1000.")
    }
    num_pat_pair_max <- min(c(num_pat_pair_max, num_rates))
    
    ## also make sure that the rates and the lists with histogram are in the 
    ## same patient order:
    for (i in seq_len(num_rates) ){
        if((names(avg_rates_m)[i] != 
        names(list_of_clon_excl_frac_trees_all_pats[[1]])[i]) ||
        (names(avg_rates_m)[i] != 
        names(list_of_clon_excl_frac_trees_all_pats[[2]])[i]) ){
            stop("[Function: generate_test_stat_hist]: The rates seem to be",
            "\nin a different patient order than the ",
            "lists with patient's histograms\nlist_of_num_trees_all_pats,",
            " list_of_clon_excl_all_pats. See for example\nthe entry number ",
            i , " of the rates.\nMake sure that patients rates and the lists",
            " of histogram is ",
            "in the same patient order!")
        }
    }
    
    ## message to user
    message("There were ", num_rates, " rates of average clonal ",
        "exclusivity provided.")
    message("The average rates will be distorted with a beta",
        " distribution and M=alpha + beta=", beta_distortion)
    message("Additionally, from the ", num_rates, " patients,",
        " the values from the data for all pairs of",
        " # clon. excl./ # occurrence across all trees is ",
        "provided.")
    message("From those, it will be sampled ", num_pairs_sim, 
        " times to obtain the test statistics and the p-values",
        " under the null.")
    message("This procedure is done for num_patients=2,...,", 
        num_pat_pair_max)
    
    lrtest_stat_list <- list()
    lrtest_stat_list[[1]] <- NULL
    ## loop over the number of patients 1, 2, ..., num_pat_pair_max
    for (i in seq(2, num_pat_pair_max)){
        message("Generate ECDF for the case where a gene pair is ",
        "mutated together in ", i, " patient(s).")
        
        ## generate samples of the test statistic under null for i patients
        this_test_stat_hist <- 
            suppressMessages(build_null_test_statistic(avg_rates_m, 
            list_of_clon_excl_frac_trees_all_pats, 
            num_pat_pair=i, num_pairs_sim=num_pairs_sim,
            beta_distortion=beta_distortion))
        stopifnot(dim(this_test_stat_hist)[2] == (2+i*4))
        
        stopifnot(is.tbl(this_test_stat_hist))
        stopifnot("test_statistic" %in% colnames(this_test_stat_hist))
        stopifnot(dim(this_test_stat_hist)[1] == num_pairs_sim)
        ## as a null reference ecdf, the pre-computed ecdf of the test 
        ## statistic under the null is chosen
        ## for that, the ecdf with num_shared_pats is taken. It is the 
        ## num_shared_pats'th ecdf in the list
        this_ecdf <- ecdf_list[[i]]
        
        ## compute a p-value given the ecdf of the test statistic ecdf(T) 
        ## from the null distribution
        ## p_value=P(T>t | H_0 true)=1-ecdf(t) #### (upper-tailed test)
        all_test_stats <- this_test_stat_hist %>% 
            pull(test_statistic)
        
        all_pvals <- vapply(all_test_stats, function(x){
            this_pval <- 1-this_ecdf(x)
        }, numeric(1))
        
        this_test_stat_hist <- this_test_stat_hist %>% 
            mutate(pval=unlist(all_pvals)) %>% 
            mutate(num_pat_pair=i)
        
        lrtest_stat_list[[i]] <- this_test_stat_hist
    }
    
    return(lrtest_stat_list)
}




#' Generate samples from the test statistic under the null distribution - 
#' here we take the average rates of clonal exclusivity
#' across trees, and also the histogram for each patient over all pairs 
#' with the values # clon. excl./#trees.
#'
#' This function simulates gene pairs for the likelihood ratio test to 
#' generate values from the test statistic under
#' the null. It draws the average rates of clonal exclusivity from the 
#' ones provided by the user. That is,
#' the average rates of clonal exclusivity have to be computed first for 
#' each patient. The number of patients the
#' simulated pairs are mutated in can be specified with \code{num_pat_pair}. 
#' This function can be used to build the
#' ecdf of the test statistic under the null hypothesis (see Examples). 
#' The patients in which the simulated
#' pairs are mutated in are randomly selected proportional to the number 
#' of pairs in a patient.
#'
#' @title Simulate pairs to generate values of the test statistic under 
#' the null distribution
#' @param avg_rates_m The average rates of clonal exclusivity to be sampled 
#' from.
#' @param list_of_clon_excl_frac_trees_all_pats The list of two lists.
#' The first one contains a list entry for each patient containing the 
#' vector with the values of the information from each pair
#' in a patient of how often it was mutated across trees. The second list 
#' entry is a list with an entry for each patient that is a vector with the 
#' values of in how many trees the pair was clonally exclusive. The patient 
#' ordering in the lists has to be the same as in \code{avg_rates_m}.
#' @param num_pat_pair The number of patients the simulated pairs are mutated 
#' in.
#' @param num_pairs_sim The number of simulated gene/pathway pairs to be 
#' generated.
#' @param beta_distortion The value \code{M=alpha + beta} for the beta 
#' distribution, with which the average rates will be distorted. The bigger 
#' the M the higher the distribution is peaked around the actual rate. 
#' Therefore, the lesser the M, the more distorted the rates will be. 
#' Default: 1000.
#' @author Ariane L. Moore
#' @return The return value is a tibble with the columns 'test_statistic', 
#' 'mle_delta', and num_pat_pair columns with the respective rates that 
#' were drawn for each of the patients, num_pat_pair columns with the 
#' respective number of mutated times across
#' trees, and num_pat_pair columns with the respective number of times 
#' of being clonally exclusive across trees, and num_pat_pair columns with
#' the rate that was distorted by the beta distribution. The 'test_statistic' 
#' is the test statistic of the likelihood ratio test. The 'mle_delta' is
#' the maximum likelihood estimate of the delta for the elevated clonal 
#' exclusivity rate in the alternative model of the likelihood ratio test.
#' @export
#' @importFrom magrittr "%>%"
#' @importFrom dplyr tibble
#' @importFrom tibble add_column
#' @importFrom stats rbeta
#' @examples
#' avg_rates_m=c(0.4, 0.3)
#' list_of_clon_excl_frac_trees_all_pats <- list(list(c(5, 4, 5), c(5, 4)), 
#'                                   list(c(4, 4, 3), c(3, 2)))
#' sim_pairs <- build_null_test_statistic(avg_rates_m,
#'              list_of_clon_excl_frac_trees_all_pats, 2, 100, 
#'              beta_distortion=100)
#' ecdf_test_stat <- 
#'   ecdf(as.numeric(as.character(sim_pairs$test_statistic)))
#' plot(ecdf_test_stat, 
#'  main="ECDF of the test statistic when num_pat_pair=2")
#' # assume the observed test statistic t=6.0, 
#' # compute a p-value given the ecdf of
#' # the test statistic ecdf(T) from the null distribution
#' # p_value=P(T>t | H_0 true)=1-ecdf(t) ## (upper-tailed test)
#' p_value <- 1-ecdf_test_stat(6.0)
build_null_test_statistic <- function(avg_rates_m, 
    list_of_clon_excl_frac_trees_all_pats,
    num_pat_pair, num_pairs_sim, 
    beta_distortion=1000){
    test_statistic <- mle_delta <- NULL
    stopifnot(is.numeric(avg_rates_m))
    stopifnot(is.list(list_of_clon_excl_frac_trees_all_pats))
    stopifnot(is.numeric(num_pat_pair))
    stopifnot(is.numeric(num_pairs_sim))
    stopifnot(is.numeric(beta_distortion))
  
    message("The rates of average clonal exclusivity will be sampled ",
        "from the ones provided.")
    num_pat_rates <- length(avg_rates_m)
    stopifnot(length(list_of_clon_excl_frac_trees_all_pats) == 2)
    stopifnot(length(list_of_clon_excl_frac_trees_all_pats[[1]]) == 
        num_pat_rates)
    stopifnot(length(list_of_clon_excl_frac_trees_all_pats[[2]]) ==
        num_pat_rates)
    message("There were ", num_pat_rates,
    " rates of average clonal exclusivity provided. From those and from",
    " the list with observed values for num_clon_excl/num_trees,",
    " it will be sampled ", num_pat_pair," times in each iteration.")
    stopifnot(num_pat_rates >= num_pat_pair)
    message("The average rates will be distorted with a beta ",
    "distribution and M=alpha+beta=", beta_distortion)
  
    ## get the number of pairs there are in each patient
    num_pairs_per_pat <- 
    unlist(lapply(list_of_clon_excl_frac_trees_all_pats[[1]], 
        function(x){length(x)}))
    stopifnot(length(num_pairs_per_pat) == num_pat_rates)
    message("On average, a patient has ", mean(num_pairs_per_pat), 
    " pairs. (max of all trees)")
  
    ## message to user
    message("Requested number of gene pairs in simulation: ", 
    num_pairs_sim)
    message("Number of patients the gene pairs are mutated in: ", 
    num_pat_pair)
  
    these_rates_m <- avg_rates_m
    num_pats_real <- length(these_rates_m)
    message("Found ", num_pats_real, " different rates.")
  
    ## each pair is mutated in num_pat_pair patients
    pairs_num_pat <- c(rep(num_pat_pair, num_pairs_sim))
  
    ## now for each gene pair, we know in how many patients it is mutated, 
    ## e.g. X patients
    ## next, we draw X indices between 1 and number of rates of average 
    ## clonal exclusivity
    ## randomly for each of the pairs
    ## these indices indicate in which patients the pairs are mutated in
    pairs_indices_pats_list <- vapply(pairs_num_pat, function(x){
        stopifnot(x == num_pat_pair)
        sample(num_pats_real, x, prob=num_pairs_per_pat) 
        ## sample without replacement and by weighting the pats 
        ## according to their num pairs
    }, numeric(num_pat_pair))
  
    ## for each pair, we know in which patients they are both mutated,
    ## for each patient, we know their rates, and we have the histogram of
    ## in how many trees a pair was clonally exclusive and in how many trees
    ## it was clonally exclusive.
    
    lr_test_res_list <- list()
    all_pairs_test_stat <- c()
    all_pairs_delta <- c()
    all_pairs_num_pat <- c()
    all_pairs_pat_rates <- list()
    all_pairs_num_trees <- list()
    all_pairs_num_clon_excl <- list()
    all_pairs_pat_rates_beta_distorted <- list()
    
    if(num_pat_pair == 1){
        lr_test_res_list <- vapply(pairs_indices_pats_list, function(x){
            this_num_pat_this_pair <- length(x) ## this means the simulated 
            ## pair is mutated in the xth patient
            these_rate_this_pair <- avg_rates_m[x]
            num_pairs_this_pat <-
                length(list_of_clon_excl_frac_trees_all_pats[[1]][[x]])
            hist_entry <- sample(num_pairs_this_pat, 1)
            num_trees_this_pair <- 
            list_of_clon_excl_frac_trees_all_pats[[1]][[x]][hist_entry]
            num_clon_excl_this_pair <- 
            list_of_clon_excl_frac_trees_all_pats[[2]][[x]][hist_entry]
            stopifnot(num_clon_excl_this_pair <= num_trees_this_pair)
            
            ## distort the actual rates with a beta distribution that is 
            ## peaked around the actual rate
            ## this makes sure that the ecdf is smooth
            stopifnot(length(these_rate_this_pair) == 1)
            this_rate <- these_rate_this_pair
            cnt <- 1
            this_alpha <- this_rate * beta_distortion
            this_beta <- (1-this_rate) * beta_distortion
            this_beta_distored_rate <- rbeta(1, shape1=this_alpha, 
                shape2=this_beta)
            beta_distorted_rates <- this_beta_distored_rate
            
            ## sanity checks
            stopifnot(this_beta_distored_rate >= 0 && 
                this_beta_distored_rate <= 1)
            if(this_beta_distored_rate == 0 && 
            num_clon_excl_this_pair[cnt] == num_trees_this_pair[cnt]){
                stop("[Function: build_null_test_statistic]:",
                " The beta distorted rate is 0,",
                "\nbut there is a pair which is clonally ",
                "exclusive across ALL trees,\ni.e. num_clon_excl ",
                "== num_trees_pair, which is == ", 
                num_trees_this_pair[cnt],
                ". The computation of the delta\nwill not be",
                " possible since the likelihood is zero. ",
                "Make sure that the rates\nare not distorted",
                " too much by the beta distribution.",
                "The original rate of\nthe current simulated pair was ",
                this_rate)
            } else if (this_beta_distored_rate == 1 && 
            num_clon_excl_this_pair[cnt] == 0){
                stop("[Function: build_null_test_statistic]: ",
                "The beta distorted rate is 1,", 
                ",\nbut there is a pair which is clonally ",
                "exclusive across ALL trees,\ni.e. num_clon_excl",
                " == num_trees_pair, which is == ",
                num_trees_this_pair[cnt],
                ". The computation of the delta\nwill not be possible",
                " since the likelihood is zero. ",
                "Make sure that the rates\nare not distorted too",
                " much by the beta distribution.",
                "The original rate of\nthe current simulated pair was ",
                this_rate)
            }
            
            this_lr_test_res <- 
            suppressMessages(compute_test_stat_avg_rate(beta_distorted_rates, 
            num_trees_this_pair, num_clon_excl_this_pair))
            this_test_statistic <- this_lr_test_res[[1]]
            this_delta <- this_lr_test_res[[2]]
            
            return(c(this_test_statistic, 
                this_delta, 
                this_num_pat_this_pair,
                these_rate_this_pair, 
                num_trees_this_pair, 
                num_clon_excl_this_pair, 
                beta_distorted_rates))
        }, numeric(7))
        
        
        stopifnot(dim(lr_test_res_list)[1] == 7)
        stopifnot(dim(lr_test_res_list)[2] == num_pairs_sim)
        
        ## extract the test statistics
        all_pairs_test_stat <- unlist(lr_test_res_list[1,])
        
        ## extract the deltas
        all_pairs_delta <- unlist(lr_test_res_list[2,])
        
        ## extract num_patients, expected_num_co, observed_num_co
        all_pairs_num_pat <- unlist(lr_test_res_list[3,])
        for (this_num in seq_len(num_pairs_sim)){
            stopifnot(all_pairs_num_pat[this_num] == num_pat_pair)
        }
        
        ## get the general rates of clonal exclusivity
        all_pairs_pat_rates[[1]] <- unlist(lr_test_res_list[4,])
        for (this_rate in seq_len(num_pairs_sim)){
            stopifnot(all_pairs_pat_rates[[1]][this_rate] >= 0)
            stopifnot(all_pairs_pat_rates[[1]][this_rate] <= 1)
        }
        
        ## get the number of times a pair was mutated across trees
        all_pairs_num_trees[[1]] <- unlist(lr_test_res_list[5,])
        
        ## get the number of times a pair was clonally exclusive across trees
        all_pairs_num_clon_excl[[1]] <- unlist(lr_test_res_list[6,])
        for (this_sim in seq_len(num_pairs_sim)){
            stopifnot(all_pairs_num_clon_excl[[1]][this_sim] <=
                all_pairs_num_trees[[1]][this_sim])
        }
        ## get the beta distorted rates of clonal exclusivity
        all_pairs_pat_rates_beta_distorted[[1]] <- unlist(lr_test_res_list[7,])
        for (this_rate in seq_len(num_pairs_sim)){
            stopifnot(all_pairs_pat_rates_beta_distorted[[1]][this_rate] >= 0)
            stopifnot(all_pairs_pat_rates_beta_distorted[[1]][this_rate] <= 1)
        }
        
        
    } else {
        lr_test_res_list <- apply(pairs_indices_pats_list, 2, function(x){
            ## this means the simulated pair
            this_num_pat_this_pair <- length(x)  
            ## is mutated in the xth patient
            these_rate_this_pair <- avg_rates_m[x]
            this_hist_sample <- vapply(x, function(y){
                num_pairs_this_pat <-
                    length(list_of_clon_excl_frac_trees_all_pats[[1]][[y]])
                hist_entry <- sample(num_pairs_this_pat, 1)
                num_trees_this_pair <- 
                list_of_clon_excl_frac_trees_all_pats[[1]][[y]][hist_entry]
                num_clon_excl_this_pair <- 
                list_of_clon_excl_frac_trees_all_pats[[2]][[y]][hist_entry]
                return(c(num_trees_this_pair, num_clon_excl_this_pair))
            }, numeric(2))
            stopifnot(is.matrix(this_hist_sample))
            stopifnot(dim(this_hist_sample)[1] == 2)
            stopifnot(dim(this_hist_sample)[2] == this_num_pat_this_pair)
            
            num_trees_this_pair <- this_hist_sample[1,]
            num_clon_excl_this_pair <- this_hist_sample[2,]
            
            ## distort the actual rates with a beta distribution that is peaked 
            ## around the actual rate
            ## this makes sure that the ecdf is smooth
            beta_distorted_rates<-vapply(these_rate_this_pair, 
                function(this_rate){
                    cnt <- which(these_rate_this_pair == this_rate)
                    this_alpha <- this_rate * beta_distortion
                    this_beta <- (1-this_rate) * beta_distortion
                    this_beta_distored_rate <- 
                        rbeta(1, shape1=this_alpha, shape2=this_beta)
                    ## sanity checks
                    stopifnot(this_beta_distored_rate >= 0 && 
                        this_beta_distored_rate <= 1)
                    if(this_beta_distored_rate == 0 && 
                        num_clon_excl_this_pair[cnt] == 
                        num_trees_this_pair[cnt]){
                        stop("[Function: build_null_test_statistic]: The beta",
                            " distorted rate is 0,",
                            "\nbut there is a pair which is clonally ",
                            "exclusive across ALL trees,\ni.e. num_clon_excl",
                            " == num_trees_pair, which is == ", 
                            num_trees_this_pair[cnt],
                            ". The computation of the delta\nwill not be",
                            " possible since the likelihood is zero. ",
                            "Make sure that the rates\nare not distorted too",
                            " much by the beta distribution.",
                            "The original rate of\nthe current simulated pair ",
                            "was ", this_rate)
                    } else if (this_beta_distored_rate == 1 && 
                        num_clon_excl_this_pair[cnt] == 0){
                        stop("[Function: build_null_test_statistic]: ",
                            " The beta distorted rate is 1,", 
                            ",\nbut there is a pair which is clonally ",
                            "exclusive across ALL trees,\ni.e. num_clon_excl",
                            " == num_trees_pair, which is == ",
                            num_trees_this_pair[cnt],
                            ". The computation of the delta\nwill not be",
                            " possible since the likelihood is zero. ",
                            "Make sure that the rates\nare not distorted too",
                            " much by the beta distribution.",
                            "The original rate of\nthe current simulated pair ",
                            "was ", this_rate)
                    }
                    return(this_beta_distored_rate)
            }, numeric(1))
            
            this_lr_test_res <- 
            suppressMessages(compute_test_stat_avg_rate(beta_distorted_rates,
            num_trees_this_pair, 
            num_clon_excl_this_pair))
            this_test_statistic <- this_lr_test_res[[1]]
            this_delta <- this_lr_test_res[[2]]
            
            return(list(this_test_statistic,this_delta,this_num_pat_this_pair,
                these_rate_this_pair,
                num_trees_this_pair, num_clon_excl_this_pair, 
                beta_distorted_rates))
        })
        
        ## extract the test statistics
        all_pairs_test_stat <- unlist(lapply(lr_test_res_list,
            function(x){x[[1]]}))
        
        ## extract the deltas
        all_pairs_delta <- unlist(lapply(lr_test_res_list, 
            function(x){x[[2]]}))
        
        ## extract num_patients, expected_num_co, observed_num_co
        all_pairs_num_pat <- unlist(lapply(lr_test_res_list,
            function(x){x[[3]]}))
        
        ## get the general rates of clonal exclusivity
        for(i in seq_len(num_pat_pair)){
            all_pairs_pat_rates[[i]] <- unlist(lapply(lr_test_res_list, 
                function(x){x[[4]][i]}))
            all_pairs_num_trees[[i]] <-  unlist(lapply(lr_test_res_list,
                function(x){x[[5]][i]}))
            all_pairs_num_clon_excl[[i]] <- unlist(lapply(lr_test_res_list,
                function(x){x[[6]][i]}))
            all_pairs_pat_rates_beta_distorted[[i]] <- 
                unlist(lapply(lr_test_res_list, function(x){x[[7]][i]}))
        }
    }
    
    stopifnot(length(all_pairs_test_stat) == num_pairs_sim)
    stopifnot(length(all_pairs_delta) == num_pairs_sim)
    stopifnot(length(all_pairs_num_pat) == num_pairs_sim)
    
    ## create final tibble to return
    ## contains the columns 'test_statistic',
    ## 'mle_delta', and num_pat_pair columns with the respective rates 
    ## that were drawn for each of the patients.
    res_all_pairs <- tibble(test_statistic=
        all_pairs_test_stat,
        mle_delta=all_pairs_delta)
    ## now add the num_pat_pair columns with the rates of average clonal 
    ## exclusivity
    for(i in seq_len(num_pat_pair)){
        col_name <- paste("pat", i, sep="")
        res_all_pairs <- res_all_pairs %>%
            add_column(new_col=all_pairs_pat_rates[[i]])
        names(res_all_pairs)[names(res_all_pairs) == "new_col"] <- col_name
    }
    ## now add the all_pairs_num_trees columns with the number of trees 
    ## the pairs were mutated
    for(i in seq_len(num_pat_pair)){
        col_name <- paste("pat", i, "_num_trees", sep="")
        res_all_pairs <- res_all_pairs %>%
          add_column(new_col=all_pairs_num_trees[[i]])
        names(res_all_pairs)[names(res_all_pairs) == "new_col"] <- col_name
    }
    ## now add the all_pairs_num_clon_excl columns with the number of 
    ## times the pair was clonally exclusive across trees
    for(i in seq_len(num_pat_pair)){
        col_name <- paste("pat", i, "_num_clon_excl", sep="")
        res_all_pairs <- res_all_pairs %>%
            add_column(new_col=all_pairs_num_clon_excl[[i]])
        names(res_all_pairs)[names(res_all_pairs) == "new_col"] <- col_name
    }
    ## now add the num_pat_pair columns with the rates of average clonal 
    ## exclusivity that were distorted by the beta distribution
    for(i in seq_len(num_pat_pair)){
        col_name <- paste("pat", i, "_beta_distorted", sep="")
        res_all_pairs <- res_all_pairs %>%
        add_column(new_col=all_pairs_pat_rates_beta_distorted[[i]])
        names(res_all_pairs)[names(res_all_pairs) == "new_col"] <- col_name
    }
    
    stopifnot(dim(res_all_pairs)[2] == (num_pat_pair*4 + 2))
    stopifnot("pat1" %in% colnames(res_all_pairs))
    stopifnot("pat1_beta_distorted" %in% colnames(res_all_pairs))
    return(res_all_pairs)
}


