
#' Compute test statistic that is based on the average rates of clonal 
#' exclusivity of a patient, and the observed
#' number of times a pair was clonally exclusive across several trees of 
#' the tree inference.
#'
#' For a given gene/pathway pair, this function takes as input: the average 
#' rates of clonal exclusivity of all patients
#' in which the pair is mutated, the number of trees among the trees in the 
#' collection of trees from each patient in which
#' the pair was occurring, and the number of times it was clonally exclusive
#' 
#' @title Compute the test statistic of the clonal exclusivity test (lrtest).
#' @param avg_rates_m The vector of average rates of clonal exclusivity of 
#' each patient the pair is mutated in. It was computed for each patient
#' separately, and is averaged over all gene pairs and all trees. Expected 
#' to be in the same (patient) order as the other inputs to this function.
#' @param num_trees_pair The vector with the number of tree inferences in 
#' which the pair was occurring in. Has to be the same order
#' as \code{avg_rates_m}.
#' @param num_clon_excl The vector with the number of times the pair was 
#' clonally exclusive in the trees in each patient.
#' Has to be the same order as \code{avg_rates_m}.
#' @author Ariane L. Moore
#' @export
#' @importFrom gtools logit inv.logit
#' @importFrom maxLik maxLik
#' @importFrom stats coef
#' @return A list with the test statistic of the clonal exclusivity test 
#' (lrtest), and the maximum likelihood estimate of delta.
#' @examples
#' compute_test_stat_avg_rate(c(0.1, 0.2), c(10, 10), c(9, 7))
#' compute_test_stat_avg_rate(c(0.05, 0.23), c(20, 20), c(8, 5))
compute_test_stat_avg_rate <- function(avg_rates_m, num_trees_pair,
    num_clon_excl){
    stopifnot(is.numeric(avg_rates_m))
    stopifnot(is.numeric(num_trees_pair))
    stopifnot(is.numeric(num_clon_excl))
    
    num_shared_pats <- length(avg_rates_m)
    stopifnot(num_shared_pats == length(num_trees_pair))
    stopifnot(num_shared_pats == length(num_clon_excl))
    for (i in seq_along(avg_rates_m)){
        stopifnot(num_clon_excl[i] <= num_trees_pair[i])
        stopifnot(num_trees_pair[i]>0)
    }
    
    ## if the rate is 0, but num_clon_excl == num_trees_pair, 
    ## this will crash this should anyways not be the case. 
    ## If the rate is 0, there should 
    ## be no pair that is clonally exclusive
    ## at all
    for(i in seq_along(avg_rates_m)){
        if(avg_rates_m[i] == 0 && 
            num_clon_excl[i] == num_trees_pair[i]){
            stop("[Function: compute_test_stat_avg_rate]: The rate is 0,",
            "\nbut there is a pair which is clonally ",
            "exclusive across ALL trees,\ni.e. num_clon_excl ==",
            " num_trees_pair, which is == ", num_trees_pair[i],
            ". This is\nnot possible. Make sure that the rates are ",
            "in the same patient\norder as the ",
            "histogram of pairs across trees, or that the rates\nare",
            " not distorted too much by the beta distribution.")
        }
    }
    
    ## compute the test statistic of the log-likelihood ratio test:
    ## H0=clonal exclusivity just by chance, explained by the general 
    ## clonal exclusivity rates in
    ##  the shared patients, likelihood: s_1_{A,B} * s_2_{A,B} * s_3_{A,B} 
    ## if the current pair {A,B} occurred in patients 1, 2, 3
    ## each of those s_i_{A,B}=(num_clon_excl/num_trees)*m_i + 
    ## ((num_trees - num_clon_excl)/num_trees)*(1-m_i)
    ##
    ## H1=clonal exclusivity at elevated rates, logit(m_star_i)=logit(m_i)
    ## + delta, likelihood: 
    ## s_star_1_{A,B} * s_star_2_{A,B} * s_star_3_{A,B}
    ##  take mle estimate for the delta
    ## LRTestStatistic ~ Chi^2 with 1 degree of freedom (under H0) OR 
    ## compare it to the empirical distribution under the null
    
    ## get likelihood of the null:
    ## helper function
    compute_weighted_p <- function(num_cl, num_trees, avg_m){
        stopifnot(num_cl <= num_trees)
        stopifnot(num_trees > 0)
        stopifnot(avg_m <= 1 && avg_m >= 0)
        res <- (num_cl/num_trees)*avg_m + 
            ((num_trees-num_cl)/num_trees) * (1-avg_m)
        stopifnot(is.numeric(res))
        stopifnot(res >= 0 && res <= 1)
        return(res)
    }
    ## compute individual factors from each patient
    l_null_factors <- apply(cbind(num_clon_excl, 
                                  num_trees_pair, 
                                  avg_rates_m),
        1, function(x){compute_weighted_p(x[1], x[2], x[3])})
    l_null <- prod(l_null_factors)
    message("l_null=", l_null)
    stopifnot(l_null >= 0 && l_null <= 1)
    
    
    ## get likelihood of alternative:
    ## get mle of delta
    ## maximize likelihood over delta
    ## to obtain an estimate for the delta (MLE)
    ## function to get the negative log-likelihood of 
    ## the data for a given delta
    
    ## get the logits of the m's: logit(m_i)=m_i/(1-m_i)
    logit_avg_rates_m <- logit(avg_rates_m) 
    ## without the gtools function, this crashes when one 
    ## of the rates is 1.0
    
    ## They should be identical: i.e. if delta is zero, m and m_star should 
    ## be the same
    if(abs((inv.logit(logit_avg_rates_m))[1]
        - avg_rates_m[1]) >= 0.01){
        stop("The sanity check with delta=0 did not have the right",
        " result. m should be identical to m_star if delta=0!!")
    }
    
    ## define the log-likelihood function
    this_LL <- function(this_delta) {
        avg_rates_m_star <- inv.logit(logit_avg_rates_m + 
            this_delta)
        l_alt_factors <- apply(cbind(num_clon_excl, num_trees_pair, 
            avg_rates_m_star), 1, 
            function(x){compute_weighted_p(x[1], x[2], x[3])})
        return(sum(log(l_alt_factors)))
    }
    
    ## get the mle estimate of the delta
    mle <- maxLik(logLik=this_LL, start=c(this_delta=0))
    mle_delta <- coef(mle)
    message("Maximum likelihood estimate of delta given the data is: ", 
        mle_delta)
    
    ## compute the alternative likelihood
    avg_rates_m_star <- inv.logit(logit_avg_rates_m + mle_delta)
    l_alt_factors <- 
        apply(cbind(num_clon_excl, num_trees_pair, avg_rates_m_star), 1, 
        function(x){compute_weighted_p(x[1], x[2], x[3])})
    l_alt <- prod(l_alt_factors)
    message("l_alt=", l_alt)
    
    ## compute the test statistic
    lr_test_statistic <- -2 * log( l_null / l_alt )
    message("lr_test_statistic=", lr_test_statistic)
    stopifnot(is.numeric(lr_test_statistic))
    
    return(list(lr_test_statistic=lr_test_statistic, 
        mle_delta=mle_delta))
}



#' Compare the likelihood ratio test statistic to its ecdf under the null 
#' for two mutated genes/pathways in clones of
#' patients.
#'
#' Tests whether the observed number of clonal exclusivities of mutated 
#' entities (genes or pathways)
#' A and B in clones of patients is significantly different from what would 
#' be expected given the average clonal exclusivity rates. The observed 
#' test statistic is compared to the ecdf of the test statistic under the 
#' null hypothesis.
#'
#' @title Compare observed likelihood ratio test statistic to its ecdf under 
#' null.
#' @param entA One gene/pathway of the pair.
#' @param entB The other gene/pathway of the pair.
#' @param clone_tbl The clone tibble as generated with 
#' \code{\link{create_tbl_tree_collection}} from several trees of
#' the tree inference, i.e. it also contains a column 'tree_id'.
#' @param avg_rates_m The average rates of clonal exclusivity for each 
#' patient. The name of each rate is the respective patient_id.
#' @param ecdf_list The list of ECDF's of the test statistic under the 
#' null distribution. Can be generated with 
#' \code{\link{generate_ecdf_test_stat}}.
#' It is important that the rates that are used for that are the same as
#'  the \code{avg_rates_m} here.
#' @param alternative The character indicating whether pairs should only be 
#' tested if delta > 0 or if all pairs should be tested.
#' Can be one of "greater" or "two.sided".
#' @author Ariane L. Moore
#' @return Returns list(p_val, num_patients, mle_delta, test_statistic), 
#' i.e. a list with the p-value, the number of patients in which both of 
#' the genes/pathways were mutated, the maximum likelihood estimate of
#' the delta, and the test statistic.
#' @export
#' @importFrom magrittr "%>%"
#' @importFrom dplyr tibble filter is.tbl pull select
#' @examples
#' clone_tbl <- dplyr::tibble("file_name"=
#'            rep(c(rep(c("fn1", "fn2"), each=3)), 2),
#'           "patient_id"=rep(c(rep(c("pat1", "pat2"), each=3)), 2),
#'           "altered_entity"=c(rep(c("geneA", "geneB", "geneC"), 4)),
#'           "clone1"=c(0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0),
#'           "clone2"=c(1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1),
#'           "tree_id"=c(rep(5, 6), rep(10, 6)))
#' clone_tbl_pat1 <- dplyr::filter(clone_tbl, patient_id == "pat1")
#' clone_tbl_pat2 <- dplyr::filter(clone_tbl, patient_id == "pat2")
#' rates_exmpl_1 <- compute_rates_clon_excl(clone_tbl_pat1)
#' rates_exmpl_2 <- compute_rates_clon_excl(clone_tbl_pat2)
#' avg_rates_m <- apply(cbind(rates_exmpl_1, rates_exmpl_2), 2, mean)
#' names(avg_rates_m) <- c(names(rates_exmpl_1)[1], 
#'     names(rates_exmpl_2)[1])
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
#'                                list_of_num_trees_all_pats, 
#'                                list_of_clon_excl_all_pats, 
#'                                num_pat_pair_max, 
#'                                num_pairs_sim)
#' ecdf_lr_test_clon_excl_avg_rate("geneA", "geneB", clone_tbl, 
#'                                   avg_rates_m,
#'                                   ecdf_list, "greater")
ecdf_lr_test_clon_excl_avg_rate <- function(entA, entB, clone_tbl, 
    avg_rates_m, 
    ecdf_list, alternative){
    altered_entity <- patient_id <- file_name <- tree_id <- NULL
    stopifnot(is.character(entA))
    stopifnot(is.character(entB))
    stopifnot(is.tbl(clone_tbl))
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
    stopifnot(entA != entB)
    
    ## get all patients
    all_pats <- unique(as.character(clone_tbl$patient_id))
    num_pats <- length(all_pats)
    stopifnot(length(avg_rates_m) == num_pats)
    message("There were ", num_pats, " different patient ids.")
    
    ## get the number of trees
    all_tree_ids <- unique(as.character(clone_tbl$tree_id))
    num_trees <- length(all_tree_ids)
    message("There were ", num_trees, " different tree inferences.")
    
    ## check in which patients the current pair is mutated, and in wich 
    ## trees
    ## this is the list of trees, where each list entry is the 
    ## vector of patient ids
    ## in which both of them are mutated
    tree_to_pats_ent_pair <- lapply(all_tree_ids, function(this_tree){
        ## select a specific tree
        clone_tbl_this_tree <- clone_tbl %>% 
            filter(tree_id == this_tree)
      
        these_pats_mutated <- lapply(c(entA, entB), function(this_ent){
            this_ent_pats <- clone_tbl_this_tree %>%
                filter(altered_entity == this_ent) %>%
                pull(patient_id)
            stopifnot(length(unique(this_ent_pats)) == length(this_ent_pats))
            this_ent_pats
        })
        ## now only those where both of them are mutated
        these_pats_both_ents_mutated <- intersect(these_pats_mutated[[1]], 
            these_pats_mutated[[2]])
    })
    
    ## the number of patients in which they are both mutated
    ## (maximum of all trees)
    shared_pats_total <- unique(unlist(tree_to_pats_ent_pair))
    num_shared_pats_total <- length(shared_pats_total)
    
    ## if there are no patients in which both of them are mutated
    if(num_shared_pats_total == 0){
        message("The current gene/pathway pair: ", entA, ", " , entB,
        " is never mutated in the same patient.\n",
        "No test will be performed.")
        return(list(".", 0, ".", "."))
    } else if (num_shared_pats_total == 1) {
        message("The current gene/pathway pair: ", entA, ", " , entB,
        " is only mutated once in the same patient.\n",
        "Currently, no test will be performed.")
        return(list(".", 1, ".", "."))
    } else {
        ## if they share enough patient(s)
        ## get the rates of clonal exclusivity from the(se) patient(s)
        these_shared_pats_rates <- vapply(shared_pats_total, function(x){
            avg_rates_m[which(names(avg_rates_m) == x)]
        }, numeric(1))
        
        ## compute for all patients: number of trees the pair is mutated, 
        ## and number of times it was also clonally exclusive
        num_trees_pair_num_clon_excl <- vapply(shared_pats_total, 
        function(x){
        ## take the clone tibble ## just for the current patient
        this_clone_tbl <- clone_tbl %>%                     
            filter(patient_id == x) %>%                      
            select(-file_name, -patient_id)
        num_trees_num_clon_excl_this_pat <- 
            suppressMessages(get_hist_clon_excl_this_pat_this_pair(entA, 
            entB, 
            this_clone_tbl))
            return(num_trees_num_clon_excl_this_pat)
        }, numeric(2))
        stopifnot(dim(num_trees_pair_num_clon_excl)[1] == 2)
        num_trees_pair <- num_trees_pair_num_clon_excl[1,]
        num_clon_excl <- num_trees_pair_num_clon_excl[2,]
        
        ## compute the test statistic and the p-value of the likelihood 
        ## ratio test
        ## by comparing the observed test statistic to the ecdf of the test 
        ## statistic under the null
        res_lr_test <- c()
        res_lr_test <- compute_test_stat_avg_rate(these_shared_pats_rates, 
            num_trees_pair, 
            num_clon_excl)
        test_stat <- res_lr_test[[1]]
        mle_delta <- res_lr_test[[2]]
        
        if(alternative == "greater" && mle_delta <= 0){
            message("The current gene/pathway pair: ", entA, ", " , 
            entB,
            " has a delta of ", mle_delta, ".\n",
            "Since the alternative is 'greater', no test will be ",
            "performed.")
            return(list(pval=".", num_patients=num_shared_pats_total, 
            mle_delta=mle_delta, test_statistic=test_stat))
        } else {
            ## as a null reference ecdf, the pre-computed ecdf of the test 
            ## statistic under the null is chosen
            ## for that, the ecdf with num_shared_pats is taken. It is the 
            ## num_shared_pats'th ecdf in the list
            if(num_shared_pats_total > length(ecdf_list)){
                stop("[Function: ecdf_lr_test_clon_excl_avg_rate] The ",
                "length of the ecdf list\nis not long enough. The current ",
                "pair ", entA, " and ", entB, " seems to be mutated ",
                "together\nin more than ",length(ecdf_list),
                " patients. More precisely, they are mutated together in ",
                num_shared_pats_total,".\n",
                "Please re-run the step for generating the ecdf and use",
                " a higher value for\nnum_pat_pair_max ",
                "(at least ", num_shared_pats_total,").")
            }
            message("Test statistic is ", test_stat)
            message("The number of patients in which both genes ",
            "are mutated: ",
            num_shared_pats_total)
            this_ecdf <- ecdf_list[[num_shared_pats_total]]
    
            ## compute a p-value given the ecdf of the test statistic
            ## ecdf(T) from the null distribution
            ## p_value=P(T>t | H_0 true)=1-ecdf(t) #### (upper-tailed test)
            message("The value of the ecdf at the test statistic: ", 
            this_ecdf(test_stat))
            
            p_val <- 1-this_ecdf(test_stat)
            message("p_val=", p_val)
            
            return(list(p_val=p_val, num_patients=num_shared_pats_total, 
            mle_delta=mle_delta, test_statistic=test_stat))
        }
    }
}


