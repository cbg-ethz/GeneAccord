
## This function takes a dataframe and returns the column means again as a df
colmeans_fn <- function(data) {
  mean_vec <- colMeans(data, na.rm = TRUE)
  mean_vec[which(is.na(mean_vec))] <- NA
  data.frame(t(mean_vec))
}

## This function takes a dataframe of mutations in each clone of a tree
## and computes whether each pair of mutations is clonally exclusive or not
cx_fn <- function(data) {
  dmat <- as.matrix(data)
  summat <- t(dmat) %*% dmat
  summat <- 1*(summat > 0) # do they ever occur in the same clones
  namemat <- matrix(paste(rep(colnames(summat), each = nrow(summat)), rownames(summat),
                          sep = "_"), nrow = nrow(summat))
  rownames(namemat) <- rownames(summat)
  colnames(namemat) <- colnames(summat)
  genes <- diag(summat)
  genes[which(genes == 0)] <- NA # do the genes occur at all
  genemat <- t(t(genes)) %*% t(genes)
  cxmat <- genemat * (1 - summat)
  cxvec <- cxmat[lower.tri(cxmat)]
  names(cxvec) <- namemat[lower.tri(namemat)]
  data.frame(t(cxvec))
}

#' Turn a data frame of clonal mutations into a data frame of pairwise
#' clonal exclusivity patterns
#'
#' @title Transform dataframes from mutations to exclusivity.
#' @param clone_df The dataframe of mutational presence in each clone.
#' The first column must be the sample id, the next column an identifier
#' of each tree obtained per sample. The third column stores the frequency of each
#' clone. The remaining columns are for the genes/positions.
#' Each row is a clone (per sample and per tree) with its mutation status
#' over the genes recorded in binary with 1 representing mutation.
#' @author Jack Kuipers
#' @export
#' @return A dataframe with one row per sample with columns for each pair
#' of genes recording the average (over trees) clonal exclusivity.
clone_df_to_cx_df <- function(clone_df) {
  clone_df %>% select(-freq) %>%
  group_by(case, tree) %>%
  nest() %>%
  mutate(cx = map(data, cx_fn)) %>%
  select(-data) %>%
  unnest(cx) %>%
  ungroup(tree) %>%
  select(-tree) %>%
  nest() %>%
  mutate(cxs = map(data, colmeans_fn)) %>%
  select(-data) %>%
  unnest(cxs) %>%
  ungroup(case) %>%
  data.frame -> cx_df
  rownames(cx_df) <- cx_df[, 1] # move rownames
  cx_df[, -1]
}

#' Filter a data frame of pairwise clonal exclusivity patterns
#' to select samples with a clonal exclusivity probability strictly
#' between 0 and 1 and gene pairs which occur in at least a given
#' number of samples.
#'
#' @title Filters a clonal exclusivity dataframe.
#' @param cx_df A dataframe with one row per sample and columns for each
#' pair of genes recording the average (over trees) clonal exclusivity.
#' @param k The minimum of samples to contain
#' @author Jack Kuipers
#' @export
#' @return A filtered dataframe of clonal exclusivity patterns with an additional
#' first column of clonal exclusivity probabilities per sample.
filter_cx_df <- function(cx_df, k) {
  cx_probs <- rowMeans(cx_df, na.rm = TRUE)
  to_keep <- which(cx_probs > 0 & cx_probs < 1)
  cx_df <- cx_df[to_keep, ]
  cx_probs <- rowMeans(cx_df, na.rm = TRUE)
  n_pairs <- colSums(!is.na(cx_df))
  to_keep <- which(n_pairs >= k)
  cx_df <- cx_df[, to_keep]
  cx_df$cx_probs <- cx_probs # add observed clonal exclusivity probabilities
  cx_df[, c(ncol(cx_df), 2:ncol(cx_df) - 1)] # reorder columns
}

# Compute the difference in log-likelihood as delta is varied

ll_delta <- function(delta, probs, obs) {
  - sum(obs*log(probs + (1 - probs)*exp(delta)) + (1 - obs)*log(1 - probs + probs*exp(-delta)))
}

# test a single gene pair with the LR statistic

test_gene_pair <- function(cx_selected, k, exact_limit = 12, filter = TRUE, MC_reps = NULL, nu = 10) {
  gene_pair_df <- cx_selected[, c(0, k) + 1] # have probabilities in first column
  to_keep <- which(!is.na(gene_pair_df[, 2])) # and gene pair in second
  gene_pair_df <- gene_pair_df[to_keep, ]
  probs <- gene_pair_df[, 1] # extract the probabilities
  obs <- gene_pair_df[, 2] # extract the observations (weighted)
  ll_max <- optimize(ll_delta, interval = c(-10, 10), probs = probs, obs = obs, maximum = TRUE)
  LR <- 2*ll_max$objective

  n_obs <- length(probs)
  if (n_obs > exact_limit && is.null(MC_reps)) { # chi-squared version
    p_value <- pchisq(LR, df = 1, lower.tail = FALSE)
  } else {
    tol <- 1e-6
    if (is.null(MC_reps)) { # exact version
      null_probs <- rep(NA, 2^(n_obs))
      null_stats <- null_probs

      for (ii in 1:2^(n_obs)) {
        obs_vec <- tail(rev(as.numeric(intToBits(ii - 1))), n_obs)
        null_probs[ii] <- prod(probs*obs_vec + (1 - probs)*(1 - obs_vec)) # probability of that vector
        null_stats[ii] <- optimize(ll_delta, interval = c(-10, 10), probs = probs, obs = obs_vec, maximum = TRUE, tol = tol)$objective
      }
    } else { # MC version
      null_probs <- rep(1/MC_reps, MC_reps)
      null_stats <- rep(NA, MC_reps)

      for (ii in 1:MC_reps) {
        probs_noisy <- rbeta(n_obs, nu*probs, nu*(1-probs))
        obs_vec <- rbinom(n_obs, 1, probs_noisy)
        null_stats[ii] <- optimize(ll_delta, interval = c(-10, 10), probs = probs_noisy, obs = obs_vec, maximum = TRUE, tol = tol)$objective
      }
    }
    more_extreme <- which(null_stats >= ll_max$objective + 2*tol)
    as_extreme <- which(null_stats >= ll_max$objective - 2*tol)
    p_more_extreme <- sum(null_probs[more_extreme]) # those clearly above
    p_as_extreme <- sum(null_probs[as_extreme]) # those the same and above
    p_value <- mean(c(p_more_extreme, p_as_extreme)) # taking the average count the same with weight half
    if (filter && is.null(MC_reps)) { # remove from testing if no chance of significance
      most_extreme <- which(null_stats >= max(null_stats) - 2*tol) # find the most extreme values possible
      p_most_extreme <- sum(null_probs[most_extreme]) / 2 # find middle point
      if (p_most_extreme > 0.05) {
        p_value <- NA
      }
    }
  }
  # output as data frame
  data.frame(gene_pair = colnames(gene_pair_df)[2], n = n_obs, n_cx = sum(obs),
             delta = ll_max$maximum, LR = LR, p = p_value)
}

#' Takes a data frame of clonal exclusivity patterns with rows representing
#' patient sample, a first column giving the overall observed exclusivity
#' probability and subsequent columns providing the rate of clonal
#' exclusivity of the tree(s) of that patient sample.
#'
#' @title Perform the geneAccord clonal exclusivity test .
#' @param clonal_exclusivity_df A data frame with one row per sample with
#' the first column containing the and following columns for each pair
#' of genes recording the average (over trees) clonal exclusivity.
#' @param exact_limit Limit above which the asymptotic chi-squared distribution
#' is employed instead of the exact test which is more computationally intensive.
#' @param filter Whether to remove untestable gene pairs (TRUE, default) to
#' reduce multiple testing
#' @param MC_reps How many repetitions to use for Monte Carlo estimation of
#' the p-values. When NULL (default) exact or chi-squared testing is used.
#' @author Jack Kuipers
#' @export
#' @return A data fram with the gene pairs, number of observed samples containing
#' the gene pair, the number of samples with clonal exclusivity, the maximum
#' likelihood estimate of delta, the LR test statistic of the clonal exclusivity
#' test and the corresponding p-value and BH corrected q-value.
geneAccord <- function(clonal_exclusivity_df, exact_limit = 12, filter = TRUE, MC_reps = NULL) {
  cx_results <- NULL
  for (kk in 2:ncol(clonal_exclusivity_df) - 1) {
    cx_results <- rbind(cx_results, test_gene_pair(clonal_exclusivity_df, kk, exact_limit = exact_limit, filter = filter, MC_reps = MC_reps))
  }
  if (filter) {
    to_keep <- which(!is.na(cx_results$p))
    cx_results <- cx_results[to_keep, ]
  }
  cx_results <- cx_results[order(cx_results$p), ]
  cx_results$q <- p.adjust(cx_results$p, method = "BH")
  row.names(cx_results) <- 1:nrow(cx_results) # make row names the rank
  cx_results
}
