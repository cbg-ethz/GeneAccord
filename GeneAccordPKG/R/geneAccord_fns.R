
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
#' The first column named "case" must be the sample id, the next column named "tree"
#' an identifier of each tree obtained per sample. The third column named "freq"
#' stores the frequency of each clone. The remaining columns are for the genes/positions.
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
#' to select gene pairs appearing in a minimum number of samples.
#'
#' @title Filters a clonal exclusivity dataframe.
#' @param cx_df A dataframe with one row per sample and columns for each
#' pair of genes recording the average (over trees) clonal exclusivity.
#' @param k The minimum of samples to contain
#' @author Jack Kuipers
#' @export
#' @return A filtered dataframe of clonal exclusivity patterns with an additional
#' first two columns of clonal exclusivity probabilities and the
#' number of gene pairs per sample.
filter_cx_df <- function(cx_df, k) {
  cx_probs <- rowMeans(cx_df, na.rm = TRUE)
  to_keep <- which(!is.na(cx_probs)) # remove trees with one node
  cx_df <- cx_df[to_keep, ]
  cx_probs <- rowMeans(cx_df, na.rm = TRUE)
  n_gp <- rowSums(!is.na(cx_df)) # count total numbers of gene pairs
  infos <- cx_probs > 0 & cx_probs < 1 # which samples are informative
  infos[which(infos == FALSE)] <- NA
  n_infos <- colSums(!is.na(infos*cx_df))
  to_keep <- which(n_infos >= k)
  cx_df <- cx_df[, to_keep]
  cx_df$cx_probs <- cx_probs # add observed clonal exclusivity probabilities
  cx_df$n_gp <- n_gp # add number of gene pairs for the sample
  cx_df[, c(ncol(cx_df)-1:0, 3:ncol(cx_df) - 2)] # reorder columns
}

# Test a single gene pair with the LR statistic
test_gene_pair <- function(cx_df, k, test_type = "placement", exact_limit = 12, filter = TRUE, MC_reps = NULL, nu = 10) {
  gene_pair_df <- cx_df[, c(-1, 0, k) + 2] # have probabilities in first column
  probs_all <- gene_pair_df[, 1] # extract the probabilities
  ns <- gene_pair_df[, 2] # and number of gene pairs
  obs_all <- which(!is.na(gene_pair_df[, 3])) # extract the observations (weighted)
  # find the informative ones
  to_keep <- obs_all[which(probs_all[obs_all] > 0 & probs_all[obs_all] < 1)]
  probs <- probs_all[to_keep]
  obs <- gene_pair_df[to_keep, 3]
  # count linear trees
  n_l <- sum(probs_all[obs_all] == 0)
  # count star trees
  n_s <- sum(probs_all[obs_all] == 1)
  n_p <- length(obs_all)
  n_obs <- length(obs)
  n_cx <- sum(obs)

  if (test_type == "occurrence") {
    ll_max <- optimize(ll_w, interval = c(-10, 10), ns = ns, probs = probs_all, obs = obs_all, maximum = TRUE)
  } else if (test_type == "combined") {
    ll_max <- optimize(ll_c, interval = c(-10, 10), ns = ns, probs = probs_all, obs = obs_all, probs_ga = probs, obs_ga = obs, maximum = TRUE)
  } else if (test_type == "placement") {
    ll_max <- optimize(ll_delta, interval = c(-10, 10), probs = probs, obs = obs, maximum = TRUE)
  }

  LR <- 2*ll_max$objective
  p_value <- pchisq(LR, df = 1, lower.tail = FALSE) # chi-squared version

  if(test_type == "combined" && exact_limit > 0) { # exact version
    null_probs <- dyn_add_2d(ns, probs_all, n_p)
    tol <- 1e-6
    null_stats <- rep(NA, n_p+1)
    for (ii in 0:n_p) {
      obs_vec <- pick_pats(probs_all, n_p, ii)
      obs_all <- which(!is.na(obs_vec))
      # find the informative ones
      to_keep <- obs_all[which(probs_all[obs_all] > 0 & probs_all[obs_all] < 1)]
      probs <- probs_all[to_keep]
      obs <- obs_vec[to_keep]
      if (length(obs_all) == n_p) {# if we could pick the right number of patients
        null_stats[ii+1] <- optimize(ll_c, interval = c(-10, 10), ns = ns, probs = probs_all, obs = obs_all, probs_ga = probs, obs_ga = obs, maximum = TRUE)$objective
      } else { # remove the case
        null_stats[ii+1] <- 0
      }
    }
    more_extreme <- which(null_stats >= ll_max$objective + 2*tol)
    as_extreme <- which(null_stats >= ll_max$objective - 2*tol)
    p_more_extreme <- sum(null_probs[more_extreme]) # those clearly above
    p_as_extreme <- sum(null_probs[as_extreme]) # those the same and above
    p_value <- mean(c(p_more_extreme, p_as_extreme)) # taking the average count the same with weight half
    if (filter) { # remove from testing if no chance of significance
      most_extreme <- which(null_stats >= max(null_stats) - 2*tol) # find the most extreme values possible
      p_most_extreme <- sum(null_probs[most_extreme]) / 2 # find middle point
      if (p_most_extreme > 0.05) {
        p_value <- NA
      }
    }
  }

  if (test_type == "placement") {
  if (n_obs < exact_limit || !is.null(MC_reps)) {
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
  }
  # output as data frame
  data.frame(gene_pair = colnames(gene_pair_df)[3], n_p = n_p,
             n = n_obs, n_cx = n_cx, n_l = n_l, n_s = n_s,
             delta = ll_max$maximum, LR = LR, p = p_value)
}

#' Takes a data frame of clonal exclusivity patterns with rows representing
#' patient sample, a first column giving the overall observed exclusivity
#' probability a second column of the number of gene pairs in that sample
#' and subsequent columns providing the rate of clonal
#' exclusivity of the tree(s) of that patient sample.
#'
#' @title Perform the geneAccord clonal exclusivity test.
#' @param clonal_exclusivity_df A data frame with one row per sample with
#' the first column containing the exclusivity rates, the second column
#' containing the number of gene pairs and following columns for each pair
#' of genes recording the average (over trees) clonal exclusivity.
#' @param test_type The type of test to perform with "placement" (default)
#' testing the location of each gene pair in the informative trees, "occurrence"
#' testing the selection of patients exhibiting the gene pair, and "combined"
#' combining both tests. The "occurrence" test uses the chi-squared approximation
#' to obtain p-values, while the "combined" test is exact unless exact_limit is
#' negative when the chi-squared approximation is used. The "placement" test
#' takes the following additional parameters:
#' @param exact_limit Limit above which the asymptotic chi-squared distribution
#' is employed instead of the exact test which is more computationally intensive.
#' @param filter Whether to remove untestable gene pairs (TRUE, default) to
#' reduce multiple testing
#' @param MC_reps How many repetitions to use for Monte Carlo estimation of
#' the p-values. When NULL (default) exact or chi-squared testing is used.
#' @author Jack Kuipers
#' @export
#' @return A data frame with the gene pairs, the total number of observed samples
#' containing the gene pair, the number of informative topologies, the number of
#' samples with clonal exclusivity, the number of linear and star topologies,
#' the maximum likelihood estimate of delta, the LR test statistic of the clonal
#' exclusivity test and the corresponding p-value and BH corrected q-value.
geneAccord <- function(clonal_exclusivity_df, test_type = "placement", exact_limit = 12, filter = TRUE, MC_reps = NULL) {
  cx_results <- NULL
  for (kk in 3:ncol(clonal_exclusivity_df) - 2) {
    cx_results <- rbind(cx_results, test_gene_pair(clonal_exclusivity_df, kk, test_type = test_type, exact_limit = exact_limit, filter = filter, MC_reps = MC_reps))
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

# Compute the difference in log-likelihood as delta is varied
ll_delta <- function(delta, probs, obs) {
  - sum(obs*log(probs + (1 - probs)*exp(delta)) + (1 - obs)*log(1 - probs + probs*exp(-delta)))
}

# Compute the difference in log-likelihood over selection of patient samples
ll_w <- function(delta, ns, probs, obs) {
  null_num <- sum(log(ns[obs]))
  null_dem <- dyn_add(ns, length(obs))
  p <- exp(delta)/(1+exp(delta))
  ws <- 2*ns*((1-probs)*p + probs*(1-p)) # new weights
  alt_num <- sum(log(ws[obs]))
  alt_dem <- dyn_add(ws, length(obs))
  alt_num - null_num - alt_dem + null_dem
}

# Compute the difference in log-likelihood over selection of patient samples
# and over the observed placement of those pairs
ll_c <- function(delta, ns, probs, obs, probs_ga, obs_ga) {
  ll_w(delta, ns, probs, obs) + ll_delta(delta, probs_ga, obs_ga)
}

# Combine weights to get normalising constant of all sets of
# size k out of the list of weights in w
dyn_add <- function(w, k, rescale = TRUE, log = TRUE) {
  if (rescale) {
    # rescale so expect around 1 at end
    scaley <- mean(w)*exp(lchoose(length(w), k)/k)
  } else {
    scaley <- 1
  }
  sums <- rep(0, k) # compute cases with dynamic programming
  for (ii in 1:length(w)) {
    sums <- sums + c(1, sums[-k])*w[ii]/scaley
  }
  if (log) { # return log value
    log(sums[k]) + k*log(scaley)
  } else { # or actual value
    sums[k]*scaley^k
  }
}

# Combine weights and probabilities to get different numbers of observed
# exclusivity patterns among sets of size k out of the list of weights in w
dyn_add_2d <- function(w, p, k, rescale = TRUE) {
  if (rescale) {
    # rescale so expect around 1 at end
    scaley <- mean(w)*exp(lchoose(length(w), k)/k)
  } else {
    scaley <- 1
  }
  sums <- matrix(0, k, k+1) # compute cases with dynamic programming
  for (ii in 1:length(w)) {
    sums <- sums + rbind(c(1, rep(0, k)), sums[-k,])*w[ii]*(1-p[ii])/scaley +
      cbind(rep(0, k), rbind(c(1, rep(0, k)), sums[-k,])[,-(k+1)])*w[ii]*p[ii]/scaley
  }
  log_probs <- log(sums[k,]) + k*log(scaley)
  unnorm_probs <- exp(log_probs - max(log_probs))
  unnorm_probs/sum(unnorm_probs)
}

# This function picks k patients with c exclusivity events
# as far as possible
pick_pats <- function(probs, k, c) {
  k_done <- 0
  c_done <- 0
  kmc <- k - c
  kmc_done <- 0
  obs <- rep(NA, length(probs))
  ii <- 0
  while (k_done < k && ii < length(probs)) {
    ii <- ii + 1
    if (probs[ii] == 0) {
      if (kmc_done < kmc) { # add linear tree
        obs[ii] <- 0
        kmc_done <- kmc_done + 1
        k_done <- k_done + 1
      }
    } else if (probs[ii] == 1) {
      if (c_done < c) { # add star tree
        obs[ii] <- 1
        c_done <- c_done + 1
        k_done <- k_done + 1
      }
    } else {
      if ((c - c_done) > (kmc - kmc_done)) {
        obs[ii] <- 1
        c_done <- c_done + 1
        k_done <- k_done + 1
      } else {
        obs[ii] <- 0
        kmc_done <- kmc_done + 1
        k_done <- k_done + 1
      }
    }
  }
  obs
}
