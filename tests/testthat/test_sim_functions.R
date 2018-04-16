

clone_tbl <- tibble::as_tibble(cbind("file_name" =
                                       rep(c(rep(c("fn1", "fn2"), each = 3)), 2),
                                     "patient_id" = rep(c(rep(c("pat1", "pat2"), each = 3)), 2),
                                     "altered_entity" = c(rep(c("geneA", "geneB", "geneC"), 4)),
                                     "clone1" = c(0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0),
                                     "clone2" = c(1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1),
                                     "tree_id" = c(rep(5, 6), rep(10, 6))))
clone_tbl_pat1 <- dplyr::filter(clone_tbl, patient_id == "pat1")
clone_tbl_pat2 <- dplyr::filter(clone_tbl, patient_id == "pat2")
rates_exmpl_1 <- compute_rates_clon_excl(clone_tbl_pat1)
rates_exmpl_2 <- compute_rates_clon_excl(clone_tbl_pat2)
avg_rates_m <- apply(cbind(rates_exmpl_1, rates_exmpl_2), 2, mean)
names(avg_rates_m) <- c(names(rates_exmpl_1)[1], names(rates_exmpl_2)[1])
values_clon_excl_num_trees_pat1 <- get_hist_clon_excl(clone_tbl_pat1)
values_clon_excl_num_trees_pat2 <- get_hist_clon_excl(clone_tbl_pat2)



test_that("generate_ecdf_test_stat", {
  list_of_num_trees_all_pats <- list(pat1=values_clon_excl_num_trees_pat1[[1]], 
                                    pat2=values_clon_excl_num_trees_pat2[[1]])
  list_of_clon_excl_all_pats <- list(pat1=values_clon_excl_num_trees_pat1[[2]],
                                    pat2=values_clon_excl_num_trees_pat2[[2]])
  num_pat_pair_max <- 2
  num_pairs_sim <- 10
  ecdf_list <- generate_ecdf_test_stat(avg_rates_m, list_of_num_trees_all_pats, list_of_clon_excl_all_pats, num_pat_pair_max, num_pairs_sim)
  assert_that(is.null(ecdf_list[[1]]))
  assert_that(length(ecdf_list) == 2)
  assert_that(is.list(ecdf_list))
  expect_error(generate_ecdf_test_stat(rev(avg_rates_m), list_of_num_trees_all_pats, list_of_clon_excl_all_pats, 2, 10))
})


test_that("generate_test_stat_hist", {
  list_of_num_trees_all_pats <- list(pat1=values_clon_excl_num_trees_pat1[[1]], 
                                    pat2=values_clon_excl_num_trees_pat2[[1]])
  list_of_clon_excl_all_pats <- list(pat1=values_clon_excl_num_trees_pat1[[2]],
                                    pat2=values_clon_excl_num_trees_pat2[[2]])
  num_pat_pair_max <- 2
  num_pairs_sim <- 10
  ecdf_list <- generate_ecdf_test_stat(avg_rates_m, list_of_num_trees_all_pats, list_of_clon_excl_all_pats, num_pat_pair_max, num_pairs_sim)
  sim_res <- generate_test_stat_hist(avg_rates_m, list_of_num_trees_all_pats, list_of_clon_excl_all_pats,
             ecdf_list, num_pat_pair_max, num_pairs_sim)
  assert_that(is.null(sim_res[[1]]))
  assert_that(dim(sim_res[[2]])[1] == num_pairs_sim)
  assert_that(dim(sim_res[[2]])[2] == 12)
  assert_that("pat1_beta_distorted" %in% colnames(sim_res[[2]]))
  assert_that("pval" %in% colnames(sim_res[[2]]))
  assert_that("pat2_beta_distorted" %in% colnames(sim_res[[2]]))
  expect_warning(generate_test_stat_hist(avg_rates_m, list_of_num_trees_all_pats, list_of_clon_excl_all_pats, ecdf_list, 2, 10, 1))
  expect_error(generate_test_stat_hist(rev(avg_rates_m), list_of_num_trees_all_pats, list_of_clon_excl_all_pats, ecdf_list, 2, 10))
})


test_that("build_null_test_statistic", {
  
  avg_rates_m = c(0.4, 0.3)
  list_of_clon_excl_frac_trees_all_pats <- list(list(c(5, 4, 5), c(5, 4)), list(c(4, 4, 3), c(3, 2)))
  sim_pairs <- build_null_test_statistic(avg_rates_m,
              list_of_clon_excl_frac_trees_all_pats, 2, 100, beta_distortion = 100)
  assert_that(dim(sim_pairs)[2] == 10)
  assert_that(dim(sim_pairs)[1] == 100)
  assert_that("test_statistic" %in% colnames(sim_pairs))
})

