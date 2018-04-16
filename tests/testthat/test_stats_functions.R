clone_tbl <- dplyr::as.tbl( as.data.frame( cbind( altered_entity = c(paste("gene", seq(1,5), sep = "")),
                                                     clone1 = c(rep(0,5)),
                                                     clone2 = c(0,1,1,0,1),
                                                     clone3 = c(1,0,1,1,0),
                                                     clone4 = c(0,1,1,0,0) ) ))


test_that("compute_test_stat_avg_rate", {
  this_res <- compute_test_stat_avg_rate(c(0.05, 0.23), c(20, 20), c(8, 5))
  assert_that(is.numeric(this_res[[1]]))
  assert_that(length(this_res) == 2)
  assert_that(is.numeric(this_res[[2]]))
  this_res <- compute_test_stat_avg_rate(c(0.1, 0.2), c(10, 10), c(9, 7))
  assert_that(this_res[[2]] > 0)
  assert_that(this_res[[1]] > 2)
  expect_error(compute_test_stat_avg_rate(c(0.05, 0.23), c(20, 20), c(21, 22)))
  expect_error(compute_test_stat_avg_rate(c(20, 20), c(21, 22), c(0.05, 0.23)))
  expect_error(compute_test_stat_avg_rate(0, 20, 20))
})


test_that("ecdf_lr_test_clon_excl_avg_rate", {
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
  list_of_num_trees_all_pats <- list(pat1=values_clon_excl_num_trees_pat1[[1]],
                                     pat2=values_clon_excl_num_trees_pat2[[1]])
  list_of_clon_excl_all_pats <- list(pat1=values_clon_excl_num_trees_pat1[[2]],
                                     pat2=values_clon_excl_num_trees_pat2[[2]])
  num_pat_pair_max <- 2
  num_pairs_sim <- 10
  ecdf_list <- generate_ecdf_test_stat(avg_rates_m, list_of_num_trees_all_pats,
                                 list_of_clon_excl_all_pats, num_pat_pair_max, num_pairs_sim)
  this_res <- ecdf_lr_test_clon_excl_avg_rate("geneA", "geneB", clone_tbl, avg_rates_m,
                                    ecdf_list, "greater")
  
  assert_that(length(this_res) == 4)
  assert_that(this_res[[1]] >= 0 && this_res[[1]] <= 1)
  assert_that(this_res[[3]] > 0)
  this_res<-ecdf_lr_test_clon_excl_avg_rate("geneA", "geneC", clone_tbl, avg_rates_m, ecdf_list, "two.sided")
  assert_that(this_res[[3]] < 0)
  small_clone_tbl <- tibble::as_tibble(cbind("file_name" = rep(c(rep(c("fn1", "fn2"), each = 3)), 2),
                                             "patient_id" = rep(c(rep(c("pat1", "pat2"), each = 3)), 2),
                                             "altered_entity" = c(rep(c("geneA", "geneB", "geneC"), 2), rep(c("geneD", "geneE", "geneF"), 2)),
                                             "clone1" = c(0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0),
                                             "clone2" = c(1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1),
                                             "tree_id" = c(rep(5, 6), rep(10, 6))))
  this_res<-ecdf_lr_test_clon_excl_avg_rate("geneA", "geneB", small_clone_tbl, avg_rates_m, ecdf_list, "two.sided")
  assert_that(this_res[[2]] == 2)
  expect_error(ecdf_lr_test_clon_excl_avg_rate("geneA", "geneC", clone_tbl, avg_rates_m, ecdf_list, "nonsensealternative"))
})

