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
ecdf_list <- generate_ecdf_test_stat(avg_rates_m, list_of_num_trees_all_pats, list_of_clon_excl_all_pats, num_pat_pair_max, num_pairs_sim)


test_that("GeneAccord", {
  alternative <- "greater"
  this_res <- GeneAccord(clone_tbl, avg_rates_m, ecdf_list, alternative)
  assert_that(dplyr::is.tbl(this_res))
  assert_that("pval" %in% colnames(this_res))
  assert_that("entity_A" %in% colnames(this_res))
  assert_that("entity_B" %in% colnames(this_res))
  assert_that(this_res$entity_A == "geneA")
  assert_that(as.numeric(as.character(this_res$pval)) <= 1.0)
  assert_that(as.numeric(dplyr::select(this_res, mle_delta)[1,]) > 0)
  
  alternative <- "two.sided"
  this_res <- GeneAccord(clone_tbl, avg_rates_m, ecdf_list, alternative)
  assert_that(dim(this_res)[2] == 7)
  assert_that(dim(this_res)[1] == 3)
  
  genes_of_interest <- c("geneB", "geneC")
  this_res <- GeneAccord(clone_tbl, avg_rates_m, ecdf_list, alternative, genes_of_interest)
  assert_that(dim(this_res)[2] == 7)
  assert_that(dim(this_res)[1] == 3)
  
  AND_OR = "AND"
  this_res <- GeneAccord(clone_tbl, avg_rates_m, ecdf_list, alternative, genes_of_interest, AND_OR)
  assert_that(dim(this_res)[1] == 1)
  assert_that(dim(this_res)[2] == 7)
  assert_that(as.numeric(dplyr::select(this_res, pval)[1,]) <= 1)
  assert_that(as.numeric(dplyr::select(this_res, pval)[1,]) >= 0)
  
  alternative <- "greater"
  genes_of_interest <- "geneC"
  AND_OR <- "AND"
  this_res <- GeneAccord(clone_tbl, avg_rates_m, ecdf_list, alternative, genes_of_interest, AND_OR)
  assert_that(dim(this_res)[1] == 0)
  assert_that(dim(this_res)[2] == 7)
  
  
  expect_error(suppressMessages(GeneAccord(clone_tbl, avg_rates_m, ecdf_list, "nonsense", genes_of_interest, AND_OR)))
  expect_warning(suppressMessages(GeneAccord(clone_tbl, avg_rates_m, ecdf_list, alternative, "aGeneNotFromTheTibble", AND_OR)))
  expect_error(suppressMessages(GeneAccord(clone_tbl, avg_rates_m, ecdf_list, alternative, genes_of_interest, "nonsense")))
 
  genes_of_interest <- c("geneA", "geneB")
  this_res <- GeneAccord(clone_tbl, avg_rates_m, ecdf_list, alternative, genes_of_interest, "AND")
  assert_that(dim(this_res)[1] == 1)
  assert_that(dim(this_res)[2] == 7)
  alternative <- "two.sided"
  this_res <- GeneAccord(clone_tbl, avg_rates_m, ecdf_list, alternative, genes_of_interest, "OR")
  assert_that(dim(this_res)[1] == 3)
  assert_that(dim(this_res)[2] == 7)
})





