data("all_genes_tbl")

test_that("plot_rates_clon_excl",{
  clone_tbl <- tibble::as_tibble(cbind("file_name" =
                                        rep(c(rep(c("fn1", "fn2"), each = 3)), 2),
                                    "patient_id" = rep(c(rep(c("pat1", "pat2"), each = 3)), 2),
                                    "altered_entity" = c(rep(c("geneA", "geneB", "geneC"), 4)),
                                    "clone1" = c(0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0),
                                    "clone2" = c(1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1),
                                    "tree_id" = c(rep(1, 6), rep(2, 6))))
  avg_rates_m <- c(pat1 = 0.014, pat2 = 0.3)
  plot_rates_clon_excl(avg_rates_m, clone_tbl, "test_plot_rates_clon_excl.pdf")
  file.remove("test_plot_rates_clon_excl.pdf")
})



test_that("plot_ecdf_test_stat", {
  avg_rates_m <- c(pat1 = 0.1, pat2 = 0.034, pat3 = 0.21, pat4 = 0.063)
  list_of_num_trees_all_pats <- list(pat1 = c(20, 20, 19), pat2 = c(20, 18, 20),
                                    pat3 = c(19, 20, 20), pat4 = c(20, 20, 20))
  list_of_clon_excl_all_pats <- list(pat1 = c(5, 0, 1), pat2 = c(10, 2, 0),
                                    pat3 = c(18, 12, 0), pat4 = c(0, 2, 0))
  num_pat_pair_max <- 2
  num_pairs_sim <- 10
  ecdf_list <- generate_ecdf_test_stat(avg_rates_m, list_of_num_trees_all_pats, 
                                      list_of_clon_excl_all_pats, num_pat_pair_max, num_pairs_sim)
  plot_ecdf_test_stat(ecdf_list, plot_idx = 2, output_pdf = "test_plot_ecdf_test_stat.pdf")
  file.remove("test_plot_ecdf_test_stat.pdf")
  expect_error(suppressMessages(plot_ecdf_test_stat(ecdf_list, 3)))
})



test_that("vis_pval_distr_num_pat",{
  out_sim_pairs <- tibble::tibble(
                              num_patients = c(rep(2,100), rep(3,100), rep(4,100)),
                              pval = c(runif(300)))
  vis_pval_distr_num_pat(out_sim_pairs, "test_vis_pval_distr_num_pat.pdf")
  file.remove("test_vis_pval_distr_num_pat.pdf")
  wrong_tbl <- out_sim_pairs %>% dplyr::select(num_patients)
  expect_error(suppressMessages(vis_pval_distr_num_pat(wrong_tbl)))
})



test_that("heatmap_clones_gene_pat", {
  sig_pairs <- tibble::tibble(entity_A = "SETD2",
                            entity_B = "BAP1")
  clone_tbl <- tibble::tibble(file_name = c("X.csv", "X.csv", "Y.csv", "Y.csv"),
                            patient_id = c("X", "X", "Y", "Y"),
                            altered_entity = c("SETD2", "BAP1", "SETD2", "BAP1"),
                            clone1 = c(0, 1, 1, 0),
                            clone2 = c(1, 0, 0, 1))
  heatmap_clones_gene_pat(sig_pairs, clone_tbl, all_genes_tbl, output_pdf = "test_heatmap_subcl_gene_pat_set.pdf")
  file.remove("test_heatmap_subcl_gene_pat_set.pdf")
})



