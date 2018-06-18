
ext_data_dir <- system.file('extdata', package = 'GeneAccord')
data("all_genes_tbl")
test_patient <- "01"
input_files_test <- paste(ext_data_dir, "/clonal_genotypes/cloe_seed", seq(5, 100, by = 5), "/", test_patient, ".csv", sep = "") 
small_clone_tbl <- dplyr::as.tbl( as.data.frame( cbind( altered_entity = c(paste("gene", seq(1,5), sep = "")),
                                                     clone1 = c(rep(0,5)),
                                                     clone2 = c(0,1,1,0,1),
                                                     clone3 = c(1,0,1,1,0),
                                                     clone4 = c(0,1,1,0,0) ) ))
test_clone_tbl <- tibble::tibble(file_name = c(paste("fn", rep(seq(1,2), each = 5), sep = ""), "fn1"),
                               patient_id = c(paste("pat", rep(seq(1,2), each = 5), sep = ""), "pat1"),
                               altered_entity = c(paste("gene", rep(seq(1,5), 2), sep = ""), "gene1"),
                               clone1 = c(1,rep(0,10)),
                               clone2 = c(0,1,1,0,1,0,1,1,0,1,1),
                               clone3 = c(1,0,1,1,0,1,1,0,1,0,1),
                               clone4 = c(0,1,1,0,0,0,1,1,0,0,0) )
test_clone_tbl2 <- dplyr::as.tbl( as.data.frame( cbind(
                          altered_entity = c(paste("gene", seq(1,10), sep = "")),
                          clone1 = c(rep(0,10)),
                          clone2 = c(sample(c(0,1), 10, replace = TRUE)),
                          clone3 = c(sample(c(0,1), 10, replace = TRUE)),
                          clone4 = c(sample(c(0,1), 10, replace = TRUE)),
                          tree_id = c(rep(5, 5), rep(10, 5)) ) ))
test_clone_tbl3 <- dplyr::as.tbl( as.data.frame( cbind(
  file_name = rep("fn1", 10),
  patient_id = rep("pat1", 10),
  altered_entity = c(paste("gene", seq(1,10), sep = "")),
  clone1 = c(rep(0,10)),
  clone2 = c(sample(c(0,1), 10, replace = TRUE)),
  clone3 = c(sample(c(0,1), 10, replace = TRUE)),
  clone4 = c(sample(c(0,1), 10, replace = TRUE)),
  tree_id = c(rep(5, 5), rep(10, 5)) ) ))

test_clone_tbl4 <- tibble::as_tibble(cbind("file_name" =
                                                rep(c(rep(c("fn1", "fn2"), each = 3)), 2),
                                    "patient_id" = rep(c(rep(c("pat1", "pat2"), each = 3)), 2),
                                    "altered_entity" = c(rep(c("geneA", "geneB", "geneC"), 4)),
                                    "clone1" = c(0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0),
                                    "clone2" = c(1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1),
                                    "tree_id" = c(rep(5, 6), rep(10, 6))))

test_that("create_tbl_ent_clones", {
  test_clone_tbl <- create_tbl_ent_clones(input_files_test[1])
  expect_that(test_clone_tbl, dplyr::is.tbl)
  assert_that("file_name" %in% colnames(test_clone_tbl))
  assert_that("patient_id" %in% colnames(test_clone_tbl))
  assert_that("altered_entity" %in% colnames(test_clone_tbl))
  assert_that("clone1" %in% colnames(test_clone_tbl))
  assert_that("clone2" %in% colnames(test_clone_tbl))
  # to make sure that the tables are read in as expected
  gene_to_num_clones_P01 <- test_clone_tbl %>% dplyr::select(-file_name, -patient_id) %>% dplyr::mutate(sum=Reduce("+",.[2:8])) %>%
     dplyr::select(altered_entity, sum)
  assert_that(gene_to_num_clones_P01 %>% dplyr::filter(altered_entity == "ENSG00000134086") %>%
                dplyr::select(sum) == 3)
  assert_that(gene_to_num_clones_P01 %>% dplyr::filter(altered_entity == "ENSG00000150471") %>%
                dplyr::select(sum) == 2)
  assert_that(gene_to_num_clones_P01 %>% dplyr::filter(altered_entity == "ENSG00000126012") %>%
                dplyr::select(sum) == 1)
})


test_that("create_tbl_tree_collection", {
  test_tbl <- create_tbl_tree_collection(input_files_test)
  expect_that(test_tbl, dplyr::is.tbl)
  assert_that("patient_id" %in% colnames(test_tbl))
  assert_that("altered_entity" %in% colnames(test_tbl))
  assert_that("clone1" %in% colnames(test_tbl))
  assert_that("clone2" %in% colnames(test_tbl))
  assert_that("clone3" %in% colnames(test_tbl))
  assert_that("tree_id" %in% colnames(test_tbl))
  assert_that(dim(test_tbl)[1] == 580)
  assert_that(dim(test_tbl)[2] == 11)
})


test_that("compute_rates_clon_excl", {  # get_all_rates_diff_branch_ent_pair
  expect_error(compute_rates_clon_excl(test_clone_tbl))
  rates <- compute_rates_clon_excl(test_clone_tbl3)
  expect_that(rates, is.numeric)
  assert_that(length(rates) == 2)
  assert_that(rates[1] >= 0)
  assert_that(rates[2] >= 0)
  assert_that(rates[1] <= 1)
  assert_that(rates[2] <= 1)
  assert_that(length(rates) == 2)
})


test_that("get_rate_diff_branch_ent_pair", {
  assert_that(suppressMessages(get_rate_diff_branch_ent_pair(small_clone_tbl)) == 0.4)
  test_clone_tbl <- tibble::tibble(altered_entity = c("testGene"), clone1 = c(1), clone2 = c(0))
  assert_that(suppressMessages(get_rate_diff_branch_ent_pair(test_clone_tbl)) == 0)
})


test_that("is_diff_branch_ent_pair", {
  expect_true(suppressMessages(is_diff_branch_ent_pair("gene1", "gene2", small_clone_tbl)))
  expect_false(suppressMessages(is_diff_branch_ent_pair("gene1", "gene3", small_clone_tbl)))
  expect_false(suppressMessages(is_diff_branch_ent_pair("gene4", "gene3", small_clone_tbl)))
  expect_true(suppressMessages(is_diff_branch_ent_pair("gene5", "gene4", small_clone_tbl)))
})

 
test_that("get_hist_clon_excl", { # get_hist_subcl_excl_this_pat
  test_tbl <- tibble::tibble(altered_entity=paste("gene", seq(1, 4)),
                             clone1=c(1,0,0,1),
                             clone2=c(1,0,1,0),
                             clone3=c(0,1,0,0),
                             clone4=c(0,1,0,0))
  expect_error(get_hist_clon_excl(test_tbl))
  test_tbl <- test_tbl %>% dplyr::mutate(tree_id=c(1,2,3,4))
  this_res <- get_hist_clon_excl(test_tbl)
  assert_that(this_res[[1]] == 0)
  assert_that(this_res[[2]] == 0)
  test_tbl$tree_id <- c(1,1,2,2)
  this_res <- get_hist_clon_excl(test_tbl)
  this_res <- get_hist_clon_excl(test_clone_tbl3)
  assert_that(is.list(this_res))
  assert_that(length(this_res) == 2)
  assert_that(length(this_res[[1]]) == 20)
})


test_that("get_hist_clon_excl_this_pat_this_pair", {
  test_tbl <- test_clone_tbl3 %>% dplyr::select(-patient_id, -file_name)
  this_res <- get_hist_clon_excl_this_pat_this_pair("gene1", "gene2", test_tbl)
  assert_that(length(this_res) == 2)
  assert_that(sum(this_res) >= 0)
})


test_that("merge_clones_identical_ents", {
  clone_tbl <- dplyr::as.tbl(as.data.frame(cbind(
                          file_name = c(rep("fn1", 4)),
                          patient_id = c(rep("pat1", 4)),
                          altered_entity = c("pw1", "pw1", "pw2", "pw3"),
                          clone1 = c(1, 0, 1, 0),
                          clone2 = c(0, 1, 0, 1),
                          clone3 = c(1, 1, 0, 1),
                          clone4 = c(0, 1, 0, 0))))
  merged_identical_tbl <- merge_clones_identical_ents(clone_tbl)
  expect_that(merged_identical_tbl, dplyr::is.tbl)
  num_entries_merged_tbl <- dim(merged_identical_tbl)[1]
  assert_that(dim(clone_tbl)[1] >= num_entries_merged_tbl)
  unique_ents_tbl <- clone_tbl %>% dplyr::select(altered_entity) %>% dplyr::distinct()
  assert_that(dim(unique_ents_tbl)[1] == num_entries_merged_tbl)
})



test_that("extract_num_clones_tbl", {
  num_clones <- suppressMessages(extract_num_clones_tbl(test_clone_tbl))
  pat_ids <- names(num_clones)
  assert_that(num_clones[which(pat_ids == "pat1")] == 4)
  assert_that(num_clones[which(pat_ids == "pat2")] == 3)
  assert_that(length(num_clones) == 2)
})


test_that("take_pairs_and_get_patients", {
  small_clone_tbl <- tibble::tibble(file_name =
                                          rep(c(rep(c("fn1", "fn2"), each = 3)), 2),
                                    patient_id = rep(c(rep(c("pat1", "pat2"), each = 3)), 2),
                                    altered_entity = c(rep(c("geneA", "geneB", "geneC"), 4)),
                                    clone1 = c(0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0),
                                    clone2 = c(1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1),
                                    tree_id = c(rep(5, 6), rep(10, 6)))
  pairs_of_interest <- tibble::tibble(entity_A = c("geneA", "geneB"), entity_B = c("geneB", "geneC"))
  this_res <- take_pairs_and_get_patients(small_clone_tbl, pairs_of_interest)
  assert_that(dplyr::is.tbl(this_res))
  assert_that("mutated_in" %in% colnames(this_res))
  assert_that("clonally_exclusive_in" %in% colnames(this_res))
  assert_that(as.character(this_res %>% dplyr::filter(entity_A == "geneB") %>% dplyr::select(clonally_exclusive_in)) == ".")
  assert_that(as.character(this_res %>% dplyr::filter(entity_A == "geneA") %>% dplyr::select(clonally_exclusive_in)) == "pat1;pat2")
})

test_that("map_pairs_to_hgnc_symbols", {
   pairs_of_interest <- tibble::tibble(entity_A = c("ENSG00000181143", "ENSG00000163939"), entity_B = c("ENSG00000141510", "ENSG00000163930"))
   this_res <- map_pairs_to_hgnc_symbols(pairs_of_interest, all_genes_tbl)
   assert_that(dplyr::is.tbl(this_res))
   assert_that(dim(this_res)[1] == 2)
   assert_that(dim(this_res)[2] == 4)
   assert_that(this_res[1,4] == "TP53")
   assert_that("hgnc_gene_A" %in% colnames(this_res))
   assert_that("hgnc_gene_B" %in% colnames(this_res))
})

test_that("write_res_pairs_to_disk", {
  sig_pairs <- tibble::tibble(hgnc_gene_A = c("VHL", "BAP1"),
                               hgnc_gene_B = c("PTEN", "KIT"),
                               pval = c(0.001, 0.002),
                               qval = c(0.01, 0.02),
                               mutated_in = c("pat1; pat2", "pat1; pat2"),
                               clonally_exclusive_in = c("pat1; pat2", "pat2"))
  avg_rates_m <- c(pat1 = 0.0034, pat2 = 0.0021)
  write_res_pairs_to_disk(sig_pairs, avg_rates_m, "test.tsv")
  file.remove("test.tsv")
})

test_that("pairs_in_patients_hist", {
  clone_tbl <- tibble::as_tibble(cbind("file_name" =
                                    rep(c(rep(c("fn1", "fn2"), each = 3)), 2),
                                    "patient_id" = rep(c(rep(c("pat1", "pat2"), each = 3)), 2),
                                    "altered_entity" = c(rep(c("geneA", "geneB", "geneC"), 4)),
                                    "clone1" = c(0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0),
                                    "clone2" = c(1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1),
                                    "tree_id" = c(rep(5, 6), rep(10, 6))))
  this_res <- pairs_in_patients_hist(clone_tbl)
  assert_that("pairs_count" %in% colnames(this_res))
  assert_that("patient_count" %in% colnames(this_res))
  expect_equal(this_res$patient_count, 2)
  expect_equal(this_res$pairs_count, 3)
  
  clone_tbl <- tibble::as_tibble(cbind("file_name" =
                                         rep(c(rep(c("fn1", "fn2"), each = 3)), 2),
                                       "patient_id" = rep(c(rep(c("pat1", "pat2"), each = 3)), 2),
                                       "altered_entity" = c("geneA", "geneB", "geneC",
                                                            "geneC", "geneB", "geneA",
                                                            "geneC", "geneB", "geneA",
                                                            "geneA", "geneB", "geneC"),
                                       "clone1" = c(0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0),
                                       "clone2" = c(1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1),
                                       "tree_id" = c(rep(5, 6), rep(10, 6))))
  this_res <- pairs_in_patients_hist(clone_tbl)
  assert_that("pairs_count" %in% colnames(this_res))
  assert_that("patient_count" %in% colnames(this_res))
  expect_equal(this_res$patient_count, 2)
  expect_equal(this_res$pairs_count, 3)
})
