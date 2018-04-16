ext_data_dir <- system.file('extdata', package = 'GeneAccord')
data("ensg_reactome_path_map")
data("all_genes_tbl")

test_that("create_ensembl_gene_tbl_hg", {
  expect_that(all_genes_tbl, dplyr::is.tbl)
})


test_that("ensembl_to_hgnc", {
  expect_equal(ensembl_to_hgnc("ENSG00000134086", all_genes_tbl), "VHL")
  expect_error(ensembl_to_hgnc("ENSG00000134086", "this_is_no_tbl"))
  expect_equal(suppressMessages(ensembl_to_hgnc("thisENSEMBLidDoesNotExist", all_genes_tbl)),
               "thisENSEMBLidDoesNotExist")
})

test_that("hgnc_to_ensembl", {
  expect_equal(hgnc_to_ensembl("VHL", all_genes_tbl), "ENSG00000134086")
  expect_equal(hgnc_to_ensembl("PBRM1", all_genes_tbl), "ENSG00000163939")
})

test_that("convert_ensembl_to_reactome_pw_tbl", {
  mutated_gene_tbl <- tibble::tibble(file_name=c("pat1.csv", "pat1.csv"),
  patient_id = c("1","1"),
  altered_entity = c("ENSG00000134086", "ENSG00000163939"),
  clone1 = c(1,0),
  clone2 = c(0,1))
  res <- convert_ensembl_to_reactome_pw_tbl(mutated_gene_tbl, ensg_reactome_path_map)
  assert_that(dplyr::is.tbl(res))
  assert_that(dim(res)[2] == 5)
})


test_that("ensembl_to_reactome", {
  ensg_gene <- "ENSG00000134086"
  res <- ensembl_to_reactome(ensg_gene, ensg_reactome_path_map)
  assert_that(is.character(res))
  res<-ensembl_to_reactome("test", ensg_reactome_path_map)
  assert_that(res == "test")
})
