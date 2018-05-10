#' Different gene id's and their chromosomal location.
#'
#' This is a tibble that contains mappings between different gene identifiers. It can be created with the
#' function \code{\link{create_ensembl_gene_tbl_hg}}. These are the human genes from the human genome
#' version \code{hg19}/\code{GRCh37}, and Ensembl Genes version 88.
#'
#' @format A tibble with 41'393 rows and seven variables:
#' \describe{
#'   \item{ensembl_gene_id}{the Ensembl gene id as a character}
#'   \item{hgnc_symbol}{the HGNC gene symbol as a character}
#'   \item{entrezgene}{the Entrez gene id as an integer}
#'   \item{uniprotswissprot}{the UniProtKB/Swiss-Prot gene id's as a character}
#'   \item{chromosome_name}{the name of the chromosome where the gene is located as a character,
#'                           e.g. "3" for chromosome three}
#'   \item{start_position}{the nucleotide start position of the gene as an integer}
#'   \item{end_position}{the nucleotide end position of the gene as an integer}
#' }
#' @source The tibble can be generated with \code{\link{create_ensembl_gene_tbl_hg}()}, which uses the R-package
#' \code{biomaRt} and the Ensembl data base \url{www.ensembl.org}.
"all_genes_tbl"



#' Ensembl gene id's and the Reactome pathways.
#'
#' This is a tibble that contains mappings between ensembl gene id's and reactome pathways. The reactome pathways
#' are from the lowest level in the hierarchy ("Lowest level pathway diagram / Subset of the pathway"), 
#' and were obtained by download from the Reactome website (\url{https://reactome.org/download-data}; "ENSEMBL to pathways").
#' The following commands were used:
#' wget https://reactome.org/download/current/Ensembl2Reactome.txt; 
#' cat Ensembl2Reactome.txt | grep "Homo sapiens" > Ensembl2Reactome_homo_sapiens.txt
#' 
#' @format A tibble with 46'141 rows and six variables:
#' \describe{
#'   \item{ensembl_gene_id}{the Ensembl gene id as a character}
#'   \item{reactome_pw_id}{the Reactome pathway stable identifier}
#'   \item{url}{The url leading to the pathway graph}
#'   \item{reactome_pw_name}{the name of the Reactome pathway}
#'   \item{evidence_code}{the evidence code}
#'   \item{species}{the species}
#' }
#' @source The tibble was created as follows: 
#' library(dplyr);
#' ensg_path_map_raw <- read.csv("Ensembl2Reactome_homo_sapiens.txt", header = F, sep = "\\t", 
#'     comment.char = "", check.names = F, skip = 0);
#' stopifnot(dim(ensg_path_map_raw)[1] == 46141);
#' stopifnot(dim(ensg_path_map_raw)[2] == 6);
#' colnames(ensg_path_map_raw) <- c("ensembl_gene_id", "reactome_pw_id", "url", "reactome_pw_name", "evidence_code", "species");
#' ensg_path_map_raw <- dplyr::as.tbl(ensg_path_map_raw);
#' ensg_reactome_path_map <- filter(filter(ensg_path_map_raw, grepl("ENSG", ensg_path_map_raw$ensembl_gene_id)), species == "Homo sapiens")
"ensg_reactome_path_map"


#' The tibble with gene-to-clone assignments from all patients and all trees
#'
#' This is a tibble that contains the information, which gene is mutated in which clone from which patient. 
#' 
#' @format A tibble containing the following columns:
#' \describe{
#'   \item{file_name}{the name of the csv-file from which the data was read}
#'   \item{patient_id}{the patient identifier}
#'   \item{altered_entity}{ensembl gene identifier of the mutated gene}
#'   \item{clone1}{the indication whether the current gene is mutated in this clone}
#'   \item{clone2}{the indication whether the current gene is mutated in this clone}
#'   \item{clone3}{the indication whether the current gene is mutated in this clone}
#'   \item{clone4}{the indication whether the current gene is mutated in this clone}
#'   \item{clone5}{the indication whether the current gene is mutated in this clone}
#'   \item{clone6}{the indication whether the current gene is mutated in this clone}
#'   \item{clone7}{the indication whether the current gene is mutated in this clone}
#'   \item{tree_id}{the identifier that tells from which tree this gene-to-clone assignment comes}
#' }
#'
#' @source The tibble can be generated for each patient separately with \code{\link{create_tbl_tree_collection}} as 
#' demonstrated in the vignette.
"clone_tbl_all_pats_all_trees"


#' The average rates of clonal exclusivity of the example data set used in the vignette
#'
#' This is a named vector that contains the average rate of clonal exclusivity for each 
#' of the 82 patients as described in the vignette. 
#'
#' @format A vector with the average rates of clonal exclusivity of each patient. The names
#' of each element is the respective patient name.
#' @source The rates can be generated for each patient separately with \code{\link{compute_rates_clon_excl}}, and then taking
#' the mean(). This is demonstrated in the vignette.
"avg_rates_m"


#' The list with the histogram of how often pairs are occurring across the collection of trees
#'
#' This is a named list whose entries for each patient are the histograms of how often pairs occurr in all 
#' trees of a patient. 
#'
#' @format A list whose entries are named after the patient, and they contain vectors with the numbers of
#' how often the pairs in this patient occur accross trees.
#' @source The histogram can be generated for each patient separately with \code{\link{get_hist_clon_excl}} as 
#' demonstrated in the vignette.
"list_of_num_trees_all_pats"


#' The list with the histogram of how often pairs are clonally exclusive across the collection of trees
#'
#' This is a named list whose entries for each patient are the histograms of how often pairs are clonally 
#' exclusive in all trees of a patient. 
#'
#' @format A list whose entries are named after the patient, and they contain vectors with the numbers of
#' how often the pairs in this patient are clonally exclusive accross trees.
#' @source The histogram can be generated for each patient separately with \code{\link{get_hist_clon_excl}} as 
#' demonstrated in the vignette.
"list_of_clon_excl_all_pats"

#' The list with the ECDF's of the test statistic under the null hypothesis
#'
#' This is a list whose entries are the empirical cumulative distribution functions for different
#' number of patients that pairs can be mutated in.  
#'
#' @format A list whose entries are the empirical cumulative distribution functions. Entry 1 is set to NULL, because
#' GeneAccord does not test pairs that occur in just one patient. Entry 2 then contains the ECDF of the test statistic
#' under the null hypothesis for the case that pairs are mutated in two patients. Entry 3 contains the ECDF for the case
#' where pairs occur in three patients.
#' @source The list was generated with the function \code{\link{generate_ecdf_test_stat}} as demonstrated in the vignette, just
#' that the following parameter was set as \code{num_pairs_sim = 100000}.
"ecdf_list"


#' Ensembl gene id's and the Reactome pathways - for mouse!
#'
#' This is a tibble that contains mappings between mouse ensembl gene id's and reactome pathways. The reactome pathways
#' are from the lowest level in the hierarchy ("Lowest level pathway diagram / Subset of the pathway"), 
#' and were obtained by download from the Reactome website (\url{https://reactome.org/download-data}; "ENSEMBL to pathways").
#' The following commands were used:
#' wget https://reactome.org/download/current/Ensembl2Reactome.txt; 
#' cat Ensembl2Reactome.txt | grep "Mus musculus" > Ensembl2Reactome_mus_musculus.txt
#' 
#' @format A tibble with 28,630 rows and six variables:
#' \describe{
#'   \item{ensembl_gene_id}{the Ensembl gene id as a character}
#'   \item{reactome_pw_id}{the Reactome pathway stable identifier}
#'   \item{url}{The url leading to the pathway graph}
#'   \item{reactome_pw_name}{the name of the Reactome pathway}
#'   \item{evidence_code}{the evidence code}
#'   \item{species}{the species}
#' }
#' @source The tibble was created as follows: 
#' library(dplyr);
#' ensmusg_path_map_raw <- read.csv("Ensembl2Reactome_mus_musculus.txt", header = F, sep = "\\t", 
#'     comment.char = "", check.names = F, skip = 0);
#' stopifnot(dim(ensmusg_path_map_raw)[1] == 28696);
#' stopifnot(dim(ensmusg_path_map_raw)[2] == 6);
#' colnames(ensmusg_path_map_raw) <- c("ensembl_gene_id", "reactome_pw_id", "url", "reactome_pw_name", "evidence_code", "species");
#' ensmusg_path_map_raw <- dplyr::as.tbl(ensmusg_path_map_raw);
#' ensmusg_reactome_path_map <- filter(filter(ensmusg_path_map_raw, grepl("ENSMUSG", ensmusg_path_map_raw$ensembl_gene_id)), species == "Mus musculus")
#' stopifnot(dim(ensmusg_reactome_path_map)[1] == 28630)
"ensmusg_reactome_path_map"


