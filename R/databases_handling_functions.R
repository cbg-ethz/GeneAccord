#' Retrieve a mapping between different gene identifiers.
#'
#' This function retrieves the ensembl gene id's from biomart 
#' together with the hgnc gene symbol, the entrez gene id, the 
#' uniprot/swissprot gene id, as well as chromosome, start and end 
#' position. This is done for the human genes from the human 
#' genome version \code{hg19}/\code{GRCh37}, and Ensembl Genes 
#' version 88. The user can also specify other human genome or 
#' ensembl versions.
#'
#' @title Get a tibble of all gene ensembl id's, gene names (hgnc), 
#' entrez gene id's, uniprot/swissprot
#' gene id's and genomic coordinates.
#' @author Ariane L. Moore
#' @return A tibble with the following columns: 
#' \code{ensembl_gene_id}, \code{hgnc_symbol},
#' \code{entrezgene}, \code{uniprotswissprot}, 
#' \code{chromosome_name}, \code{start_position},
#' \code{end_position}. The entrez gene id, as well as the start 
#' and end positions are numeric, and
#' the other columns are characters. The chromosome is specified 
#' without "chr", i.e. the chromosome
#' 13 for example, would be specified with "13".
#' @param GRCh The human genome version. Default: 37.
#' @param ensembl_version The version of the ensembl data base. Default: 88.
#' @export
#' @importFrom dplyr as.tbl
#' @importFrom biomaRt getBM useEnsembl
#' @examples
#' \dontrun{
#' create_ensembl_gene_tbl_hg()
#' }
create_ensembl_gene_tbl_hg <- function(GRCh=37, ensembl_version=88) {
    message("GRCh version: ", 
        GRCh, ", 
        and Ensembl Genes version: ", 
        ensembl_version)
    message("Obtain a tibble with various gene id's",
    " and their genomic coordinates...")
    all_genes_tbl <- tryCatch(
    {
        ensembl <- 
        useEnsembl(biomart="ensembl", 
            dataset="hsapiens_gene_ensembl",
            GRCh=GRCh, version=ensembl_version)
        message("Retrieve ensembl gene",
            " id's and hgnc symbols from biomart...")
        all_genes <- 
            getBM(attributes=c('ensembl_gene_id', 
            'hgnc_symbol', 
            'entrezgene', 
            'uniprotswissprot',
            'chromosome_name', 
            'start_position',
            'end_position'), 
            mart=ensembl)
        stopifnot(is.data.frame(all_genes))
        stopifnot(dim(all_genes)[2] == 7)
        stopifnot(is.character(all_genes[,1]))
        stopifnot(is.character(all_genes[,2]))
        stopifnot(is.numeric(all_genes[,3]))
        stopifnot(is.character(all_genes[,4]))
        stopifnot(is.character(all_genes[,5]))
        stopifnot(is.numeric(all_genes[,6]))
        stopifnot(is.numeric(all_genes[,7]))
        all_genes_tbl <- as.tbl(all_genes)
        stopifnot(dim(all_genes_tbl)[1] == dim(all_genes)[1])
        stopifnot(dim(all_genes_tbl)[2] == dim(all_genes)[2])
        all_genes_tbl
    },
    error=function(cond) {
        message("BiomaRt function call did not work. Possibly",
        " biomart or ensembl are currently down.")
        message("Here's the original error message:")
        message(cond)
        message("\nLoad 'all_genes_tbl' from the GeneAccord saved data:")
        all_genes_tbl <- GeneAccord::all_genes_tbl
        return(all_genes_tbl)
    },
    warning=function(cond) {
        message("BiomaRt function call caused a warning. Possibly",
        " biomart or ensembl are currently down.")
        message("Here's the original warning message:")
        message(cond)
        message("\nLoad 'all_genes_tbl' from the GeneAccord saved data:")
        all_genes_tbl <- GeneAccord::all_genes_tbl
        return(all_genes_tbl)
    },
    finally={
        message("...done")
    })
    return(all_genes_tbl)
}



#' Map a given ensembl gene id to the hgnc gene symbol.
#'
#' For an ensembl id and a tibble with all genes as input, this 
#' function returns the matching
#' hgnc gene symbol. The tibble with all genes can be generated 
#' with \code{\link{create_ensembl_gene_tbl_hg}}.
#'
#' @title Get the hgnc gene symbol for an ensembl gene id.
#' @param this_ensembl The ensembl id of a gene.
#' @param all_genes_tbl A tibble with all genes ensembl id's 
#' and hgnc symbols.
#' @author Ariane L. Moore
#' @return The matching hgnc gene symbol.
#' @export
#' @importFrom magrittr "%>%"
#' @importFrom dplyr is.tbl filter select distinct
#' pull 
#' @examples
#' \dontrun{
#' all_genes_tbl <- create_ensembl_gene_tbl_hg()
#' ensembl_to_hgnc("ENSG00000134086", all_genes_tbl)
#' ensembl_to_hgnc("ENSG00000141510", all_genes_tbl)
#' }
ensembl_to_hgnc <- function(this_ensembl, all_genes_tbl){
    ensembl_gene_id <- hgnc_symbol <- NULL
    stopifnot(is.character(this_ensembl))
    stopifnot(length(this_ensembl) == 1)
    stopifnot(is.tbl(all_genes_tbl))
    stopifnot("ensembl_gene_id" %in% colnames(all_genes_tbl))
    stopifnot("hgnc_symbol" %in% colnames(all_genes_tbl))

    ## translate ensembl into gene id
    this_hgnc_character <- all_genes_tbl %>% 
        filter(ensembl_gene_id == this_ensembl) %>%
        select(hgnc_symbol) %>%
        distinct() %>% 
        pull(hgnc_symbol)
    
    num_gene_symbols=length(this_hgnc_character)
    
    if (num_gene_symbols == 0){  ## i.e., if there was no such ensembl id
        message("Ensembl id ", this_ensembl," 
        cannot be found in the tibble!")
        return(this_ensembl)
    } else {
        if (num_gene_symbols > 1){
            warning("Ensembl id ", this_ensembl,
            " was mapped to several HGNC symbols:")
            for (i in seq_len(num_gene_symbols)){
                warning(this_hgnc_character[i])
            }
            warning("Which one to choose? Will return ensembl id instead.")
            return(this_ensembl)
        } else {  ## mapped to exactly one gene symbol
            if (this_hgnc_character == ""){
                message("The HGNC symbol for this Ensembl id ", 
                this_ensembl," is empty!")
                return(this_ensembl)
            } else{
                return(this_hgnc_character)
            }
        }
    }
}



#' Map a given hgnc gene symbol to the ensembl gene id.
#'
#' For a hgnc gene symbol and a tibble with all genes as 
#' input, this function returns the matching ensembl gene 
#' id. The tibble with all genes can be generated with 
#' \code{\link{create_ensembl_gene_tbl_hg}}.
#'
#' @title Get the ensembl gene id for a hgnc gene symbol.
#' @param this_hgnc The hgnc gene symbol of a gene.
#' @param all_genes_tbl A tibble with all genes ensembl id's 
#' and hgnc gene symbols.
#' @author Ariane L. Moore
#' @return The matching ensembl gene id. In case several
#' ensembl gene id's were found, they
#' are all returned with ";" as a separator.
#' @export
#' @importFrom magrittr "%>%"
#' @importFrom dplyr is.tbl filter select distinct pull
#' @examples
#' \dontrun{
#' all_genes_tbl <- create_ensembl_gene_tbl_hg()
#' hgnc_to_ensembl("VHL", all_genes_tbl)
#' hgnc_to_ensembl("PBRM1", all_genes_tbl)
#' }
hgnc_to_ensembl <- function(this_hgnc, all_genes_tbl){
    ensembl_gene_id <- hgnc_symbol <- NULL
    stopifnot(is.character(this_hgnc))
    stopifnot(length(this_hgnc) == 1)
    stopifnot(is.tbl(all_genes_tbl))
    stopifnot("ensembl_gene_id" %in% colnames(all_genes_tbl))
    stopifnot("hgnc_symbol" %in% colnames(all_genes_tbl))

    ## translate gene symbol into ensembl
    this_ensembl_character <- all_genes_tbl %>% 
        filter(hgnc_symbol == this_hgnc) %>%
        select(ensembl_gene_id) %>%
        distinct() %>%
        filter(ensembl_gene_id != "") %>%
        pull(ensembl_gene_id)
    
    num_ensembl_ids=length(this_ensembl_character)
    
    if (num_ensembl_ids == 0){  ## i.e., if there was no 
        ## such ensembl id or the entrez id
        ## was empty
        message("hgnc symbol ", 
        this_hgnc,
        " cannot be found in the tibble or there was",
        " no matching ensembl id!")
        return(this_hgnc)
    } else {
        if (num_ensembl_ids > 1){
            warning("hgnc symbol ", 
            this_hgnc,
            " was mapped to several ensembl id's.")
            toReturn <- 
                paste(this_ensembl_character, collapse=";")
            return(toReturn)
        } else {  ## mapped to exactly one gene id
            return(this_ensembl_character)
        }
    }
}



#' For a tibble that contains the information which ensembl 
#' gene id is mutated in which clone, map the ensembl
#' gene id to the reactome pathways that contain this gene.
#'
#' Such a tibble can be generated with e.g. the function 
#' \code{\link{create_tbl_tree_collection}}. If the altered entities
#' in the lists were the ensembl gene id's, this function 
#' can convert the tibble into a tibble with the altered
#' reactome pathways. It has the columns 'file_name', 
#' 'patient_id', 'altered_entity', 'clone1', 'clone2', ... up 
#' to the maximal number of clones (Default: until 'clone7'). 
#' If the mutated entities are ensembl gene id's, they can 
#' be mapped with this function to the pathways from 
#' 'reactome'. The pathways are from the lowest level of hierarchy.
#'
#' @title Map ensembl gene id clone tibble to reactome pathway 
#' clone tibble.
#' @param mutated_gene_tbl The tibble containing the information
#' of which ensembl gene id is altered in which patient
#'  and clone. Can be created with e.g. 
#'  \code{\link{create_tbl_ent_clones}}.
#' @param ensg_reactome_path_map A tibble with all ensembl id's 
#' and their reactome pathways. Can be loaded with 
#' \code{data("ensg_reactome_path_map")}.
#' @author Ariane L. Moore
#' @return The tibble containing the information of which 
#' pathway is altered in which clone.
#' @export
#' @importFrom magrittr "%>%"
#' @importFrom dplyr tibble is.tbl filter select as.tbl distinct
#' @examples
#' data("ensg_reactome_path_map")
#' mutated_gene_tbl <- 
#'   dplyr::tibble(file_name=c("pat1.csv", "pat1.csv"),
#' patient_id=c("1","1"),
#' altered_entity=c("ENSG00000134086", 
#' "ENSG00000141510"),
#' clone1=c(1,0),
#' clone2=c(0,1))
#' convert_ensembl_to_reactome_pw_tbl(mutated_gene_tbl, 
#'     ensg_reactome_path_map)
convert_ensembl_to_reactome_pw_tbl <- function(mutated_gene_tbl, 
    ensg_reactome_path_map){
    file_name <- patient_id <- altered_entity <- NULL
    stopifnot(is.tbl(mutated_gene_tbl))
    stopifnot(is.tbl(ensg_reactome_path_map))
    stopifnot("file_name" %in% colnames(mutated_gene_tbl))
    stopifnot("patient_id" %in% colnames(mutated_gene_tbl))
    stopifnot("altered_entity" %in% colnames(mutated_gene_tbl))
    stopifnot("clone1" %in% colnames(mutated_gene_tbl))
    stopifnot("clone2" %in% colnames(mutated_gene_tbl))
    stopifnot("ensembl_gene_id" %in% colnames(ensg_reactome_path_map))
    stopifnot("reactome_pw_name" %in% colnames(ensg_reactome_path_map))
    if("tree_id" %in% colnames(mutated_gene_tbl)){
        stop("Please do the gene-to-pathway mapping separately"," 
        for each tree inference.")
    }
    
    ## sanity check: to see whether the entities really are 
    ## ensembl gene id's and whether there are no ensembl gene 
    ## id's double in there
    all_unique_ents <- unique(mutated_gene_tbl$altered_entity)
    which_ones_ENS <- sum(vapply(all_unique_ents, 
        function(x){grep("ENS", x)}, 
        integer(length=1)))
    num_all_unique_ents <- length(all_unique_ents)
    all_ents <- mutated_gene_tbl$altered_entity
    stopifnot(which_ones_ENS == num_all_unique_ents)
    stopifnot(length(all_ents) == num_all_unique_ents)
    
    ## make sure that this is just from one patient
    this_pat <- unique(as.character(mutated_gene_tbl$patient_id))
    stopifnot(length(this_pat) == 1)
    
    ## message to user
    message("Found ", num_all_unique_ents, 
    " different Ensembl gene id's in the provided alteration tibble.")
    
    ## check how many clones are in the tibble
    ## minus the columns for file_name, patient_id, and altered_entity
    num_clones <- dim(mutated_gene_tbl)[2] - 3
    ## check how many colnames are named "clone*"
    stopifnot(sum(grepl("clone",colnames(mutated_gene_tbl))) == num_clones)
  
    ## map all ensembl gene id's to pathways
    ## if a gene id is unmapped, the ensembl gene id remains
    ensembl_to_pw_list <- lapply(all_unique_ents, function(x){
        suppressMessages(ensembl_to_reactome(x, 
        ensg_reactome_path_map))
    })
    names(ensembl_to_pw_list) <- all_unique_ents
    ## if there is no pathway to be found, the ensembl id is retained
    num_unmapped <- sum(grepl("ENS", unlist(ensembl_to_pw_list)))
    message(num_unmapped, 
    " could not be mapped to a reactome pathway.")
    
    ## convert each row of the gene tibble to a tibble itself where the
    ## gene id is replaced by the pathways
    ## this assumes that each altered entity only occurs once in the
    ## tibble
    converted_tbl_list <- lapply(all_unique_ents, function(this_ent){
        these_pws <- ensembl_to_pw_list[[this_ent]]
        ## extract the tibble row with just this gene and
        ## create two tibbles, one with the info of 
        ## file_name and patient_id, and one with
        ## the clones
        mutated_gene_tbl_this_ent_pat_file_info <- 
            mutated_gene_tbl %>%
            filter(altered_entity == this_ent) %>%
            select(file_name, patient_id)
        mutated_gene_tbl_this_ent_only_clones <- 
            mutated_gene_tbl %>%
            filter(altered_entity == this_ent) %>%
            select(-file_name, -patient_id, -altered_entity)
        
        ## this gene is just once in the tibble
        stopifnot(dim(mutated_gene_tbl_this_ent_pat_file_info)[1] == 1)
        stopifnot(dim(mutated_gene_tbl_this_ent_only_clones)[1] == 1)
        
        ## convert the pathways into a data frame
        these_pws_df <- data.frame(altered_entity=these_pws)
        
        ## now we create a new tibble with the column 
        ## 'altered_entity' again, but with the pathways
        ## of this gene
        this_converted_tbl <- as.data.frame(cbind(
            mutated_gene_tbl_this_ent_pat_file_info,
            altered_entity=these_pws_df,
            mutated_gene_tbl_this_ent_only_clones))
        
        ## in order to use bind_rows, we have to make sure that 
        ## the altered ent is character and not factor
        this_converted_tbl$altered_entity <-
            as.character(this_converted_tbl$altered_entity)
        
        stopifnot(dim(this_converted_tbl)[1] == length(these_pws))
        return(this_converted_tbl)
    })
    converted_tbl <- as.tbl(bind_rows(converted_tbl_list))
    
    ## now it could be that several genes were mapped to the same pathway
    ## i.e. we have in the same patient the same mutated entity 
    ## (here pathway) possibly with different clone assignments
    ## this will be merged here:
    merged_converted_tbl <- 
        suppressMessages(merge_clones_identical_ents(converted_tbl))
    
    ## sanity check how many ensembl id's could not be mapped
    num_unmapped_ids <- dim(merged_converted_tbl %>% 
        filter(grepl("ENS", altered_entity)) %>%
        select(altered_entity) %>% 
        distinct())[1]
    stopifnot(num_unmapped_ids == num_unmapped)
    
    ## message to user
    num_pws_unique <- dim(merged_converted_tbl %>% 
        filter(!grepl("ENS", altered_entity)) %>%
        select(altered_entity) %>% 
        distinct())[1]
    num_pws <- dim(merged_converted_tbl %>% 
        filter(!grepl("ENS", altered_entity)) %>%
        select(altered_entity))[1]
    message("The remaining ", 
        num_all_unique_ents-num_unmapped_ids, 
        " different Ensembl gene id's were",
        " mapped to ", num_pws_unique, 
        " different reactome pathways.",
        "\nThere are ", num_pws, 
        " clone-pathway associations in the converted tibble.")
        
    return(merged_converted_tbl)
}




#' Map a given ensembl gene id to the reactome pathways that 
#' contain this gene.
#'
#' As input, an ensembl gene id is given as well as the tibble 
#' 'ensg_reactome_path_map'. It can be loaded with 
#' \code{data("ensg_reactome_path_map")}, and contains the 
#' ensembl gene id to reactome pathway mappings. The reactome
#' pathways are from the lowest level of the hierarchy. 
#' This function returns the reactome pathways for the input gene.
#'
#' @title Get the reactome pathways for an ensembl gene id.
#' @param this_ensembl The ensembl id of a gene.
#' @param ensg_reactome_path_map A tibble with all ensembl 
#' id's and their reactome pathways. Can be loaded with 
#' \code{data("ensg_reactome_path_map")}.
#' @author Ariane L. Moore
#' @return The pathways that contain this gene as a character 
#' vector.
#' @export
#' @importFrom magrittr "%>%"
#' @importFrom dplyr is.tbl filter select
#' @examples
#' data("ensg_reactome_path_map")
#' ensg_gene <- "ENSG00000134086"
#' ensembl_to_reactome(ensg_gene, ensg_reactome_path_map)
ensembl_to_reactome <- function(this_ensembl, 
    ensg_reactome_path_map) {
    reactome_pw_name <- ensembl_gene_id <- NULL
    stopifnot(is.character(this_ensembl))
    stopifnot(length(this_ensembl) == 1)
    stopifnot(is.tbl(ensg_reactome_path_map))
    stopifnot("ensembl_gene_id" %in% colnames(ensg_reactome_path_map))
    stopifnot("reactome_pw_name" %in% colnames(ensg_reactome_path_map))
    
    
    ## retrieve the corresponding reactome pathways
    these_pws_tbl <- ensg_reactome_path_map %>% 
        filter(ensembl_gene_id == this_ensembl) %>%
        select(reactome_pw_name)
    ## convert to character vector
    these_pws<-
        unique(as.character(as.vector(these_pws_tbl$reactome_pw_name)))
    
    ## check how many pathways were found
    num_pws=length(these_pws)
    
    if (num_pws == 0){  ## i.e., if there was no such 
        ## ensembl id in the tibble
        message("The ensembl ID ", this_ensembl,
        " could not be mapped to a reactome pathway.")
        return(this_ensembl)
    } else {
        message("The ensembl ID ", 
        this_ensembl,
        " was mapped to ", num_pws, 
        " reactome pathways.")
        return(these_pws)
    }
}
