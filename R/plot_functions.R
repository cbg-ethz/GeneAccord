#' This function plots the average rates of clonal exclusivity for each 
#' patient.
#'
#' In addition to the average rate of clonal exclusivity, it also visualizes
#' the average number of clones of each patient. 
#'
#' @title Barplot of rates of clonal exclusivity and number of clones.
#' @param avg_rates_m A named vector with the average rates of clonal 
#' exclusivity for each patient. The name of each
#' element is the patient id to be used in the barplot.
#' @param clone_tbl The tibble containing the gene-to-clone assignments 
#' from all patients and all trees from the collection
#' of trees.
#' @param output_pdf The name of the pdf to be generated. Or if output_pdf 
#' is "direct", then the plot is
#' generated directly and not to a pdf. Default: "direct"
#' @author Ariane L. Moore
#' @return None, the function plots the average rates of clonal exclusivity.
#' @import 
#' dplyr
#' ggplot2
#' tibble
#' @examples
#' clone_tbl <- tibble::as_tibble(cbind(
#'             "file_name"=rep(c(rep(c("fn1", "fn2"), each=3)), 2),
#'             "patient_id"=rep(c(rep(c("pat1", "pat2"), each=3)), 2),
#'             "altered_entity"=c(rep(c("geneA", "geneB", "geneC"), 4)),
#'             "clone1"=c(0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0),
#'             "clone2"=c(1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1),
#'             "tree_id"=c(rep(1, 6), rep(2, 6))))
#' avg_rates_m <- c(pat1=0.014, pat2=0.3)
#' plot_rates_clon_excl(avg_rates_m, clone_tbl)
plot_rates_clon_excl <- function(avg_rates_m, clone_tbl, 
    output_pdf="direct") {
    tree_id <- pat_ids_c <- rates_m <- num_c <- pat_ids_r <- NULL 
    stopifnot(is.numeric(avg_rates_m))
    stopifnot(dplyr::is.tbl(clone_tbl))
    stopifnot("tree_id" %in% colnames(clone_tbl))
    stopifnot(is.character(output_pdf))
    
    ## find out the average number of clones
    all_tree_ids <- unique(as.character(clone_tbl$tree_id))
    num_trees <- length(all_tree_ids)
    num_pats <- length(unique(as.character(clone_tbl$patient_id)))
    stopifnot(num_pats == length(avg_rates_m))
    ## message to user
    message(paste("There are rates from ", num_pats, " patients.", 
        sep=""))
    
    all_trees_num_clones <- c(rep(0, num_pats))
    for (this_tree in all_tree_ids){
        clone_tbl_this_tree <- clone_tbl %>% 
        dplyr::filter(tree_id == this_tree) %>% 
        dplyr::select(-tree_id)
        this_tree_num_clones <- 
            suppressMessages(extract_num_clones_tbl(clone_tbl_this_tree))
        ## this counts the number of clones that have at least one
        ## gene assigned to it
        all_trees_num_clones <- all_trees_num_clones + this_tree_num_clones
    }
    all_trees_avg_num_clones <- all_trees_num_clones/num_trees
   
    ## put it into the right ordering
    avgNumClones <- 
        all_trees_avg_num_clones[order(match(names(all_trees_avg_num_clones), 
        names(avg_rates_m)))]
    stopifnot(is.numeric(avgNumClones))
  
    ## extract patient id's and make sure that avg_rates_m and avgNumClones 
    ## are in the same order
    pat_ids_rates <- names(avg_rates_m)
    pat_ids_clones <- names(avgNumClones)
    stopifnot(is.character(pat_ids_rates))
    stopifnot(length(setdiff(pat_ids_rates, pat_ids_clones)) == 0)
    avgNumClones <- avgNumClones[match(pat_ids_rates, pat_ids_clones)]
    pat_ids_clones <- pat_ids_clones[match(pat_ids_rates, pat_ids_clones)]
  
    myOrder <- order(nchar(pat_ids_rates), pat_ids_rates)
    rates_and_clones_tbl <- tibble::tibble(pat_ids_r=pat_ids_rates[myOrder],
        pat_ids_c=pat_ids_clones[myOrder],
        rates_m=avg_rates_m[myOrder],
        num_c=avgNumClones[myOrder]) %>% 
        dplyr::filter(pat_ids_r == pat_ids_c)
    stopifnot(dim(rates_and_clones_tbl)[1] == num_pats)
    ## message to user
    message(paste("The average rate of clonal exclusivity is between ", 
    round(min(avg_rates_m), digits=2), 
    "-", round(max(avg_rates_m), digits=2), sep=""))
    
    ## this is to set the order among the ggplot bars
    rates_and_clones_tbl$pat_ids_r <- factor(rates_and_clones_tbl$pat_ids_r, 
        levels=rev(rates_and_clones_tbl$pat_ids_r))
    rates_and_clones_tbl$pat_ids_c <- factor(rates_and_clones_tbl$pat_ids_c, 
        levels=rev(rates_and_clones_tbl$pat_ids_c))
  
    if(output_pdf != "direct")
        pdf(output_pdf, height=10, width=5)
    
    this_plot <- ggplot2::ggplot(rates_and_clones_tbl, 
            ggplot2::aes(x=pat_ids_r, 
            y=rates_m, fill=num_c)) +
        ggplot2::geom_bar(stat="identity") +
        ggplot2::ggtitle("Mean rates of clonal exclusivity in each patient") +
        ggplot2::ylab("Rate m") +
        ggplot2::xlab("Patients") +
        ggplot2::scale_fill_gradient(low="lightblue", high="darkblue", 
            guide="colourbar") +
        ggplot2::guides(fill=ggplot2::guide_colourbar(title=paste0("Average",
            "number of clones"))) +
        ggplot2::coord_flip() + 
        ggplot2::theme_gray()
    print(this_plot)
    if(output_pdf != "direct"){
        dev.off()
        stopifnot(file.exists(output_pdf))
    }
}


#' This function plots the ECDFs of the test statistic under the null 
#' hypothesis.
#'
#' The ECDF's of the test statistic under the null for a data set can be 
#' generated with \code{\link{generate_ecdf_test_stat}}. 
#' Afterwards, they can be visualized with this function. It is assumed 
#' that the first ECDF in the ecdf_list is
#' the ECDF for the case where pairs are mutated in two patients.
#'
#' @title Plot empirical cumulative distribution functions of the test 
#' statistic under the null.
#' @param ecdf_list The list of ECDF's as generated with 
#' \code{\link{generate_ecdf_test_stat}}.
#' @param plot_idx The index of which of the list entries of the ecdf_list 
#' to plot. Default: c(2,3).
#' @param num_panel_rows The ECDF's will be plotted altogether, hence 
#' \code{par(mfrow=c(x, y))} is used.
#' Here, \code{x} is the number of panel rows, which has to be set with 
#' this parameter, and \code{y} will be taken as
#' ceil(#ECDF's/x). E.g., if you have 20 ECDF's in total, you can set 
#' \code{num_panel_rows=4}, and then your 20
#' ECDF's will be plotted in panels with four rows, and five columns. 
#' Default=1.
#' @param output_pdf The name of the pdf to be generated. Or if output_pdf 
#' is "direct", then the plot is
#' generated directly and not to a pdf. Default: "direct".
#' @author Ariane L. Moore
#' @return None, the function plots ecdf curves.
#' @import 
#' graphics
#' @examples
#' avg_rates_m <- c(pat1=0.1, pat2=0.034, pat3=0.21, pat4=0.063)
#' list_of_num_trees_all_pats <- list(pat1=c(20, 20, 19), 
#'                               pat2=c(20, 18, 20),
#'                               pat3=c(19, 20, 20), 
#'                               pat4=c(20, 20, 20))
#' list_of_clon_excl_all_pats <- list(pat1=c(5, 0, 1), 
#'                               pat2=c(10, 2, 0),
#'                               pat3=c(18, 12, 0), 
#'                               pat4=c(0, 2, 0))
#' num_pat_pair_max <- 2
#' num_pairs_sim <- 10
#' ecdf_list <- generate_ecdf_test_stat(avg_rates_m, 
#'                                      list_of_num_trees_all_pats, 
#'                                      list_of_clon_excl_all_pats, 
#'                                      num_pat_pair_max, 
#'                                      num_pairs_sim)
#' plot_ecdf_test_stat(ecdf_list, plot_idx=2)
plot_ecdf_test_stat <- function(ecdf_list, plot_idx=c(2,3), 
    num_panel_rows=1, output_pdf="direct"){
    stopifnot(is.list(ecdf_list))
    stopifnot(is.numeric(plot_idx))
    stopifnot(is.numeric(num_panel_rows))
    stopifnot(is.character(output_pdf))
  
    num_ecdfs <- length(ecdf_list)
    stopifnot(max(plot_idx) <= num_ecdfs)
    num_ecdfs_to_plot <- length(plot_idx)
    ## the first list entry is just set to NULL!
    real_num_ecdfs <- num_ecdfs-1
    message(paste("Found ", num_ecdfs_to_plot, " ECDF's to plot.", 
        sep=""))
  
    stopifnot(num_panel_rows <= num_ecdfs_to_plot)
    stopifnot(num_panel_rows > 0)
  
    ## plot the ecdf's
    num_panel_cols <- ceiling(num_ecdfs_to_plot/num_panel_rows)
    ## make sure it is an integer
    if( num_panel_cols != round(num_panel_cols, digits=0) || 
    num_panel_rows != round(num_panel_rows, digits=0) ||
    (num_panel_rows*num_panel_cols < num_ecdfs_to_plot)){
        stop(paste0("The number of rows and columns for the panels needs to",
        "be integers and need to multiply to a number which is greater or",
        " equal to the number of ECDF's to plot!"))
    }
  
    if(output_pdf != "direct"){
        this_height <- num_panel_rows*4
        this_width <- (num_panel_cols/num_panel_rows)*this_height
        pdf(output_pdf, height=this_height, width=this_width)
    }
  
    ## set the number of panel rows and columns
    par(mfrow=c(num_panel_rows, num_panel_cols))
  
    ## plot ecdf's
    for (i in plot_idx){
        if(i == 1){ ## the first ecdf is just NULL
            if(is.null(ecdf_list[[i]])){
                stop(paste0("Cannot plot the first ECDF of the list because by",
                " default, it is set to NULL."))
            }
        }
        this_ecdf <- ecdf_list[[i]]
        graphics::plot(this_ecdf, 
            main=paste("ECDF (Number of patients a pair",
            " is mutated in=", i,")", 
            sep=""))
        grid()
    }
  
    if(output_pdf != "direct"){
        dev.off()
        stopifnot(file.exists(output_pdf))
    } 
}




#' This function visualizes the distribution of p-values.
#'
#' It is especially useful, when exploring the results with simulated 
#' data under the null hypothesis, i.e. when delta is zero.
#' In that scenario, the p-values are expected to be uniformly distributed. 
#' This function can take the p-values from
#' \code{\link{generate_test_stat_hist}} where the concatenated tibble
#'  contains different values for 'num_pat_pair', 
#' i.e. the number of patients the simulated pairs are mutated in. The 
#' input tibble is expected to have the two columns 
#' 'pval', and 'num_patients'. Left panel: histogram of all p-values from 
#' the whole tibble.
#' Right panel: ecdf of the p-values with different colors for different
#'  numbers of patients that the pairs were mutated in.
#'
#' @title Plot histogram and empirical cumulative distribution function of 
#' p-values.
#' @param res_sim tibble containing the simulated pairs of genes/pathways. 
#' It contains the columns 'num_patients',
#' and 'pval', and can be generated with 
#' \code{\link{generate_test_stat_hist}}
#' and then concatenating the tibbles.
#' @param output_pdf The name of the pdf to be generated. Or if output_pdf 
#' is "direct", then the plot is
#' generated directly and not to a pdf. Default: "direct"
#' @author Ariane L. Moore
#' @return None, the function plots a p-value histogram.
#' @import 
#' dplyr
#' RColorBrewer
#' @examples
#' res_sim <- tibble::tibble(num_patients=c(rep(2,100), 
#'                           rep(3,100), rep(4,100)),
#'                           pval=c(runif(300)))
#' vis_pval_distr_num_pat(res_sim)
vis_pval_distr_num_pat <- function(res_sim, output_pdf="direct"){
    pval <- num_patients <- NULL
    stopifnot(is.character(output_pdf))
    stopifnot(dplyr::is.tbl(res_sim))
    stopifnot("pval" %in% colnames(res_sim))
    stopifnot("num_patients" %in% colnames(res_sim))
  
    ## extract the p-values
    p_values <- as.numeric(res_sim$pval)
    num_p_values <- length(p_values)
  
    ## plot the histogram and qqplot
    if(output_pdf != "direct")
        pdf(output_pdf, height=6, width=10)
    par(mfrow=c(1,2), oma=c(3, 3, 3, 3)) 
    ##  make the outer margin at the bottom of the plot large
  
    ## histogram
    hist(p_values, main="Histogram of p-values", xlab="P-values", 
        breaks=25, xlim=c(0,1))
  
    ####
    ## qq-plot
    ## check how many observations there are for the different num_patients
    num_patients_tally <- res_sim %>% 
        dplyr::select(num_patients) %>% 
        dplyr::group_by(num_patients) %>% 
        dplyr::tally()
    number_different_patients <- dim(num_patients_tally)[1]
    num_patients_colors <- c()
    if(number_different_patients > 8){
        num_patients_colors <- rainbow(number_different_patients)
    } else {
        num_patients_colors <- c(RColorBrewer::brewer.pal(8, "Dark2"))
    }
  
  
    ## first num_patients
    this_num_patients <- as.numeric(num_patients_tally[1,1])
    this_num_pat_pvals <- 
        as.numeric(as.vector(as.data.frame(res_sim %>% 
        dplyr::filter(num_patients == this_num_patients) %>%
        dplyr::select(pval))[,1]))
    this_num_p_values <- length(this_num_pat_pvals)
  
    p_vals_idx <- seq(1, this_num_p_values)
    max_points <- 5000
    if (this_num_p_values > max_points){
        p_vals_idx <- sample(p_vals_idx, max_points)
    }
    plot(cbind(sort(runif(length(p_vals_idx))), 
        sort(this_num_pat_pvals[p_vals_idx])), main="Uniform Q-Q Plot", 
        xlab="Theoretical Quantiles",
        ylab="Sample Quantiles", xlim=c(0,1), ylim=c(0,1), 
        col=num_patients_colors[1], pch=".")
  
    these_num_patients <- c(this_num_patients)
  
    ## Potential other num_patients
    if(number_different_patients > 1){
        for(i in seq(2, number_different_patients)){
            this_num_patients <- as.numeric(num_patients_tally[i,1])
            this_num_pat_pvals <- 
            as.numeric(as.vector(as.data.frame(res_sim %>% 
                dplyr::filter(num_patients == this_num_patients) %>%
                dplyr::select(pval))[,1]))
            this_num_p_values <- length(this_num_pat_pvals)
            p_vals_idx <- seq(1, this_num_p_values)
            max_points <- 5000
            if (this_num_p_values > max_points){
                p_vals_idx <- sample(p_vals_idx, max_points)
            }
            points(sort(runif(length(p_vals_idx))), 
            sort(this_num_pat_pvals[p_vals_idx]), 
            col=num_patients_colors[i], pch=".")
  
            these_num_patients <- c(these_num_patients, this_num_patients)
        }
    }
    abline(0, 1, col="lightgrey")
    grid()
  
  
    ## overlay the entire figure region with a new, single plot. Then call 
    ## legend with the location ("bottom", ...)
    par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), mar=c(0, 0, 0, 0), 
        new=TRUE)
    plot(0, 0, type="n", bty="n", xaxt="n", yaxt="n")
    legend("bottom", inset=c(0, 0), xpd=TRUE, horiz=TRUE,
        legend=these_num_patients, 
        col=num_patients_colors, pch=19, pt.cex=0.8, 
        title="Number of patients")
    ## xpd=TRUE tells R that it is OK to plot outside the region
    ## horiz=TRUE tells R that I want a horizontal legend
    ## inset=c(x,y) tells R how to move the legend relative to the 
    ## 'bottom' location
  
  
    if(output_pdf != "direct"){
        dev.off()
        stopifnot(file.exists(output_pdf))
    }
}





#' This function plots the heatmaps of final gene clone matrices.
#'
#' After running the \code{\link{GeneAccord}}, one may want to visualize 
#' the gene clone
#' heatmap for significant gene pairs.
#'
#' @title Heatmaps of gene pairs of interest
#' @param pairs_of_interest The tibble containing the pairs of 
#' genes/pathways that should be visualized in the heatmap.
#' This may be, e.g. the gene pairs were mle_delta > 0, qval < 0.1, 
#' and num_patients > 1. It contains the columns 'entity_A',
#' and 'entity_B', and can be generated with \code{\link{GeneAccord}}. 
#' For the plot, the function will attempt
#' to map the gene ID's from ensembl ID to gene name. However, if the 
#' input genes are not ensembl IDs, it does not matter.
#' @param clone_tbl The tibble containing the information of which 
#' gene/pathway is mutated in which
#' clone from which patient. Here, it is assumed that only one tree 
#' from the collection of trees was chosen per
#' patient.
#' @param all_genes_tbl A tibble with all genes ensembl id's and 
#' hgnc symbols. Can be created with \code{\link{create_ensembl_gene_tbl_hg}}.
#' @param first_clone_is_N Logical indicating whether the first 
#' clone column is actually representing the normal or germline, 
#' and is not a tumor
#' clone. In that case, it will have the name 'N', and all other 
#' columns will be one clone number smaller, e.g. 'clone2' is then actually
#' 'clone1' etc. Default: FALSE.
#' @param output_pdf The name of the pdf to be generated. Or if 
#' output_pdf is "direct", then the plot is
#' generated directly and not to a pdf. Default: "direct"
#' @author Ariane L. Moore
#' @return None, the function plots a gene-to-clone assignment heatmap.
#' @import 
#' dplyr
#' ggplot2
#' ggpubr
#' reshape2
#' @examples
#' pairs_of_interest <- tibble::tibble(entity_A="SETD2",
#'                            entity_B="BAP1")
#' clone_tbl <- tibble::tibble(
#'                file_name=c("X.csv", "X.csv", "Y.csv", "Y.csv"),
#'                patient_id=c("X", "X", "Y", "Y"),
#'                altered_entity=c("SETD2", "BAP1", "SETD2", "BAP1"),
#'                clone1=c(0, 1, 1, 0),
#'                clone2=c(1, 0, 0, 1))
#' \dontrun{all_genes_tbl <- create_ensembl_gene_tbl_hg()}
#' all_genes_tbl_example <- tibble::tibble(
#'                  ensembl_gene_id=c("ENSG00000181555", 
#'                  "ENSG00000163930"),
#'                  hgnc_symbol=c("SETD2", "BAP1"))
#' heatmap_clones_gene_pat(pairs_of_interest, clone_tbl, 
#' all_genes_tbl_example)
heatmap_clones_gene_pat <- function(pairs_of_interest, clone_tbl, 
    all_genes_tbl, 
    first_clone_is_N=FALSE, 
    output_pdf="direct"){
    file_name <- patient_id <- entity_A <- entity_B <- variable <- 
        altered_entity <- value <- n <- NULL
    stopifnot(dplyr::is.tbl(pairs_of_interest))
    stopifnot(dplyr::is.tbl(clone_tbl))
    stopifnot(dplyr::is.tbl(all_genes_tbl))
    stopifnot(is.character(output_pdf))
    stopifnot("entity_A" %in% colnames(pairs_of_interest))
    stopifnot("entity_B" %in% colnames(pairs_of_interest))
    stopifnot("file_name" %in% colnames(clone_tbl))
    stopifnot("patient_id" %in% colnames(clone_tbl))
    stopifnot("altered_entity" %in% colnames(clone_tbl))
    stopifnot(is.logical(first_clone_is_N))
    if("tree_id" %in% colnames(clone_tbl)){
        stop(paste("Can only plot for one tree at a time. Make sure",
        " that the clone tbl is just from one tree id, and does not contain the column 'tree_id'", sep=""))
    }
  
    ## extract the genes
    all_A <- as.character(pairs_of_interest$entity_A)
    all_B <- as.character(pairs_of_interest$entity_B)
    genes_of_interest <- unique(c(all_A, all_B))
  
    ## clone columns are renamed
    ## if the column labeled 'clone1' is actually the normal/germline, 
    ## the cancer clone numbers are also reduced by one
    clone_col_names <- colnames(clone_tbl)[grepl("clone", 
        colnames(clone_tbl))]
    new_clone_col_names <- c()
    for (this_col in clone_col_names){
        if(first_clone_is_N){
            if (this_col == "clone1"){
                this_new_col <- "N"
            } else {
                this_new_col <- paste0("Clone", as.numeric(sub("clone",
                "", this_col))-1)
            }
        } else {
            this_new_col <- sub("clone", "Clone", this_col)
        }
        new_clone_col_names <- c(new_clone_col_names, this_new_col)
    }
    colnames(clone_tbl)[grepl("clone", colnames(clone_tbl))] <- 
        new_clone_col_names
  
    ## here we remove columns that are just clone zero in all patients,
    ## for plotting
    filterd_clone_tbl_just_clones <- clone_tbl %>% 
       dplyr::select(-file_name, -patient_id, -altered_entity)
    filterd_clone_tbl_nonzero_clones <- cbind(clone_tbl['patient_id'],
       clone_tbl['altered_entity'],
    filterd_clone_tbl_just_clones[, colSums(filterd_clone_tbl_just_clones) > 0])
    ## extract from the clone tbl just the entries with the genes of interest
    filterd_clone_tbl_goi <- filterd_clone_tbl_nonzero_clones %>% 
        dplyr::filter(altered_entity %in% genes_of_interest)
  
  
    ## map the ensembl gene ID's to the gene names
    these_ens_ids <- as.vector(as.data.frame(filterd_clone_tbl_goi %>% 
        dplyr::select(altered_entity))[,1])
    hgnc_symbols <- c()
    for(this_ens in these_ens_ids){
        this_hgnc <- suppressMessages(ensembl_to_hgnc(this_ens, all_genes_tbl))
        hgnc_symbols <- c(hgnc_symbols, this_hgnc)
    }
    filterd_clone_tbl_goi_gene_names <- 
        cbind(as.data.frame(filterd_clone_tbl_goi), hgnc_symbols)
  
    ## plot the heatmap of gene pairs of interest
    
    ## do the heatmap for each patient separate, and then arrange them together to one plot
    all_pats <- 
        unique(as.character(as.vector(as.data.frame(filterd_clone_tbl_goi_gene_names %>%
        dplyr::select(patient_id) %>%
        dplyr::group_by(patient_id) %>%
        dplyr::tally() %>%
        dplyr::filter(n >= 2) %>%
        dplyr::select(patient_id))[,1])))
    num_pats_to_plot <- length(all_pats)
    
    
    ## define the pdf output
    if(output_pdf != "direct")
        pdf(output_pdf, width=num_pats_to_plot*5)
  
    plot_list <- list()
    cnt <- 0
  
    for (this_pat in all_pats){
        tbl_to_plot <- filterd_clone_tbl_goi_gene_names %>%
        dplyr::filter(patient_id == this_pat) %>%
        dplyr::select(-altered_entity, -patient_id)
  
        ## here we remove columns from the end that are just clone zero in this current patient, for plotting
        tbl_to_plot_just_clones <- tbl_to_plot %>% 
            dplyr::select(-hgnc_symbols)
        col_sums_clones <- colSums(tbl_to_plot_just_clones)
        idx_clones_this_pat <- seq(1, max( which(col_sums_clones != 0 )))
        tbl_to_plot_nonzero_clones <- cbind(tbl_to_plot['hgnc_symbols'],
            tbl_to_plot_just_clones[, idx_clones_this_pat])
      
  
        cnt <- cnt + 1
  
        tbl_to_plot_melted <- 
            suppressMessages(reshape2::melt(tbl_to_plot_nonzero_clones))
  
        colors <- rev(c("darkblue", "lightblue"))
        this_plot <- ggplot2::ggplot() +
            ggplot2::geom_tile(data=tbl_to_plot_melted, 
                ggplot2::aes(x=hgnc_symbols, 
                y=variable, 
                fill=factor(value))) +
            ggplot2::geom_tile(data=tbl_to_plot_melted, 
               ggplot2::aes(x=hgnc_symbols, 
               y=variable), 
               size=1, fill=NA, color="black") +
            ggplot2::scale_fill_manual(values=colors, name="Mutation status") +
            ggplot2::theme_minimal() +
            ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90, 
                vjust=0.5, hjust=1)) + 
        ## Vertical text on x axis
        ggplot2::labs(x="") + ggplot2::labs(y="") +
        ggplot2::ggtitle(paste("Patient ", this_pat, sep="")) +
        ggplot2::theme(plot.title=ggplot2::element_text(size=14,
            face="bold"),
            text=ggplot2::element_text(size=14),
            axis.title=ggplot2::element_text(size=14),
            axis.text.x=ggplot2::element_text(size=14),
            axis.text.y=ggplot2::element_text(size=14))
        plot_list[[cnt]] <- this_plot
    }
          
    if (cnt > 156) {
        my_letters <- NULL
    } else if(cnt > 26){
        my_letters <- c(LETTERS,
            paste(LETTERS, ".1", sep=""),
            paste(LETTERS, ".2", sep=""),
            paste(LETTERS, ".3", sep=""),
            paste(LETTERS, ".4", sep=""),
            paste(LETTERS, ".5", sep=""))
    } else {
        my_letters <- LETTERS[seq(1, cnt)]
    }
    plot_altogether <- ggpubr::ggarrange(plotlist=plot_list, 
        ncol=cnt, labels=my_letters)
    print(plot_altogether)
    if(output_pdf != "direct"){
        dev.off()
        stopifnot(file.exists(output_pdf))
    }
}





