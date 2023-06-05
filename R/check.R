# Copyright: Priit Paluoja and Priit Palta; University of Tartu and University of Helsinki.
# All rights reserved, unauthorised usage and distribution are prohibited.
# Contact: priit.palta@gmail.com


#' Check reference data frame for consistency
#' 
#' This function checks the reference data frame for consistency, ensuring that it has the
#' expected columns, the sample size is not zero, and the number of PCA components is appropriate.
#'
#' @param reference A data frame holding reference set
#' @param use_pca  Should PCA normalization be used?
#' @param num_comp The number of PCA components to use?
check_reference <- function(reference, use_pca, num_comp) {
  ref_expected_cols <- c("chr", "start", "end", "focus", "reads", "sample")
  
  for (col in ref_expected_cols) {
    if (!col %in% colnames(reference)) {
      stop(paste0("Reference file has a missing column: '", col, "'."))
    }
  }
  
  ref_size <- nrow(reference |>
                     dplyr::select(sample) |>
                     dplyr::distinct(sample))
  
  if (nrow(reference) == 0) {
    stop("Reference file cannot be empty")
  }
  
  if (ref_size < 10) {
    warning("Reference group has only ", ref_size, " samples.")
  } else {
    message("Reference group has ", ref_size, " samples.")
  }
  
  
  if (use_pca & (!is.null(num_comp))) {
    if(ref_size < num_comp){
      stop("Number of PCA components (", num_comp, ") cannot be greater than reference sample size (", ref_size, ").")
    }
  }
}


#' Check that analyzable sample name is not present in reference set.
#'
#' This check is important as otherwise using group_by in the metrics
#' calculations will lead to incorrect calculations.
#'
#' @param binned_reads binned BAM
#' @param reference A data frame holding reference set
check_uniq <- function(binned_reads, reference) {
  similar_names_in_ref <- nrow(
    binned_reads |>
      dplyr::select(sample) |>
      dplyr::distinct(sample) |>
      dplyr::inner_join(reference |> dplyr::select(sample))
  )
  
  if (similar_names_in_ref != 0) {
    stop("BAM name is present in the reference group. Stopping.")
  }
}
