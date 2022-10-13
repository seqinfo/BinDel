# Copyright: Priit Paluoja and Priit Palta; University of Tartu and University of Helsinki.
# All rights reserved, unauthorised usage and distribution are prohibited.
# Contact: priit.palta@gmail.com


#' Sanity check for the reference data frame (after initial read in).
#'
#' Checks if:
#' 1. required column names are present.
#' 2. reference sample size is not 0.
#' 3. number of PCA components suit with reference sample size.
#'
#' @importFrom magrittr %>%
#' @param reference A data frame holding reference set
#' @param use_pca Is PCA-normalization switched on?
#' @param nComp How many PCA components are used?
check_reference <- function(reference, use_pca, nComp) {
  ref_expected_cols <-
    c("chr", "start", "end", "focus", "reads", "sample")
  
  for (col in ref_expected_cols) {
    if (!col %in% colnames(reference)) {
      stop(paste0("Reference file has a missing column: '", col, "'."))
    }
  }
  
  ref_size <- nrow(reference %>%
                     dplyr::select(sample) %>%
                     dplyr::distinct(sample))
  
  if (nrow(reference) == 0) {
    stop("Reference file cannot be empty")
  }
  
  if (ref_size < 10) {
    warning("Reference group has only ", ref_size, " samples.")
  } else {
    message("Reference group has ", ref_size, " samples.")
  }
  
  
  if (use_pca & (ref_size < nComp)) {
    stop("nComp (", nComp, ") > reference sample size (", ref_size, ").")
  }
}


#' Check that analyzable sample name is not present in reference set.
#'
#' This check is important as otherwise using group_by in the metrics
#' calculations will lead to incorrect calculations.
#'
#' @importFrom magrittr %>%
#' @param binned_reads binned BAM
#' @param reference A data frame holding reference set
check_uniq <- function(binned_reads, reference) {
  similar_names_in_ref <- nrow(
    binned_reads %>%
      dplyr::select(sample) %>%
      dplyr::distinct(sample) %>%
      dplyr::inner_join(reference %>% dplyr::select(sample))
  )
  
  if (similar_names_in_ref != 0) {
    stop("BAM name is present in the reference group. Stopping.")
  }
}
