# Copyright: Priit Paluoja and Priit Palta; University of Tartu and University of Helsinki.
# All rights reserved, unauthorised usage and distribution are prohibited.
# Contact: priit.palta@gmail.com


#' Calculate bin statistics (microduplication and microdeletion risk bin stat).
#'
#' Calculate statistics for each genomic bin compared to a reference group.
#' This function should be called after normalization.
#'
#' Required columns in the input data frame are: \emph{chr}, \emph{start},
#' \emph{reference}, \emph{gc_corrected}.
#'
#' @param samples A normalized data frame containing the required columns.
#' @return A data frame with additional columns \emph{dup_score} and \emph{md_score}.
calculate_bin_stat <- function(samples) {
  # Calculate each reference bin i mean
  reference <- samples |>
    dplyr::filter(reference) |>
    dplyr::group_by(chr, start, pcar) |>
    dplyr::summarise(mean_ref_bin = mean(gc_corrected),
                     mean_ref_sd = sd(gc_corrected)) |>
    dplyr::ungroup()
  

  samples <- samples |>
    dplyr::right_join(reference) |>
    dplyr::mutate(
      z_score = (gc_corrected - mean_ref_bin) / mean_ref_sd,
      over_mean = as.integer(gc_corrected >= mean_ref_bin),
      under_mean = as.integer(gc_corrected <= mean_ref_bin)
    ) |>
    dplyr::filter(!is.na(z_score)) |>
    dplyr::group_by(sample, chr, focus, reference, pcar) |>
    dplyr::mutate(
      n = dplyr::n(),
      # Normalize with under_mean/over_mean + Laplace smoothing
      dup_score = ((z_score / sqrt(n)) + 1 / n) / (sum(under_mean) + 2 / n),
      md_score = ((z_score / sqrt(n)) + 1 / n) / (sum(over_mean) + 2 / n)
    ) |>
    dplyr::ungroup()
  
  return(samples)
}


#' Calculate summary statistics including Mahalanobis distances.
#'
#' @param samples A binned statistics data frame.
#' @return A data frame with probabilities
calculate_summary <- function(samples) {
  
  # Summarize duplication and deletion scores
  samples <- samples |>
    dplyr::group_by(sample, chr, focus, reference, pcar) |>
    dplyr::summarise(dup_score = sum(dup_score),
                     md_score = sum(md_score)) |>
    dplyr::group_by(reference, chr, focus, pcar) |>
    dplyr::mutate(mean_dup = mean(dup_score),
                  mean_md = mean(md_score)) |>
    dplyr::ungroup()
  
  # Split data by chromosome and focus (region)
  samples <- samples |>
    dplyr::group_by(chr, focus, pcar) |>
    dplyr::group_split() |>
    purrr::map_dfr( ~ {
      # Calculate covariance and Mahalanobis distance for duplication score
      cov <- stats::cov(.x |>
                          dplyr::filter(reference) |>
                          dplyr::select(dup_score))
      center <- .x |>
        dplyr::filter(reference) |>
        dplyr::select(mean_dup) |>
        dplyr::distinct()
      
      dist_dup <- tryCatch({
        stats::mahalanobis(.x |>
                             dplyr::select(dup_score),
                           c(center$mean_dup),
                           cov,
                           tol = 1e-20)
      }, error = function(e) {
        warning("Error in Mahalanobis distance calculation for dup_score: ", e)
        NA
      })
      
      
      # Calculate covariance and Mahalanobis distance for deletion score
      cov <- stats::cov(.x |>
                          dplyr::filter(reference) |>
                          dplyr::select(md_score))
      
      center <- .x |>
        dplyr::filter(reference) |>
        dplyr::select(mean_md) |>
        dplyr::distinct()
      
      dist_md <- tryCatch({
        stats::mahalanobis(.x |>
                             dplyr::select(md_score),
                           c(center$mean_md),
                           cov,
                           tol = 1e-20)
      }, error = function(e) {
        warning("Error in Mahalanobis distance calculation for md_score: ", e)
        NA
      })
      
      
      dplyr::bind_cols(
        sample = .x$sample,
        chr = .x$chr,
        subregion = .x$focus,
        reference = .x$reference,
        pcar = .x$pcar,
        dist_md = dist_md,
        md_prob = -log10(stats::pchisq(dist_md, df = 1, lower.tail = F) + 1e-100),
        dist_dup = dist_dup,
        dup_prob = -log10(stats::pchisq(dist_dup, df = 1, lower.tail = F) + 1e-100)
      )
    }) 
}
