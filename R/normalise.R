# Copyright: Priit Paluoja and Priit Palta; University of Tartu and University of Helsinki.
# All rights reserved, unauthorised usage and distribution are prohibited.
# Contact: priit.palta@gmail.com


#' Correct GC% bias in bins
#'
#' Sample wise, based on PMID: 28500333 and PMID: 20454671. Expects columns
#' \emph{chr}, \emph{start}, \emph{end}, \emph{focus}, \emph{reads},
#' \emph{samples}. Creates column \emph{gc_correct}.
#' 
#' @importFrom magrittr %>%
#' @param samples A data frame to GC% correct.
#' @return A GC% corrected data frame.
gc_correct <- function(samples) {
  message("Applying GC% correct.")
  return(
    samples %>% 
      dplyr::left_join(find_gc(
        dplyr::select(.data = ., chr, start, end, focus) %>% 
          dplyr::distinct(chr, start, end, focus)
      )) %>% 
      # Do GC correct
      dplyr::group_by(sample, gc) %>% 
      dplyr::mutate(reads_gc_interval = mean(reads))%>% 
      dplyr::ungroup() %>% 
      dplyr::group_by(sample) %>% 
      dplyr::mutate(gc_corrected = reads * mean(reads) / reads_gc_interval) %>% 
      dplyr::filter(!is.na(gc_corrected)) %>% 
      dplyr::ungroup()
  )
}


#' Normalize bin read count by sample total read count and bin length.
#'
#' @param samples \code{\link{gc_correct}} output.
#' @return A normalized data frame.
normalize_reads <- function(samples) {
  message("Normalizing by read count and bin length.")
  return(
    samples |>
      # Sample read count correct
      dplyr::group_by(sample) |>
      dplyr::mutate(gc_corrected = gc_corrected / sum(gc_corrected)) |>
      dplyr::ungroup() |>
      # Sample bin length correct
      dplyr::mutate(gc_corrected = gc_corrected / (end - start)) |> 
      dplyr::ungroup()
  )
}


#' PCA-normalize
#'
#' This function performs PCA-based normalization on a normalized data frame that has undergone previous normalization steps such as read count normalization and GC-correction.
#'
#' @param samples A data frame resulting from the \code{\link{normalize_reads}} function.
#' @param num_comp The number of PCA components to use in the normalization. If not specified, the function will use the cumulative variance to determine the number of components.
#' @param cumulative_variance The percentage of cumulative variance to be explained by the number of PCA components. If not specified, the function will use the \code{num_comp} parameter to determine the number of components.
#' @return A normalized data frame.
#' 
#' @references
#' This function is based on a \href{https://stats.stackexchange.com/questions/229092/how-to-reverse-pca-and-reconstruct-original-variables-from-several-principal-com}{StackExchange answer}.
#' 
#' @seealso
#' \code{\link{normalize_reads}}, \code{\link{gc_correct}}
#' 
pca_correct <- function(samples, num_comp, cumulative_variance) {

  if (!is.null(num_comp) & is.null(cumulative_variance)) {
    message("Applying PCA normalization with ", num_comp, " components.")
  } else if (is.null(num_comp) & !is.null(cumulative_variance)) {
    message("Applying PCA normalization with cumulative variance of ", cumulative_variance,".")
  } else if (!is.null(num_comp) & !is.null(cumulative_variance)) {
    stop("num_comp and cumulative_variance cannot both be specified. Please choose one or the other.")
  } else {
    stop("Either num_comp or cumulative_variance must be specified.")
  }
  
  # For PCA sort ()
  samples <- samples |>
    dplyr::arrange(reference)
  
  # Pivot wide for PCA normalization
  message(head(samples %>%  dplyr::select(focus) %>% dplyr::distinct()))
  
  wider <- samples |>
    dplyr::select(focus, start, sample, reference, gc_corrected) |>
    tidyr::pivot_wider(
      names_from = c(focus, start),
      id_cols = c(sample, reference),
      values_from = gc_corrected,
      names_sep = ":"
    )
  
  # Train PCA
  ref <- wider |>
    dplyr::filter(reference) |>
    dplyr::select(-reference,-sample)
  
  mu <- colMeans(ref, na.rm = T)
  
  refPca <- stats::prcomp(ref)

  if(is.null(num_comp)){
    num_comp <- (factoextra::get_eig(refPca) |> 
                dplyr::mutate(PCA = dplyr::row_number()) |> 
                dplyr::filter(cumulative.variance.percent >= cumulative_variance) |> 
                dplyr::select(PCA))[1,1]
    
    message("Setting PCA automatically to ", num_comp, " based on the target ", cumulative_variance, "% of cumulative variance calculated from the reference set.")
  }
  
  Xhat <- refPca$x[, 1:num_comp] %*% t(refPca$rotation[, 1:num_comp])
  Xhat <- scale(Xhat, center = -mu, scale = FALSE)
  
  # Use trained PCA on other samples
  pred <- wider |>
    dplyr::filter(!reference) |>
    dplyr::select(-reference,-sample)
  
  rm(wider)
  
  Yhat <-
    stats::predict(refPca, pred)[, 1:num_comp] %*% t(refPca$rotation[, 1:num_comp])
  Yhat <- scale(Yhat, center = -mu, scale = FALSE)
  
  # Actual PCA normalization and conversion back to long:
  normalized <-
    dplyr::bind_rows(as.data.frame(as.matrix(pred) / as.matrix(Yhat)),
                     as.data.frame(as.matrix(ref) / as.matrix(Xhat))) |>
    tidyr::pivot_longer(
      names_sep = ":",
      names_to = c("focus", "start"),
      cols = dplyr::everything(),
      values_to = "gc_corrected"
    )
  
  normalized$sample <- samples$sample
  normalized$reference <- samples$reference
  normalized$chr <- samples$chr
  
  # Replace NAs with 0.
  normalized$gc_corrected[is.na(normalized$gc_corrected)] <- 0
  return(normalized)
}
