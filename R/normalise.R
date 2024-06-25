# Copyright: Priit Paluoja and Priit Palta; University of Tartu and University of Helsinki.
# All rights reserved, unauthorised usage and distribution are prohibited.
# Contact: priit.palta@gmail.com


#' Apply GC% correction to samples.
#'
#' This function applies GC% correction to the read counts in the given samples data frame.
#' It assumes the presence of columns: \emph{chr}, \emph{start}, \emph{end}, \emph{focus}, \emph{reads}, and \emph{sample}.
#'
#' @param samples A data frame containing the required columns.
#' @return A data frame with an added column \emph{gc_corrected}.
gc_correct <- function(samples) {
  message("Applying GC% correct.")
  
  # Ensure required columns are present
  required_cols <- c("chr", "start", "end", "focus", "reads", "sample", "reference")
  if (!all(required_cols %in% colnames(samples))) {
    stop("The samples data frame must contain the following columns: ", paste(required_cols, collapse = ", "))
  }
  
  return(
    samples %>% 
      dplyr::left_join(find_gc(
        dplyr::distinct(.data = ., chr, start, end, focus))) %>% 
      # Do GC correct
      dplyr::group_by(sample, gc, reference) %>% 
      dplyr::mutate(reads_gc_interval = mean(reads))%>% 
      dplyr::ungroup() %>% 
      dplyr::group_by(sample, reference) %>% 
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
      dplyr::group_by(sample, reference) |>
      dplyr::mutate(gc_corrected = gc_corrected / sum(gc_corrected)) |>
      dplyr::ungroup() |>
      # Sample bin length correct
      dplyr::mutate(gc_corrected = gc_corrected / (end - start)) |> 
      dplyr::ungroup()
  )
}


#' Perform PCA-based normalization on a normalized data frame.
#'
#' This function applies PCA-based normalization to a data frame that has undergone
#' previous normalization steps such as read count normalization and GC-correction.
#'
#' @param samples A data frame resulting from the \code{\link{normalize_reads}} function.
#' @param cumulative_variance The cumulative variance threshold (%) to determine the number of PCA components.
#' @return A normalized data frame.
#'
#' @references
#' This function is based on a \href{https://stats.stackexchange.com/questions/229092/how-to-reverse-pca-and-reconstruct-original-variables-from-several-principal-com}{StackExchange answer}.
#'
pca_correct <- function(samples, cumulative_variance) {
  # Arrange samples by reference for consistency
  samples <- samples |>
    dplyr::arrange(reference)
  
  # Pivot wide for PCA normalization
  wider <- samples |>
    tidyr::pivot_wider(
      names_from = c(focus, start),
      id_cols = c(sample, reference),
      values_from = gc_corrected,
      names_sep = ":"
    )
  
  # Extract reference data for PCA training
  ref <- wider |>
    dplyr::filter(reference) |>
    dplyr::select(-reference,-sample)
  
  mu <- colMeans(ref, na.rm = T)
  
  refPca <- stats::prcomp(ref)
  
  # Determine the number of principal components required to reach the specified cumulative variance
  # Calculate eigenvalues from the standard deviations of the principal components
  eig_values <- refPca$sdev^2
  # Calculate the cumulative variance percentage explained by each principal component
  cumulative_variance_percent <- cumsum(eig_values / sum(eig_values)) * 100
  
  # Create a data frame with PCA component numbers and their cumulative variance percentages
  num_comp <- data.frame(
    PCA = seq_along(cumulative_variance_percent),  # PCA component numbers
    cumulative.variance.percent = cumulative_variance_percent) |>  # Cumulative variance percentages
    # Filter to retain rows where the cumulative variance percentage is greater than or equal to the specified threshold
    dplyr::filter(cumulative.variance.percent >= cumulative_variance) |>
    # Extract the PCA component numbers that meet the criteria
    dplyr::pull(PCA) |>
    # Select the first PCA component that meets the criteria
    dplyr::first()
    
  message("Using ", num_comp, " PCA components based on the target ", cumulative_variance, "% of cumulative variance calculated from the reference set.")
  
  
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
  normalized$pcar <- cumulative_variance
  
  # Replace NAs with 0.
  normalized$gc_corrected[is.na(normalized$gc_corrected)] <- 0
  return(normalized)
}
