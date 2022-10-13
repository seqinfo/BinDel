# Copyright: Priit Paluoja and Priit Palta; University of Tartu and University of Helsinki.
# All rights reserved, unauthorised usage and distribution are prohibited.
# Contact: priit.palta@gmail.com


#' GC% bin-correct
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
      dplyr::mutate(reads_gc_interval = mean(reads)) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(sample) %>%
      dplyr::mutate(gc_corrected = reads * mean(reads) / reads_gc_interval) %>%
      dplyr::filter(!is.na(gc_corrected)) %>%
      dplyr::ungroup()
  )
}


#' Normalize bin read count by sample total read count and bin length.
#'
#' @importFrom magrittr %>%
#' @param samples \code{\link{gc_correct}} output.
#' @return A normalized data frame.
normalize_reads <- function(samples) {
  message("Normalizing by read count and bin length.")
  return(
    samples %>%
      # Sample read count correct
      dplyr::group_by(sample) %>%
      dplyr::mutate(gc_corrected = gc_corrected / sum(gc_corrected)) %>%
      dplyr::ungroup() %>%
      # Sample bin length correct
      dplyr::mutate(gc_corrected = gc_corrected / (end - start)) %>% 
      dplyr::ungroup()
  )
}


#' PCA-normalize.
#'
#' Expects other normalization steps (\code{\link{normalize_reads}} and
#' \code{\link{gc_correct}}) to be completed.
#' Based on \href{https://stats.stackexchange.com/questions/229092/how-to-reverse-pca-and-reconstruct-original-variables-from-several-principal-com}{StackExchange answer}.
#'
#' @importFrom magrittr %>%
#' @param samples \code{\link{normalize_reads}} output.
#' @param nComp How many PCA components to use in the normalization.
#' @return A normalized data frame.
pca_correct <- function(samples, nComp) {
  message("Applying PCA normalization with ", nComp, " components.")
  # For PCA sort ()
  samples <- samples %>%
    dplyr::arrange(reference)
  
  # Pivot wide for PCA normalization
  wider <- samples %>%
    dplyr::select(focus, start, sample, reference, gc_corrected) %>%
    tidyr::pivot_wider(
      names_from = c(focus, start),
      id_cols = c(sample, reference),
      values_from = gc_corrected,
      names_sep = ":"
    )
  
  # Train PCA
  ref <- wider %>%
    dplyr::filter(reference) %>%
    dplyr::select(-reference,-sample)
  
  mu <- colMeans(ref, na.rm = T)
  refPca <- stats::prcomp(ref)
  
  
  Xhat <- refPca$x[, 1:nComp] %*% t(refPca$rotation[, 1:nComp])
  Xhat <- scale(Xhat, center = -mu, scale = FALSE)
  
  # Use trained PCA on other samples
  pred <- wider %>%
    dplyr::filter(!reference) %>%
    dplyr::select(-reference,-sample)
  
  rm(wider)
  
  Yhat <-
    stats::predict(refPca, pred)[, 1:nComp] %*% t(refPca$rotation[, 1:nComp])
  Yhat <- scale(Yhat, center = -mu, scale = FALSE)
  
  # Actual PCA normalization and conversion back to long:
  normalized <-
    dplyr::bind_rows(as.data.frame(as.matrix(pred) / as.matrix(Yhat)),
                     as.data.frame(as.matrix(ref) / as.matrix(Xhat))) %>%
    tidyr::pivot_longer(
      names_sep = ":",
      names_to = c("focus", "start"),
      cols = dplyr::everything(),
      values_to = "gc_corrected"
    )
  
  normalized$sample <- samples$sample
  normalized$reference <- samples$reference
  normalized$chr <- samples$chr
  return(normalized)
}
