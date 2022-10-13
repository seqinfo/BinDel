# Copyright: Priit Paluoja and Priit Palta; University of Tartu and University of Helsinki.
# All rights reserved, unauthorised usage and distribution are prohibited.
# Contact: priit.palta@gmail.com

#' Infer sample high-risk probability.
#'
#' This function does the following:
#' \enumerate{
#' \item Apply bin-based GC% correct.
#' \item Normalizing by read count.
#' \item Normalize by bin length.
#' \item Apply PCA-based normalization
#' \item Calculate reference group statistics per focus region.
#' \item Calculate sample Mahalanobis distance from the reference group.
#' \item Calculate -log10 chi-squared distribution probabilities.
#' }
#'
#' Outputs results to different files.
#'
#' @importFrom magrittr %>%
#' @param bam_location Path to the input file.
#' @param reference_location Path to the reference file.
#' @param use_pca Use PCA based normalization? Increases sensitivity by reducing variation between samples.
#' @param nComp How many components to use in PCA-based normalization? Cannot be higher than the number of samples in the reference set. It affects the minimum read count and the fetal fraction that the tool can reliably operate on.
#' @param bin_plot Create and save detailed bin plots?
#' @param save_bins Save bins?
#' @export
#' @examples
#'
#' bam <- "sample.bam"
#' reference <- "reference.tsv"
#' infer_normality(bam, reference)
#' head("sample.bam.tsv")
infer_normality <- function(bam_location,
                            reference_location,
                            use_pca = TRUE,
                            nComp = 80,
                            bin_plot = TRUE,
                            save_bins = FALSE)  {
  message_package_version()
  
  
  message("Reading reference set from: ", reference_location)
  reference <-
    readr::read_tsv(reference_location) %>%
    dplyr::mutate(reference = TRUE) %>%
    dplyr::group_by(sample) %>%
    # de-identify reference.
    dplyr::mutate(sample = as.character(dplyr::cur_group_id())) %>%
    dplyr::ungroup()
  
  check_reference(reference, use_pca, nComp)
  
  message("Reading and binning: ", bam_location)
  binned_reads <- bin_bam(
    bam_location,
    reference %>%
      dplyr::select(chr, start, end, focus) %>%
      dplyr::distinct(chr, start, end, focus)
  ) %>%
    dplyr::mutate(reference = FALSE)
  
  check_uniq(binned_reads, reference)
  
  
  message("Merging BAM and reference for calculations.")
  samples <- reference %>%
    dplyr::bind_rows(binned_reads)
  
  rm(binned_reads, reference)
  
  
  samples <- samples %>%
    gc_correct() %>%
    normalize_reads() %>%
    # Optimize memory
    dplyr::select(chr, focus, start, sample, reference, gc_corrected)
  
  if (use_pca) {
    samples <- pca_correct(samples, nComp)
  }
  
  
  samples <- calculate_bin_stat(samples)
  
  sample_name <- basename(bam_location)
  
  
  if (bin_plot) {
    message("Creating and saving region plots.")
    save_bin_plot(samples, sample_name)
  }
  
  
  if (save_bins) {
    message("Writing bins to file.")
    readr::write_tsv(
      samples %>%
        dplyr::ungroup() %>%
        dplyr::filter(sample == sample_name) %>%
        dplyr::select(chr, focus, start, PPDX_norm) %>%
        dplyr::rename(`Normalized Z-Score` = PPDX_norm),
      paste0(sample_name, ".bins", ".tsv")
    )
  }
  
  
  samples <- calculate_summary(samples)
  
  
  message("Saving metrics.")
  samples %>%
    dplyr::filter(!reference) %>%
    dplyr::select(-reference) %>%
    dplyr::rename(subregion = focus) %>%
    dplyr::rename(`z_score` = PPDX) %>%
    dplyr::rename(`normalized_z_score` = PPDX_norm) %>%
    dplyr::mutate_if(is.numeric, round, 3) %>%
    readr::write_tsv(paste0(sample_name, ".tsv"))
}
