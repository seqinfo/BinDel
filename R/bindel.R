# Copyright: Priit Paluoja and Priit Palta; University of Tartu and University of Helsinki.
# All rights reserved, unauthorised usage and distribution are prohibited.
# Contact: priit.palta@gmail.com

#' Infer sample high-risk probability based on a reference set.
#'
#' Note: The function assumes that the reference set has already been processed and is available for use in the calculations.
#'
#' This function does the following:
#' \enumerate{
#' \item Apply bin-based GC% correct.
#' \item Normalizing by read count.
#' \item Normalize by bin length.
#' \item Apply PCA-based normalization
#' \item Calculate reference group statistics per focus region.
#' \item Calculate sample Mahalanobis distance from the reference group.
#' \item Calculate \code{-log10} chi-squared distribution probabilities.
#' }
#' @param bam_file_path A character string specifying the path to the BAM file containing the samples to be analyzed.
#' @param ref_file_path A character string specifying the path to the reference file.
#' @param use_pca A logical value indicating whether to use PCA-based normalization. If \code{TRUE}, PCA-based normalization will be used, otherwise not.
#' @param use_pca_per_region A logical value indicating whether to use PCA-based normalization per region.
#' @param cumulative_variance A numeric value specifying the cumulative percentage of variance to be retained by the number of PCA components in the normalization (based on the reference set). The default value is \code{95.0}.
#' @param cumulative_variance_per_region A numeric value specifying the cumulative percentage of variance to be retained by the number of PCA components in the per region normalization. The default value is \code{50.0}.
#' @param num_comp An integer specifying the number of PCA components to use for normalization. If \code{cumulative_variance} is not specified, this value will be used instead. The default value is \code{NULL}.
#' @param num_comp_per_region An integer specifying the number of PCA components to use for normalization per region. If \code{cumulative_variance_per_region} is not specified, this value will be used instead. The default value is \code{NULL}.
#' @param create_bin_plot A logical value indicating whether to create and save detailed bin plot. If \code{TRUE}, detailed bin plots will be created and saved, otherwise not. The default value is \code{FALSE}.
#' @param save_bins A logical value indicating whether to save the bins. If \code{TRUE}, the bins will be saved to file, otherwise not. The default value is \code{FALSE}.
#' @param output_file_path Set file name (path) for output.
#' @param internally_deidentify_reference  A logical value indicating whether to rename the reference samples to avoid naming conflicts with the analyzable samples during computations. If \code{TRUE}, reference samples will be internally renamed, otherwise not. The default value is \code{TRUE}
#' @param output_reference_scores A logical value indicating whether to output the high-risk probabilities for the reference set during computations. If \code{TRUE}, the high-risk probabilities for the reference set will be outputted, otherwise not. The default value is \code{FALSE}
#' @param write_probabilities_dataframe A logical value indicating whether to also write the high-risk probabilities to a file. If \code{TRUE}, the high-risk probabilities will be written to file. The default value is \code{TRUE}
#' @param output_intermediate_scores A logical value indicating whether to output all intermediate scores used for calculating probabilities. If \code{TRUE}, all intermediate scores will be outputted, otherwise not. The default value is \code{FALSE}
#' @param score_file_extension Specify the desired file extension for the output score file.
#' @param bin_file_extension Specify the desired file extension for the output bin file.
#' @export
#' @examples
#' infer_normality("sample.bam", "reference.gz")
#' 
infer_normality <- function(bam_file_path,
                            ref_file_path,
                            use_pca = TRUE,
                            use_pca_per_region = FALSE,
                            cumulative_variance = 95.0,
                            cumulative_variance_per_region = 50.0,
                            num_comp = NULL,
                            num_comp_per_region = NULL,
                            create_bin_plot = FALSE,
                            save_bins = FALSE,
                            output_file_path = NULL,
                            internally_deidentify_reference = TRUE,
                            output_reference_scores = FALSE,
                            write_probabilities_dataframe = TRUE,
                            output_intermediate_scores = FALSE,
                            score_file_extension = ".scores",
                            bin_file_extension = ".bins")  {
  message_package_version()
  
  
  message("Reading reference set from: ", ref_file_path)
  reference <-
    readr::read_tsv(ref_file_path) |>
    dplyr::mutate(reference = TRUE)
  
  if (internally_deidentify_reference) {
    reference <- reference |>
      dplyr::group_by(sample) |>
      # de-identify reference.
      dplyr::mutate(sample = as.character(dplyr::cur_group_id())) |>
      dplyr::ungroup()
  }
  
  
  check_reference(reference, use_pca, num_comp)
  
  message("Reading and binning: ", bam_file_path)
  binned_reads <- bin_bam(
    bam_file_path,
    reference |>
      dplyr::select(chr, start, end, focus) |>
      dplyr::distinct(chr, start, end, focus)
  ) |>
    dplyr::mutate(reference = FALSE)
  
  check_uniq(binned_reads, reference)
  
  
  message("Merging BAM and reference for calculations.")
  samples <- reference |>
    dplyr::bind_rows(binned_reads)
  
  rm(binned_reads, reference)
  
  
  samples <- samples |>
    gc_correct() |>
    normalize_reads() |>
    # Optimize memory
    dplyr::select(chr, focus, start, sample, reference, gc_corrected)
  
  if (use_pca) {
    samples <- pca_correct(samples, num_comp, cumulative_variance)
  }
  
  if(use_pca_per_region){
    # Normalize per target/focus region
    # Split the data by the "focus" variable
    sample_groups <- split(samples, samples$focus)
    sample_groups <- lapply(sample_groups, 
                            pca_correct, 
                            num_comp = num_comp_per_region,
                            cumulative_variance = cumulative_variance_per_region)
    
    # Combine the results back into a single data frame
    samples <- do.call(rbind, sample_groups)
  }
  
  
  samples <- calculate_bin_stat(samples)
  
  sample_name <- basename(bam_file_path)
  
  if (is.null(output_file_path)) {
    output_file_path = sample_name
  }
  
  if (create_bin_plot) {
    message("Creating and saving region plots.")
    save_bin_plot(samples, sample_name, output_file_path)
  }
  
  
  if (save_bins) {
    message("Writing bins to file.")
    
    if (!output_reference_scores) {
      bins <- samples |>
        dplyr::filter(!reference) |>
        dplyr::select(-reference)
    }else{
      bins <- samples
    }
    
    readr::write_tsv(
      bins |>
        dplyr::select(sample, chr, focus, start, PPDX_norm, PPDX) |>
        dplyr::ungroup() |>
        dplyr::rename(`z_score` = PPDX) |>
        dplyr::rename(`Normalized Z-Score` = PPDX_norm) |>
        dplyr::mutate_if(is.numeric, function(x) format(x, scientific = FALSE)),
      paste0(output_file_path, bin_file_extension)
    )
    
    rm(bins)
  }
  
  
  samples <- calculate_summary(samples)
  
  
  if (!output_reference_scores) {
    samples <- samples |>
      dplyr::filter(!reference)
  }
  
  message("Saving metrics.")
  output <- samples  |>
    dplyr::rename(subregion = focus) |>
    dplyr::rename(`z_score` = PPDX) |>
    dplyr::rename(`normalized_z_score` = PPDX_norm) |>
    dplyr::mutate_if(is.numeric, function(x) format(x, scientific = FALSE))
    
  
  
  if(!output_intermediate_scores){
    output <- output |> 
      dplyr::select(sample, chr, subregion, greedy_prob, conservative_prob)
  }
  
  
  if (write_probabilities_dataframe) {
    readr::write_tsv(output, paste0(output_file_path, score_file_extension))
  } else{
    output
  }
}
