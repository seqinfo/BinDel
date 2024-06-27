# Copyright: Priit Paluoja and Priit Palta; University of Tartu and University of Helsinki.
# All rights reserved, unauthorised usage and distribution are prohibited.
# Contact: priit.palta@gmail.com

#' Infer sample high-risk probability based on a reference set.
#'
#' @param bam_file_path character. The path to the BAM file to be processed.
#' @param ref_file_path character. The path to the reference file in TSV format.
#' @param cumulative_variance numeric, optional. The cumulative variance threshold (%) for general PCA correction over all the genome. Default is 95 (%).
#' @param cumulative_variance_per_region numeric vector, optional. A vector of cumulative variance thresholds (%) for PCA correction to apply per region. Default is NULL.
#' @param output_file_path character, optional. The path to the output file where results will be saved. If not provided, results will be returned as a data frame. Default is NULL.
#' @param bin_file_path character, optional. The path to the file where bin data will be saved. Default is NULL.
#' @return If `output_file_path` is not provided, returns a data frame with the processed data. Otherwise, results are saved to the specified output file.
#' @export
#' @examples
#' infer_sample_high_risk_probability(
#'   bam_file_path = "sample.bam",
#'   ref_file_path = "reference.tsv",
#'   cumulative_variance = 95,
#'   cumulative_variance_per_region = c(0, 10),
#'   output_file_path = "output.tsv"
#' )
infer_normality <- function(bam_file_path,
                            ref_file_path,
                            cumulative_variance = 95.0,
                            cumulative_variance_per_region = NULL,
                            output_file_path = NULL,
                            bin_file_path = NULL)  {
  start.time <- Sys.time()
  message_package_version()
  
  message("Reading reference set from: ", ref_file_path)
  reference <-
    readr::read_tsv(ref_file_path) |>
    dplyr::mutate(reference = TRUE)
  
  message("Reading and binning: ", bam_file_path)
  binned_reads <- reference |>
    dplyr::distinct(chr, start, end, focus) |>
    (function(.) bin_bam(bam_file_path, .))() |>
    dplyr::mutate(reference = FALSE)
  
  
  message("Merging BAM and reference for calculations.")
  samples <- reference |>
    dplyr::bind_rows(binned_reads)
  
  rm(binned_reads, reference)
  
  samples <- samples |>
    gc_correct() |>
    normalize_reads()
  

  if (!is.null(cumulative_variance) ) {
    message("Applying PCA normalization.")
    samples <- pca_correct(samples, cumulative_variance)
  }
  samples$pcar <- 0
  
  if(is.null(cumulative_variance_per_region)){
    cumulative_variance_per_region <- c(0)
  }
  
  
  if(!is.null(cumulative_variance_per_region) && length(cumulative_variance_per_region) > 0) {
    message("Applying PCA normalization per subregion.")
    
    # Create an empty list to store all normalized samples
    all_normalized_samples <- list()
    
    # Iterate over each value in the cumulative_variance_per_region vector
    for (pcar in cumulative_variance_per_region) {
      message(paste("Processing with cumulative variance:", pcar))
      
      # Split the data by the "focus" variable
      sample_groups <- split(samples, samples$focus)
      
      if (pcar != 0) {
        # Normalize each sample group using the current variance value
        normalized_groups <- lapply(sample_groups, pca_correct, cumulative_variance = pcar)
        
        # Combine the results back into a single data frame
        normalized_samples <- do.call(rbind, normalized_groups)
      } else {
        # If variance is 0, skip PCA correction and keep original samples
        normalized_samples <- do.call(rbind, sample_groups)
      }
      
      # Append the normalized samples to the list
      all_normalized_samples[[length(all_normalized_samples) + 1]] <- normalized_samples
    }
    
    # Combine all normalized samples into a single data frame
    samples <- do.call(rbind, all_normalized_samples)
  }
  

  message("Analysing and comparing bins")
  samples <- calculate_bin_stat(samples)
  
  if (!is.null(bin_file_path)) {
    message("Writing bins to file.")
    readr::write_tsv(samples |>
                       dplyr::filter(!reference) |>
                       dplyr::select(-reference),
                     bin_file_path)
  }
  
  message("Calculating risk scores")
  samples <- calculate_summary(samples) |> 
    dplyr::filter(!reference) |>
    dplyr::select(-reference)
  

  end.time <- Sys.time()
  time.taken <- end.time - start.time
  message(paste0("Time elapsed:", time.taken))
  
  if (!is.null(output_file_path)) {
    readr::write_tsv(samples, output_file_path)
  } else{
    samples
  }
}
