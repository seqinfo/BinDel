# Copyright: Priit Paluoja and Priit Palta; University of Tartu and University of Helsinki.
# All rights reserved, unauthorised usage and distribution are prohibited.
# Contact: priit.palta@gmail.com


#' Create a BinDel reference
#'
#' This function creates a BinDel reference from a set of input BAM files and a coordinates file.
#'
#' Note: The function assumes that the input BAM files and the coordinates file have already been prepared and are available for use in the function.
#'
#' @param bam_locations A vector of file paths to BAM files that will be used to generate the reference.
#' @param coordinates_file A file path to the coordinates file. The file should contain the following columns: \emph{chr}, \emph{start}, \emph{end}, \emph{focus}, \emph{length}.
#' @param output_name The name of the output reference file. If the file extension \code{.gz} is appended to the name, the resulting file will be compressed.
#' @param col_names A logical value that specifies whether to include column names in the output file.
#' @param anonymise A logical value that specifies whether to include the original sample name in the reference file. If \code{TRUE}, the sample names are replaced with a \code{prefix + integer}.
#' @param prefix A character string that is used as a prefix when anonymising the sample names.
#'
#' @export
#' @examples
#' write_reference(c("sample1.bam", "sample2.bam", "sample3.bam"), "locations.info.tsv", "reference.gz")
#' write_reference(c("sample.bam"), "locations.info.tsv", "reference.tsv", col_names = FALSE)
write_reference <-
  function(bam_locations,
           coordinates_file,
           output_name,
           col_names = TRUE,
           anonymise = TRUE,
           prefix = "S") {
    message("Creating BinDel reference.")
    
    if (length(bam_locations) == 0) {
      stop("No .bam files found in '", bam_locations, "'.")
    }
    
    if (!file.exists(coordinates_file)) {
      stop("Coordinates file not found in '", coordinates_file, "'.")
    }
    
    if (file.exists(output_name)) {
      stop("Reference already exists, aborting.")
    }
    
    bed <- divide_bins(coordinates_file)
    
    if (col_names) {
      df_cols <- c("chr", "start", "end", "focus", "reads", "sample")
      empty_df <-
        data.frame(matrix(ncol = length(df_cols), nrow = 0))
      colnames(empty_df) <- df_cols
      readr::write_tsv(x = empty_df, output_name)
    }
    
    
    i <- 1
    lapply(bam_locations, function(x) {
      message("Processing '", x, "'.")
      
      binned <- bin_bam(x, bed)
      
      if (anonymise) {
        binned <- binned |>
          dplyr::mutate(sample = paste0(prefix, i))
        i <<- i + 1
      }
      
      readr::write_tsv(
        file = output_name,
        x = binned,
        append = TRUE,
        col_names = F
      )
    })
  }
