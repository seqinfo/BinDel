# Copyright: Priit Paluoja and Priit Palta; University of Tartu and University of Helsinki.
# All rights reserved, unauthorised usage and distribution are prohibited.
# Contact: priit.palta@gmail.com


#' Create a reference file.
#'
#' Takes a location to the folder with .bam files and .bed file and bins the
#' .bam files with GRCh38. 
#' 
#' Outputs .tsv reference file with columns:
#' \enumerate{
#' \item \emph{chr}
#' \item \emph{start}
#' \item \emph{end}
#' \item \emph{focus}
#' \item \emph{reads}
#' \item \emph{sample}
#'}
#' @param bam_locations The path to the folder, where the reference files exists.
#' @param bed_location A location to the .bed file with columns: \emph{chr}, \emph{start}, \emph{end}, \emph{focus}.
#' @param reference_name The name of the output file. File must not exist.
#'
#' @export
create_reference <-
  function(bam_locations,
           bed_location,
           reference_name) {
    files <- list.files(
      path = bam_locations,
      pattern = "*.bam",
      full.names = TRUE,
      recursive = FALSE
    )
    
    if (length(files) == 0) {
      stop("No .bam files found in '", bam_locations, "'.")
    }
    
    if (!file.exists(bed_location)) {
      stop(".bed file not found in '", bed_location, "'.")
    }
    
    if (file.exists(reference_name)) {
      stop("Output '", reference_name, "' already exists.")
    }
    
    bed <- readr::read_tsv(bed_location)
    
    
    df_cols <- c("chr", "start", "end", "focus", "reads", "sample")
    empty_df <- data.frame(matrix(ncol = length(df_cols), nrow = 0))
    colnames(empty_df) <- df_cols
    readr::write_tsv(x = empty_df, reference_name)
    
    
    lapply(files, function(x) {
      message("Processing '", x, "'.")
      readr::write_tsv(
        file = reference_name,
        x = bin_bam(x, bed),
        append = TRUE,
        col_names = FALSE
      )
    })
  }
