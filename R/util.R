# Copyright: Priit Paluoja and Priit Palta; University of Tartu and University of Helsinki.
# All rights reserved, unauthorised usage and distribution are prohibited.
# Contact: priit.palta@gmail.com


#' Bin aligned sequences (from \emph{.bam}) into genomic bins based on the \emph{.bed} file.
#'
#' @importFrom magrittr %>%
#' @param bam_location A location to the BAM-file to bin.
#' @param bed A data frame in .bed format with columns: \emph{chr}, \emph{start}, \emph{end}, \emph{focus}.
#' @return A data frame in bed format with GC%.
#' @export
bin_bam <- function(bam_location, bed) {
  bam <-
    GenomicAlignments::readGAlignments(bam_location,
                                       param = Rsamtools::ScanBamParam(
                                         flag = Rsamtools::scanBamFlag(isDuplicate = FALSE,
                                                                       isSecondaryAlignment = FALSE),
                                         what = c("pos")
                                       ))
  
  binned_counts <- bed %>% 
    dplyr::mutate(reads = SummarizedExperiment::assay(
      GenomicAlignments::summarizeOverlaps(GenomicRanges::makeGRangesFromDataFrame(.), bam, inter.feature=FALSE, mode = "IntersectionStrict"))[,1]) %>% 
    dplyr::mutate(sample = basename(bam_location))
  
  return(binned_counts)
}


#' Find GC% for GRCh38 based .bed data frame.
#'
#' @param bed A data frame in .bed format with columns: \emph{chr}, \emph{start}, \emph{end}, \emph{focus}.
#' @return A data frame in bed format with GC%.
find_gc <- function(bed) {
  reads <- bed |>
    dplyr::mutate(gc = Biostrings::letterFrequency(
      Biostrings::getSeq(
        BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
        GenomicRanges::GRanges(chr, IRanges::IRanges(start = start, end = end))
      ),
      "GC",
      as.prob = T
    )) |>
    dplyr::mutate(gc = round(gc, 1))
}


#' Message current package name and version.
#'
#' @export
message_package_version <- function() {
  message(paste(
    utils::packageName(),
    "version:",
    utils::packageVersion(utils::packageName())
  ))
}


#' Create based on the location file.
#'
#' @param locations_file A tab separated text file describing how to create bins with columns: \emph{chr}, \emph{start}, \emph{end}, \emph{focus}, \emph{length}.
#' @return A data frame containing bins.
#' @export
divide_bins <- function(locations_file) {
  message("Dividing '", locations_file, "' to bins.")
  
  df <- readr::read_tsv(locations_file)
  
  header <- c('chr', 'start', 'end', 'focus') 
  
  df_out <- data.frame(matrix(ncol = length(header), nrow = 0))
  names(df_out) <- header
  
  for (i in 1:nrow(df)) {
    line <- df[i, ]
    
    chromosome <- line[["chr"]]
    start <- as.integer(line[["start"]])
    end <- as.integer(line[["end"]])
    focus <- line[["focus"]]
    bin_width <- as.integer(line[["length"]])
    
    while (start + bin_width < end) {
      output <- c(chromosome, as.integer(start), as.integer(min(start + bin_width - 1, end)), focus)
      output <- data.frame(t(output))
      names(output) <- header
      df_out <- rbind(df_out, output)
      
      start <- start + bin_width
    }
  }
  return(df_out)
}
