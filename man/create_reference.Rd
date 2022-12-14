% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reference.R
\name{create_reference}
\alias{create_reference}
\title{Create a reference file.}
\usage{
create_reference(bam_locations, bed_location, reference_name)
}
\arguments{
\item{bam_locations}{The path to the folder, where the reference files exists.}

\item{bed_location}{A location to the .bed file with columns: \emph{chr}, \emph{start}, \emph{end}, \emph{focus}.}

\item{reference_name}{The name of the output file. File must not exist.}
}
\description{
Takes a location to the folder with .bam files and .bed file and bins the
.bam files with GRCh38.
}
\details{
Outputs .tsv reference file with columns:
\enumerate{
\item \emph{chr}
\item \emph{start}
\item \emph{end}
\item \emph{focus}
\item \emph{reads}
\item \emph{sample}
}
}
