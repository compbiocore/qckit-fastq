% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/overrep_reads.R
\name{overrep_reads}
\alias{overrep_reads}
\title{Sort all sequences per read by count.}
\usage{
overrep_reads(infile, output_file = NA)
}
\arguments{
\item{infile}{Path to gzippped FASTQ file.}

\item{output_file}{File to save data frame to. Default NA.}
}
\value{
Table of sequences sorted by count.
}
\description{
Sort all sequences per read by count.
}
\examples{
infile <- system.file("extdata", "10^5_reads_test.fq.gz",
    package = "qckitfastq")
overrep_reads(infile)[1:5,]
}
