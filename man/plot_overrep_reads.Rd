% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_overrep_reads.R
\name{plot_overrep_reads}
\alias{plot_overrep_reads}
\title{Plot the top 5 seqeunces}
\usage{
plot_overrep_reads(overrep_reads, output_file = NA)
}
\arguments{
\item{overrep_reads}{the table that sorts the sequence content and corresponding counts in descending order}

\item{output_file}{File to save plot to. Will not write to file if NA. Default NA.}
}
\value{
plot of the top 5 overrepresented sequences
}
\description{
Plot the top 5 seqeunces
}
\examples{

infile <- system.file("extdata", "10^5_reads_test.fq.gz", package = "qckitfastq")
overrep_df <- overrep_reads(infile)
plot_overrep_reads(overrep_df)

}
