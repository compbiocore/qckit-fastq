% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kmer_count.R
\name{kmer_count}
\alias{kmer_count}
\title{Return kmer count per sequence for the length of kmer desired}
\usage{
kmer_count(infile, k, output_file = NA)
}
\arguments{
\item{infile}{the object that is the path to gzippped FASTQ file}

\item{k}{the length of kmer}

\item{output_file}{File to save plot to. Default NA.}
}
\value{
kmers counts per sequence
}
\description{
Return kmer count per sequence for the length of kmer desired
}
\examples{

infile <- system.file("extdata", "10^5_reads_test.fq.gz",
    package = "qckitfastq")
km<-kmer_count(infile,k=4)
km[1:20,1:10]

}
