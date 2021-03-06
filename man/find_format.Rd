% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{find_format}
\alias{find_format}
\title{Gets quality score encoding format from the FASTQ file. Return possibilities are Sanger(/Illumina1.8),
Solexa(/Illumina1.0), Illumina1.3, and Illumina1.5. This encoding is heuristic based and may not be 100% accurate
since there is overlap in the encodings used, so it is best if you already know the format.}
\usage{
find_format(infile, reads_used)
}
\arguments{
\item{infile}{A string giving the path for the fastq file}

\item{reads_used}{int, the number of reads to use to determine the encoding format.}
}
\value{
A string denoting the read format. Possibilities are Sanger, Solexa, Illumina1.3, and Illumina1.5.
}
\description{
Gets quality score encoding format from the FASTQ file. Return possibilities are Sanger(/Illumina1.8),
Solexa(/Illumina1.0), Illumina1.3, and Illumina1.5. This encoding is heuristic based and may not be 100% accurate
since there is overlap in the encodings used, so it is best if you already know the format.
}
\examples{
infile <- system.file("extdata", "10^5_reads_test.fq.gz", package = "qckitfastq")
find_format(infile,100)
}
