% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util.R
\name{coord_to_arm}
\alias{coord_to_arm}
\title{Genomic coordinate to chromosome arm}
\usage{
coord_to_arm(
  chromosome,
  position,
  assembly = "hg19",
  full = FALSE,
  mergesmallarms = FALSE
)
}
\arguments{
\item{chromosome}{Character or numeric vector, with chromosome of genomic coordinate}

\item{position}{Numeric vector, with genomic position within chromosome}

\item{assembly}{a string specifying which genome assembly version should be applied
to determine chromosome arms. Allowed options are "hg38", hg19", "hg18", "hg17"
and "hg16" (corresponding to the five latest human genome annotations in the
UCSC genome browser).}

\item{full}{}

\item{mergesmallarms}{for very small acrocentric arms just use chromosome}
}
\value{
Character vector, with choromosome arm of given genomic coordinates
}
\description{
Returns chromosome arms for given chromosome and genomic position.
}
