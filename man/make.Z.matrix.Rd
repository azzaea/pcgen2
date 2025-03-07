% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make.Z.matrix.R
\name{make.Z.matrix}
\alias{make.Z.matrix}
\title{Create the incidence matrix}
\usage{
make.Z.matrix(genotype)
}
\arguments{
\item{genotype}{a vector (character or factor) with genotype labels}
}
\description{
The incidence (design) matrix links observations (replicates) to unique genotypes. In
expriements without replicates, this reduces to the Identity matrix
}
