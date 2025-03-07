% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{gemmaMVLMM}
\alias{gemmaMVLMM}
\title{R interface to 'Genome-wide Efficient Mixed Model Association' (GEMMA)}
\usage{
gemmaMVLMM(
  genoinputs,
  kfile,
  colnums,
  k_mode = 1L,
  miss = 0.05,
  maf = 0.01,
  r2 = 0.9999,
  hwe = 0,
  notsnp = FALSE,
  lmmMode = 1L,
  gxe = "",
  outprefix = "result",
  outdir = "output",
  license = FALSE
)
}
\arguments{
\item{genoinputs}{Character vector of input files. It should be either the
prefix name of PLINK binary ped file, or the 3 BIMBAM files' names (the
menan genotypes, the phenotyes and SNP annotaiton files)}

\item{kfile}{Character variable for the relatedness matrix file name (can be
in gzip compressed format)}

\item{colnums}{Numeric vector for specifying column number from the
phenotypes file used for association testing}

\item{k_mode}{Integer variable for the type of the kinship/relatedness
matrix type}

\item{miss, maf, r2, hwe}{Floating point variables for filtering SNPs at
missingness cutt-off (default 0.05), minor allele frequency (default
0.01), r2 threshold (default 0.9999) and HWE test p value threshold
(default 0; no test)}

\item{notsnp}{Boolean variable to use real values as covariates (minor
allele frequency cutoff is not used)}

\item{lmmMode}{Integer specifying the frequenist test to use: 1 for the Wald
test, 2 for the likelihood ratio test, 3 for score test; and 4 for
performing all the three tests}

\item{gxe}{Character variable for the file containing a column of
environmental variables. If provided, GEMMA fits a LMM controlling for
both the SNP main effect and the environmental main effect while testing
for interaction effect.}

\item{outprefix}{Character variable for output file prefix}

\item{outdir}{Character variable for output directory path}

\item{license}{Boolean variable for printing GEMMA's license information}
}
\description{
GEMMA is an efficient association testing software of Linear Mixed Models
(LMM)s and related GWAS models. This interface is for the Multivariate
linear mixed model (mvLMM) functionality which also corrects for population
structure and sample (non)exchangability
}
