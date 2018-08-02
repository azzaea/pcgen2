\name{getResiduals}
\alias{getResiduals}
\title{Residuals from the GBLUP}
\usage{ getResiduals(suffStat, covariates = NULL, cov.method = "uni", K = NULL) }

\arguments{
\item{suffStat}{A data.frame, of which the first column is the factor genotype,  and subsequent columns contain the traits, and optionally some QTLs. The name of the first column should be genotype. }

\item{covariates}{A data.frame containing covariates, that should always be used in each conditional independence test. Should be either \code{NULL} (default) or a data.frame with the same number of rows as suffStat. An intercept is already included for each trait in suffStat; covariates should not contain a column of ones.}

\item{cov.method}{(A string, specifying which method should be used to compute the GBLUP. Options are \code{"us"} (unstructured multi-trait model fitted using sommer) and \code{"uni"} (based on univariate GBLUPs). Default is \code{"uni"}). }

\item{K}{A genetic relatedness matrix. If \code{NULL} (default), independent genetic effects are assumed.}
}

\value{ A data-frame with the residuals.}

\description{Residuals from the best linear unbiased predictor of the genetic effects (GBLUP), which is computed given REML-estimates of the variance components.
}

\details{
If \code{cov.method = "uni"}, the GBLUP and the residuals are computed separately for each trait in suffStat. The covariance of each trait is then assumed to be \deqn{\sigma_G^2 Z K Z^t + \sigma_E^2 I_n} where \eqn{Z} is a binary incidence matrix, assigning plants or plots to genotypes. \eqn{Z} is based on the first column in \code{suffStat}. If there is a single observation per genotype (typically a genotypic mean), \eqn{Z} is the identity matrix, and the relatedness matrix \eqn{K} should be specified. If there are replicates for at least some of the genotypes, and no \eqn{K} is provided, independent genetic effects are assumed (\eqn{K} will be the identity matrix). It is also possible to have replicates and specify a non-diagonal \eqn{K}.
Whenever \eqn{K} is specified, sommer (mmer2) will be used; otherwise lmer (lme4). The mmer2 is also used when \code{cov.method = "us"}, in which case the multivariate GBLUP is computed, for all traits in \code{suffStat} simultaneously. This is only possible for a limited number of traits.
}

\examples{
%\donttest{
data(simdata)
rs <- getResiduals(suffStat= simdata)
%}
}

\references{
Covarrubias-Pazaran, G., 2016. Genome-assisted prediction of quantitative traits using the R package sommer. PloS one, 11(6), p.e0156744. \cr
}

\author{
Willem Kruijer and Pariya Behrouzi.
        Maintainers: Willem Kruijer \email{willem.kruijer@wur.nl} and
       Pariya Behrouzi \email{pariya.behrouzi@gmail.com}
}