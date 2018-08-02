\name{pcRes}
\alias{pcRes}
\title{The pc algorithm applied to residuals}
\usage{
pcRes(suffStat, alpha= 0.01, K = NULL, m.max = Inf, verbose = FALSE,
  covariates = NULL, QTLs = integer(), cov.method = "uni",
  use.GBLUP = FALSE)
}

\arguments{
\item{suffStat}{A data.frame, of which the first column is the factor genotype,  and subsequent columns contain the traits, and optionally some QTLs. The name of the first column should be genotype. }

\item{alpha}{The significance level used in the test. Default is 0.01.}

\item{K}{A genetic relatedness matrix. If \code{NULL} (the default), independent genetic effects are assumed.}

\item{m.max}{Maximum size of the conditioning set, in the pc-algorithm on the residuals. }

\item{verbose}{If \code{TRUE}, p-values for the conditional independence tests are printed. }

\item{covariates}{A data.frame containing covariates, that should always be used in each conditional independence test. Should be either \code{NULL} (default) or a data.frame with the same number of rows as \code{suffStat}. An intercept is already included for each trait in \code{suffStat}; covariates should not contain a column of ones.}

\item{QTLs}{Column numbers in suffStat that correspond to QTLs. }

\item{cov.method}{ A string, specifying which method should be used to compute the GBLUP. Options are \code{'us'} (unstructured multi-trait model fitted using sommer) and \code{'uni'} (based on univariate GBLUPs).}

\item{use.GBLUP}{Use the GBLUP itself, instead of the residuals}
}

\value{
A graph (an object with S3 class \code{"pcgen"})
}

\description{
The standard pc algorithm applied to GBLUP residuals, or to the GBLUP itself.
}

\details{
If \code{use.GBLUP = FALSE}, GBLUP residuals are used as input for the pc-stable algorithm of Colombo and Maathuis (2014). This closely resembles the residual networks of Valente et al., (2010) and Topner et al., (2017) (who used different ways to predict the genetic effects, and applied other causal inference algorithms to the residuals). When \code{use.GBLUP = TRUE}, pc-stable is applied to the GBLUP itself, which resembles the genomic networks of Topner et al., (2017). If \code{cov.method = "uni"}, the GBLUP and the residuals are computed separately for each trait in suffStat. The covariance of each trait is assumed to be
\deqn{\sigma_G^2 Z K Z^t + \sigma_E^2 I_n}
where \eqn{Z} is a binary incidence matrix, assigning plants or plots to genotypes. \eqn{Z} is based on the first column in suffStat. If there is a single observation per genotype (typically a genotypic mean), \eqn{Z} is the identity matrix, and the relatedness matrix \eqn{K} should be specified. If there are replicates
for at least some of the genotypes, and no \eqn{K} is provided, independent genetic effects are assumed (\eqn{K} will be the identity matrix). It is also possible to have replicates and specify a non-diagonal \eqn{K}. Whenever  \eqn{K} is specified, sommer (mmer2) will be used; otherwise lmer (lme4).
mmer2 is also used when \code{cov.method = "us"}, in which case the multivariate GBLUP is computed, for all traits in suffStat simultaneously. This is only possible for a limited number of traits.
}

\examples{
%\donttest{
data(simdata)
out <- pcRes(suffStat = simdata, alpha = 0.01, verbose= FALSE)
%}
}

\references{
1. Colombo, D. and Maathuis, M.H., 2014. Order-independent constraint-based causal structure learning. \emph{The Journal of Machine Learning Research}, 15(1), pp.3741-3782. \cr
2. Kruijer, W., Behrouzi, P., Rodriguez-Alvarez, M. X., Wit, E. C., Mahmoudi, S. M., Yandell, B., Van Eeuwijk, F., (2018, in preparation), Reconstruction of networks with direct and indirect genetic effects. \cr
3. Topner, K., Rosa, G.J., Gianola, D. and Schon, C.C., 2017. Bayesian Networks Illustrate Genomic and Residual Trait Connections in Maize (Zea mays L.). \emph{G3: Genes, Genomes, Genetics}, pp.g3-117. \cr
4. Valente, B.D., Rosa, G.J., Gustavo, A., Gianola, D. and Silva, M.A., 2010. Searching for recursive causal structures in multivariate quantitative genetics mixed models. \emph{Genetics}.
}
\author{
Willem Kruijer and Pariya Behrouzi.
        Maintainers: Willem Kruijer \email{willem.kruijer@wur.nl} and
       Pariya Behrouzi \email{pariya.behrouzi@gmail.com}
}
