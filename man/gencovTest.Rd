\name{gencovTest}
\alias{gencovTest}
\title{Estimate genetic covariances between all pairs of traits, and test their significance}
\usage{
  gencovTest(suffStat, max.iter = 200, out.cor = TRUE)
}

\arguments{
  \item{suffStat}{A data.frame with (p + 1) columns, of which the first column is the factor G (genotype),  and subsequent p columns contain traits. It should not contain covariates or QTLs.}

  \item{max.iter}{Maximum number of iterations in the EM-algorithm, used to fit the bivariate mixed model}.
  
  \item{out.cor}{If \code{TRUE}, the output will contain estimates of genetic correlations; otherwise covariances. The pvalues are always for genetic covariance.
  }
  
}

\value{
  A list with elements pvalues and out.cor, which are both p x p matrices
}

\description{
  For each pair of traits in suffStat, we fit a bivariate mixed model, and 
  perform a likelihood ratio test for the null-hypothesis of zero genetic covariance.
}


\references{
  Kruijer, W., Behrouzi, P., Rodriguez-Alvarez, M. X., Wit, E. C., Mahmoudi, S. M., Yandell, B., Van Eeuwijk, F., (2018, in preparation), Reconstruction of networks with direct and indirect genetic effects. \cr
}

\author{
  Willem Kruijer and Pariya Behrouzi.
  Maintainers: Willem Kruijer \email{willem.kruijer@wur.nl} and
  Pariya Behrouzi \email{pariya.behrouzi@gmail.com}
}

\examples{
  \donttest{
    data(simdata)
    test <- gencovTest(suffStat= simdata, max.iter = 200, out.cor= TRUE )
  }
}
