\name{simdata}
\alias{simdata}
\docType{data}
\title{Simulated data}

\description{
Simulated data, for two replicates of genotypes g1,...,g200.
Three traits were simulated (Y1, Y2 and Y3), using a structural equation model
defined by Y1 -> Y2 -> Y3, and direct genetic effects on Y1 and Y3.
}

\usage{data(simdata)}
\format{
  A data frame of dimension \eqn{ 4 \times 400 }. The first column is the factor G (genotype); the subsequent columns contain \eqn{y_1, y_2} and \eqn{y_3}.
}

\examples{
\donttest{
data(simdata)
out <- pcgen(simdata)
out2 <- pcRes(suffStat = simdata, alpha = 0.01, verbose= FALSE)
}
}

\keyword{datasets}
