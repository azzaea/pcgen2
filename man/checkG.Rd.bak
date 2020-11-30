\name{checkG}
\alias{checkG}
\title{Check for consistency in genetic effects}
\usage{
checkG(pcgen.output, suffStat, alpha = 0.01, covariates = NULL)
}

\arguments{
\item{pcgen.output}{A graph with nodes G (genotype) and a number of traits. Typically output from \code{pcgen} or \code{pcgenFast}.}

\item{suffStat}{A data.frame, of which the first column is the factor G (genotype),  and subsequent columns contain the traits, and optionally some QTLs. The name of the first column should be G.}

\item{alpha}{The significance level used in each conditional independence test. Default is 0.01. }

\item{covariates}{A data.frame containing covariates, to be used in each conditional independence test.
Cannot contain factors. Should be either \code{NULL} (default) or a data.frame with the same number of rows as suffStat. An intercept is already included for each trait in suffStat; covariates should not contain a column of ones.}
}

\value{
A logical matrix of dimension \eqn{(p+1) \times (p+1)}, \eqn{p} being the number of traits. Most entries are \code{FALSE}, except those in the first row and column for which there are conflicts.  Entries \eqn{[1, j]} and \eqn{[j, 1]} are \code{TRUE} if the \eqn{j}th trait has significant genetic variance, but there is no         partially directed path from G towards that trait. The matrix can then be used in a subsequent run of \code{pcgen} or \code{pcgenFast}, in the \code{fixedEdges} argument. The arguments \code{suffStat}, \code{alpha} and \code{covariates} should stay the same throughout (first run of \code{pcgen}, \code{checkG}, second run of \code{pcgen}).
}

\description{
Given output from pcgen or pcgenFast, this function checks whether the estimated graph is consistent with the set of traits having significant genetic variance. The function detects traits that have significant genetic variance but for which there is no partially directed path from G.
}

\references{
Kruijer, W., Behrouzi, P., Rodriguez-Alvarez, M. X., Wit, E. C., Mahmoudi, S. M., Yandell, B., Van Eeuwijk, F., (2018, in preparation), Reconstruction of networks with direct and indirect genetic effects. \cr
}

\author{
Willem Kruijer and Pariya Behrouzi.
        Maintainers: Willem Kruijer \email{willem.kruijer@wur.nl} and
       Pariya Behrouzi \email{pariya.behrouzi@gmail.com}
}
