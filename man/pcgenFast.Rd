\name{pcgenFast}
\alias{pcgenFast}
\title{pcgen with residual-based screening}
\usage{
pcgenFast(suffStat, alpha = 0.01, m.max = Inf, res.m.max = Inf, verbose = FALSE, 
          covariates = NULL, fixedEdges = NULL, QTLs = integer(), max.iter = 50,
          stop.if.significant = TRUE, cov.method = 'uni', use.res = FALSE)
}
\arguments{
\item{suffStat}{A data.frame, of which the first column is the factor  genotype,  and subsequent columns contain the traits, and optionally some QTLs. The name of the first column should be genotype.}

\item{alpha}{The significance level used in each conditional independence test. Default is 0.01. }

\item{m.max}{Maximum size of the conditioning set, in the pcgen algorithm.}

\item{res.m.max}{Maximum size of the conditioning set, in the pc-algorithm on the residuals (used for prior screening).}

\item{verbose}{If \code{TRUE}, p-values for the conditional independence tests are printed. }

\item{covariates}{A data.frame containing covariates, to be used in each conditional independence test. Cannot contain factors. Should be either \code{NULL} (default) or a data.frame with the same number of rows as suffStat. An intercept is already included for each trait in suffStat; covariates should not contain a column of ones.}

\item{fixedEdges}{A logical matrix of dimension \eqn{(p+1) \times (p+1)}, where \eqn{p} is the number of traits. The first row and column refer to the node genotype, and subsequent rows and columns to the traits. As in the \code{pcalg} package, the edge \eqn{i - j} is never considered for removal if the entry \eqn{[i, j]} or \eqn{[j, i]} (or both) are \code{TRUE}. In that case, the edge is guaranteed to be present in the resulting graph.}

\item{QTLs}{Column numbers in \code{suffStat} that correspond to QTLs.}

\item{max.iter}{Maximum number of iterations in the EM-algorithm, used to fit the bivariate mixed model (when \code{use.res = FALSE)}.}

\item{stop.if.significant}{If \code{TRUE}, the EM-algorithm used in some of the conditional independence tests (when \code{use.res = FALSE}) will be stopped whenever the p-value becomes significant, i.e. below \code{alpha}. This will speed up calculations, and can be done because (1) the PC algorithm only needs an accept/reject decision (2) In EM the likelihood is nondecreasing. Should be put to \code{FALSE} if the precise p-values are of interest.}

\item{cov.method}{A string, specifying which method should be used to compute the GBLUP. Options are \code{'us'} (unstructured multi-trait model fitted using sommer) and \code{'uni'} (based on univariate GBLUPs). Default is \code{'uni'}.}

\item{use.res}{If \code{FALSE}, residuals from GBLUP are only used for screening with the standard pc algotihm. After that, the standard pcgen algorithm is run on the remaining edges; the test for conditional independence of 2 traits given a set of other traits and G is based on bivariate mixed models. If \code{TRUE}, this test is
based on the residuals. In this case, no further edges between traits are removed after screening and \code{pcgen} will only infer the orientation, and the direct genetic effects. }
}

\value{
A graph (an object with S3 class \code{"pcgen"})
}

\description{
The pcgen algorithm starting with a skeleton estimated using the standard pc-algorithm, based on residuals from the GBLUP.
}

\references{
1. Kruijer, W., Behrouzi, P., Rodriguez-Alvarez, M. X., Wit, E. C., Mahmoudi, S. M., Yandell, B., Van Eeuwijk, F., (2018, in preparation), Reconstruction of networks with direct and indirect genetic effects. \cr
2. Colombo, D. and Maathuis, M.H., 2014. Order-independent constraint-based causal structure learning. \emph{The Journal of Machine Learning Research}, 15(1), pp.3741-3782. \cr
}

\author{
Willem Kruijer and Pariya Behrouzi.
        Maintainers: Willem Kruijer \email{willem.kruijer@wur.nl} and
       Pariya Behrouzi \email{pariya.behrouzi@gmail.com}
}

\seealso{\code{\link{getResiduals}}}

\examples{
%\donttest{
data(simdata)
out <- pcgenFast(suffStat = simdata, alpha = 0.01, verbose= FALSE, use.res = TRUE)
}
%}

