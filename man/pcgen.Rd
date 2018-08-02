\name{pcgen}
\alias{pcgen}
\title{Causal inference with genetic effects}
\usage{
pcgen(suffStat, covariates = NULL, QTLs = integer(), alpha = 0.01, m.max = Inf,
fixedEdges = NULL, fixedGaps = NULL, verbose = FALSE, use.res = FALSE, 
res.cor = NULL, max.iter = 50, stop.if.significant = TRUE)
}

\arguments{
\item{suffStat}{A data.frame, of which the first column is the factor genotype,  and subsequent columns contain the traits, and optionally some QTLs. The name of the first column should be genotype. Should not contain covariates.}

\item{covariates}{A data.frame containing covariates, that should always be used in each conditional independence test Should be either \code{NULL} (default) or a data.frame with the same number of rows as \code{suffStat}. An intercept is already included for each trait in suffStat; covariates should not contain a column of ones.}

\item{QTLs}{Column numbers in \code{suffStat} that correspond to QTLs.}

\item{alpha}{The significance level used in each conditional independence test. Default is 0.01}

\item{m.max}{Maximum size of the conditioning sets}

\item{fixedEdges}{A logical matrix of dimension \eqn{(p+1) \times (p+1)}, where \eqn{p} is the number of traits. The first row and column refer to the node genotype, and subsequent rows and columns to the traits. As in the pcalg package, the edge \eqn{i - j} is never considered for removal if the entry \eqn{[i, j]} or \eqn{[j, i]} (or both) are \code{TRUE}. In that case, the edge is guaranteed to be present in the resulting graph.}

\item{fixedGaps}{A logical matrix of dimension \eqn{(p+1) \times (p+1)}, where \eqn{p} is the number of traits. The first row and column refer to the node genotype, and subsequent rows and columns to the traits. As in the pcalg package, the edge \eqn{i - j} is removed before starting the algorithm if the entry \eqn{[i, j]} or \eqn{[j, i]} (or both) are \code{TRUE}. In that case, the edge is guaranteed to be absent in the resulting graph.}

\item{verbose}{If \code{TRUE}, p-values for the conditional independence tests are printed}

\item{use.res}{If \code{TRUE}, the test for conditional independence of 2 traits given a set of other traits and \eqn{G} is based on residuals from GBLUP. If \code{FALSE} (the default), it is based on bivariate mixed models.
}

% SHOULD BE REMOVED !!
\item{res.cor}{If \code{use.res = TRUE}, res.cor should be the correlation matrix of the residuals from the GBLUP. These can be obtained with the getResiduals function. See the example below.}

\item{max.iter}{Maximum number of iterations in the EM-algorithm, used to fit the bivariate mixed model (when \code{use.res = FALSE)}. }

\item{stop.if.significant}{If \code{TRUE}, the EM-algorithm used in some of the conditional independence tests (when \code{use.res = FALSE}) will be stopped whenever the p-value becomes significant, i.e. below \code{alpha}. This will speed up calculations, and can be done because (1) the PC algorithm only needs an accept/reject decision (2) In EM the likelihood is nondecreasing. Should be put to \code{FALSE} if the precise p-values are of interest.}
}

\value{
a graph (an object with S3 class \code{"pcgen"})
}
\description{
Reconstruction of directed networks with random genetic effects, based on phenotypic observations. The pcgen algorithm is a modification of the pc-stable algorithm of Colombo & Maathuis (2014) . It is assumed that there are replicates, and independent genetic effects.
}

\details{
The \code{pcgen} function is based on the \code{pc} function from the pcalg package (Kalisch et al. (2012) and Hauser and Buhlmann (2012)). 
}
\references{
1. Kruijer, W., Behrouzi, P., Rodriguez-Alvarez, M. X., Wit, E. C., Mahmoudi, S. M., Yandell, B., Van Eeuwijk, F., (2018, in preparation), Reconstruction of networks with direct and indirect genetic effects. \cr
2. Colombo, D. and Maathuis, M.H., 2014. Order-independent constraint-based causal structure learning. \emph{The Journal of Machine Learning Research}, 15(1), pp.3741-3782. \cr
3. Kalisch, M., Machler, M., Colombo, D., Maathuis, M.H. and Buhlmann, P., 2012. Causal inference using graphical models with the R package pcalg. \emph{Journal of Statistical Software}, 47(11), pp.1-26. \cr
4. Hauser, A. and Buhlmann, P., 2012. Characterization and greedy learning of interventional Markov equivalence classes of directed acyclic graphs. Journal of Machine Learning Research, 13(Aug), pp.2409-2464.\cr
}
\author{
Willem Kruijer and Pariya Behrouzi.
        Maintainers: Willem Kruijer \email{willem.kruijer@wur.nl} and
       Pariya Behrouzi \email{pariya.behrouzi@gmail.com}
}

\seealso{\code{\link{getResiduals}}}

\examples{
\donttest{
data(simdata)
out <- pcgen(simdata) 
}
data(simdata)
rs <- getResiduals(suffStat = simdata)
pc.fit1 <- pcgen(suffStat = simdata, alpha = 0.01, verbose = TRUE, 
           use.res = TRUE, res.cor = cor(rs))
}