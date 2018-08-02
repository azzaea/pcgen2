\name{pcgenTest}
\alias{pcgenTest}
\title{The conditional independence test in pcgen}
\usage{
pcgenTest(x, y, S, suffStat, QTLs = integer(), covariates = NULL, alpha = 0.01, 
          max.iter = 50, stop.if.significant = TRUE, use.res = FALSE, res.cor = NULL)
}

\arguments{
\item{x, y}{Column numbers in \code{suffStat} that should be tested for conditional independence given the variables in S.}

\item{S}{vector of integers defining the conditioning set, where the integers refer to column numbers in \code{suffStat}. May be numeric(), i.e. the empty set.}

\item{suffStat}{A data.frame, of which the first column is the factor genotype,  and subsequent columns contain the traits, and optionally some QTLs. The name of the first column should be genotype. It should not contain covariates.}

\item{QTLs}{Column numbers in \code{suffStat} that correspond to QTLs. These may be partly in \eqn{S} and \eqn{x} and \eqn{y}, but \eqn{x} and \eqn{y} cannot be both QTLs.}

\item{covariates}{A data.frame containing covariates. It should be either \code{NULL} (default) or a data.frame with the same number of rows as \code{suffStat}. An intercept is already included for each trait in suffStat; covariates should not contain a column of ones.}

\item{alpha}{The significance level used in the test. The test itself of course does not depend on this, but it is used in the EM-algorithm to speed up calculations. When \code{ stop.if.significant = TRUE}, the EM-algorithm is stopped once the p-value is below the significance level. Default is 0.01.}

\item{max.iter}{Maximum number of iterations in the EM-algorithm, used to fit the bivariate mixed model (when \code{use.res = FALSE)}.}

\item{stop.if.significant}{If \code{TRUE}, the EM-algorithm used in some of the conditional independence tests (when \code{use.res = FALSE}) will be stopped whenever the p-value becomes significant, i.e. below \code{alpha}. This will speed up calculations, and can be done because (1) the PC algorithm only needs an accept/reject decision (2) In EM the likelihood is nondecreasing. It should be put to \code{FALSE} if the precise p-value is of interest.}

\item{use.res}{If \code{TRUE}, the test for conditional independence of 2 traits given a set of other traits and G is based on residuals from GBLUP. If \code{FALSE} (the default), it is based on bivariate mixed models.
}

% SHOULD BE REMOVED !!
\item{res.cor}{If \code{use.res = TRUE}, res.cor should be the correlation matrix of the residuals from the GBLUP. These can be obtained with the getResiduals function. See the example below.}
}

\value{
A p-value
}

\description{
This performs the conditional independence test used in the pcgen algorithm, assuming there are replicates, and independent genetic effects.
}

\details{\code{pcgenTest} tests for conditional independence between \eqn{x} and \eqn{y} given \eqn{S}.
It distinguishes 2 situations:
(i) if one of \eqn{x} and \eqn{y} (say \eqn{x}) is the factor genotype, \code{pcgenTest} will test if the genetic variance in \eqn{y} is zero, given the traits in S. (ii) if \eqn{x} and \eqn{y} are both traits, \code{pcgenTest} tests if the residual covariance between them is zero, given the traits in \eqn{S} and the factor genotype. The factor genotype is automatically included in the conditioning set \eqn{S} (\eqn{S} does not need to contain the integer 1). This test is either based on a bivariate mixed model (when \code{use.res=FALSE}), or on residuals from GBLUP (\code{use.res=T}), obtained with the getResiduals function. In the latter case, \code{res.cor} must be provided.
}

\references{
Kruijer, W., Behrouzi, P., Rodriguez-Alvarez, M. X., Wit, E. C., Mahmoudi, S. M., Yandell, B., Van Eeuwijk, F., (2018, in preparation), Reconstruction of networks with direct and indirect genetic effects. \cr
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
rs <- getResiduals(suffStat= simdata)
pcgenTest(suffStat= simdata, x= 2, y= 3, S= 4)
pcgenTest(suffStat= simdata, x= 2, y= 3, S= c(1,4))
pcgenTest(suffStat= simdata, x= 2, y= 3, S= 4, use.res= TRUE, res.cor= cor(rs))
pcgenTest(suffStat= simdata, x= 2, y= 1, S= 4)
%}
}
