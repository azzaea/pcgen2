#' Residuals from GBLUP
#'
#' Residuals from the best linear unbiased predictor of the genetic effects
#' (GBLUP), which is computed given REML-estimates of the variance components.
#' The residuals are used in \code{pcgenFast} and \code{pcRes}
#'
#' If \code{cov.method = "uni"}, the GBLUP and the residuals are computed
#' separately for each trait in suffStat. The covariance of each trait is then
#' assumed to be \deqn{\sigma_G^2 Z K Z^t + \sigma_E^2 I_n} where \eqn{Z} is a
#' binary incidence matrix, assigning plants or plots to genotypes. \eqn{Z} is
#' based on the first column in \code{suffStat}. If there is a single
#' observation per genotype (typically a genotypic mean), \eqn{Z} is the
#' identity matrix, and the relatedness matrix \eqn{K} should be specified. If
#' there are replicates for at least some of the genotypes, and no \eqn{K} is
#' provided, independent genetic effects are assumed (\eqn{K} will be the
#' identity matrix). It is also possible to have replicates and specify a
#' non-diagonal \eqn{K}. Whenever \eqn{K} is specified, sommer (mmer2) will be
#' used; otherwise lmer (lme4). The mmer2 is also used when \code{cov.method =
#' "us"}, in which case the multivariate GBLUP is computed, for all traits in
#' \code{suffStat} simultaneously. This is only possible for a limited number of
#' traits.
#'
#' @references Covarrubias-Pazaran, G., 2016. Genome-assisted prediction of
#'   quantitative traits using the R package sommer. \emph{PloS one}, 11(6),
#'   p.e0156744.
#'
#' @author Willem Kruijer and Pariya Behrouzi. Maintainers: Willem Kruijer
#'   \email{willem.kruijer@wur.nl} and Pariya Behrouzi
#'   \email{pariya.behrouzi@gmail.com}
#'
#' @inheritParams pcgen
#'
#' @param cov.method A string, specifying which method should be used to compute
#'   the G-BLUP. Options are \code{"us"} (unstructured multi-trait model fitted
#'   using \code{sommer}) and \code{"uni"} (based on univariate GBLUPs).
#'   (Default \code{"uni"})
#'
#' @return A data-frame with the residuals
#'
#' @examples
#' \dontrun{
#' data(simdata)
#' rs <- getResiduals(suffStat= simdata)
#' }
#'
#'
#' @export
#'
#' @importFrom sommer mmer vs unsm
#' @importFrom lme4 lmer VarCorr
#' @importFrom Matrix Matrix

getResiduals <- function(suffStat, covariates=NULL, cov.method = 'uni',
                         K = NULL) {

  #suffStat = dr[, 1:10]; covariates = dr[, 37:38]; cov.method = 'uni'; K = NULL

# cov.method :
#           'us' (fit an unstructured multi-trait model using sommer)
#           and 'uni' (univariate : single trait models with lme4)

# to do: * further (i.e. nongenetic) covariates
#        * incorporate K

# suffStat = d; cov.method = 'us'; covariates = NULL

    stopifnot(cov.method %in% c('uni','us'))

    names(suffStat)[1] <- 'G'

    max.rep <- max(table(suffStat$G))

    if (max.rep==1 & is.null(K)) {
      if (!is.null(covariates)) {
        covariates <- NULL
        cat('covariates is set to NULL, as the data seem to contain genotypic means', '\n')
      }
    }

    #if (max(as.numeric(unlist(lapply(cbind(suffStat, covariates), function(x){sum(is.na(x))})))) > 0.001) {
    if (max(as.numeric(unlist(lapply(suffStat, function(x){sum(is.na(x))})))) > 0.001) {
      if (cov.method == 'us') {stop('Missing values detected. Choose cov.method = uni')}
    }

    p           <- ncol(suffStat)

    res         <- matrix(NA, nrow(suffStat), p-1)

    if (p==2 & cov.method!='uni') {
      cov.method <- 'uni'
      cat('Warning: there seems to be only 1 trait; cov.method is put to uni')
    }

    if (!is.null(covariates)) {
      covariates        <- as.data.frame(covariates)
      names(covariates) <- paste0('CoVaRiaTe_',1:ncol(covariates))
      stopifnot(nrow(covariates)==nrow(suffStat))
      suffStat          <- cbind(suffStat, covariates)
      cv                <- paste('~', paste(names(covariates), collapse ='+'), '+')
      cv2               <- cv
      cv3               <- paste0('cbind(',paste(names(suffStat)[2:p], collapse=','),')~', paste(names(covariates), collapse ='+'))
    } else {
      cv                <- '~'
      cv2               <- '~ 1'
      cv3               <- paste0('cbind(',paste(names(suffStat)[2:p], collapse=','),')~ 1')
    }

    if (cov.method != 'uni' | !is.null(K)) {

      Z.t <- as.matrix(make.Z.matrix(suffStat$G))

      # define Ks, the kinship to be used by sommer
      if (is.null(K)) {
        Ks <- round(Z.t %*% t(Z.t))
      } else {
        Ks <- Z.t %*% as.matrix(K) %*% t(Z.t)
      }

      rownames(Ks) <- colnames(Ks) <- rownames(suffStat) <- paste0(suffStat$G,'_', 1:nrow(suffStat))

      suffStat$id <- rownames(suffStat)

    }



    if (cov.method=='us') {

      #out <- mmer2(fixed = as.formula(cv3), random=~us(trait):g(id),
      #                     rcov=~us(trait):units, data=suffStat, draw = F,
      #                     silent = T, G = list(id = Ks))

      out <- mmer(fixed  = as.formula(cv3),
                  random = ~ vs(id, Gu = Ks, Gtc = unsm(ncol(suffStat) -2) ),
                  rcov   = ~ vs(units, Gtc = unsm(ncol(suffStat) -2)),
                  data   =  suffStat)

      if (out$convergence == TRUE) {

        gb <- matrix(NA, nrow(suffStat), p-1)
        for (j in 2:p) {
          gb[,j-1] <- out$U[[1]][[j-1]][suffStat$id]
        }

        #res <- out$residuals - gb
        res <- suffStat[,2:p] - gb
      } else {
        cov.method <- 'uni'
        cat('Warning: no convergence. cov.method is set to uni','\n')
      }
    }

    if (cov.method=='uni') {

      if (is.null(K)) {

        lmer.cv <- paste(cv, '(1|G)')

        for (j in 2:p) {

          out      <- lmer(as.formula(paste0(names(suffStat)[j], lmer.cv)),
                           data = suffStat)

          res[which(!is.na(suffStat[,j])), j - 1] <- as.numeric(residuals(out))

        }

      } else {

        for (j in 2:p) {

          #out <- mmer2(fixed = as.formula(paste0(names(suffStat)[j], cv2)),
          #             random = ~ g(id),
          #             data = suffStat, draw = F,
          #             silent = T, G = list(id = Ks))

          out <- mmer(fixed  = as.formula(paste0(names(suffStat)[j], cv2)),
                      random = ~vs(id, Gu = Ks),
                      rcov   = ~vs(units), data=suffStat)

          #res[which(!is.na(suffStat[,j])), j - 1] <- as.numeric(out$residuals) - out$U[[1]][[1]][suffStat$id]
          res[which(!is.na(suffStat[,j])), j - 1] <-  suffStat[,j]- out$U[[1]][[1]][suffStat$id]
        }
      }
      rownames(res) <- rownames(suffStat)
    }

    colnames(res) <- names(suffStat)[2:p]

return(as.data.frame(res))
}

# if (cov.method=='fe') {
#
#   Vg <- Ve <- matrix(0, p - 1, p - 1)
#
#   if (is.null(covariates)) {
#     cv <- '~ (1|G)'
#     X.t <- Matrix(rep(1, nrow(suffStat)))
#   } else {
#     stopifnot(nrow(suffStat)==nrow(covariates))
#     covariates        <- as.data.frame(covariates)
#     names(covariates) <- paste0('covariate_', 1:ncol(covariates))
#     suffStat          <- cbind(suffStat, covariates)
#     cv <- paste('~',paste(names(covariates), collapse ='+'), '+ (1|G)')
#     X.t <- Matrix(cbind(rep(1, nrow(suffStat)), as.matrix(covariates)))
#   }
#
#   for (j in 1:(p-1)) {
#
#     out      <- lmer(as.formula(paste0(names(suffStat)[j+1], cv)), data = suffStat)
#
#     temp <- as.data.frame(VarCorr(out))[,'vcov']
#
#     Vg[j,j] <- temp[1]
#
#     Ve[j,j] <- temp[2]
#
#   }
#
#   for (j in 1:(p-2)) {
#
#     for (k in (j+1):(p-1)) {
#
#       em.vec        <- c(suffStat[,j+1], suffStat[,k+1])
#       names(em.vec) <- rep(as.character(suffStat$G), 2)
#
#       out <- fitEM(y = em.vec, X.t = X.t, Z.t = Z.t, max.iter = 300)
#
#       Vg[j,k] <- Vg[k,j] <- out$variances$Vg[3]
#
#       Ve[j,k] <- Ve[k,j] <- out$variances$Ve[3]
#
#     }
#   }
#
#   Vg <- as.matrix(nearPD(Vg, keepDiag=T)$mat)
#
#   Ve <- as.matrix(nearPD(Ve, keepDiag=T, maxit=500)$mat)
#
#   ###
#
#   if (!is.null(covariates)) {
#     cv.sommer <- paste0('cbind(',paste(names(suffStat)[2:p], collapse=','),')~', paste(names(covariates), collapse ='+'))
#   } else {
#     cv.sommer <- paste0('cbind(',paste(names(suffStat)[2:p], collapse=','),')~ 1')
#   }
#
#   out <- mmer2(fixed = as.formula(cv.sommer), random = ~ G, data=suffStat,
#                draw = F, silent = T, forced = list(`g(id)` = Vg, units = Ve))
#
#   res <- out$cond.residuals
# }
