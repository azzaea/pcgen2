#' @importFrom Matrix determinant
#'
deviance2 <- function(C, Gi, ng, Ri, N, Rinv, res, edf) {
	log_det_C <- determinant(C)$modulus
	log_det_G <- ng*determinant(Gi)$modulus
	log_det_R <- N*determinant(Ri)$modulus
	deviance  <- log_det_C + log_det_G + log_det_R + t(res)%*%Rinv%*%res + edf
	deviance
}
