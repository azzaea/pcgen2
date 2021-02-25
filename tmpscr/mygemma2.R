mygemma2<- function (max_iter = 10000, max_prec = 1/1e+03, eval, X, Y, V_g,
                     V_e, verbose_output = FALSE, reduced) {
  n_size <- length(eval)
  c_size <- nrow(X)
  d_size <- nrow(Y)
  dc_size <- c_size * d_size
  XXt <- X %*% t(X)
  logl_const <- -(n_size - c_size) * d_size * log(2 * pi)/2 +
    d_size * log(det(XXt))/2
  out <- list()
  for (t in 1:max_iter) {
    ep_out <- eigen_proc(V_g, V_e)
    logdet_Ve <- ep_out[[1]]
    UltVeh <- ep_out[[2]]
    UltVehi <- ep_out[[3]]
    D_l <- ep_out[[4]]
    cq_out <- calc_qi(eval, D_l, X)
    Qi <- cq_out[[1]]
    lndetQ <- cq_out[[2]]
    UltVehiY <- UltVehi %*% Y
    xHiy <- calc_XHiY(eval, D_l, X, UltVehiY)
    logl_new <- logl_const + gemma2:::MphCalcLogL(eval = eval, xHiy = xHiy,
                                         D_l = D_l, UltVehiY = UltVehiY, Qi = Qi) - 0.5 *
      n_size * logdet_Ve
    logl_new <- logl_new - 0.5 * (lndetQ - c_size * logdet_Ve)
    if (t > 1) {
      if (logl_new - logl_old < max_prec) {
        break
      }
    }
    logl_old <- logl_new
    co_out <- calc_omega(eval, D_l)
    OmegaU <- co_out[[1]]
    OmegaE <- co_out[[2]]
    UltVehiB <- gemma2:::UpdateRL_B(xHiy, Qi, d_size = nrow(Y))
    UltVehiBX <- UltVehiB %*% X
    UltVehiU <- update_u(OmegaE, UltVehiY, UltVehiBX)
    UltVehiE <- gemma2:::update_e(UltVehiY, UltVehiBX, UltVehiU)
    U_hat <- t(UltVeh) %*% UltVehiU
    E_hat <- t(UltVeh) %*% UltVehiE
    B <- t(UltVeh) %*% UltVehiB
    cs_out <- gemma2:::calc_sigma(eval = eval, D_l = D_l, X = X, OmegaU = OmegaU,
                         OmegaE = OmegaE, UltVeh = UltVeh, Qi = Qi)
    Sigma_ee <- cs_out[[1]]
    Sigma_uu <- cs_out[[2]]
    uv_out <- gemma2:::update_v(eval, U_hat, E_hat, Sigma_uu, Sigma_ee)
    V_g <- uv_out[[1]];
    V_e <- uv_out[[2]];
    if (reduced) {
      V_g <- diag(diag(V_g))
      V_e <- diag(diag(V_e))
    }
    if (verbose_output) {
      out[[t]] <- list(logl_new = logl_new, Vg = V_g, Ve = V_e,
                       logl_const = logl_const, it = t
                       )
    }
    else {
      out[[1]] <- list(logl_new = logl_new, Vg = V_g, Ve = V_e,
                       Sigma_uu = Sigma_uu, Sigma_ee = Sigma_ee, B = B,
                       U_hat = U_hat, E_hat = E_hat, OmegaU = OmegaU,
                       OmegaE = OmegaE, logdet_Ve = logdet_Ve, UltVeh = UltVeh,
                       UltVehi = UltVehi, Dl = D_l, xHiy = xHiy, logl_const = logl_const,
                       UltVehiU = UltVehiU, it = t)
    }
  }
  if (length(out) == max_iter) {
    warning("gemma2 EM algorithm didn't converge.")
  }
  return(out)
}
