#' Method for only linear effects Gibbs sampler.
#'
#' @inheritParams sample_all_params
#'
#' @export
sample_all_params.auto_lin <- function(pe, mm) {
  order_p <- pe$order_p
  for (d in 1:pe$DD2) {
    id_phi_tmp  <- (pe$id_phi[d] + 1):pe$id_phi[d + 1]
    id_betz_tmp <- (pe$id_zet[d] + 1):pe$id_zet[d + 1]
    id_zet_tmp  <- (pe$id_zet[d] + 1):pe$id_zet[d + 1]

    # id_regs_z_tmp <- (pe$id_reg_z[d] + 1 + (d - 1) * order_p):(pe$id_reg_z[d + 1] + order_p * d)
    id_regs_z_tmp <- (pe$id_reg_z[d] + 1):(pe$id_reg_z[d + 1])

    dd_range_nn <- pe$dd_list_nn[[d]]
    Xtmp <- as.matrix(pe$X[, d, mm - 1, ])
    Ztmp <- pe$Z[(1 + order_p):pe$TT, id_zet_tmp, , drop = FALSE]

    pe$sig_sq_x[d, mm] <- sample_sig_sq_x_al(phi_x = pe$phi_x[id_phi_tmp,
                                                              mm - 1],
                                             bet_z = pe$bet_z[id_betz_tmp,
                                                              mm - 1],
                                             X = Xtmp,
                                             regs_z = Ztmp,
                                             prior_ig = c(pe$prior_ig_a,
                                                          pe$prior_ig_b),
                                             iter_range_NN = dd_range_nn,
                                             TT = pe$TT,
                                             order_p)
    beta_sampled <- sample_bet_z_al(sig_sq_x = pe$sig_sq_x[d, mm],
                                    X = Xtmp,
                                    regs_z =  pe$regs_z[, id_regs_z_tmp, ,
                                                        drop = FALSE],
                                    TT = pe$TT,
                                    dim_bet_z = pe$dim_bet_z[d],
                                    prior_vcm_bet_z = pe$prior_vcm_bet_z[[d]],
                                    iter_range_NN = dd_range_nn,
                                    order_p)

    pe$phi_x[id_phi_tmp, mm]  <- beta_sampled[1:order_p]
    pe$bet_z[id_betz_tmp, mm] <- beta_sampled[-(1:order_p)]

    pe$Regs_beta[, d, ] <- get_regs_beta_l(Z  = pe$Z[, id_zet_tmp, ],
                                           TT = pe$TT,
                                           pe$bet_z[id_betz_tmp, mm, drop = FALSE],
                                           iter_range_NN = 1:pe$NN)
  }
  cat("MCMC iteration number:", mm, "\n")
}
#' Draws a (particle) Gibbs sample of the std. deviation parameter
#'
#' The standard deviation is per component \code{d} of the DD-dimensional latent
#' state process and drawn from the inverse Gamma. This version is adjusted to
#' incorporate linear regressors and autoregressive states processes.
#'
#' @inheritParams sample_sig_sq_x_alr
#'
#' @return one sample from the inverse gamma distribution for the standard
#'   deviation parameter
#' @export
sample_sig_sq_x_al <- function(phi_x,
                               bet_z,
                               X,
                               regs_z,
                               prior_ig,
                               iter_range_NN,
                               TT,
                               order_p) {
  err_sig_sq_x_all <- 0

  x_lhs <- X[(order_p + 1):TT, , drop = FALSE]
  x_rhs <- get_x_rhs(X, order_p, TT)

  for(n in iter_range_NN) {
    x_rhs_tmp    <- matrix(x_rhs[, , n, drop = FALSE],
                           nrow = TT - order_p, ncol = order_p)
    err_sig_sq_x <- x_lhs[, n] - (x_rhs_tmp %*% phi_x + regs_z[, , n, drop = FALSE] %*% bet_z)

    err_sig_sq_x_all <- err_sig_sq_x_all + sum(err_sig_sq_x ^ 2)
  }
  out <- 1/stats::rgamma(n = 1,
                         prior_ig[1],
                         prior_ig[2] + err_sig_sq_x_all/2)
  return(out)
}
#' Draws a (particle) Gibbs sample of the covariance matrix of random effects
#'
#' This is the autoregressive- and linear effects regressor-type sampler.
#'
#' @inheritParams sample_bet_z_alr
#'
#' @export
sample_bet_z_al <- function(sig_sq_x,
                            X,
                            regs_z,
                            TT,
                            dim_bet_z,
                            prior_vcm_bet_z,
                            iter_range_NN,
                            order_p) {
  omega_tmp_all <- 0
  mu_tmp_all    <- 0

  x_rhs_all <- get_x_rhs(X, order_p, TT)
  regs_z[, 1:order_p, ] <- x_rhs_all
  for (n in iter_range_NN) {
    regs_tmp      <- matrix(regs_z[, , n, drop = FALSE], nrow = TT - order_p)

    omega_tmp     <-  crossprod(regs_tmp, regs_tmp) / sig_sq_x
    omega_tmp_all <- omega_tmp_all + omega_tmp

    mu_tmp     <- crossprod(regs_tmp, X[(order_p + 1):TT,
                                        n, drop = FALSE]) / sig_sq_x
    mu_tmp_all <- mu_tmp_all + mu_tmp
  }
  Omega_bet <- solveme(omega_tmp_all + prior_vcm_bet_z)
  mu_bet    <- Omega_bet %*% mu_tmp_all
  out <- rnorm_fast_n1(mu = mu_bet, Sigma = Omega_bet, dim_bet_z + order_p)
  # out   <- MASS::mvrnorm(n = 1, mu = mu_bet, Sigma = Omega_bet)
  while (check_stationarity(out[1:order_p], order_p)) {
    # out <- MASS::mvrnorm(n = 1, mu = mu_bet, Sigma = Omega_bet)
    out <- rnorm_fast_n1(mu = mu_bet, Sigma = Omega_bet, dim_bet_z + order_p)
  }
  return(out)
}
