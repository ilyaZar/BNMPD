#' Method for only linear effects Gibbs sampler.
#'
#' @inheritParams sample_all_params
#'
#' @export
sample_all_params.auto_re <- function(pe, mm) {
  order_p <- pe$order_p
  for (d in 1:pe$DD) {
    id_phi_tmp  <- (pe$id_phi[d] + 1):pe$id_phi[d + 1]
    id_betu_tmp <- (pe$id_bet_u[d] + 1):pe$id_bet_u[d + 1]
    id_uet_tmp  <- (pe$id_uet[d] + 1):pe$id_uet[d + 1]

    dd_range_nn <- pe$dd_list_nn[[d]]
    Xtmp <- as.matrix(pe$X[, d, mm - 1, ])
    Utmp <- pe$U[(1 + order_p):pe$TT, id_uet_tmp, , drop = FALSE]

    pe$sig_sq_x[d, mm] <- sample_sig_sq_x_ar(phi_x = pe$phi_x[id_phi_tmp,
                                                              mm - 1],
                                             bet_u = pe$bet_u[id_betu_tmp,
                                                              mm - 1,,
                                                              drop = FALSE],
                                             X = Xtmp,
                                             regs_u = Utmp,
                                             prior_ig = c(pe$prior_ig_a,
                                                          pe$prior_ig_b),
                                             iter_range_NN = dd_range_nn,
                                             TT = pe$TT,
                                             order_p)
    pe$vcm_bet_u[[d]][, , mm] <- sample_vcm_bet_u(pe$bet_u[id_betu_tmp,
                                                           mm - 1, ,
                                                           drop = FALSE],
                                                  pe$dim_bet_u[d],
                                                  pe$dof_vcm_bet_u[d],
                                                  pe$prior_vcm_bet_u2[[d]],
                                                  dd_range_nn)
    pe$bet_u[id_betu_tmp, mm,
             dd_range_nn] <- sample_bet_u_ar(pe$sig_sq_x[d, mm],
                                             pe$phi_x[id_phi_tmp, mm - 1],
                                             pe$vcm_bet_u[[d]][, , mm],
                                             pe$dim_bet_u[d],
                                             Xtmp, Utmp,
                                             dd_range_nn,
                                             pe$TT,
                                             order_p)
    beta_sampled <- sample_bet_z_ar(sig_sq_x = pe$sig_sq_x[d, mm],
                                    vcm_bet_u = pe$vcm_bet_u[[d]][, , mm],
                                    X = Xtmp,
                                    U = Utmp,
                                    TT = pe$TT,
                                    prior_vcm_bet_z = pe$prior_vcm_bet_z[[d]],
                                    iter_range_NN = dd_range_nn,
                                    order_p)

    pe$phi_x[id_phi_tmp, mm]  <- beta_sampled[1:order_p]

    pe$Regs_beta[, d, ] <- get_regs_beta_r(U = pe$U,
                                           id_uet = c(pe$id_uet[d],
                                                      pe$id_uet[d + 1]),
                                           TT = pe$TT,
                                           bet_u = pe$bet_u[, mm, ,
                                                            drop = FALSE],
                                           id_bet_u = c(pe$id_bet_u[d],
                                                        pe$id_bet_u[d + 1]),
                                           iter_range_NN = 1:pe$NN)
  }
  cat("MCMC iteration number:", mm, "\n")
}
#' Draws a (particle) Gibbs sample of the std. deviation parameter
#'
#' The standard deviation is per component \code{d} of the DD-dimensional latent
#' state process and drawn from the inverse Gamma. This version is adjusted to
#' incorporate linear random effect regressors and autoregressive states
#' processes.
#'
#' @inheritParams sample_sig_sq_x_alr
#'
#' @return one sample from the inverse gamma distribution for the standard
#'   deviation parameter
#' @export
sample_sig_sq_x_ar <- function(phi_x,
                               bet_u,
                               X,
                               regs_u,
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
    err_sig_sq_x <- x_lhs[, n] - (x_rhs_tmp %*% phi_x +
                                    regs_u[, , n] %*%  bet_u[ , , n])
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
sample_bet_z_ar <- function(sig_sq_x,
                            vcm_bet_u,
                            X,
                            U,
                            TT,
                            prior_vcm_bet_z,
                            iter_range_NN,
                            order_p) {
  omega_tmp_all <- 0
  mu_tmp_all    <- 0

  x_rhs_all <- get_x_rhs(X, order_p, TT)
  for (n in iter_range_NN) {
    regs_tmp <- matrix(x_rhs_all[, , n, drop = FALSE], nrow = TT - order_p)
    omega_tmp <-  compute_vcm_WBinv(sig_sq_x = sig_sq_x,
                                    matLHS = regs_tmp,
                                    U = matrix(U[, , n, drop = FALSE],
                                               nrow = TT - order_p),
                                    C = vcm_bet_u,
                                    matRHS =  regs_tmp)
    omega_tmp_all <- omega_tmp_all + omega_tmp
    mu_tmp <- compute_vcm_WBinv(sig_sq_x = sig_sq_x,
                                matLHS =  regs_tmp,
                                U = matrix(U[, , n, drop = FALSE],
                                           nrow = TT - order_p),
                                C = vcm_bet_u,
                                matRHS = X[(order_p + 1):TT, n, drop = FALSE])
    mu_tmp_all <- mu_tmp_all + mu_tmp
  }
  Omega_bet <- solveme(omega_tmp_all + prior_vcm_bet_z)
  mu_bet    <- Omega_bet %*% mu_tmp_all
  out <- rnorm_fast_n1(mu = mu_bet, Sigma = Omega_bet, order_p)
  # out   <- MASS::mvrnorm(n = 1, mu = mu_bet, Sigma = Omega_bet)
  while (check_stationarity(out[1:order_p], order_p)) {
    # out <- MASS::mvrnorm(n = 1, mu = mu_bet, Sigma = Omega_bet)
    out <- rnorm_fast_n1(mu = mu_bet, Sigma = Omega_bet, order_p)
  }
  return(out)
}
#' Sampling random effects under autoregressive only latent processes.
#'
#' This version neglects linear z-type regressors and only considers random
#' effects and autoregressive states of order p per component d.
#'
#' @inheritParams sample_bet_u_alr
#'
#' @return a sample
sample_bet_u_ar <- function(sig_sq_x,
                            phi_x,
                            vcm_bet_u,
                            dim_bet_u,
                            X,
                            U,
                            iter_range_NN,
                            TT,
                            order_p) {
  out_mat          <- matrix(0, nrow = dim_bet_u, ncol = length(iter_range_NN))
  vcm_bet_u_inv    <- solveme(vcm_bet_u)

  nn <- 1
  x_lhs <- X[(order_p + 1):TT, , drop = FALSE]
  x_rhs <- get_x_rhs(X, order_p, TT)
  for (n in iter_range_NN) {
    x_rhs_tmp <- matrix(x_rhs[, , n], nrow = TT - order_p, ncol = order_p)
    Omega_bet_u <- matrix(0, nrow = dim_bet_u, ncol = dim_bet_u)
    mu_bet_u    <- matrix(0, nrow = dim_bet_u, ncol = 1)

    x_n   <- x_lhs[, n] - (x_rhs_tmp %*% phi_x)
    Umat <- matrix(U[, , n, drop = FALSE], nrow = TT - order_p)

    Omega_bet_u <- solveme(crossprod(Umat, Umat)/sig_sq_x + vcm_bet_u_inv)
    mu_bet_u    <- Omega_bet_u %*% (crossprod(Umat, x_n) / sig_sq_x)

    out_mat[, nn] <- rnorm_fast_n1(mu = mu_bet_u,
                                   Sigma = Omega_bet_u,
                                   dim_bet_u)
    nn <- nn + 1
  }
  return(out_mat)
}
