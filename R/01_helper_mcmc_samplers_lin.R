#' Method for only linear effects Gibbs sampler.
#'
#' @inheritParams sample_all_params
#'
#' @export
sample_all_params.lin <- function(pe, mm) {
  for (d in 1:pe$DD) {
    id_betz_tmp <- (pe$id_bet_z[d] + 1):pe$id_bet_z[d + 1]
    id_zet_tmp  <- (pe$id_zet[d] + 1):pe$id_zet[d + 1]

    dd_range_nn <- pe$dd_list_nn[[d]]
    Xtmp <- as.matrix(pe$X[, d, mm - 1, ])
    Ztmp <- pe$Z[, id_zet_tmp, , drop = FALSE]

    pe$sig_sq_x[d, mm] <- sample_sig_sq_x_lin(bet_z = pe$bet_z[id_betz_tmp,
                                                               mm - 1],
                                              X = Xtmp,
                                              regs_z = Ztmp,
                                              prior_ig = c(pe$prior_ig_a,
                                                           pe$prior_ig_b),
                                              iter_range_NN = dd_range_nn,
                                              TT = pe$TT)
    beta_sampled <- sample_bet_z_lin(sig_sq_x = pe$sig_sq_x[d, mm],
                                     X = Xtmp,
                                     regs_z = Ztmp,
                                     TT = pe$TT,
                                     id_reg_z = c(pe$id_reg_z[d],
                                                  pe$id_reg_z[d + 1]),
                                     dim_bet_z = pe$dim_bet_z[d],
                                     prior_vcm_bet_z = pe$prior_vcm_bet_z[[d]],
                                     iter_range_NN = dd_range_nn)

    pe$bet_z[id_betz_tmp, mm] <- beta_sampled

    pe$Regs_beta[, d, ] <- get_regs_beta_lin(Z  = pe$Z[, id_zet_tmp, ,
                                                       drop = FALSE],
                                             TT = pe$TT,
                                             pe$bet_z[id_betz_tmp, mm],
                                             iter_range_NN = 1:pe$NN)
  }
  cat("MCMC iteration number:", mm, "\n")
}
#' Draws a (particle) Gibbs sample of the std. deviation parameter
#'
#' The standard deviation is per component \code{d} of the DD-dimensional
#' latent state process and drawn from the inverse Gamma. This version is
#' adjusted to only incorporate linear regressors and not random effects.
#'
#' @inheritParams sample_sig_sq_x_alr
#' @export
sample_sig_sq_x_lin <- function(bet_z,
                                X,
                                regs_z,
                                prior_ig,
                                iter_range_NN,
                                TT) {
  err_sig_sq_x_all <- 0
  for(n in iter_range_NN) {
    tmp_regs <- regs_z[, , n]
    err_sig_sq_x <- X[, n]- f(x_tt = rep(0, times = TT),
                              regs  = tmp_regs,
                              phi_x = 0,
                              bet_reg = bet_z)
    err_sig_sq_x_all <- err_sig_sq_x_all + crossprod(err_sig_sq_x)
  }
  out <- 1/stats::rgamma(n = 1,
                         prior_ig[1],
                         prior_ig[2] + err_sig_sq_x_all/2)
  return(out)
}
#' Draws a (particle) Gibbs sample of the covariance matrix of random effects
#'
#' This is the linear-only reg-type sampler.
#'
#' @param inheritParams sample_bet_all
#'
#' @return a sample
#' @export
sample_bet_z_lin <- function(sig_sq_x,
                             X,
                             regs_z,
                             TT,
                             id_reg_z,
                             dim_bet_z,
                             prior_vcm_bet_z,
                             iter_range_NN,
                             order_p = 1) {
  browser()
  vcm_x_errors     <- diag(rep(sig_sq_x^(-1), times = TT))
  # vcm_x_errors2    <- compute_vcm_x_errors(vcm_x_errors_rhs, NULL,
  #                                          NULL, TT= TT, type = "lin")
  omega_tmp_all <- 0
  mu_tmp_all    <- 0

  for (n in iter_range_NN) {
    regs_tmp  <- regs_z[, , n]
    omega_tmp <- crossprod(regs_tmp,
                           regs_tmp) / sig_sq_x
    # omega_tmp <- crossprod(regs_tmp, vcm_x_errors) %*% regs_tmp
    omega_tmp_all <- omega_tmp_all + omega_tmp
    mu_tmp <- crossprod(regs_tmp,  X[, n, drop = FALSE]) / sig_sq_x
    # mu_tmp2 <- crossprod(regs_tmp, vcm_x_errors) %*% X[, n, drop = FALSE]
    mu_tmp_all <- mu_tmp_all + mu_tmp
  }
  Omega_bet <- solveme(omega_tmp_all + prior_vcm_bet_z)
  mu_bet    <- Omega_bet %*% mu_tmp_all
  out <- rnorm_fast_n1(mu = mu_bet, Sigma = Omega_bet, dim_bet_z)
  # out   <- MASS::mvrnorm(n = 1, mu = mu_bet, Sigma = Omega_bet)
  return(out)
}
