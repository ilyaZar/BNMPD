#' Method for linear and random effects Gibbs sampler.
#'
#' @inheritParams sample_all_params
#'
#' @export
sample_all_params.auto_lin_re <- function(pe, mm) {
  for (d in 1:pe$DD) {
    browser()
    id_betz_tmp <- (pe$id_bet_z[d] + 1):pe$id_bet_z[d + 1]
    id_betu_tmp <- (pe$id_bet_u[d] + 1):pe$id_bet_u[d + 1]
    id_zet_tmp  <- (pe$id_zet[d] + 1):pe$id_zet[d + 1]
    id_uet_tmp  <- (pe$id_uet[d] + 1):pe$id_uet[d + 1]

    dd_range_nn <- pe$dd_list_nn[[d]]
    Xtmp <- as.matrix(pe$X[, d, mm - 1, ])
    Ztmp <- pe$Z[2:pe$TT, id_zet_tmp, , drop = FALSE]
    Utmp <- pe$U[2:pe$TT, id_uet_tmp, , drop = FALSE]

    pe$sig_sq_x[d, mm] <- sample_sig_sq_x(phi_x = pe$phi_x[d, mm - 1],
                                          bet_z = pe$bet_z[id_betz_tmp, mm - 1],
                                          bet_u = pe$bet_u[id_betu_tmp, mm - 1,,
                                                           drop = FALSE],
                                          X = Xtmp,
                                          regs_z = Ztmp, regs_u = Utmp,
                                          prior_ig = c(pe$prior_ig_a,
                                                       pe$prior_ig_b),
                                          iter_range_NN = dd_range_nn,
                                          TT = pe$TT)
    # pe$sig_sq_x[d, mm] <- c(0.21, 0.32, 0.43)[d]

    pe$vcm_bet_u[[d]][, , mm] <- sample_vcm_bet_u(pe$bet_u[id_betu_tmp, mm, ,
                                                           drop = FALSE],
                                                  pe$dim_bet_u[d],
                                                  pe$dof_vcm_bet_u[d],
                                                  pe$prior_vcm_bet_u2[[d]],
                                                  dd_range_nn)
    # pe$vcm_bet_u[[d]][, , mm] <- c(0.2771419, 0.02868432, 0.07917361)[d]
    pe$bet_u[id_betu_tmp, mm,
             dd_range_nn] <- sample_bet_u_auto_lin_re(pe$sig_sq_x[d, mm],
                                          pe$phi_x[d, mm - 1],
                                          pe$bet_z[id_betz_tmp,
                                                   mm - 1],
                                          pe$vcm_bet_u[[d]][, , mm],
                                          pe$dim_bet_u[d],
                                          Xtmp, Ztmp, Utmp,
                                          dd_range_nn,
                                          pe$TT)
    beta_sampled <- sample_bet_z(sig_sq_x = pe$sig_sq_x[d, mm],
                                 vcm_bet_u = pe$vcm_bet_u[[d]][, , mm],
                                 X = Xtmp,
                                 regs_z = pe$regs_z,
                                 U = Utmp,
                                 TT = pe$TT,
                                 id_reg_z = c(pe$id_reg_z[d],
                                              pe$id_reg_z[d + 1]),
                                 dim_bet_z = pe$dim_bet_z[d],
                                 prior_vcm_bet_z = pe$prior_vcm_bet_z[[d]],
                                 iter_range_NN = dd_range_nn)

    pe$phi_x[d, mm] <- beta_sampled[1]
    pe$bet_z[id_betz_tmp, mm] <- beta_sampled[-1]

    # pe$phi_x[d, mm] <- c(0.35, 0.55, 0.75)[d]
    # pe$phi_x[d, mm] <- c(0.0, 0.0, 0.0)[d]
    # pe$bet_z[id_betz_tmp, mm] <- list(c(-2.5, 3.0),
    #                                   c(2.0, -4.0),
    #                                   c(0.4, -0.7))[[d]]

    pe$Regs_beta[, d, ] <- get_regs_beta(Z  = pe$Z[, id_zet_tmp, ],
                                         U = pe$U,
                                         id_uet = c(pe$id_uet[d],
                                                    pe$id_uet[d + 1]),
                                         TT = pe$TT,
                                         pe$bet_z[id_betz_tmp, mm],
                                         bet_u = pe$bet_u[, mm, ,
                                                          drop = FALSE],
                                         id_bet_u = c(pe$id_bet_u[d],
                                                      pe$id_bet_u[d + 1]),
                                         iter_range_NN = 1:pe$NN)
    # Z_beta[,d, ]    <- regs_bet[[1]]
    # U_beta[,d, ]    <- regs_bet[[2]]
    # Regs_beta[,d, ] <- regs_bet[[3]]
  }
  # pe$phi_x[, mm] <- phi_x[, mm]
  # pe$bet_z[, mm] <- bet_z[, mm]
  # pe$bet_u[, mm, ] <- bet_u[, mm, ]
  # pe$sig_sq_x[, mm] <- sig_sq_x[, mm]
  # for (d in 1:DD) {
  #   pe$vcm_bet_u[[d]][, , mm] <- [[d]][, , mm]
  # }
  # pe$Regs_beta <- Regs_beta
  cat("MCMC iteration number:", mm, "\n")
}
#' Draws a (particle) Gibbs sample of the covariance matrix of random effects
#'
#' This is the autoregressive-linear and random effects regressor-type sampler
#' so some details: The covariance matrix of random effects is per component
#' \code{d} of the DD-dimensional latent state process and drawn from the
#' inverse Wishart.
#'
#' @param sig_sq_x m'th sample of standard deviation parameter in a component
#'   \code{d} of latent state process
#' @param vcm_bet_u covariance matrix of the random effects
#' @param X state matrix sliced over PGAS component \code{m} and state component
#'   \code{d}
#' @param regs_z Z regressors sliced over corresponding \code{d} component
#' @param U regressor matrix of random effects
#' @param TT number of timer periods
#' @param id_reg_z a vector of 2; Z regressor id for components \code{d} and
#'   \code{d + 1}
#' @param prior_vcm_bet_z covariance matrix prior for the variance of the
#'   \code{beta_z} coefficient
#' @param iter_range_NN iteration range i.e. the cross sectional components
#'   that are actually contributing to \code{d}
#' @param order_p integer; autoregressive lag order where \code{order_p = 0}
#'   reduces the sampling to a standard (non-dynamic) linear regression model
#'
#' @return a sample
#' @export
sample_bet_z_auto_lin_re <- function(sig_sq_x,
                                     vcm_bet_u,
                                     X,
                                     regs_z,
                                     U,
                                     TT,
                                     id_reg_z,
                                     dim_bet_z,
                                     prior_vcm_bet_z,
                                     iter_range_NN,
                                     order_p = 1) {
  browser()
  # vcm_x_errors_rhs     <- diag(rep(sig_sq_x, times = TT - 1))
  # vmc_x_errors_rhs_inv <- solveme(vcm_x_errors_rhs)
  # vcm_bet_u_inv        <- solveme(vcm_bet_u)
  omega_tmp_all <- 0
  mu_tmp_all    <- 0
  # if (m %in% (c(2,30))) browser()
  x_lhs     <- X[1:(TT - order_p), , drop = FALSE]
  x_rhs_all <- X[(1 + order_p):TT, , drop = FALSE]
  browser()
  if (order_p == 1) regs_z[, id_reg_z[1] + 1,] <- x_lhs
  if (order_p == 0) regs_z <- regs_z[, -c(id_reg_z[1] + 1), ]
  for (n in iter_range_NN) {
    # regs_z[, id_reg_z[1] + 1, n] <- x_lhs[, n]

    # Umat             <- matrix(U[,, n, drop = FALSE], nrow = TT - 1)
    # vcm_x_errors_lhs <- Umat %*% vcm_bet_u %*% t(Umat)
    # vcm_x_errors2     <- vcm_x_errors_lhs + vcm_x_errors_rhs
    # vcm_x_errors2     <- solveme(vcm_x_errors2)
    vcm_x_errors    <- compute_vcm_x_errors_inv(sig_sq_x,
                                                matrix(U[,, n, drop = FALSE],
                                                       nrow = TT - order_p),
                                                vcm_bet_u,
                                                TT - order_p,
                                                type = "lin_re")

    regs_tmp  <- regs_z[, (id_reg_z[1] + 1):id_reg_z[2], n]
    # omega_tmp <- crossprod(regs_tmp,
    #                        regs_tmp)/sig_sq_x[d, m]
    omega_tmp <- crossprod(regs_tmp,
                           vcm_x_errors) %*% regs_tmp
    omega_tmp_all <- omega_tmp_all + omega_tmp
    # mu_tmp <- crossprod(regs_tmp, x_rhs)/sig_sq_x[d, m]
    mu_tmp <- crossprod(regs_tmp, vcm_x_errors) %*% x_rhs_all[, n]
    mu_tmp_all <- mu_tmp_all + mu_tmp
  }
  # browser()
  Omega_bet <- solveme(omega_tmp_all + prior_vcm_bet_z)
  mu_bet    <- Omega_bet %*% mu_tmp_all
  out <- rnorm_fast_n1(mu = mu_bet, Sigma = Omega_bet, dim_bet_z + 1)
  # out   <- MASS::mvrnorm(n = 1, mu = mu_bet, Sigma = Omega_bet)
  if (order_p >= 1) {
    while ((abs(abs(out[1]) - 1) < 0.01) | abs(out[1]) > 1) {
      # out <- MASS::mvrnorm(n = 1, mu = mu_bet, Sigma = Omega_bet)
      out <- rnorm_fast_n1(mu = mu_bet, Sigma = Omega_bet, dim_bet_z + 1)
    }
  }
  return(out)
}
#' Draws a (particle) Gibbs sample of the std. deviation parameter
#'
#' The standard deviation is per component \code{d} of the DD-dimensional
#' latent state process and drawn from the inverse Gamma. This version is
#' adjusted to incorporate linear regressors as well as random effects.
#'
#' @param phi_x m'th sample of auto-regressive parameter of the state component
#'   \code{d}
#' @param bet_z m'th sample of beta coefficient of the state component \code{d}
#' @param bet_u m'th sample of beta coefficient of the state component \code{d}
#' @param X state matrix sliced over PGAS component \code{m} and state component
#'   \code{d}
#' @param regs_z Z regressors sliced over corresponding \code{d} component
#' @param regs_u U regressors sliced over corresponding \code{d} component
#' @param prior_ig a numeric vector of two: \code{prior_ig_a} and
#'   \code{prior_ig_b}
#' @param iter_range_NN iteration range i.e. the cross sectional components
#'   that are actually contributing to \code{d}
#' @param TT time series length
#'
#' @return one sample from the inverse gamma distribution for the standard
#'   deviation parameter
#' @export
sample_sig_sq_x_auto_lin_re <- function(phi_x,
                                        bet_z,
                                        bet_u,
                                        X,
                                        regs_z,
                                        regs_u,
                                        prior_ig,
                                        iter_range_NN,
                                        TT) {
  err_sig_sq_x_all <- 0
  x_lhs <- X[2:TT, , drop = FALSE]
  x_rhs <- X[1:(TT - 1), , drop = FALSE]
  for(n in iter_range_NN) {
    bet_u_n <- bet_u[ , , n]
    tmp_regs <- cbind(regs_z[, , n], regs_u[, , n])
    err_sig_sq_x <- x_lhs[, n]- f(x_tt = x_rhs[, n],
                                  regs  = tmp_regs,
                                  phi_x = phi_x,
                                  bet_reg = c(bet_z, bet_u_n))
    err_sig_sq_x_all <- err_sig_sq_x_all + crossprod(err_sig_sq_x)
  }
  out <- 1/stats::rgamma(n = 1,
                         prior_ig[1],
                         prior_ig[2] + err_sig_sq_x_all/2)
  return(out)
}
