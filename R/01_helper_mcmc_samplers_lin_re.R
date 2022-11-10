#' Method for linear and random effects Gibbs sampler.
#'
#' @inheritParams sample_all_params
#'
#' @export
sample_all_params.lin_re <- function(pe, mm) {
  if (is.null(pe$phi_x)) pe$phi_x <- 0
  for (d in 1:pe$DD) {
    id_betz_tmp <- (pe$id_bet_z[d] + 1):pe$id_bet_z[d + 1]
    id_betu_tmp <- (pe$id_bet_u[d] + 1):pe$id_bet_u[d + 1]
    id_zet_tmp  <- (pe$id_zet[d] + 1):pe$id_zet[d + 1]
    id_uet_tmp  <- (pe$id_uet[d] + 1):pe$id_uet[d + 1]

    dd_range_nn <- pe$dd_list_nn[[d]]
    Xtmp <- as.matrix(pe$X[, d, mm - 1, ])
    Ztmp <- pe$Z[2:pe$TT, id_zet_tmp, , drop = FALSE]
    Utmp <- pe$U[2:pe$TT, id_uet_tmp, , drop = FALSE]

    pe$sig_sq_x[d, mm] <- sample_sig_sq_x_lin_re(phi_x = pe$phi_x,
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

    pe$vcm_bet_u[[d]][, , mm] <- sample_vcm_bet_u(pe$bet_u[id_betu_tmp, mm - 1, ,
                                                           drop = FALSE],
                                                  pe$dim_bet_u[d],
                                                  pe$dof_vcm_bet_u[d],
                                                  pe$prior_vcm_bet_u2[[d]],
                                                  dd_range_nn)
    # pe$vcm_bet_u[[d]][, , mm] <- c(27.71419)[d]
    pe$bet_u[id_betu_tmp, mm,
             dd_range_nn] <- sample_bet_u(pe$sig_sq_x[d, mm],
                                          pe$phi_x,
                                          pe$bet_z[id_betz_tmp,
                                                   mm - 1],
                                          pe$vcm_bet_u[[d]][, , mm],
                                          pe$dim_bet_u[d],
                                          Xtmp, Ztmp, Utmp,
                                          dd_range_nn,
                                          pe$TT)
    beta_sampled <- sample_bet_z_lin_re(sig_sq_x = pe$sig_sq_x[d, mm],
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

    # pe$phi_x[d, mm] <- beta_sampled[1]
    # pe$bet_z[id_betz_tmp, mm] <- beta_sampled[-1]
    pe$bet_z[id_betz_tmp, mm] <- beta_sampled

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
  }
  cat("MCMC iteration number:", mm, "\n")
}
#' Draws a (particle) Gibbs sample of the std. deviation parameter
#'
#' The standard deviation is per component \code{d} of the DD-dimensional
#' latent state process and drawn from the inverse Gamma. This version is
#' adjusted to incorporate linear regressors as well as random effects.
#'
#' @inheritParams sample_sig_sq_x_auto_lin_re
#'
#' @return one sample from the inverse gamma distribution for the standard
#'   deviation parameter
#' @export
sample_sig_sq_x_lin_re <- function(phi_x,
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
#' Draws a (particle) Gibbs sample of the covariance matrix of random effects
#'
#' This is the linear and random effects regressor-type sampler so some details:
#' The covariance matrix of random effects is per component \code{d} of the
#' DD-dimensional latent state process and drawn from the inverse Wishart.
#'
#' @param inheritParams sample_bet_z_lin_re
#'
#' @return a sample
#' @export
sample_bet_z_lin_re <- function(sig_sq_x,
                                vcm_bet_u,
                                X,
                                regs_z,
                                U,
                                TT,
                                id_reg_z,
                                dim_bet_z,
                                prior_vcm_bet_z,
                                iter_range_NN) {
  browser()
  vcm_x_errors_rhs     <- diag(rep(sig_sq_x, times = TT - 1))
  # vmc_x_errors_rhs_inv <- solveme(vcm_x_errors_rhs)
  # vcm_bet_u_inv        <- solveme(vcm_bet_u)
  omega_tmp_all <- 0
  mu_tmp_all    <- 0
  x_lhs     <- X[1:(TT - 1), , drop = FALSE]
  x_rhs_all <- X[2:TT, , drop = FALSE]
  for (n in iter_range_NN) {
    # regs_z[, id_reg_z[1] + 1, n] <- x_lhs[, n]

    # Umat             <- matrix(U[,, n, drop = FALSE], nrow = TT - 1)
    # vcm_x_errors_lhs <- Umat %*% vcm_bet_u %*% t(Umat)
    # vcm_x_errors2     <- vcm_x_errors_lhs + vcm_x_errors_rhs
    # vcm_x_errors2     <- solveme(vcm_x_errors2)
    vcm_x_errors    <- compute_vcm_x_errors(vcm_x_errors_rhs,
                                            matrix(U[,, n, drop = FALSE],
                                                   nrow = TT - 1),
                                            vcm_bet_u, type = "lin_re")

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
  browser()
  Omega_bet <- solveme(omega_tmp_all + prior_vcm_bet_z)
  mu_bet    <- Omega_bet %*% mu_tmp_all
  out <- rnorm_fast_n1(mu = mu_bet, Sigma = Omega_bet, dim_bet_z + 1)
  # out   <- MASS::mvrnorm(n = 1, mu = mu_bet, Sigma = Omega_bet)
  return(out)
}