#' Draws a (particle) Gibbs sample of the std. deviation parameter
#'
#' The standard deviation is per component \code{d} of the DD-dimensional
#' latent state process and drawn from the inverse Gamma.
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
sample_sig_sq_x <- function(phi_x,
                            bet_z,
                            bet_u,
                            X,
                            regs_z,
                            regs_u,
                            prior_ig,
                            iter_range_NN,
                            TT) {
  err_sig_sq_x_all <- 0
  x_lhs <- X[2:TT, ]
  x_rhs <- X[1:(TT - 1), ]
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
#' The covariance matrix of random effects is per component \code{d} of the
#' DD-dimensional latent state process and drawn from the inverse Wishart.
#'
#' @param bet_u beta coefficient of of the state component \code{d}
#' @param dim_bet_u dimension of the beta_u component
#' @param dof_vcm_bet_u degrees of freedom for the Wishart distribution
#' @param prior_vcm_bet_scl prior of the covariance matrix of the random effects
#' @param iter_range_NN iteration range i.e. the cross sectional components
#'   that are actually contributing to \code{d}
#'
#' @return a sample from the inverse Wishart distribution for the VCM of
#'   \code{beta_u}
#' @export
sample_vcm_bet_u <- function(bet_u,
                             dim_bet_u,
                             dof_vcm_bet_u,
                             prior_vcm_bet_scl,
                             iter_range_NN) {
  # browser()
  scale_mat_vcm_bet_u <- matrix(0, nrow = dim_bet_u, ncol = dim_bet_u)
  for (n in iter_range_NN) {
    bet_u_n <- bet_u[, , n]
    scale_mat_vcm_bet_u <- scale_mat_vcm_bet_u + tcrossprod(bet_u_n)
  }
  # browser()
  scale_mat_vcm_bet_u   <- solve(scale_mat_vcm_bet_u + prior_vcm_bet_scl)
  out <- solve(stats::rWishart(1, dof_vcm_bet_u,
                               scale_mat_vcm_bet_u)[, , 1])
  return(out)
}
#' Draws a (particle) Gibbs sample of the random effects
#'
#' The covariance matrix of random effects is per component \code{d} of the
#' DD-dimensional latent state process and drawn from the inverse Wishart.
#'
#' @param sig_sq_x m'th sample of standard deviation parameter in a component
#'   \code{d} of latent state process
#' @param phi_x m'th sample of auto-regressive parameter of the state component
#'   \code{d}
#' @param bet_z m'th sample of beta coefficient of the state component \code{d}
#' @param vcm_bet_u covariance matrix of the random effects
#' @param dim_bet_u dimension of the beta_u regressors
#' @param X state matrix sliced over PGAS component \code{m} and state component
#'   \code{d}
#' @param regs_z Z regressors sliced over corresponding \code{d} component
#' @param U random effects regressors
#' @param iter_range_NN iteration range i.e. the cross sectional components
#'   that are actually contributing to \code{d}
#' @param TT time series length
#'
#' @return a list of \code{length(iter_range_NN)} containing one sample of the
#'   \code{beta_u} coefficients
#' @export
sample_bet_u <- function(sig_sq_x,
                         phi_x,
                         bet_z,
                         vcm_bet_u,
                         dim_bet_u,
                         X,
                         regs_z,
                         U,
                         iter_range_NN,
                         TT) {

  out_mat         <- matrix(0, nrow = dim_bet_u, ncol = length(iter_range_NN))
  vcm_x_errors     <- diag(rep(sig_sq_x, times = TT - 1))
  vmc_x_errors_inv <- solve(vcm_x_errors)
  vcm_bet_u_inv    <- solve(vcm_bet_u)

  nn <- 1
  for (n in iter_range_NN) {
    Omega_bet_u <- matrix(0, nrow = dim_bet_u, ncol = dim_bet_u)
    mu_bet_u    <- matrix(0, nrow = dim_bet_u, ncol = 1)

    x_lhs <- X[2:TT, n]
    x_rhs <- X[1:(TT - 1), n]
    x_n   <- x_lhs - f(x_tt = x_rhs,
                       regs  = regs_z[, , n],
                       phi_x = phi_x,
                       bet_reg = bet_z)
    Umat <- matrix(U[, , n, drop = FALSE], nrow = TT - 1)
    Omega_bet_u <- crossprod(Umat, vmc_x_errors_inv) %*% Umat + vcm_bet_u_inv
    Omega_bet_u <- solve(Omega_bet_u)
    mu_bet_u    <- Omega_bet_u %*% (crossprod(Umat, vmc_x_errors_inv) %*% x_n)

    out_mat[, nn] <- rnorm_fast_n1(mu = mu_bet_u,
                                   Sigma = Omega_bet_u,
                                   dim_bet_u)
    nn <- nn + 1
  }
  return(out_mat)
}
#' Draws a (particle) Gibbs sample of the covariance matrix of random effects
#'
#' The covariance matrix of random effects is per component \code{d} of the
#' DD-dimensional latent state process and drawn from the inverse Wishart.
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
#'   \code{beta_z} coefficinet
#' @param iter_range_NN iteration range i.e. the cross sectional components
#'   that are actually contributing to \code{d}
#'
#' @return a sample
#' @export
sample_beta_all <- function(sig_sq_x,
                            vcm_bet_u,
                            X,
                            regs_z,
                            U,
                            TT,
                            id_reg_z,
                            dim_bet_z,
                            prior_vcm_bet_z,
                            iter_range_NN) {
  vcm_x_errors_rhs     <- diag(rep(sig_sq_x, times = TT - 1))
  vmc_x_errors_rhs_inv <- solve(vcm_x_errors_rhs)
  vcm_bet_u_inv        <- solve(vcm_bet_u)
  omega_tmp_all <- 0
  mu_tmp_all    <- 0
  # if (m %in% (c(2,30))) browser()
  x_lhs     <- X[1:(TT - 1), ]
  x_rhs_all <- X[2:TT, ]
  for (n in iter_range_NN) {
    regs_z[, id_reg_z[1] + 1, n] <- x_lhs[, n]
    x_rhs                     <- x_rhs_all[, n]

    Umat             <- matrix(U[,, n, drop = FALSE], nrow = TT - 1)
    vcm_x_errors_lhs <- Umat %*% vcm_bet_u %*% t(Umat)
    vcm_x_errors     <- vcm_x_errors_lhs + vcm_x_errors_rhs
    vcm_x_errors     <- solve(vcm_x_errors)

    regs_tmp  <- regs_z[, (id_reg_z[1] + 1):id_reg_z[2], n]
    # omega_tmp <- crossprod(regs_tmp,
    #                        regs_tmp)/sig_sq_x[d, m]
    omega_tmp <- crossprod(regs_tmp,
                           vcm_x_errors) %*% regs_tmp
    omega_tmp_all <- omega_tmp_all + omega_tmp
    # mu_tmp <- crossprod(regs_tmp, x_rhs)/sig_sq_x[d, m]
    mu_tmp <- crossprod(regs_tmp, vcm_x_errors) %*% x_rhs
    mu_tmp_all <- mu_tmp_all + mu_tmp
  }
  Omega_bet <- solve(omega_tmp_all + prior_vcm_bet_z)
  mu_bet    <- Omega_bet %*% mu_tmp_all
  out <- rnorm_fast_n1(mu = mu_bet, Sigma = Omega_bet, dim_bet_z + 1)
  # out   <- MASS::mvrnorm(n = 1, mu = mu_bet, Sigma = Omega_bet)
  while ((abs(abs(out[1]) - 1) < 0.01) | abs(out[1]) > 1) {
    # out <- MASS::mvrnorm(n = 1, mu = mu_bet, Sigma = Omega_bet)
    out <- rnorm_fast_n1(mu = mu_bet, Sigma = Omega_bet, dim_bet_z + 1)
  }
  return(out)
}
sample_all <- function(pe, mm) {
  for (d in 1:pe$DD) {
    id_betz_tmp <- (pe$id_bet_z[d] + 1):pe$id_bet_z[d + 1]
    id_betu_tmp <- (pe$id_bet_u[d] + 1):pe$id_bet_u[d + 1]
    id_zet_tmp  <- (pe$id_zet[d] + 1):pe$id_zet[d + 1]
    id_uet_tmp  <- (pe$id_uet[d] + 1):pe$id_uet[d + 1]

    dd_range_nn <- pe$dd_list_nn[[d]]

    pe$sig_sq_x[d, mm] <- sample_sig_sq_x(phi_x = pe$phi_x[d, mm - 1],
                                          bet_z = pe$bet_z[id_betz_tmp, mm - 1],
                                          bet_u = pe$bet_u[id_betu_tmp, mm - 1,,
                                                           drop = FALSE],
                                          X = pe$X[, d, mm - 1, ],
                                          regs_z = pe$Z[2:pe$TT, id_zet_tmp, ],
                                          regs_u = pe$U[2:pe$TT, id_uet_tmp, ,
                                                        drop = FALSE],
                                          prior_ig = c(pe$prior_ig_a,
                                                       pe$prior_ig_b),
                                          iter_range_NN = dd_range_nn,
                                          TT = pe$TT)

    pe$vcm_bet_u[[d]][, , mm] <- sample_vcm_bet_u(pe$bet_u[id_betu_tmp, mm, ,
                                                           drop = FALSE],
                                                  pe$dim_bet_u[d],
                                                  pe$dof_vcm_bet_u[d],
                                                  pe$prior_vcm_bet_u2[[d]],
                                                  dd_range_nn)
    pe$bet_u[id_betu_tmp, mm,
             dd_range_nn] <- sample_bet_u(pe$sig_sq_x[d, mm],
                                          pe$phi_x[d, mm - 1],
                                          pe$bet_z[id_betz_tmp,
                                                   mm - 1],
                                          pe$vcm_bet_u[[d]][, , mm],
                                          pe$dim_bet_u[d],
                                          pe$X[, d, mm - 1, ],
                                          pe$Z[2:pe$TT, id_zet_tmp, ],
                                          pe$U[2:pe$TT, id_uet_tmp, ,
                                               drop = FALSE],
                                          dd_range_nn,
                                          pe$TT)
    beta_sampled <- sample_beta_all(sig_sq_x = pe$sig_sq_x[d, mm],
                                    vcm_bet_u = pe$vcm_bet_u[[d]][, , mm],
                                    X = pe$X[, d, mm - 1, ],
                                    regs_z = pe$regs_z,
                                    U = pe$U[2:pe$TT, id_uet_tmp, ,
                                             drop = FALSE],
                                    TT = pe$TT,
                                    id_reg_z = c(pe$id_reg_z[d],
                                                 pe$id_reg_z[d + 1]),
                                    dim_bet_z = pe$dim_bet_z[d],
                                    prior_vcm_bet_z = pe$prior_vcm_bet_z[[d]],
                                    iter_range_NN = dd_range_nn)

    pe$phi_x[d, mm] <- beta_sampled[1]
    pe$bet_z[id_betz_tmp, mm] <- beta_sampled[-1]

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
}
#
#
#
#
#
# phi_x[d, m] <- beta_sampled[1]
# bet_z[(id_bet_z[d] + 1):id_bet_z[d + 1], m] <- beta_sampled[-1]
# [[d]][, , m] <-
#  sample_vcm_bet_u(bet_u = bet_u[(id_bet_u[d] + 1):id_bet_u[d + 1], m, ],
#                                dim_bet_u = dim_bet_u[d],
#                                dof_vcm_bet_u = dof_vcm_bet_u[d],
#                                prior_vcm_bet_scl = prior_vcm_bet_u2[[d]],
#                                 iter_range_NN = 1:NN)
# monitor_mcmc <- function(states_true, states_drawn) {
# }
# tmp_scale_mat_vcm_bet_u <- matrix(0, nrow = dim_bet_u[d], ncol = dim_bet_u[d])
# for (n in 1:NN) {
#   tmp_scale_mat_vcm_bet_u <- tmp_scale_mat_vcm_bet_u + tcrossprod(bet_u[id_betu_tmp, m, n])
# }
# tmp_scale_mat_vcm_bet_u <- solve(tmp_scale_mat_vcm_bet_u + prior_vcm_bet_u2[[d]])
# set.seed(42)
# vcm_bet_u[[d]][, , m]   <- solve(stats::rWishart(1, dof_vcm_bet_u[d],
#                                                  tmp_scale_mat_vcm_bet_u)[, , 1])
# browser()
# set.seed(42)
#
# #
# #
# #
# vcm_x_errors_rhs[[d]] <- diag(rep(sig_sq_x[d, m], times = TT - 1))
# vmc_x_errors_rhs_inv  <- solve(vcm_x_errors_rhs[[d]])
# vcm_bet_u_inv         <- solve(vcm_bet_u[[d]][, , m])
# #
# #
# #
# #
# #
# #
# # set.seed(42) bet_z[id_betz_tmp, m] <- beta_sampled[-1]
# for (n in 1:NN) {
#   Omega_bet_u <- matrix(0, nrow = dim_bet_u[d], ncol = dim_bet_u[d])
#   mu_bet_u    <- matrix(0, nrow = dim_bet_u[d], ncol = 1)
#   x_tilde_n   <- X[2:TT, d, m - 1, n] - f(x_tt = X[1:(TT - 1), d, m - 1, n],
#                                           regs  = Z[2:TT, id_zet_tmp, n],
#                                           phi_x = phi_x[d, m - 1],
#                                           bet_reg = bet_z[id_betz_tmp, m - 1])
#   Umat <- matrix(U[2:TT, id_uet_tmp, n, drop = FALSE], nrow = TT - 1)
#   Omega_bet_u <- crossprod(Umat, vmc_x_errors_rhs_inv) %*% Umat + vcm_bet_u_inv
#   Omega_bet_u <- solve(Omega_bet_u)
#   mu_bet_u    <- Omega_bet_u %*% (crossprod(Umat, vmc_x_errors_rhs_inv) %*% x_tilde_n)
#
#   bet_u[id_betu_tmp, m, n] <- rnorm_fast_n1(mu = mu_bet_u, Sigma = Omega_bet_u, dim_bet_u[d])
# }
#
#
#
#
#
# monitor_pgas_states(states_drawn = cbind(exp(X1[1, ]), exp(X2[1, ]),
#                                          exp(X3[1, ]), exp(X4[1, ]),
#                                          exp(X5[1, ]), exp(X6[1, ])),
#                     states_comp = cbind(exp(states_init_1), exp(states_init_2),
#                                          exp(states_init_3), exp(states_init_4),
#                                          exp(states_init_5), exp(states_init_6)),
#                       # NULL,
#                       # cbind(xa1_t, xa2_t,
#                       #                    xa3_t, xa4_t,
#                       #                    xa5_5, xa6_t),
#                     current = 1, total = 1, num_prints = 1)
# monitor_pgas_states(states_drawn = cbind(exp(X1[1, ]), exp(X2[1, ]),
#                                          exp(X3[1, ]), exp(X4[1, ]),
#                                          exp(X5[1, ]), exp(X6[1, ])),
#                     states_comp  = cbind(xa1_t, xa2_t,
#                                          xa3_t, xa4_t,
#                                          xa5_t, xa6_t),
#                     #            cbind(exp(states_init_1), exp(states_init_2),
#                     #                  exp(states_init_3), exp(states_init_4),
#                     #                  exp(states_init_5), exp(states_init_6)),
#                     current = 1, total = 1, num_prints = 1)
# monitor_pgas_states(states_drawn = cbind(exp(X1[m, ]), exp(X2[m, ]),
#                                          exp(X3[m, ]), exp(X4[m, ]),
#                                          exp(X5[m, ]), exp(X6[m, ])),
#                     # states_comp = cbind(exp(states_init_1), exp(states_init_2),
#                     #                      exp(states_init_3), exp(states_init_4),
#                     #                      exp(states_init_5), exp(states_init_6)),
#                     states_comp = cbind(exp(X1[m - 1, ]), exp(X2[m - 1, ]),
#                                         exp(X3[m - 1, ]), exp(X4[m - 1, ]),
#                                         exp(X5[m - 1, ]), exp(X6[m - 1, ])),
#                     # NULL,
#                     # cbind(xa1_t, xa2_t, xa3_t,
#                     #                    xa4_t, xa5_t, xa6_t),
#                     current = m, total = MM,
#                     num_prints = num_plots_states)
# monitor_pgas_time(m, MM, len = MM)
# monitor_pgas_mcmc2(m, MM, len = MM,
#                    val_init = par_init,
#                    current_pars = cbind(sig_sq_x1[1:m], phi_x1[1:m],
#                                         t(bet_z1)[1:m,],
#                                         sig_sq_x2[1:m], phi_x2[1:m],
#                                         t(bet_z2)[1:m,],
#                                         sig_sq_x3[1:m], phi_x3[1:m],
#                                         t(bet_z3)[1:m,],
#                                         sig_sq_x4[1:m], phi_x4[1:m],
#                                         t(bet_z4)[1:m,],
#                                         sig_sq_x5[1:m], phi_x5[1:m],
#                                         t(bet_z5)[1:m,],
#                                         sig_sq_x6[1:m], phi_x6[1:m],
#                                         t(bet_z6)[1:m,]),
#                    dim_all = dim_all)
# monitor_pgas_mcmc2(m, MM, len = MM,
#                    val_true = par_true,
#                    val_init = par_init,
#                    current_pars = cbind(sig_sq_x1[1:m], phi_x1[1:m],
#                                         t(bet_z1)[1:m,],
#                                         sig_sq_x2[1:m], phi_x2[1:m],
#                                         t(bet_z2)[1:m,],
#                                         sig_sq_x3[1:m], phi_x3[1:m],
#                                         t(bet_z3)[1:m,],
#                                         sig_sq_x4[1:m], phi_x4[1:m],
#                                         t(bet_z4)[1:m,],
#                                         sig_sq_x5[1:m], phi_x5[1:m],
#                                         t(bet_z5)[1:m,],
#                                         sig_sq_x6[1:m], phi_x6[1:m],
#                                         t(bet_z6)[1:m,]),
#                    dim_all = dim_all)
# vcm_x_errors_rhs[[d]] <- diag(rep(sig_sq_x[d, m], times = TT - 1))
#       vmc_x_errors_rhs_inv  <- solve(vcm_x_errors_rhs[[d]])
#       vcm_bet_u_inv         <- solve(vcm_bet_u[[d]][, , m])
#       omega_tmp_all <- 0
#       mu_tmp_all <- 0
#       # if (m %in% (c(2,30))) browser()
#       for (n in 1:NN) {
#         regs_z[, id_reg_z[d] + 1, n]  <- X[1:(TT - 1), d, m - 1, n]
#         x_lhs        <- X[2:TT, d, m - 1, n]
#
#         Umat <- matrix(U[2:TT, id_uet_tmp, n, drop = FALSE], nrow = TT - 1)
#         vcm_x_errors_lhs[[d]][, , n] <- Umat %*% vcm_bet_u[[d]][, , m] %*% t(Umat)
#         vcm_x_errors          <- vcm_x_errors_lhs[[d]][, , n] + vcm_x_errors_rhs[[d]]
#         vcm_x_errors          <- solve(vcm_x_errors)
#
#         regs_tmp              <- regs_z[, (id_reg_z[d] + 1):id_reg_z[d + 1], n]
#         # omega_tmp <- crossprod(regs_tmp,
#         #                        regs_tmp)/sig_sq_x[d, m]
#         omega_tmp <- crossprod(regs_tmp,
#                                vcm_x_errors) %*% regs_tmp
#         omega_tmp_all <- omega_tmp_all + omega_tmp
#         # mu_tmp <- crossprod(regs_tmp, x_lhs)/sig_sq_x[d, m]
#         mu_tmp <- crossprod(regs_tmp, vcm_x_errors) %*% x_lhs
#         mu_tmp_all <- mu_tmp_all + mu_tmp
#       }
#       Omega_bet    <- solve(omega_tmp_all + prior_vcm_bet_z[[d]])
#       mu_bet       <- Omega_bet %*% mu_tmp_all
#       #
#       #
#       #
#       #
#       #
#       # beta_sampled   <- MASS::mvrnorm(n = 1, mu = mu_bet, Sigma = Omega_bet)
#       beta_sampled <- rnorm_fast_n1(mu = mu_bet, Sigma = Omega_bet, dim_bet_z[d] + 1)
#       while ((abs(abs(beta_sampled[1]) - 1) < 0.01) | abs(beta_sampled[1]) > 1) {
#         # beta_sampled <- MASS::mvrnorm(n = 1, mu = mu_bet, Sigma = Omega_bet)
#         beta_sampled <- rnorm_fast_n1(mu = mu_bet, Sigma = Omega_bet, dim_bet_z[d] + 1)
#       }
#       phi_x[d, m] <- beta_sampled[1]
#       bet_z[id_betz_tmp, m] <- beta_sampled[-1]
