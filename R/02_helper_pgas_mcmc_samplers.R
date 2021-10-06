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
                            iter_range_NN) {
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
#' Draws a (particle) Gibbs sample of the covariance matrix of random effects
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
                         iter_range_NN) {

  out_mat         <- matrix(0, nrow = dim_bet_u, ncol = length(iter_range_NN))
  vcm_x_errors     <- diag(rep(sig_sq_x, times = TT - 1))
  vmc_x_errors_inv <- solve(vcm_x_errors)
  vcm_bet_u_inv    <- solve(vcm_bet_u)

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

    out_mat[, n] <- rnorm_fast_n1(mu = mu_bet_u, Sigma = Omega_bet_u, dim_bet_u)
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
#
#
#
#
#
# phi_x[d, m] <- beta_sampled[1]
# bet_z[(id_bet_z[d] + 1):id_bet_z[d + 1], m] <- beta_sampled[-1]
# vcm_bet_u[[d]][, , m] <-
#  sample_vcm_bet_u(bet_u = bet_u[(id_bet_u[d] + 1):id_bet_u[d + 1], m, ],
#                                dim_bet_u = dim_bet_u[d],
#                                dof_vcm_bet_u = dof_vcm_bet_u[d],
#                                prior_vcm_bet_scl = prior_vcm_bet_u2[[d]],
#                                 iter_range_NN = 1:NN)
# monitor_mcmc <- function(states_true, states_drawn) {
# }
