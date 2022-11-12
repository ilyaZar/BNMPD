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
  scale_mat_vcm_bet_u <- matrix(0, nrow = dim_bet_u, ncol = dim_bet_u)
  for (n in iter_range_NN) {
    bet_u_n <- bet_u[, , n]
    scale_mat_vcm_bet_u <- scale_mat_vcm_bet_u + tcrossprod(bet_u_n)
  }
  # browser()
  scale_mat_vcm_bet_u   <- solveme(scale_mat_vcm_bet_u + prior_vcm_bet_scl)
  out <- solveme(stats::rWishart(1, dof_vcm_bet_u,
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
  # browser()
  out_mat          <- matrix(0, nrow = dim_bet_u, ncol = length(iter_range_NN))
  vmc_x_errors_inv <- diag(rep(sig_sq_x^(-1), times = TT))
  vcm_bet_u_inv    <- solveme(vcm_bet_u)

  nn <- 1
  for (n in iter_range_NN) {
    Omega_bet_u <- matrix(0, nrow = dim_bet_u, ncol = dim_bet_u)
    mu_bet_u    <- matrix(0, nrow = dim_bet_u, ncol = 1)

    x_n   <- X[, n] - f(x_tt = rep(0, times = TT),
                        regs  = regs_z[, , n],
                        phi_x = 0,
                        bet_reg = bet_z)
    Umat <- matrix(U[, , n, drop = FALSE], nrow = TT)
    Omega_bet_u <- crossprod(Umat, vmc_x_errors_inv) %*% Umat + vcm_bet_u_inv
    # Omega_bet_u2 <- crossprod(Umat, Umat)/sig_sq_x + vcm_bet_u_inv
    Omega_bet_u <- solveme(Omega_bet_u)
    mu_bet_u   <- Omega_bet_u %*% (crossprod(Umat, x_n) / sig_sq_x)
    # mu_bet_u2    <- Omega_bet_u %*% (crossprod(Umat, vmc_x_errors_inv) %*% x_n)

    out_mat[, nn] <- rnorm_fast_n1(mu = mu_bet_u,
                                   Sigma = Omega_bet_u,
                                   dim_bet_u)
    nn <- nn + 1
  }
  # browser()
  return(out_mat)
}
