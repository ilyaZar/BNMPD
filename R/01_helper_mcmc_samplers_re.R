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
  scale_mat_vcm_bet_u   <- solveme(scale_mat_vcm_bet_u + prior_vcm_bet_scl)
  out <- solveme(stats::rWishart(1, dof_vcm_bet_u,
                                 scale_mat_vcm_bet_u)[, , 1])
  return(out)
}
