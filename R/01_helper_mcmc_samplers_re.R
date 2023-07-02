#' Method for only linear effects Gibbs sampler.
#'
#' @inheritParams sample_all_params
#'
#' @export
sample_all_params.re <- function(pe, mm) {
  for (d in 1:pe$DD) {
    id_betu_tmp <- (pe$id_uet[d] + 1):pe$id_uet[d + 1]
    id_uet_tmp  <- (pe$id_uet[d] + 1):pe$id_uet[d + 1]

    dd_range_nn <- pe$dd_list_nn[[d]]
    Xtmp <- as.matrix(pe$X[, d, mm - 1, ])
    Utmp <- pe$U[, id_uet_tmp, , drop = FALSE]

    pe$sig_sq_x[d, mm] <- sample_sig_sq_x_r(bet_u = pe$bet_u[id_betu_tmp,
                                                              mm - 1, ,
                                                              drop = FALSE],
                                             X = Xtmp,
                                             regs_u = Utmp,
                                             prior_ig = c(pe$prior_ig_a,
                                                          pe$prior_ig_b),
                                             iter_range_NN = dd_range_nn,
                                             TT = pe$TT)
    pe$vcm_bet_u[[d]][, , mm] <- sample_vcm_bet_u(pe$bet_u[id_betu_tmp,
                                                           mm - 1, ,
                                                           drop = FALSE],
                                                  pe$dim_bet_u[d],
                                                  pe$dof_vcm_bet_u[d],
                                                  pe$prior_vcm_bet_u2[[d]],
                                                  dd_range_nn)
    pe$bet_u[id_betu_tmp, mm,
             dd_range_nn] <- sample_bet_u_r(pe$sig_sq_x[d, mm],
                                            pe$vcm_bet_u[[d]][, , mm],
                                            pe$dim_bet_u[d],
                                            Xtmp, Utmp,
                                            dd_range_nn,
                                            pe$TT)
    pe$Regs_beta[, d, ] <- get_regs_beta_r(U = pe$U,
                                           id_uet = c(pe$id_uet[d],
                                                      pe$id_uet[d + 1]),
                                           TT = pe$TT,
                                           bet_u = pe$bet_u[, mm, ,
                                                            drop = FALSE],
                                           id_uet = c(pe$id_uet[d],
                                                        pe$id_uet[d + 1]),
                                           iter_range_NN = 1:pe$NN)
  }
  cat("MCMC iteration number:", mm, "\n")
}
#' Draws a (particle) Gibbs sample of the std. deviation parameter
#'
#' The standard deviation is per component \code{d} of the DD-dimensional latent
#' state process and drawn from the inverse Gamma. This version is adjusted to
#' incorporate linear random effect regressors only.
#'
#' @inheritParams sample_sig_sq_x_alr
#'
#' @return one sample from the inverse gamma distribution for the standard
#'   deviation parameter
#' @export
sample_sig_sq_x_r <- function(bet_u,
                              X,
                              regs_u,
                              prior_ig,
                              iter_range_NN,
                              TT) {
  err_sig_sq_x_all <- 0
  for (n in iter_range_NN) {
    tmp_regs <- matrix(regs_u[, , n], nrow = TT)
    tmp_bets <- matrix(bet_u[, , n], nrow = ncol(tmp_regs))
    err_sig_sq_x <- X[, n, drop = FALSE] - (tmp_regs %*%  tmp_bets)
    err_sig_sq_x_all <- err_sig_sq_x_all + sum(err_sig_sq_x ^ 2)
  }
  out <- 1 / stats::rgamma(n = 1,
                           prior_ig[1],
                           prior_ig[2] + err_sig_sq_x_all / 2)
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
sample_bet_u_r <- function(sig_sq_x,
                           vcm_bet_u,
                           dim_bet_u,
                           X,
                           U,
                           iter_range_NN,
                           TT) {
  out_mat          <- matrix(0, nrow = dim_bet_u, ncol = length(iter_range_NN))
  vcm_bet_u_inv    <- solveme(vcm_bet_u)

  nn <- 1
  for (n in iter_range_NN) {
    Omega_bet_u <- matrix(0, nrow = dim_bet_u, ncol = dim_bet_u)
    mu_bet_u    <- matrix(0, nrow = dim_bet_u, ncol = 1)

    x_n   <- X[, n, drop = FALSE]
    Umat <- matrix(U[, , n, drop = FALSE], nrow = TT)

    Omega_bet_u <- solveme(crossprod(Umat, Umat) / sig_sq_x + vcm_bet_u_inv)
    mu_bet_u    <- Omega_bet_u %*% (crossprod(Umat, x_n) / sig_sq_x)

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
