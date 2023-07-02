#' Method for linear and random effects Gibbs sampler.
#'
#' @inheritParams sample_all_params
#'
#' @export
sample_all_params.lin_re <- function(pe, mm) {
  for (d in 1:pe$DD) {
    id_betz_tmp <- (pe$id_zet[d] + 1):pe$id_zet[d + 1]
    id_betu_tmp <- (pe$id_uet[d] + 1):pe$id_uet[d + 1]
    id_zet_tmp  <- (pe$id_zet[d] + 1):pe$id_zet[d + 1]
    id_uet_tmp  <- (pe$id_uet[d] + 1):pe$id_uet[d + 1]

    dd_range_nn <- pe$dd_list_nn[[d]]
    Xtmp <- as.matrix(pe$X[, d, mm - 1, ])
    Ztmp <- pe$Z[1:pe$TT, id_zet_tmp, , drop = FALSE]
    Utmp <- pe$U[1:pe$TT, id_uet_tmp, , drop = FALSE]

    pe$sig_sq_x[d, mm] <- sample_sig_sq_x_lr(bet_z = pe$bet_z[id_betz_tmp,
                                                              mm - 1],
                                             bet_u = pe$bet_u[id_betu_tmp,
                                                              mm - 1,,
                                                              drop = FALSE],
                                             X = Xtmp,
                                             regs_z = Ztmp, regs_u = Utmp,
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
             dd_range_nn] <- sample_bet_u_lr(pe$sig_sq_x[d, mm],
                                             pe$bet_z[id_betz_tmp,
                                                      mm - 1],
                                             pe$vcm_bet_u[[d]][, , mm],
                                             pe$dim_bet_u[d],
                                             Xtmp, Ztmp, Utmp,
                                             dd_range_nn,
                                             pe$TT)
    beta_sampled <- sample_bet_z_lr(sig_sq_x = pe$sig_sq_x[d, mm],
                                    vcm_bet_u = pe$vcm_bet_u[[d]][, , mm],
                                    X = Xtmp,
                                    regs_z = pe$regs_z,
                                    U = Utmp,
                                    TT = pe$TT,
                                    id_reg_z = c(pe$id_reg_z[d],
                                                 pe$id_reg_z[d + 1]),
                                    dim_bet_z = pe$dim_bet_z[d],
                                    prior_vcm_bet_z=pe$prior_vcm_bet_z[[d]],
                                    iter_range_NN = dd_range_nn)
    # browser()
    pe$bet_z[id_betz_tmp, mm] <- beta_sampled

    pe$Regs_beta[, d, ] <- get_regs_beta(Z  = pe$Z[, id_zet_tmp, ],
                                         U = pe$U,
                                         id_uet = c(pe$id_uet[d],
                                                    pe$id_uet[d + 1]),
                                         TT = pe$TT,
                                         pe$bet_z[id_betz_tmp, mm],
                                         bet_u = pe$bet_u[, mm, ,
                                                          drop = FALSE],
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
#' @inheritParams sample_sig_sq_x_alr
#'
#' @return one sample from the inverse gamma distribution for the standard
#'   deviation parameter
#' @export
sample_sig_sq_x_lr <- function(bet_z,
                               bet_u,
                               X,
                               regs_z,
                               regs_u,
                               prior_ig,
                               iter_range_NN,
                               TT) {
  err_sig_sq_x_all <- 0
  for(n in iter_range_NN) {
    bet_u_n <- bet_u[ , , n]
    tmp_regs <- cbind(regs_z[, , n], regs_u[, , n])
    err_sig_sq_x <- X[, n, drop = FALSE] - regs_z[, , n] %*% bet_z -
      regs_u[, , n] %*% bet_u_n
    err_sig_sq_x_all <- err_sig_sq_x_all + sum(err_sig_sq_x^2)
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
#' @inheritParams sample_bet_z_alr
#' @param id_reg_z a vector of 2; Z regressor id for components \code{d} and
#'   \code{d + 1}
#'
#' @return a sample
#' @export
sample_bet_z_lr <- function(sig_sq_x,
                            vcm_bet_u,
                            X,
                            regs_z,
                            U,
                            TT,
                            id_reg_z,
                            dim_bet_z,
                            prior_vcm_bet_z,
                            iter_range_NN) {
  omega_tmp_all <- 0
  mu_tmp_all    <- 0
  for (n in iter_range_NN) {
    regs_tmp  <- regs_z[, (id_reg_z[1] + 1):id_reg_z[2], n]
    omega_tmp <- compute_vcm_WBinv(sig_sq_x = sig_sq_x,
                                   matLHS = regs_tmp,
                                   U = matrix(U[,,n, drop = FALSE], nrow = TT),
                                   C = vcm_bet_u,
                                   matRHS = regs_tmp)
    omega_tmp_all <- omega_tmp_all + omega_tmp
    mu_tmp <- compute_vcm_WBinv(sig_sq_x = sig_sq_x,
                                matLHS = regs_tmp,
                                U = matrix(U[,,n, drop = FALSE], nrow = TT),
                                C = vcm_bet_u,
                                matRHS = X[, n, drop = FALSE])
    mu_tmp_all <- mu_tmp_all + mu_tmp
  }
  Omega_bet <- solveme(omega_tmp_all + prior_vcm_bet_z)
  mu_bet    <- Omega_bet %*% mu_tmp_all
  out <- rnorm_fast_n1(mu = mu_bet, Sigma = Omega_bet, dim_bet_z)
  return(out)
}
#' Draws a (particle) Gibbs sample of the random effects
#'
#' The covariance matrix of random effects is per component \code{d} of the
#' DD-dimensional latent state process and drawn from the inverse Wishart.
#'
#' @inheritParams sample_bet_u_alr
#'
#' @return a list of \code{length(iter_range_NN)} containing one sample of the
#'   \code{beta_u} coefficients
#' @export
sample_bet_u_lr <- function(sig_sq_x,
                            bet_z,
                            vcm_bet_u,
                            dim_bet_u,
                            X,
                            regs_z,
                            U,
                            iter_range_NN,
                            TT) {
  out_mat          <- matrix(0, nrow = dim_bet_u, ncol = length(iter_range_NN))
  vcm_bet_u_inv    <- solveme(vcm_bet_u)

  nn <- 1
  for (n in iter_range_NN) {
    Omega_bet_u <- matrix(0, nrow = dim_bet_u, ncol = dim_bet_u)
    mu_bet_u    <- matrix(0, nrow = dim_bet_u, ncol = 1)

    x_n   <- X[,n] - regs_z[,,n] %*% bet_z
    Umat <- matrix(U[, , n, drop = FALSE], nrow = TT)

    Omega_bet_u <- solveme(crossprod(Umat, Umat)/sig_sq_x + vcm_bet_u_inv)
    mu_bet_u    <- Omega_bet_u %*% (crossprod(Umat, x_n) / sig_sq_x)

    out_mat[, nn] <- rnorm_fast_n1(mu = mu_bet_u,
                                   Sigma = Omega_bet_u,
                                   dim_bet_u)
    nn <- nn + 1
  }
  return(out_mat)
}
