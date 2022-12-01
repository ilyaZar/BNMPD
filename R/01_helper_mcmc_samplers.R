#' Multiplication of regressor matrix with coefficient vector
#'
#' Multiplication of regressor matrices Z and U with corresponding coefficient
#' vectors and adding up. The result needs to be passed to the SMC sampler. This
#' is the version for linear and random effect type regressors.
#'
#' @param Z Z regressor matrix sliced at corresponding component \code{d}
#' @param U U regressor matrix sliced at corresponding component \code{d}
#' @param id_uet id of \code{d}-component random effects regressor matrix
#' @param bet_z m'th sample of beta coefficient of the state component \code{d}
#' @param bet_u m'th sample of beta coefficient of the state component \code{d}
#' @param id_bet_u id of \code{d}-component random effects
#' @param iter_range_NN iteration range i.e. the cross sectional components
#'   that are actually contributing to \code{d}
#'
#' @return a vector of dimension \code{TT x NN}, containing the corresponding
#'   multiplication result for the d'th component
#' @export
get_regs_beta <- function(Z, U, id_uet, TT,
                          bet_z, bet_u, id_bet_u,
                          iter_range_NN) {
  # browser()
  NN <- length(iter_range_NN)
  Z_beta    <- matrix(0, nrow = TT, ncol = NN)
  U_beta    <- matrix(0, nrow = TT, ncol = NN)
  Regs_beta <- matrix(0, nrow = TT, ncol = NN)
  # browser()
  for (n in iter_range_NN) {
    Z_tmp <- Z[, , n]
    U_tmp <- matrix(U[, (id_uet[1] + 1):id_uet[2], n, drop = FALSE], nrow = TT)
    bet_u_tmp <- matrix(bet_u[(id_bet_u[1] + 1):id_bet_u[2], , n, drop = FALSE])

    Z_beta[, n]    <- Z_tmp %*% bet_z
    U_beta[, n]    <- U_tmp %*% bet_u_tmp
    Regs_beta[, n] <- Z_beta[, n] + U_beta[, n]
  }
  return(Regs_beta)
}
#' Multiplication of regressor matrix with coefficient vector
#'
#' Multiplication of regressor matrices Z with corresponding coefficient
#' vectors and adding up (this is the linear-regressor only type version). The
#' result needs to be passed to the SMC sampler.
#'
#' @inheritParams get_regs_beta
#'
#' @return a vector of dimension \code{TT x NN}, containing the corresponding
#'   multiplication result for the d'th component
#' @export
get_regs_beta_l <- function(Z, TT, bet_z, iter_range_NN) {
  NN <- length(iter_range_NN)
  Z_beta    <- matrix(0, nrow = TT, ncol = NN)
  Regs_beta <- matrix(0, nrow = TT, ncol = NN)
  for (n in iter_range_NN) {
    Z_tmp <- Z[, , n]

    Regs_beta[, n] <- Z_tmp %*% bet_z
  }
  return(Regs_beta)
}
#' Multiplication of regressor matrix with coefficient vector
#'
#' Multiplication of regressor matrices U with corresponding coefficient
#' vectors and adding up (this is the random-effects-regressor only type
#' version). The result needs to be passed to the SMC sampler.
#'
#' @inheritParams get_regs_beta
#'
#' @return a vector of dimension \code{TT x NN}, containing the corresponding
#'   multiplication result for the d'th component
#' @export
get_regs_beta_r <- function(U, id_uet, TT,
                            bet_u, id_bet_u,
                            iter_range_NN) {
  # browser()
  NN <- length(iter_range_NN)
  # Z_beta    <- matrix(0, nrow = TT, ncol = NN)
  U_beta    <- matrix(0, nrow = TT, ncol = NN)
  Regs_beta <- matrix(0, nrow = TT, ncol = NN)
  # browser()
  for (n in iter_range_NN) {
    # Z_tmp <- Z[, , n]
    U_tmp <- matrix(U[, (id_uet[1] + 1):id_uet[2], n, drop = FALSE], nrow = TT)
    bet_u_tmp <- matrix(bet_u[(id_bet_u[1] + 1):id_bet_u[2], , n, drop = FALSE])

    # Z_beta[, n]    <- Z_tmp %*% bet_z
    # U_beta[, n]    <- U_tmp %*% bet_u_tmp
    # Regs_beta[, n] <- Z_beta[, n] + U_beta[, n]
    Regs_beta[, n] <- U_tmp %*% bet_u_tmp
  }
  return(Regs_beta)
}
get_x_rhs <- function(X, order_p, TT) {
  # if (order_p == 1) {
  # x_rhs <- X[1:(TT - order_p), , drop = FALSE]
  # tmp_dim <- unname(dim(x_rhs))
  # dim(x_rhs) <- c(TT = tmp_dim[1], NN = tmp_dim[2])
  # } else if (order_p > 1) {
  x_rhs <- array(0, dim = c(TT - order_p, order_p, dim(X)[2]))
  for (p in 1:order_p) {
    x_rhs[, p, ] <- X[(order_p - p + 1):(TT - p), , drop = FALSE]
  }
  tmp_dim <- unname(dim(x_rhs))
  dim(x_rhs) <- c(TT = tmp_dim[1], order_p = tmp_dim[2], NN = tmp_dim[3])
  # }
  return(x_rhs)
}
check_stationarity <- function(vals, order) {
  if (order == 1) {
    check <- all((1 - abs(vals[1])) < 0.01 || abs(vals[1]) > 1)
  } else {
    check <- all((1 - sum(abs(vals[1:order]))) < 0.01 || sum(abs(vals[1:order])) > 1)
  }
  return(check)
}
