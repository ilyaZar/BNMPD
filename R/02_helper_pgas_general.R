#' Multiplication of regressor matrix with coefficient vector
#'
#' Multiplication of regressor matrices Z and U with corresponding coefficient
#' vectors and adding up. The result needs to be passed to the SMC sampler.
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
get_regs_beta <- function(Z, U, id_uet, TT, NN,
                          bet_z, bet_u, id_bet_u,
                          iter_range_NN) {
  Z_beta    <- matrix(0, nrow = TT, ncol = NN)
  U_beta    <- matrix(0, nrow = TT, ncol = NN)
  Regs_beta <- matrix(0, nrow = TT, ncol = NN)
  # browser()
  for (n in iter_range_NN) {
    Z_tmp <- Z[, , n]
    U_tmp <- matrix(U[, (id_uet[1] + 1):id_uet[2], n, drop = FALSE], nrow = TT)
    bet_u_tmp <-matrix(bet_u[(id_bet_u[1] + 1):id_bet_u[2], , n, drop = FALSE])

    Z_beta[, n]    <- Z_tmp %*% bet_z
    U_beta[, n]    <- U_tmp %*% bet_u_tmp
    Regs_beta[, n] <- Z_beta[, n] + U_beta[, n]
  }
  return(Regs_beta)
}
