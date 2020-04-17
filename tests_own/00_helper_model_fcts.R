#' State transition
#'
#' Helper function computing the deterministic state transition, or, to put
#' differently the one-period ahead conditional mean of the latent state
#' process.
#'
#' @param x_tt state value in t-1 i.e. x_{t-1}
#' @param z regressor values i.e. z_{t} (matrix)
#' @param phi_x autoregressive parameter phi
#' @param bet_x regressor parameters/coefficients at z_{t} (matrix)
#'
#' @return
#' @export
f <- function(x_tt, z, phi_x, bet_x) {
  # xt <- phi_x*xtt
  x_t <- phi_x*x_tt + z %*% bet_x
  # xt <- phi_x*xtt + 8*cos(1.2*t)
  # xt <- phi_x*xtt + 25*xtt/(1 + xtt^2)
  # xt <- phi_x*xtt + 25*xtt/(1 + xtt^2) + 8*cos(1.2*t)
  return(x_t)
}
