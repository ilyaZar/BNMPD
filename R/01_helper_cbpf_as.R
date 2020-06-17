#' helps to compute as weights
#'
#' helps to compute as weights
#'
#' @param M lhs
#' @param x rhs
#'
#' @return some return value needed within the computations of \code{cbpf_as_R()}
helper_as <- function(M, x) {
  apply(X = x,
        MARGIN = 1,
        function(x) {drop(crossprod(crossprod(M, x), x))})
}
#' Computes log=particle weights for the Multinomial model
#'
#' Computes normalized particle weights for particles in the
#' \code{cbpf_as_m_R()} function.
#'
#' @param y measurements: multinomial counts; matrix of dimension \code{DDxTT}
#' @param N number of particles
#' @param xa state trajectories: matrix of DDxTT
#' @param num_counts measurement: dirichlet-multinomial total counts per time
#'   period; vector of dimension \code{TT}
#' @param D number of categories/classes
#'
#' @return \code{N}-dimensional vector of normalized weights
w_cbpf_m_R <- function(y, N, xa, num_counts, D = DD) {
  log_num_counts_fac <- lfactorial(num_counts)

  alphas <- matrix(cbind(exp(xa[, 1:(D - 1)]), rep(1, times = N)), nrow = N, ncol = D)

  ys <- matrix(rep(as.vector(y), times = N), ncol = D, nrow = N, byrow = TRUE)

  log_lhs <- log_num_counts_fac - sum(lfactorial(y))
  log_rhs <- .rowSums(log(alphas/.rowSums(alphas, m = N, n = D)) * ys, m = N, n = D)
  w <- log_lhs + log_rhs

  if (any(is.nan(w)) || any(is.na(w))) {
    browser()
    stop("NAN or NA values in weight computation!")
  }
  return(w)
}
#' Computes log particle weights for the Dirichlet Multinomial model
#'
#' Computes normalized particle weights for particles in the
#' \code{cbpf_as_dm_R()} function.
#'
#' @param y measurements: dirichlet-multinomial fractions/shares; matrix of
#'   dimension \code{DDxTT}
#' @param N number of particles
#' @param xa state trajectories: matrix of DDxTT
#' @param num_counts measurement: dirichlet-multinomial total counts per time
#'   period; vector of dimension \code{TT}
#' @param D number of categories/classes
#'
#' @return \code{N}-dimensional vector of normalized weights
w_cbpf_dm_R <- function(y, N, xa, num_counts, D = DD) {
  alphas <- matrix(exp(xa), nrow = N, ncol = D)
  # alphas[alphas == 0] <- 1e-300
  # log_Balpha <- rowSums(lgamma(alphas)) - lgamma(rowSums(alphas))
  # log_denom  <- (alphas - 1) %*% t(log(y))
  # w <- log_denom - log_Balpha
  # browser()
  ys <- matrix(rep(as.vector(y), times = N), ncol = D, nrow = N, byrow = TRUE)
  log_lhs <- (lgamma(.rowSums(x = alphas, m = N, n = D)) -
                lgamma(.rowSums(x = alphas, m = N, n = D) + num_counts))
  log_rhs <- .rowSums(lgamma(alphas + ys) - lgamma(alphas),
                      m = N, n = D)
  w <- log_lhs + log_rhs

  if (any(is.nan(w)) || any(is.na(w))) {
    # browser()
    # stop("NAN or NA values in weight computation!")
  }
  return(w)
  # list(.rowSums(x = alphas, m = N, n = D))
  # list(-lgamma(.rowSums(x = alphas, m = N, n = D) + num_counts))
  # lgamma(.rowSums(x = alphas, m = N, n = D)),
  # -lgamma(.rowSums(x = alphas, m = N, n = D) + num_counts),
  # lgamma(.rowSums(x = alphas, m = N, n = D)) - lgamma(.rowSums(x = alphas, m = N, n = D) + num_counts)
}
w_normalize <- function(w) {
  w_max   <- max(w)
  w_tilde <- exp(w - w_max)
  w  <- w_tilde/sum(w_tilde)
  if (any(is.nan(w)) || any(is.na(w)) || any(w < 0)) {
    browser()
    stop("NAN or NA values in weight computation!")
  }
  return(w)
}
#' State transition
#'
#' Helper function computing the deterministic state transition, or, to put
#' differently the one-period ahead conditional mean of the latent state
#' process.
#'
#' @param x_tt state value in t-1 i.e. x_{t-1}
#' @param phi_x autoregressive parameter phi
#' @param z_beta result of regressor values i.e. z_{t} (matrix) multiplied by
#'   parameters/coefficients (vector)
#'
#' @return \code{T}-dimensional vector of (deterministically computed) state
#'   transitions (conditional means)
f_cbpf <- function(x_tt, phi_x, z_beta) {
  # xt <- phi_x*xtt
  x_t <- phi_x*x_tt + z_beta
  # xt <- phi_x*xtt + 8*cos(1.2*t)
  # xt <- phi_x*xtt + 25*xtt/(1 + xtt^2)
  # xt <- phi_x*xtt + 25*xtt/(1 + xtt^2) + 8*cos(1.2*t)
  return(x_t)
}
