#' helps to compute as weights
#'
#' helps to compute as weights
#'
#' @param M lhs
#' @param x rhs
#'
#' @return some return value needed within the computations of
#'  \code{cbpf_as_r()}
helper_as <- function(M, x) {
  apply(X = x,
        MARGIN = 1,
        function(x) {
          drop(crossprod(crossprod(M, x), x))})
}
#' Computes log particle weights for the Multinomial model
#'
#' Computes normalized particle weights for particles in the
#' \code{cbpf_as_m_r()} function.
#'
#' @param y measurements: multinomial counts; matrix of dimension \code{DDxTT}
#' @param NN number of particles
#' @param xa state trajectories: matrix of DDxTT
#' @param DD number of categories/classes
#'
#' @return \code{NN}-dimensional vector of normalized weights
#' @export
w_cbpf_m_r <- function(y, NN, xa, DD) {
  xs <- xa[, 1:(DD - 1)]
  ps <- exp(xs)
  ys <- matrix(rep(as.vector(y[-DD]), times = NN),
               ncol = DD - 1, nrow = NN, byrow = TRUE)

  out_tmp <- xs - log(.rowSums(x = ps, m = NN, n = (DD - 1), na.rm = TRUE) + 1)
  out_tmp <- out_tmp * ys

  w       <- .rowSums(out_tmp, m = NN, n = (DD - 1), na.rm = TRUE)

  check_w_computation(w)
  return(w)
}
#' Computes log particle weights for the Dirichlet-multinomial model
#'
#' Computes normalized particle weights for particles in the
#' \code{cbpf_as_dm_r()} function.
#'
#' @param y measurements: Dirichlet-Multinomial fractions/shares; matrix of
#'   dimension \code{DDxTT}
#' @param NN number of particles
#' @param xa state trajectories: matrix of \code{DDxTT}
#' @param num_counts measurement: Dirichlet-multinomial total counts per time
#'   period; vector of dimension \code{TT}
#' @param DD number of categories/classes
#'
#' @return \code{NN}-dimensional vector of normalized weights
#' @export
w_cbpf_dm_r <- function(y, NN, xa, num_counts, DD) {
  alphas <- matrix(exp(xa), nrow = NN, ncol = DD)
  # alphas[alphas == 0] <- 1e-300
  # log_Balpha <- rowSums(lgamma(alphas)) - lgamma(rowSums(alphas))
  # log_denom  <- (alphas - 1) %*% t(log(y))
  # w <- log_denom - log_Balpha
  # browser()
  ys <- matrix(rep(as.vector(y), times = NN),
               ncol = DD,
               nrow = NN,
               byrow = TRUE)
  log_lhs <- (lgamma(.rowSums(x = alphas, m = NN, n = DD)) -
                lgamma(.rowSums(x = alphas, m = NN, n = DD) + num_counts))
  log_rhs <- .rowSums(lgamma(alphas + ys) - lgamma(alphas),
                      m = NN, n = DD)
  w <- log_lhs + log_rhs

  check_w_computation(w)
  return(w)
  # list(.rowSums(x = alphas, m = NN, n = DD))
  # list(-lgamma(.rowSums(x = alphas, m = NN, n = DD) + num_counts))
  # lgamma(.rowSums(x = alphas, m = NN, n = DD)),
  # -lgamma(.rowSums(x = alphas, m = NN, n = DD) + num_counts),
  # lgamma(.rowSums(x = alphas, m = NN, n = DD)) -
  #   lgamma(.rowSums(x = alphas, m = NN, n = DD) + num_counts)
}
#' Computes log particle weights for the Dirichlet model
#'
#' Computes normalized particle weights for particles in the
#'   \code{cbpf_as_d_r()} function.
#'
#' @param y measurements: Dirichlet fractions/shares; matrix of
#'   dimension \code{DDxTT}
#' @param xa state trajectories: matrix of \code{DDxTT}
#' @param NN number of particles
#' @param DD number of categories/classes
#'
#' @return \code{NN}-dimensional vector of normalized weights
#' @export
w_cbpf_d_r <- function(y, xa, NN, DD) {
  alphas <- matrix(exp(xa), nrow = NN, ncol = DD)
  ys <- matrix(rep(as.vector(y), times = NN),
               ncol = DD,
               nrow = NN,
               byrow = TRUE)
  log_lhs <- (lgamma(.rowSums(x = alphas, m = NN, n = DD)) -
                .rowSums(x = lgamma(alphas), m = NN, n = DD))
  # EVALUATE HERE INVERSE OF THE MULTIVARIATE BETA FUNCTION
  # check if the matrix multiplication actually performs pointwise
  log_rhs <- .rowSums((alphas - 1) * log(ys),
                      m = NN, n = DD)
  # print(log_lhs)
  # print(log_rhs)
  # browser()
  w <- log_lhs + log_rhs

  check_w_computation(w)
  return(w)
}
w_normalize <- function(w) {
  w_max   <- max(w)
  w_tilde <- exp(w - w_max)
  w  <- w_tilde / sum(w_tilde)

  check_w_computation(w)
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
  x_t <- phi_x * x_tt + z_beta
  # xt <- phi_x*xtt + 8*cos(1.2*t)
  # xt <- phi_x*xtt + 25*xtt/(1 + xtt^2)
  # xt <- phi_x*xtt + 25*xtt/(1 + xtt^2) + 8*cos(1.2*t)
  return(x_t)
}
#' Checks for numeric anomalies in the weight computations.
#'
#' An error is thrown if any of the weights passed to the argument \code{weight}
#' are 'NA' or 'NaN'.
#'
#' @param weights a numeric vector of weights to check
#'
#' @return either an error is thrown or returning invisible NULL
#' @export
check_w_computation <- function(weights) {
  check_me <- FALSE
  check_me <- any(is.nan(weights)) || any(is.na(weights))
  if (check_me) {
    stop("NAN or NA values in weight computation!")
  }
  check_me <- any(is.infinite(weights))
  if (check_me) {
    stop("INFINITE values in weight computation!")
  }
  check_me <- any(weights == 0)
  if (check_me) {
    warning("ZERO values in weight computation!")
  } else {
    return(invisible(NULL))
  }
}
