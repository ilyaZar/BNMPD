helper_as <- function(M, x) {
  apply(X = x,
        MARGIN = 1,
        function(x) {drop(crossprod(crossprod(M, x), x))})
}
w_BPF <- function(y, N, xa1, xa2, xa3, xa4, xa5, xa6, num_counts, D = 6) {
  alphas <- matrix(c(exp(xa1), exp(xa2), exp(xa3),
                     exp(xa4), exp(xa5), exp(xa6)),
                   nrow = N,
                   ncol = D)
  alphas[alphas == 0] <- 1e-300
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
  # if (sum(is.nan(w) | is.na(w))) {
  #   stop("NAN or NA values in weight computation!")
  # }
  # w
  # list(.rowSums(x = alphas, m = N, n = D))
  # list(-lgamma(.rowSums(x = alphas, m = N, n = D) + num_counts))
  # lgamma(.rowSums(x = alphas, m = N, n = D)),
  # -lgamma(.rowSums(x = alphas, m = N, n = D) + num_counts),
  # lgamma(.rowSums(x = alphas, m = N, n = D)) - lgamma(.rowSums(x = alphas, m = N, n = D) + num_counts)
}
w_cBPF <- function(y, N, xa, num_counts, D = 6) {
  alphas <- matrix(exp(xa), nrow = N, ncol = D)
  alphas[alphas == 0] <- 1e-300
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
  w_max   <- max(w)
  w_tilde <- exp(w - w_max)
  w  <- w_tilde/sum(w_tilde)
  # if (sum(is.nan(w) | is.na(w))) {
  #   stop("NAN or NA values in weight computation!")
  # }
  # w
  # list(.rowSums(x = alphas, m = N, n = D))
  # list(-lgamma(.rowSums(x = alphas, m = N, n = D) + num_counts))
  # lgamma(.rowSums(x = alphas, m = N, n = D)),
  # -lgamma(.rowSums(x = alphas, m = N, n = D) + num_counts),
  # lgamma(.rowSums(x = alphas, m = N, n = D)) - lgamma(.rowSums(x = alphas, m = N, n = D) + num_counts)
}
f_cbpf <- function(x_tt, phi_x, z_beta) {
  # xt <- phi_x*xtt
  x_t <- phi_x*x_tt + z_beta
  # xt <- phi_x*xtt + 8*cos(1.2*t)
  # xt <- phi_x*xtt + 25*xtt/(1 + xtt^2)
  # xt <- phi_x*xtt + 25*xtt/(1 + xtt^2) + 8*cos(1.2*t)
  return(x_t)
}
