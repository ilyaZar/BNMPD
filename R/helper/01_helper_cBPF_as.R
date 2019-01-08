helper_as <- function(M, x) {
  apply(X = x,
        MARGIN = 1,
        function(x) {drop(crossprod(crossprod(M, x), x))})
}
w_BPF <- function(y, yz, KK, xa1, xa2, xa3, xa4) {
  N   <- length(xa1)
  zks <- matrix(rep(yz, times = N), nrow = N, byrow = TRUE)
  d <- zks/exp(xa2)
  d <- d^exp(xa1)
  d <- d/(1 + d)

  if (any(is.nan(d))) {d[is.nan(d)] <- 1}

  F_gb2 <- pbeta(q = d, shape1 = exp(xa3), shape2 = exp(xa4))
  F_gb2 <- cbind(F_gb2, rep(1, times = N))

  pi_prob <- F_gb2[, 2:(KK + 1)] - F_gb2[, 1:KK]

  if (any(pi_prob < 0)) {pi_prob[pi_prob < 0] <- 0}

  pi_prob <- log(pi_prob)

  pi_prob <- t(pi_prob)*y
  w <- .colSums(pi_prob, m = KK, n = N)

  if (sum(is.nan(pi_prob))) {
    stop("NA values in weight computation!")
  }
  return(w)
}
