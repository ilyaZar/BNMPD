helper_as <- function(M, x) {
  apply(X = x,
        MARGIN = 1,
        function(x) {drop(crossprod(crossprod(M, x), x))})
}
w_BPF <- function(y, N, xa1, xa2, xa3, xa4, xa5) {
  alphas <- matrix(c(exp(xa1), exp(xa2), exp(xa3), exp(xa4), exp(xa5)), nrow = N)
  log_Balpha <- rowSums(lgamma(alphas)) - lgamma(rowSums(alphas))
  log_denom  <- (alphas - 1) %*% t(log(y))
  w <- log_denom - log_Balpha

  if (sum(is.nan(w) | is.na(w))) {
    stop("NA values in weight computation!")
  }
  return(w)
}
