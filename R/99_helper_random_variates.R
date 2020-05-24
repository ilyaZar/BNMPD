#' My rnorm function for the case n = 1
#'
#' This case is needed for the mcmc sampling part; it's basically a simplified
#' wrapper around \code{MASS:mvrnorm()}
#'
#' @param mu mean
#' @param Sigma variance-covariance matrix
#' @param len length of the mean vector (or nrows/ncols of \code{Sigma})
#'
#' @return a numeric vector of length \code{len} giving the sampled values
#' @export
rnorm_fast_n1 <- function(mu, Sigma, len) {

  tol <- sqrt(.Machine$double.eps)
  if (!isSymmetric(Sigma, tol = tol, check.attributes = FALSE)) {
    warning("Sigma must be a symmetric matrix")
  }
  eS <- eigen(Sigma, symmetric = TRUE)
  evals <- eS$values
  if (any(evals < -tol * abs(evals[1L]))) {
    warning("'Sigma' is not positive definite")
  }
  # return(mu + eS$vectors %*% (diag(sqrt(pmax(evals, 0)), len) %*% rnorm(len)))
  return(mu + eS$vectors %*% matrix((sqrt(pmax(evals, 0))) * rnorm(len), ncol = 1))
  #
  #
  #
  #
  #
  ##############################################################################
  ##############################################################################
  ##############################################################################
  # n = 1, tol = 1e-06, method = c("eigen", "svd", "chol")
  # p <- length(mu)
  # nm <- names(mu)
  # if (is.null(nm) && !is.null(dn <- dimnames(Sigma)))
  #   nm <- dn[[1L]]
  # dimnames(X) <- list(nm, NULL)
  # empirical = FALSE
  # if (empirical) {
  #   X <- scale(X, TRUE, FALSE)
  #   X <- X %*% svd(X, nu = 0)$v
  #   X <- scale(X, FALSE, TRUE)
  # }
  # if (!isSymmetric(Sigma, tol = sqrt(.Machine$double.eps), check.attributes = FALSE)) {
  #   stop("Sigma must be a symmetric matrix")
  # }
  # if (length(mu) != nrow(Sigma))
  #   stop("mu and Sigma have non-conforming size")
  # if (!all(dim(Sigma) == c(p, p)))
  #   stop("incompatible arguments")
  # eS <- eigen(Sigma, symmetric = TRUE)
  # ev <- eS$values
  # if (!all(ev >= -tol * abs(ev[1L])))
  #   stop("'Sigma' is not positive definite")
  # X <- matrix(rnorm(p * n), n)
  # X <- drop(mu) + eS$vectors %*% diag(sqrt(pmax(ev, 0)), p) %*% t(X)
  #
  # if (n == 1)
  #   drop(X)
  # else t(X)
  #
  #
  #
  #
  ##############################################################################
  ##############################################################################
  ##############################################################################
  # method <- match.arg(method)
  # R <- if (method == "eigen") {
  #   ev <- eigen(Sigma, symmetric = TRUE)
  #   if (!all(ev$values >= -sqrt(.Machine$double.eps) * abs(ev$values[1]))) {
  #     warning("Sigma is numerically not positive semidefinite")
  #   }
  #   t(ev$vectors %*% (t(ev$vectors) * sqrt(pmax(ev$values, 0))))
  # }
  # else if (method == "svd") {
  #   s. <- svd(Sigma)
  #   if (!all(s.$d >= -sqrt(.Machine$double.eps) * abs(s.$d[1]))) {
  #     warning("Sigma is numerically not positive semidefinite")
  #   }
  #   t(s.$v %*% (t(s.$u) * sqrt(pmax(s.$d, 0))))
  # }
  # else if (method == "chol") {
  #   R <- chol(Sigma, pivot = TRUE)
  #   R[, order(attr(R, "pivot"))]
  # }
  # retval <- matrix(rnorm(n * ncol(Sigma)), nrow = n, byrow = !pre0.9_9994) %*%
  #   R
  # retval <- sweep(retval, 2, mean, "+")
  # colnames(retval) <- names(mean)
  # retval
}

