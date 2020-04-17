#' Generates random samples from dirichlet distribution
#'
#' Generates random samples from dirichlet distribution; the dimension of the
#' Dirichlet distribution (i.e. the number of shares of fractions) is taken as
#' the number of columns in the alpha-matrix; the number of samples are taken
#' as the rows of the alpha matrix; hence, each row of \code{alpha} corresponds
#' to \code{D} alpha parameters for which a \code{D}-dimensional random draw
#' is generated, and all these \code{n} draws are returned in matrix form.
#'
#' @param alpha alpha parameters of dirichlet distribution given as a matrix
#'
#' @return a nxD dimensional matrix of Dirichlet draws
my_rdirichlet <- function(alpha) {

  n <- nrow(alpha)
  l <- ncol(alpha)
  n <- n*l
  x <- matrix(stats::rgamma(n = n, shape = t(alpha)), ncol = l, byrow = TRUE)
  x_colsums <- as.vector(x %*% rep(1, l))
  x/x_colsums
}
#' Generates random samples from a mulinomial distribution
#'
#' Generates random samples from multinomial distribution; the dimension of the
#' multinomial distribution (i.e. the number of counts that fall into the
#' categories) is taken as the number of columns in the \code{probs}-matrix; the
#' number of samples are taken as the rows of the \code{probs} matrix; hence,
#' each row of \code{robs} corresponds to the \code{D} parameters of the
#' mulitinomial distribution for which a \code{D}-dimensional random draw is
#' generated, and all these \code{n} draws are returned in matrix form.
#'
#' @param probs probabilities for the different categories of the multinomial
#'   distribution given as a matrix
#' @param num_counts counts for the different categories of the multinomial
#'
#' @return a nxD dimensional matrix of multinomial draws
my_rmultinomial <- function(probs, num_counts) {
  n <- nrow(probs)
  l <- ncol(probs)
  x <- matrix(0, ncol = l, nrow = n)
  for (t in 1:n) {
    x[t, ] <- stats::rmultinom(n = 1, size = num_counts[t], prob = probs[t, ])
  }
  x
}
#' Generates random samples from a Dirichlet-multinomial distribution
#'
#' Generates random samples from Dirichlet-dirichlet distribution; the dimension
#' of the Dirichlet-multinomial distribution (i.e. the number of shares of
#' fractions) is taken as the number of columns in the \code{alpha} matrix; the
#' number of samples are taken as the rows of the \code{alpha} matrix; hence,
#' each row of \code{alpha} corresponds to \code{D} \code{alpha} parameters for
#' which a \code{D}-dimensional random draw is generated, and all these \code{n}
#' draws are returned in matrix form.
#'
#' @param alpha alpha parameters of dirichlet distribution given as a matrix
#' @param num_counts number of counts for the dirichlet shares/fractions
#'
#' @return a nxD dimensional matrix of Dirichlet-multinomial draws
my_rmult_diri <- function(alpha, num_counts) {
  n <- nrow(alpha)
  l <- ncol(alpha)
  num_probs <- n*l
  probs <- matrix(stats::rgamma(n = num_probs, shape = t(alpha)), ncol = l, byrow = TRUE)
  probs_colsums <- as.vector(probs %*% rep(1, l))
  probs <- probs/probs_colsums

  x <- matrix(0, ncol = l, nrow = n)
  for (t in 1:n) {
    x[t, ] <- stats::rmultinom(n = 1, size = num_counts[t], prob = probs[t, ])
  }
  x
}
# Testing my dirichelet vs dirichelt from gtools-package
# a1 <- 1:5
# a2 <- 21:25
# a_all <- matrix(c(a1, a2), nrow = 2, byrow = T)
# set.seed(123)
# rdirichlet(1, a1)
# rdirichlet(1, a2)
# set.seed(123)
# my_rdirichlet(2, a_all)
# parameter_fct_log_norm <- function(exp_mu, exp_sd) {
#   log_mu  <- log(exp_mu/sqrt( 1 + (exp_sd^2/exp_mu^2) ))
#   log_var <- log(1 + exp_sd^2/exp_mu^2)
#   return(list(log_mu, log_var))
# }
# parameter_fct_log_norm_test <- function(log_mu, log_sd) {
#   exp_mu  <- exp(log_mu + log_sd^2/2)
#   exp_var <- (exp(log_sd^2) - 1)*(exp(2*log_mu + log_sd^2))
#   return(list(exp_mu, exp_var))
# }
