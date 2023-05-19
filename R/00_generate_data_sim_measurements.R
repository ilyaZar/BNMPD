generate_measurements <- function(x_states, X_LOG_SCALE, distribution) {
  tmp_dim <- dim(x_states)
  TT <- tmp_dim[1]
  DD <- tmp_dim[2]
  NN <- tmp_dim[3]

  data_part1 <- generate_y_x_containter(distribution, NN = NN, TT = TT, DD = DD)

  if (distribution == "multinomial" || distribution == "dirichlet_mult") {
    data_part2 <- matrix(0, nrow = TT, ncol = NN)
    rownames(data_part2) <- paste0("t_", seq_len(TT))
    colnames(data_part2) <- paste0("n_", seq_len(NN))
  } else {
    data_part2 <- NULL
  }
  out_data <-  list(part1 = data_part1,
                    part2 = data_part2)
  attr(out_data, which = "distribution") <- distribution

  if (sum(any(x_states <= 0)) & X_LOG_SCALE == FALSE) {
    stop("some state process (x1_t, x2_t, ... or xD_t) not positive!")
  }
  if (X_LOG_SCALE) {
    x <- exp(x_states)
  }
  if (distribution == "normal") {
    out_data[["part1"]] <- x_states
    return(out_data) # early return with y=x for a Gaussian linear model spec.
  }
  for (n in 1:NN) {
    if (distribution == "dirichlet") {
      # yraw <- 0
      # while (any(yraw <= 0.01)) {
      yraw <- my_rdirichlet(alpha = x[, , n])
      if (sum(rowSums(yraw)) != TT || any(yraw == 0)) {
        msg <- paste0("Bad Dirichelet simulation: ",
                      "fractions don't sum to 1 and/or zero componenent!")
        stop(msg)
      }
      print(paste0("Simulatiing Dirichlet data at cross section: ",
                   n))
      out_data[["part1"]][, , n] <- yraw
      # }
    }
    if (distribution == "multinomial") {
      num_counts <- sample(x = 80000:120000, size = TT)
      tmp_x <- x[, , n]
      tmp_x / rowSums(tmp_x)
      yraw <- my_rmultinomial(probs = tmp_x, num_counts = num_counts)
      out_data[["part1"]][, , n] <- yraw
      out_data[["part2"]][, n]   <- num_counts
    }
    if (distribution == "dirichlet_mult") {
      num_counts <- sample(x = 80000:120000, size = TT)
      yraw <- my_rmult_diri(alpha =  x[, , n],
                            num_counts = num_counts)
      out_data[["part1"]][, , n] <- yraw
      out_data[["part2"]][, n]   <- num_counts
    }
  }
  return(out_data)
}
#' Generates random samples from Dirichlet distribution
#'
#' Generates random samples from Dirichlet distribution; the dimension of the
#' Dirichlet distribution (i.e. the number of shares or fractions) is taken as
#' the number of columns in the alpha-matrix; the number of samples are taken as
#' the number of rows of the alpha matrix; hence, each row of \code{alpha}
#' corresponds to \code{D} alpha parameters for which a \code{D}-dimensional
#' random draw from the Dirichlet is generated.
#'
#' @param alpha a matrix of alpha parameters of a Dirichlet distribution
#'
#' @return a \code{n x D} dimensional matrix of Dirichlet draws
my_rdirichlet <- function(alpha) {

  n <- nrow(alpha)
  l <- ncol(alpha)
  n <- n*l
  x <- matrix(stats::rgamma(n = n, shape = t(alpha)), ncol = l, byrow = TRUE)
  x_colsums <- as.vector(x %*% rep(1, l))
  x/x_colsums
}
#' Generates random samples from generalized Dirichlet distribution
#'
#' Generates random samples from generalized Dirichlet distribution; the
#' dimension of the Dirichlet distribution (i.e. the number of shares or
#' fractions) is taken as the number of columns in the alpha-matrix; the number
#' of samples are taken as the number of rows of the alpha matrix; hence, each
#' row of \code{alpha} corresponds to \code{D} alpha parameters for which a
#' \code{D}-dimensional random draw from the generalized Dirichlet is generated.
#'
#' @param alpha a matrix of alpha parameters of a generalized Dirichlet
#'   distribution
#' @param beta a matrix of beta parameters of a generalized Dirichlet
#'   distribution
#'
#' @return a nxD dimensional matrix of generalized Dirichlet draws
my_r_generalized_dirichlet <- function(alpha, beta) {
  browser()
  if(!(nrow(alpha) == nrow(beta)) || !(ncol(alpha) == ncol(beta))) {
    stop("Arguments 'alpha' and 'beta' must have the same number of rows/cols!")
  }
  n <- nrow(alpha)
  l <- ncol(alpha)
  n <- n*l
  x <- matrix(stats::rbeta(n = n,
                           shape1 = t(alpha),
                           shape2 = t(beta)), ncol = l, byrow = TRUE)
  x_cumsums <- t(apply(x, MARGIN = 1, cumsum))
  x[, 2:l] <- x[, 2:l]*(1-x_cumsums[1:l])
  stop("THIS FUNCTION STILL NEEDS TESTING!")
  return(x)
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
#' Generates random samples from Dirichlet-multinomial distribution; the
#' dimension of the Dirichlet-multinomial distribution (i.e. the number of
#' shares of fractions) is taken as the number of columns in the \code{alpha}
#' matrix; the number of samples are taken as the rows of the \code{alpha}
#' matrix; hence, each row of \code{alpha} corresponds to \code{D} \code{alpha}
#' parameters for which a \code{D}-dimensional random draw is generated, and all
#' these \code{n} draws are returned in matrix form.
#'
#' @param alpha alpha parameters of dirichlet distribution given as a matrix
#' @param num_counts number of counts for the dirichlet shares/fractions
#'
#' @return a nxD dimensional matrix of Dirichlet-multinomial draws
my_rmult_diri <- function(alpha, num_counts) {
  n <- nrow(alpha)
  l <- ncol(alpha)
  num_probs <- n*l
  probs <- matrix(stats::rgamma(n = num_probs, shape = t(alpha)),
                  ncol = l, byrow = TRUE)
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
