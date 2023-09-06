#' Generating observations/measurements
#'
#' @param x_states state values as generated via [generate_data_t()] for each
#'    \code{n=1,...,NN}
#' @inheritParams new_dataSim
#' @param dims dimensions as returned from [get_dimension()] with \code{dim}
#'    argument set to 'all' applied to an object of \code{class} 'trueParams'
#'
#' @return an array of size \code{TT x DD x NN} with measurements from the
#'    distribution passed via the \code{distribution} argument
#' @export
generate_measurements <- function(x_states, X_LOG_SCALE, distribution, dims, options_include) {
  check_distribution(distribution)
  check_x_states(x_states, X_LOG_SCALE)
  ### EXTREMELY AKWARD BUT MUST BE DONE SINCE OTHERWISE RANDOM NUMBERS ARE NOT
  ### THE SAME: INTEGER VS DOUBLE AND NAMING OF DIM CHANGE THE EXACT RESULTS!
  dimns <- names(dims)
  dims  <- as.integer(dims)
  names(dims) <- dimns
  TT  <- dims["TT"]
  NN  <- dims["NN"]
  DD  <- dims["DD"]

  out_data <- generate_data_container(distribution, NN = NN, TT = TT, DD = DD)
  if (X_LOG_SCALE) {
    for (n in 1:NN) {
      zero_cols <- options_include[[n]][["zeros"]]
      if (is.logical(zero_cols)) {
        x_states[, zero_cols, n] <- -Inf
      } else {
        next
      }
    }
    x <- exp(x_states)
  } else {
    x <- x_states
  }

  out_data <- switch(
    distribution,
    "normal" = generate_normal_obs(x, out_data),
    "dirichlet" = generate_dirichlet_obs(x, NN, TT, out_data),
    "gen_dirichlet" = generate_gen_dirichlet_obs(x, NN, TT, DD, out_data),
    "dirichlet_mult" = generate_dirichlet_mult_obs(x, NN, TT, out_data, options_include),
    "gen_dirichlet_mult" = generate_gen_dirichlet_mult_obs(x, NN, TT, DD, out_data)
  )
  return(out_data)
}
generate_normal_obs <- function(x, out_data) {
  if (any(x == 0)) {browser(); out_data[["part1"]]; return(out_data);}
  out_data[["part1"]] <- x
  return(out_data) # early return with y=x for a Gaussian linear model spec.
}
generate_multinomial_obs <- function(x, NN, TT, DD, out_data) {
  for (n in 1:NN) {
    if (any(x[, , n] == 0)) {
      browser()
      out_data[["part1"]][, , n] <- 0
      out_data[["part2"]][, n]   <- 0
      return(out_data)
    }
    num_counts <- sample(x = 80000:120000, size = TT)
    tmp_x <- x[, , n]
    tmp_x / rowSums(tmp_x)
    yraw <- my_rmultinomial(probs = tmp_x, num_counts = num_counts)
    print(paste0("Simulatiing Multinomial data at cross section: ", n))
    out_data[["part1"]][, , n] <- yraw
    out_data[["part2"]][, n]   <- num_counts
  }
  return(out_data)
}
generate_dirichlet_obs <- function(x, NN, TT, out_data) {
  for (n in 1:NN) {
    if (any(x[, , n] == 0)) {
      browser()
      out_data[["part1"]][, , n] <- 0
      return(out_data)
    }
    yraw <- my_rdirichlet(alpha = x[, , n])
    if (sum(rowSums(yraw)) != TT || any(yraw == 0)) {
      msg <- paste0("Bad Dirichelet simulation: ",
                    "fractions don't sum to 1 and/or zero componenent!")
      stop(msg)
    }
    print(paste0("Simulatiing Dirichlet data at cross section: ", n))
    out_data[["part1"]][, , n] <- yraw
  }
  return(out_data)
}
generate_gen_dirichlet_obs <- function(x, NN, TT, DD, out_data) {
  for (n in 1:NN) {
    xa <- x[, grepl("A", colnames(x)), , drop = FALSE]
    xb <- x[, grepl("B", colnames(x)), , drop = FALSE]
    if (any(xa[, , n] == 0) && all(xb[, , n] == 0)) {
      browser()
      out_data[["part1"]][, , n] <- 0
      return(out_data)
    }
    yraw <- my_r_generalized_dirichlet(alpha = xa[, , n, drop = FALSE],
                                       beta = xb[, , n, drop = FALSE],
                                       DD)
    if (sum(rowSums(yraw)) != TT || any(yraw == 0)) {
      msg <- paste0("Bad Dirichelet simulation: ",
                    "fractions don't sum to 1 and/or zero componenent!")
      stop(msg)
    }
    print(paste0("Simulatiing Gen. Dirichlet data at cross section: ", n))
    out_data[["part1"]][, , n] <- yraw
  }
  return(out_data)
}
generate_dirichlet_mult_obs <- function(x, NN, TT, out_data, options_include) {
  for (n in 1:NN) {
    if (is.logical(options_include[[n]]$zeros)) {
      no_zero_cols <- !options_include[[n]]$zeros
    } else if (is.null(options_include[[n]]$zeros)) {
      no_zero_cols <- seq_len(ncol(x))
    }
    num_counts <- sample(x = 80000:120000, size = TT)
    yraw <- my_rmult_diri(alpha =  x[, no_zero_cols, n],
                          num_counts = num_counts)
    out_data[["part1"]][, no_zero_cols, n] <- yraw
    out_data[["part2"]][, n]   <- num_counts
  }
  return(out_data)
}
generate_gen_dirichlet_mult_obs <- function(x, NN, TT, DD, out_data) {
  for (n in 1:NN) {
    xa <- x[, grepl("A", colnames(x)), , drop = FALSE]
    xb <- x[, grepl("B", colnames(x)), , drop = FALSE]
    if (any(xa[, , n] == 0) && all(xb[, , n] == 0)) {
      browser()
      out_data[["part1"]][, , n] <- 0
      out_data[["part2"]][, n]   <- 0
      return(out_data)
    }
    num_counts <- sample(x = 80000:120000, size = TT)

    yraw <- my_r_generalized_dirichlet_mult(alpha = xa[, , n, drop = FALSE],
                                            beta = xb[, , n, drop = FALSE],
                                            DD,
                                            num_counts)
    print(paste0("Simulatiing Gen. Dirichlet Mult. data at cross section: ", n))
    out_data[["part1"]][, , n] <- yraw
    out_data[["part2"]][, n]   <- num_counts
  }
  return(out_data)
}
generate_data_container <- function(distribution, NN = NN, TT = TT, DD = DD) {
  data_part1 <- generate_y_x_containter(distribution,
                                        NN = NN, TT = TT, DD = DD,
                                        type = "y")
  data_part2 <- generate_data_part2_containter(distribution, TT, NN)
  out_data <-  list(part1 = data_part1, part2 = data_part2)
  attr(out_data, which = "distribution") <- distribution
  return(out_data)
}
generate_data_part2_containter <- function(distribution, TT, NN) {
  if (distribution %in% c("multinomial",
                          "dirichlet_mult",
                          "gen_dirichlet_mult")) {
    data_part2 <- matrix(0, nrow = TT, ncol = NN)
    rownames(data_part2) <- paste0("t_", seq_len(TT))
    colnames(data_part2) <- paste0("n_", seq_len(NN))
  } else {
    data_part2 <- NULL
  }
  return(data_part2)
}
check_x_states <- function(x_states, X_LOG_SCALE) {
  if (sum(any(x_states <= 0)) & X_LOG_SCALE == FALSE) {
    stop("some state process (x1_t, x2_t, ... or xD_t) not positive!")
  }
  if (isTRUE(X_LOG_SCALE)) {
    if (sum(any(exp(x_states) == Inf))) {
      stop("some state process (x1_t, x2_t, ... or xD_t) are infinite!")
    }
  } else {
    if (sum(any(abs(x_states) == Inf))) {
      stop("some state process (x1_t, x2_t, ... or xD_t) are infinite!")
    }
  }
  return(invisible(x_states))
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
#' @param DD multivariate dimension of the space to sample
#'
#' @return a nxD dimensional matrix of generalized Dirichlet draws
my_r_generalized_dirichlet <- function(alpha, beta, DD) {
  if(!(nrow(alpha) == nrow(beta)) || !(ncol(alpha) == ncol(beta))) {
    stop("Arguments 'alpha' and 'beta' must have the same number of rows/cols!")
  }
  TT <- nrow(alpha)
  x  <- matrix(0, nrow = TT, ncol = DD)
  for (t in 1:TT) {
    tmp_sum <- 0
    for (k in seq_len(DD - 1)) {
      x[t, k] <- stats::rbeta(n = 1,
                              shape1 = alpha[t, k, 1],
                              shape2 = beta[t, k, 1]) * (1 - tmp_sum)
      tmp_sum <- tmp_sum + x[t, k]
    }

  }
  x[, DD] <- 1 - rowSums(x[, 1:(DD - 1), drop = FALSE])
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
my_r_generalized_dirichlet_mult <- function(alpha, beta, DD, num_counts) {
  if(!(nrow(alpha) == nrow(beta)) || !(ncol(alpha) == ncol(beta))) {
    stop("Arguments 'alpha' and 'beta' must have the same number of rows/cols!")
  }
  TT <- nrow(alpha)
  x  <- matrix(0, nrow = TT, ncol = DD)
  for (t in 1:TT) {
    x[t, ] <- MGLM::rgdirmn(n = 1, size = num_counts[t],
                            alpha[t, , 1], beta[t, , 1])
  }
  stopifnot(`Total number must equal rowsums in gen. dir. mult simulation` =
              rowSums(x) == num_counts)
  return(x)
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
