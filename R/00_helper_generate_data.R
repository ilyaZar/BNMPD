#' Generates panel data for various models.
#'
#' The data is a panel data from a Dirichlet, Multinomial or
#' Dirichlet-multinomial distribution. The mulivariate draws can vary along the
#' time and cross-sectional dimensions. This is because the parameter of the
#' distributions that generate the draws are modelled as a function of
#' regressors and latent states where the regressors and latent states can vary
#' over time and cross section.
#'
#' @param distribution specifies the distribution; "dirichlet", "multinomial" or
#'   "dirichlet-multinomial"
#' @param NN number of cross sectional units
#' @param TT number of time periods
#' @param DD number of shares/fractions (for dirichlet or dirichlet-multinomial)
#'   or the number of categories for a multinomial distribution
#' @param par_true list of true parameters that describe the latent state
#'   process
#' @param x_levels target "mean" levels of the states around which they
#'   fluctuate
#' @param x_log_scale logical; if \code{TRUE}, x_levels are taken as logs and
#'   the random number generation for the states is performed in logs
#' @param intercept_include logical; if \code{TRUE} include an intercept at the
#'   cross sectional unit
#'
#' @return NN-dimensinal list of lists of two: \code{[[1]]} -> regressors and
#'   \code{[[2]]} -> latent states
#' @export
generate_data_t_n <- function(distribution = "dirichlet",
                              NN, TT, DD,
                              par_true,
                              x_levels,
                              x_log_scale,
                              intercept_include) {
  out_data <- list()
  for (n in 1:NN) {
    par_true_current <- list(par_true[[1]][, n],
                             par_true[[2]][, n],
                             par_true[[3]][[n]])
    out_data[[n]] <- generate_data_t(distribution = distribution,
                                     TT = TT, DD = DD,
                                     par_true = par_true_current,
                                     x_levels = x_levels[, n],
                                     x_log_scale = x_log_scale,
                                     intercept_include = intercept_include[, n])
  }
  return(out_data)
}
#' Generates time series data for various models.
#'
#' The data is a time series, for a given cross sectional unit, from a
#' Dirichlet, Multinomial or Dirichlet-multinomial distribution. The mulivariate
#' draws can vary along the time dimension. This is because the parameter of the
#' distributions that generate the draws are modelled as a function of
#' regressors and latent states where the regressors and latent states can vary
#' over time.
#'
#' @param distribution specifies the distribution; "dirichlet", "multinomial" or
#'   "dirichlet-multinomial"
#' @param TT number of time periods
#' @param DD number of shares/fractions (for dirichlet or dirichlet-multinomial)
#'   or the number of categories for a multinomial distribution
#' @param par_true list of true parameters that describe the latent state
#'   process
#' @param x_levels vector of target "means"/"levels" of the states around which
#'   they fluctuate
#' @param x_log_scale logical; if \code{TRUE}, then \code{x_levels[i]} ist taken
#'   as log and the random number generation for the states is performed in logs
#' @param intercept_include logical vector; if component is \code{TRUE},
#'   includes an intercept at the cross sectional unit
#'
#' @return a list of two: \code{[[1]]} -> regressors and \code{[[2]]} -> latent
#'   states
generate_data_t <- function(distribution = "dirichlet",
                            TT, DD,
                            par_true,
                            x_levels,
                            x_log_scale,
                            intercept_include) {
  xa <- matrix(nrow = TT, ncol = DD, 0)
  za <- list()
  sig_sq_xa <- par_true[[1]]
  phi_xa <- par_true[[2]]
  bet_xa <- list()

  for (d in 1:DD) {
    bet_xa[[d]] <- par_true[[3]][[d]]
  }

  if (distribution %in% c("multinomial", "mult-diri", "mult-gen-diri")) {
    num_counts <- sample(x = 80000:120000, size = TT)
  }

  x_sd_level <- c(0.0125, 0.1, 0.025, 0.1, 0.1, 0.1)

  for (d in 1:DD) {
    res <- generate_x_z(phi_x = phi_xa[d], sig_sq_x = sig_sq_xa[d], bet_x = bet_xa[[d]],
                        x_level = x_levels[d],
                        x_sd = x_sd_level[d],
                        process_log_scale = x_log_scale,
                        intercept   = intercept_include[d],
                        TT = TT)
    xa[, d] <- res[[1]]
    za[[d]] <- res[[2]]
  }
  za <- Reduce(cbind, za)
  xalphas <- xa

  if (distribution == "dirichlet") {
    yraw <- my_rdirichlet(alpha = xalphas)
    if (sum(rowSums(yraw)) != TT) {
      stop("Something is wrong with Dirichelet: fractions don't sum up to 1!")
    }
    return(list(yraw, xa, za))
  }
  if (distribution == "multinomial") {
    xalphas <- xalphas/rowSums(xalphas)
    yraw <- my_rmultinomial(probs = xalphas, num_counts = num_counts)
    return(list(yraw, xa, za, num_counts = num_counts))
  }
  if (distribution == "mult-diri") {
    yraw <- my_rmult_diri(alpha = xalphas, num_counts = num_counts)
    return(list(yraw, xa, za, num_counts = num_counts))
  }
}
#' Simulates regressors and latent states.
#'
#' The core function that simulates the regressors and latent states given
#' certain parameters.
#'
#' @param phi_x autoregressive parameter of latent states
#' @param sig_sq_x standard deviation of error term in latent state transition
#'   equation
#' @param bet_x regressor coefficients/parameters of latent state process
#' @param x_level target levels of latent states around which they fluctuate (a
#'   tuning parameter that ensures that states at each multivariate component of
#'   the response fluctuate around a particular level)
#' @param x_sd standard deviation allowed for \code{x_level} target (a tuning
#'   parameter that ensures that deviations from the \code{x_level} target are
#'   not too severe)
#' @param process_log_scale logical; if \code{TRUE} process is simulated on the
#'   log-scale
#' @param intercept logical; if \code{TRUE}, includes an intercept term i.e. a
#'   constant level regressor
#' @param policy_dummy logical; if \code{TRUE} includes a dummy that jumps from
#'   zero to one (e.g. a policy or other jump effect to be modelled)
#' @param zero_pattern TO-BE-EXPLAINED-LATER
#' @param drift  TO-BE-EXPLAINED-LATER
#' @param TT number of time periods
#'
#' @return a \code{length(beta_x)}x\code{TT} matrix of regressors and a \code{TT}-dimensional vector of latent states
generate_x_z <- function(phi_x, sig_sq_x, bet_x,
                         x_level,
                         x_sd,
                         process_log_scale,
                         intercept,
                         policy_dummy = FALSE,
                         zero_pattern = 0,
                         drift = FALSE,
                         TT) {

  if (process_log_scale) {
    x_level <- log(x_level)
  }
  dim_reg <- length(bet_x)
  x <- rep(0, TT)
# BEGINNING OF REGRESSOR SIMULATION: --------------------------------------
  if (dim_reg == 1) {
    if (intercept) {
      const_level <- x_level * (1 - phi_x)/bet_x
      z <- matrix(const_level, nrow = TT, ncol = dim_reg)
    } else {
      const_mean <- x_level * (1 - phi_x)/bet_x
      z          <- matrix(stats::rnorm(T*dim_reg, mean = const_mean, sd = x_sd),
                           nrow = TT,
                           ncol = dim_reg,
                           byrow = TRUE)
    }
  } else {
    if (dim_reg == 2) {
      if (intercept) {
        zmeans     <- 1
      } else {
        zmeans <- stats::rnorm(dim_reg - 1, mean = 0, sd = 3)
      }
    } else if (dim_reg > 2) {
      if (intercept) {
        zmeans <- stats::rnorm(dim_reg - 2, mean = 0, sd = 3)
        zmeans <- c(1, zmeans)
      } else {
        zmeans <- stats::rnorm(dim_reg - 1, mean = 0, sd = 3)
      }
    }
    last_zmean <- x_level * (1 - phi_x) - sum(zmeans * bet_x[-dim_reg])
    last_zmean <- last_zmean/bet_x[dim_reg]
    zmeans     <- c(zmeans, last_zmean)
    z          <- matrix(stats::rnorm(TT*dim_reg, mean = zmeans, sd = x_sd),
                         nrow = TT,
                         ncol = dim_reg,
                         byrow = TRUE)
    if (intercept) {
      z[, 1] <- 1
    }
  }
# END OF REGRESSOR SIMULATION: --------------------------------------------
  xinit <- x_level # 0

  # browser()
  # set.seed(123)
  if (policy_dummy == TRUE) {
    z <- cbind(c(rep(0, times = round(0.5*TT, digits = 0)),
                 rep(1, times = T - round(0.5*TT, digits = 0))),
               z)
    bet_x <- c(-2, bet_x)
  }

  x[1] <- f(x_tt = xinit, z = z[1, ], phi_x = phi_x, bet_x = bet_x)
  x[1] <- x[1] + sqrt(sig_sq_x)*stats::rnorm(n = 1)

  for (t in 1:TT) {
    if (t < TT) {
      x[t + 1] <- f(x_tt = x[t], z = z[t + 1, ],
                    phi_x = phi_x, bet_x = bet_x)
      x[t + 1] <- x[t + 1] + sqrt(sig_sq_x)*stats::rnorm(n = 1)
    }
  }
  if (process_log_scale) {
    x <- exp(x)
  }
  # plot(x, type = "l")
# START TEST: -------------------------------------------------------------
  # set.seed(123)
  # z2 <- cbind(c(rep(0, times = round(0.3*T, digits = 0)),
  #               rep(1, times = T - round(0.3*T, digits = 0))),
  #             z)
  # bet2 <- c(-1, bet_x)
  #
  # x[1] <- f(x_tt = xinit, z = z2[1, ], phi_x = phi_x, bet_x = bet2)
  # x[1] <- x[1] + sqrt(sig_sq_x)*stats::rnorm(n = 1)
  #
  # for (t in 1:T) {
  #   if (t < T) {
  #     x[t + 1] <- f(x_tt = x[t], z = z2[t + 1, ],
  #                   phi_x = phi_x, bet_x = bet2)
  #     x[t + 1] <- x[t + 1] + sqrt(sig_sq_x)*stats::rnorm(n = 1)
  #   }
  # }
  # if (process_log_scale) {
  #   x <- exp(x)
  # }
  # plot(x, type = "l")
# END TEST ----------------------------------------------------------------

  if (sum(any(x <= 0)) & process_log_scale == FALSE) {
    stop("some state process (xa1_t, xa2_t, ... or xa6_t) not positive!")
  }
  return(list(x, z))
}
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
#' @return \code{T}-dimensional vector of (deterministically computed) state
#'   transitions (conditional means)
f <- function(x_tt, z, phi_x, bet_x) {
  # xt <- phi_x*xtt
  x_t <- phi_x*x_tt + z %*% bet_x
  # xt <- phi_x*xtt + 8*cos(1.2*t)
  # xt <- phi_x*xtt + 25*xtt/(1 + xtt^2)
  # xt <- phi_x*xtt + 25*xtt/(1 + xtt^2) + 8*cos(1.2*t)
  return(x_t)
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
# if (plot_measurements) {
#   names_title <- "Measurement components"
#   names_ylab  <- "measurements: y_t's"
#   names_xlab <- paste0("ya1_t (black),", " ya2_t (red),",
#                        " ya3_t (green),", " ya4_t (blue)",
#                        " ya5_t (turq.),", " and", " ya6_t (pink)")
#
#   par(mfrow = c(1,1))
#   all_measurms <- cbind(yraw[, 1], yraw[, 2],
#                         yraw[, 3], yraw[, 4],
#                         yraw[, 5], yraw[, 6])
#   matplot(all_measurms,
#           type = "l",
#           main = names_title,
#           ylab = names_ylab,
#           xlab = names_xlab
#   )
#   matplot(all_measurms/rowSums(all_measurms),
#           type = "l",
#           main = names_title,
#           ylab = names_ylab,
#           xlab = names_xlab
#   )
# }
# if (plot_states) {
#   names_title <- "True States"
#   names_ylab  <- "states: xt's"
#   names_xlab  <- paste0("xa1_t (black),", " xa2_t (red),",
#                         " xa3_t (green),",  "xa4_t (blue)",
#                         " xa5_t (turq.)", " and", " xa6_t (pink)")
#
#   par(mfrow = c(1,1))
#   all_states <- cbind(xa1, xa2, xa3, xa4, xa5, xa6)
#   matplot(all_states,
#           type = "l",
#           main = names_title,
#           ylab = names_ylab,
#           xlab = names_xlab
#   )
#   matplot(all_states/rowSums(all_states),
#           type = "l",
#           main = names_title,
#           ylab = names_ylab,
#           xlab = names_xlab
#   )
# }
