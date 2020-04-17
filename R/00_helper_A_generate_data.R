#' Generates data for various models.
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
#' @param TT number of time periods
#' @param DD number of shares/fractions (for dirichlet or dirichlet-multinomial)
#'   or the number of categories for a multinomial distribution
#' @param par_true list of true parameters that describe the latent state
#'   process
#' @param x_levels target "mean" levels of the states around which they fluctuate
#' @param x_log_scale logical; if \code{TRUE}, x_levels are taken as logs and
#'   the random number generation for the states is performed in logs
#' @param intercept_include logical; if \code{TRUE} include an intercept at the
#'   cross sectional unit
#'
#' @return a list of two: \code{[[1]]} -> regressors and \code{[[2]]} -> latent states
generate_data <- function(distribution = "dirichlet",
                          TT, DD,
                          par_true,
                          x_levels,
                          x_log_scale,
                          intercept_include) {
  xa1 <- rep(0, TT)
  xa2 <- rep(0, TT)
  xa3 <- rep(0, TT)
  xa4 <- rep(0, TT)
  xa5 <- rep(0, TT)
  xa6 <- rep(0, TT)

  sig_sq_xa1 <- par_true[[1]][[1]]
  phi_xa1    <- par_true[[1]][[2]]
  bet_xa1    <- par_true[[1]][[3]]
  sig_sq_xa2 <- par_true[[2]][[1]]
  phi_xa2    <- par_true[[2]][[2]]
  bet_xa2    <- par_true[[2]][[3]]
  sig_sq_xa3 <- par_true[[3]][[1]]
  phi_xa3    <- par_true[[3]][[2]]
  bet_xa3    <- par_true[[3]][[3]]
  sig_sq_xa4 <- par_true[[4]][[1]]
  phi_xa4    <- par_true[[4]][[2]]
  bet_xa4    <- par_true[[4]][[3]]
  sig_sq_xa5 <- par_true[[5]][[1]]
  phi_xa5    <- par_true[[5]][[2]]
  bet_xa5    <- par_true[[5]][[3]]
  sig_sq_xa6 <- par_true[[6]][[1]]
  phi_xa6    <- par_true[[6]][[2]]
  bet_xa6    <- par_true[[6]][[3]]

  if (distribution %in% c("multinomial", "mult-diri", "mult-gen-diri")) {
    num_counts <- sample(x = 80000:120000, size = TT)
  }

  res_a1 <- generate_x_z(phi_x = phi_xa1, sig_sq_x = sig_sq_xa1, bet_x = bet_xa1,
                         x_level = x_levels[1],
                         x_sd = 0.0125,
                         process_log_scale = x_log_scale[1],
                         intercept   = intercept_include[1],
                         TT = TT)
  xa1    <- res_a1[[1]]
  za1    <- res_a1[[2]]
  res_a2 <- generate_x_z(phi_x = phi_xa2, sig_sq_x = sig_sq_xa2, bet_x = bet_xa2,
                         x_level = x_levels[2],
                         x_sd = 0.1,
                         process_log_scale = x_log_scale[2],
                         intercept   = intercept_include[2],
                         TT = TT)
  xa2    <- res_a2[[1]]
  za2    <- res_a2[[2]]
  res_a3 <- generate_x_z(phi_x = phi_xa3, sig_sq_x = sig_sq_xa3, bet_x = bet_xa3,
                         x_level = x_levels[3],
                         x_sd = 0.025,
                         process_log_scale = x_log_scale[3],
                         intercept   = intercept_include[3],
                         policy_dummy = FALSE,
                         TT = TT)
  xa3    <- res_a3[[1]]
  za3    <- res_a3[[2]]
  res_a4 <- generate_x_z(phi_x = phi_xa4, sig_sq_x = sig_sq_xa4, bet_x = bet_xa4,
                         x_level = x_levels[4],
                         x_sd = 0.1,
                         process_log_scale = x_log_scale[4],
                         intercept   = intercept_include[4],
                         TT = TT)
  xa4 <- res_a4[[1]]
  za4 <- res_a4[[2]]
  res_a5 <- generate_x_z(phi_x = phi_xa5, sig_sq_x = sig_sq_xa5, bet_x = bet_xa5,
                         x_level = x_levels[5],
                         x_sd = 0.1,
                         process_log_scale = x_log_scale[5],
                         intercept   = intercept_include[5],
                         TT = TT)
  xa5 <- res_a5[[1]]
  za5 <- res_a5[[2]]
  res_a6 <- generate_x_z(phi_x = phi_xa6, sig_sq_x = sig_sq_xa6, bet_x = bet_xa6,
                         x_level = x_levels[6],
                         x_sd = 0.1,
                         process_log_scale = x_log_scale[6],
                         intercept   = intercept_include[6],
                         TT = TT)
  xa6 <- res_a6[[1]]
  za6 <- res_a6[[2]]

  xalphas <- cbind(xa1, xa2, xa3, xa4, xa5, xa6)

  if (distribution == "dirichlet") {
    yraw <- my_rdirichlet(alpha = xalphas)
  }
  if (distribution == "multinomial") {
    xalphas <- xalphas/rowSums(xalphas)
    yraw <- my_rmultinomial(probs = xalphas, num_counts = num_counts)
  }
  if (distribution == "mult-diri") {
    yraw <- my_rmult_diri(alpha = xalphas, num_counts = num_counts)
  }

  if (distribution == "dirichlet") {
    if (sum(rowSums(yraw)) != TT) {
      stop("Something is wrong with the Dirichelet: y-fractions don't sum up to
           1!")
    }
    return(list(yraw,
                list(xa1, xa2, xa3, xa4, xa5, xa6),
                list(za1, za2, za3, za4, za5, za6)))
  }
  if (distribution == "multinomial") {
    return(list(yraw,
                list(xa1, xa2, xa3, xa4, xa5, xa6),
                list(za1, za2, za3, za4, za5, za6),
                num_counts = num_counts))
  }
  if (distribution == "mult-diri") {
    return(list(yraw,
                list(xa1, xa2, xa3, xa4, xa5, xa6),
                list(za1, za2, za3, za4, za5, za6),
                num_counts = num_counts))
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
#' @export
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
