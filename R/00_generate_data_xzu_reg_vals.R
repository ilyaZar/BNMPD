#' Simulates regressors and latent states for fixed NN/DD but all TT.
#'
#' Core function to simulate regressors and latent states given certain
#' parameters along the time dimension (i.e. for a fixed cross sectional unit
#' and multivariate component).
#'
#' @param TT number of time periods
#' @param phi_x autoregressive parameter of latent states
#' @param sig_sq_x standard deviation of error term in latent state transition
#'   equation
#' @param bet_z regressor coefficients/parameters of latent state process
#'   referring to all cross sectional units
#' @param bet_u regressor coefficients/parameters of latent state process;
#'   random effects that refer to each individual cross sectional unit
#' @param bet_z_spl regressor coefficients/parameters of latent state process
#'   referring to all cross sectional units and are non-linear i.e. splines
#' @param bet_u_spl regressor coefficients/parameters of latent state process;
#'   random effects that refer to each individual cross sectional unit and are
#'   non-linear i.e. splines
#' @param options_reg_simul a list of 3 elements giving the options for
#'   regressor simulation:
#'   \itemize{
#'   \item{\code{x_level:}}{target levels i.e. stationary mean/level of latent
#'   states (ensures that states at each multivariate component of the response
#'   fluctuates around that particular level too)}
#'   \item{\code{reg_var_within:}}{variance of the regressor values within each
#'   component (tunable)}
#'   \item{\code{reg_var_among:}}{variance of random coefficient value draws
#'   (tunable)}
#'   }
#' @param modelling_reg_types named logical vector of dimension 4 where each
#'   component indicates, if \code{TRUE}, that the corresponding regressor type
#'   is to be generated in the following order:
#'   \itemize{
#'     \item{modelling_reg_types, component 1: }{includes z regressors}
#'     \item{modelling_reg_types, component 2: }{includes u regressors}
#'     \item{modelling_reg_types, component 3: }{includes z spline regressors}
#'     \item{modelling_reg_types, component 4: }{includes u spline regressors}
#'   }
#' @param intercept_z logical; if \code{TRUE} includes an intercept term i.e. a
#'   constant level regressor
#' @param intercept_u logical; if \code{TRUE} includes an intercept term i.e. a
#'   constant level regressor
#' @param policy_dummy logical; if \code{TRUE} includes a dummy that jumps from
#'   zero to one (e.g. a policy or other jump effect to be modelled)
#' @param zero_pattern double with possible values 1, 2, 3 or 4:
#'   \itemize{
#'     \item{1: }{A dummy pattern that starts at the beginning with zeros and
#'     jumps after half of the overall time period}
#'     \item{2: }{A dummy pattern that starts at the beginning with ones and
#'     plummets to zeros after half of the overall time period}
#'     \item{3: }{A dummy pattern that starts at the beginning with one, then
#'     plummets to zeros after a third of the overall time periods, and then
#'     reverts back to ones for the last third of the time}
#'     \item{4: }{A dummy pattern that starts at the beginning with zeros,
#'     then jumps to ones after a third of the overall time periods, and then
#'     reverts back to zeros for the last third of the time}
#'   }
#' @param drift  TO-BE-EXPLAINED-LATER
#'
#' @return a \code{length(beta_x)}x\code{TT} matrix of regressors and a
#'   \code{TT}-dimensional vector of latent states
generate_x_z_u <- function(TT,
                           phi_x,
                           sig_sq_x,
                           bet_z = NULL,
                           bet_u = NULL,
                           bet_z_spl = NULL,
                           bet_u_spl = NULL,
                           options_reg_simul,
                           modelling_reg_types,
                           intercept_z,
                           intercept_u,
                           policy_dummy = FALSE,
                           zero_pattern = NULL,
                           drift = FALSE) {
  order_p <- length(phi_x)
  if (!modelling_reg_types[["autoregression"]]) phi_x <- 0
  x_level <- options_reg_simul[["x_level"]]
  x_sd_within  <- sqrt(options_reg_simul[["reg_var_within"]])
  x_sd_among   <- sqrt(options_reg_simul[["reg_var_among"]])
  spl_lvl <- c(0.7, 0.3)
  dim_z <- length(bet_z)
  dim_u <- length(bet_u)
  # BEGINNING OF REGRESSOR SIMULATION: --------------------------------------
  if (modelling_reg_types[["z-linear-regressors"]]) {
    x_level_split_z <- ifelse(modelling_reg_types[["z-linear-regressors"]],
                              x_level * spl_lvl[1],
                              x_level)
    policy_dummy_z  <- ifelse(modelling_reg_types[["u-linear-regressors"]],
                              FALSE,
                              policy_dummy)
    z <- generate_reg_vals(TT = TT,
                           bet_reg = bet_z,
                           dim_reg = dim_z,
                           phi_x = phi_x,
                           x_level = x_level_split_z,
                           x_sd_within = x_sd_within,
                           x_sd_among = x_sd_among,
                           intercept = intercept_z,
                           policy_dummy = policy_dummy_z,
                           zero_pattern = NULL)
  } else {
    z <- NULL
  }
  if (modelling_reg_types[["u-linear-regressors"]]) {
    x_level_split_u <- ifelse(modelling_reg_types[["z-linear-regressors"]],
                              x_level * spl_lvl[2],
                              x_level)
    u <- generate_reg_vals(TT = TT,
                           bet_reg = bet_u,
                           dim_reg = dim_u,
                           phi_x = phi_x,
                           x_level = x_level_split_u,
                           x_sd_within = x_sd_within,
                           x_sd_among = x_sd_among,
                           intercept = intercept_u,
                           policy_dummy = policy_dummy,
                           zero_pattern = NULL)
  } else {
    u <- NULL
  }
  regs_all <- cbind(z, u)
  # END OF REGRESSOR SIMULATION: --------------------------------------------
  out <- list()
  # browser()
  out$x <- simulate_x(x_level, regs = regs_all,
                      phi_x, sig_sq_x,
                      bet_reg = c(bet_z, bet_u),
                      TT, order_p)
  out$z <- z
  out$u <- u
  return(out)
}
simulate_x <- function(x_level, regs, phi_x, sig_sq_x, bet_reg, TT, order_p) {
  if (is.null(regs) && is.null(bet_reg)) {
    regs    <- matrix(0, nrow = TT + 1, ncol = 1)
    bet_reg <- matrix(0, nrow = 1, ncol = 1)
  }
  x <- rep(0, TT)
  xinit <- x_level
  iter  <- 1
  x[1] <- f(x_tt = xinit,
            regs = regs[1, ],
            phi_x = phi_x[1],
            bet_reg = bet_reg)
  x[1] <- x[1] + sqrt(sig_sq_x) * stats::rnorm(n = 1)
  while (iter < order_p) {
    iter    <- iter + 1
    x[iter] <- f(x_tt = c(x[(iter - 1):1], xinit),
                 regs = regs[iter, ],
                 phi_x = phi_x[1:iter],
                 bet_reg = bet_reg)
    x[iter] <- x[iter] + sqrt(sig_sq_x) * stats::rnorm(n = 1)
  }

  for (t in iter:TT) {
    if (t < TT) {
      x[t + 1] <- f(x_tt = x[t:(t - order_p + 1)],
                    regs = regs[t + 1, ],
                    phi_x = phi_x, bet_reg = bet_reg)
      x[t + 1] <- x[t + 1] + sqrt(sig_sq_x) * stats::rnorm(n = 1)
    }
  }
  return(x)
}
#' Generates the actual regressor values
#'
#' Generates all types of regressors values including \code{z} (standard
#' "overall regressors"), \code{u} i.e. random effect regressors affecting only
#' the cross sectional unit and the corresponding spline version \code{z_spl}
#' and \code{u_spl}
#'
#' @param TT number of time periods
#' @param bet_reg merged betas including those for \code{z} (\code{bet_z}),
#'   \code{u} (\code{bet_u}), \code{z_spl} (\code{bet_z_spl}) and
#'   \code{u_spl} (\code{bet_u_spl})
#' @param dim_reg overall dimension including all regressors (i.e. length of
#'   \code{bet_reg})
#' @param phi_x autoregressive parameter
#' @param x_level target level of the states to simulate
#' @param x_sd_within standard deviations for the regressor value simulation
#'   within each time series
#' @param x_sd_among standard deviations for the regressor value simulation
#'   among regressor components (for a given cross section)
#' @param intercept logical; if \code{TRUE}, then an intercept is included
#' @param policy_dummy logical; if \code{TRUE}, then a policy dummy is included
#' @param zero_pattern double with possible values 1, 2, 3 or 4:
#'   \itemize{
#'     \item{1: }{A dummy pattern that starts at the beginning with zeros and
#'     jumps after half of the overall time period}
#'     \item{2: }{A dummy pattern that starts at the beginning with ones and
#'     plummets to zeros after half of the overall time period}
#'     \item{3: }{A dummy pattern that starts at the beginning with one, then
#'     plummets to zeros after a third of the overall time periods, and then
#'     reverts back to ones for the last third of the time}
#'     \item{4: }{A dummy pattern that starts at the beginning with zeros,
#'     then jumps to ones after a third of the overall time periods, and then
#'     reverts back to zeros for the last third of the time}
#'   }
#'
#' @return a matrix of regressors of dimension \code{TT} \eqn{x} \code{dim_bet}
generate_reg_vals <- function(TT, bet_reg, dim_reg, phi_x,
                              x_level, x_sd_within, x_sd_among,
                              intercept, policy_dummy, zero_pattern = NULL) {
  order_p <- length(phi_x)
  dummy_to_use <- get_pattern_policy_zeros(policy_dummy, zero_pattern, TT)

  reg_means <- stats::rnorm(max(dim_reg - 3, 0), mean = 0, sd = x_sd_among)
  if (!intercept && !policy_dummy) {
    num_add_sims <- max(min(2, dim_reg - 1), 0)
    reg_means <- c(stats::rnorm(num_add_sims, mean = 0, sd = x_sd_among),
                   reg_means)
  } else if (intercept && !policy_dummy) {
    num_add_sims <- max(min(2, dim_reg - 2), 0)
    reg_means <- c(stats::rnorm(num_add_sims, mean = 0, sd = x_sd_among),
                   reg_means)
    if (dim_reg > 1) reg_means <- c(1, reg_means); message("Appr. simul. type1")
  } else if (!intercept && policy_dummy) {
    num_add_sims <- max(min(2, dim_reg - 2), 0)
    reg_means <- c(stats::rnorm(num_add_sims, mean = 0, sd = x_sd_among),
                   reg_means)
    if (dim_reg > 1) reg_means <- c(0, reg_means); message("Appr. simul. type2")
  } else if (intercept && policy_dummy) {
    if (dim_reg < 2) stop("Can't simulate intercept&policy_dummy for dim < 2!")
    num_add_sims <- max(min(2, dim_reg - 3), 0)
    reg_means <- c(stats::rnorm(num_add_sims, mean = 0, sd = x_sd_among),
                   reg_means)
    if (dim_reg > 2) reg_means <- c(1, 0, reg_means); message("Approx. type3")
    if (dim_reg > 1) reg_means <- c(1, reg_means); message("Approx. type4")
  }
  last_reg_mean <- x_level * (1 - sum(phi_x))
  last_reg_mean <- last_reg_mean - sum(reg_means * bet_reg[-dim_reg])
  last_reg_mean <- last_reg_mean/bet_reg[dim_reg]
  reg_means     <- c(reg_means, last_reg_mean)

  reg_len       <- length(reg_means)
  regs          <- matrix(stats::rnorm(TT*reg_len,
                                       mean = reg_means,
                                       sd = x_sd_within),
                          nrow = TT,
                          ncol = reg_len,
                          byrow = TRUE)

  if (intercept && !policy_dummy) regs[, 1]   <- 1
  if (!intercept && policy_dummy) regs[, 1]   <- policy_dummy
  if (intercept && policy_dummy)  regs[, 1:2] <- c(1, policy_dummy)

  return(regs)
}
get_pattern_policy_zeros <- function(policy_dummy, zero_pattern, TT) {
  if (policy_dummy) {
    if (zero_pattern == 1) {
      dummy_to_use <- c(rep(0, times = round(0.5*TT, digits = 0)),
                        rep(1, times = TT - round(0.5*TT, digits = 0)))
    } else if (zero_pattern == 2) {
      dummy_to_use <- c(rep(1, times = TT - round(0.5*TT, digits = 0)),
                        rep(0, times = round(0.5*TT, digits = 0)))
    } else if (zero_pattern == 3) {
      my_third <-  round(1/3*TT, digits = 0)
      dummy_to_use <- c(rep(1, times = my_third),
                        rep(0, times = my_third),
                        rep(1, times = TT - 2*my_third))
    } else if (zero_pattern == 4) {
      my_third <-  round(1/3*TT, digits = 0)
      dummy_to_use <- c(rep(0, times = my_third),
                        rep(1, times = my_third),
                        rep(0, times = TT - 2*my_third))
    } else {
      stop(paste0("Undefined zero patterns: please use 1L to 4L for different",
                  "patterns and check the doc/help for their meaning!"))
    }
    return(dummy_to_use)
  } else {
    return(NULL)
  }
}
#' State transition
#'
#' Helper function computing the deterministic state transition, or, to put
#' differently the one-period ahead conditional mean of the latent state
#' process.
#'
#' @param x_tt state value in t-1 i.e. x_{t-1}
#' @param regs regressor values e.g. z_{t}, or u_{t} or both (matrix)
#' @param phi_x autoregressive parameter phi
#' @param bet_reg regressor parameters/coefficients e.g. at z_{t}, or u_{t} or
#'   both (matrix)
#'
#' @return \code{T}-dimensional vector of (deterministically computed) state
#'   transitions (conditional means)
f <- function(x_tt, regs, phi_x, bet_reg) {
  # xt <- phi_x*xtt
  x_t <- sum(phi_x*x_tt) + regs %*% bet_reg
  # x_t <- phi_x*x_tt + z %*% bet_z + u %*% bet_u

  # xt <- phi_x*xtt + 8*cos(1.2*t)
  # xt <- phi_x*xtt + 25*xtt/(1 + xtt^2)
  # xt <- phi_x*xtt + 25*xtt/(1 + xtt^2) + 8*cos(1.2*t)
  return(x_t)
}
