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
                           # policy_dummy = FALSE,
                           zero_pattern = NULL,
                           drift = FALSE) {
  if (!modelling_reg_types[["autoregression"]]) phi_x <- 0;
  order_p <- length(phi_x)
  x_level <- options_reg_simul[["x_level"]]
  x_sd_within  <- sqrt(options_reg_simul[["reg_var_within"]])
  x_sd_among   <- sqrt(options_reg_simul[["reg_var_among"]])
  spl_lvl <- c(0.6, 0.4)
  dim_z <- length(bet_z)
  dim_u <- length(bet_u)
  # BEGINNING OF REGRESSOR SIMULATION: --------------------------------------
  if (modelling_reg_types[["z-linear-regressors"]]) {
    # x_level_split_z <- ifelse(modelling_reg_types[["u-linear-regressors"]],
    #                           x_level * spl_lvl[1],
    #                           x_level)
    # policy_dummy_z  <- ifelse(modelling_reg_types[["u-linear-regressors"]],
    #                           FALSE,
    #                           policy_dummy)
    lvl_split_z <- get_sub_level_x(modelling_reg_types[["u-linear-regressors"]],
                                   level_split = spl_lvl[1],
                                   level_target = x_level,
                                   phi_x)
    z <- generate_reg_vals(TT = TT,
                           bet_reg = bet_z,
                           dim_reg = dim_z,
                           x_level = lvl_split_z,
                           x_sd_within = x_sd_within,
                           x_sd_among = x_sd_among,
                           intercept = intercept_z,
                           # policy_dummy = policy_dummy,
                           zero_pattern = zero_pattern)
  } else {
    z <- NULL
  }
  if (modelling_reg_types[["u-linear-regressors"]]) {
    # x_level_split_u <- ifelse(modelling_reg_types[["z-linear-regressors"]],
    #                           x_level * spl_lvl[2],
    #                           x_level)
    lvl_split_u <- get_sub_level_x(modelling_reg_types[["z-linear-regressors"]],
                                   level_split = spl_lvl[2],
                                   level_target = x_level,
                                   phi_x)
    u <- generate_reg_vals(TT = TT,
                           bet_reg = bet_u,
                           dim_reg = dim_u,
                           x_level = lvl_split_u,
                           x_sd_within = x_sd_within,
                           x_sd_among = x_sd_among,
                           intercept = intercept_u,
                           # policy_dummy = policy_dummy,
                           zero_pattern = zero_pattern)
  } else {
    u <- NULL
  }
  # kappa_u <- kappa(t(u) %*% u)
  # kappa_z <- kappa(t(z) %*% z)
  # if (kappa_z > 10^5) stop("Z regs might be ill-conditioned.")
  # if (kappa_u > 10^5) stop("U regs might be ill-conditioned.")
  regs_all <- cbind(z, u)
  # END OF REGRESSOR SIMULATION: --------------------------------------------
  out <- list()
  out$x <- simulate_x(x_level, regs = regs_all,
                      phi_x, sig_sq_x,
                      bet_reg = c(bet_z, bet_u),
                      TT, order_p)
  out$z <- z
  out$u <- u
  return(out)
}

get_sub_level_x <- function(other_regtype, level_split,
                            level_target, phi) {
  x_level_tmp <- ifelse(other_regtype,
                        level_target * level_split,
                        level_target)
  x_level_tmp * (1-sum(phi))
}
simulate_x <- function(x_level, regs, phi_x, sig_sq_x, bet_reg, TT, order_p) {
  if (all(regs == 0)) return(rep(0, TT))
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
#' @param x_level target level of the states to simulate
#' @param x_sd_within standard deviations for the regressor value simulation
#'   within each time series
#' @param x_sd_among standard deviations for the regressor value simulation
#'   among regressor components (for a given cross section)
#' @param intercept logical; if \code{TRUE}, then an intercept is included
#' @param zero_pattern double with possible values 1, 2, 3 or 4:
#'   \itemize{
#'     \item{1: }{A dummy pattern representing structural zeros i.e. only zeros
#'     present}
#'     \item{2: }{A dummy pattern that starts at the beginning with zeros and
#'     jumps after half of the overall time period}
#'     \item{3: }{A dummy pattern that starts at the beginning with ones and
#'     plummets to zeros after half of the overall time period}
#'     \item{4: }{A dummy pattern that starts at the beginning with one, then
#'     plummets to zeros after a third of the overall time periods, and then
#'     reverts back to ones for the last third of the time}
#'     \item{5: }{A dummy pattern that starts at the beginning with zeros,
#'     then jumps to ones after a third of the overall time periods, and then
#'     reverts back to zeros for the last third of the time}
#'   }
#'
#' @return a matrix of regressors of dimension \code{TT} \eqn{x} \code{dim_bet}
generate_reg_vals <- function(TT, bet_reg, dim_reg,
                              x_level, x_sd_within, x_sd_among,
                              intercept, zero_pattern = NULL) {
  # order_p <- length(phi_x)
  zero_dummy_to_use <- get_pattern_policy_zeros(zero_pattern, TT)
  ZERO_PATTERN      <- get_zero_pattern_flag(zero_dummy_to_use)
  # reg_means <- stats::rnorm(max(dim_reg - 3, 0), mean = 0, sd = x_sd_among)
  if (!intercept && isFALSE(ZERO_PATTERN)) {
    weights <- rep(1/dim_reg, times = dim_reg) #(abs(bet_reg)/sum(abs(bet_reg)))
    reg_means <- x_level * weights
    reg_means <- reg_means / bet_reg
    # num_add_sims <- max(min(2, dim_reg - 1), 0)
    # reg_means <- c(stats::rnorm(num_add_sims, mean = 0, sd = x_sd_among),
    #                reg_means)
  } else if (intercept && isFALSE(ZERO_PATTERN)) {
    bet_tmp <- bet_reg[-1]
    dim_tmp <- length(bet_tmp)
    x_lvl_tmp <- x_level - bet_reg[1]
    weights <- rep(1/dim_tmp, times = dim_tmp) #(abs(bet_tmp)/sum(abs(bet_tmp)))
    reg_means <- x_lvl_tmp * weights
    reg_means <- reg_means / bet_reg[-1]
    # num_add_sims <- max(min(2, dim_reg - 2), 0)
    # reg_means <- c(stats::rnorm(num_add_sims, mean = 0, sd = x_sd_among),
    #                reg_means)
    # if (dim_reg > 1) reg_means <- c(1, reg_means); message("Appr. simul. type1")
  } else if (!intercept && isTRUE(ZERO_PATTERN)) {
    if (all(zero_dummy_to_use == 0)) {
      reg_means   <- rep(0, times = dim_reg)
      x_sd_within <- 0
    } else {
      stop("Not yet implemented.")
    }
    # num_add_sims <- max(min(2, dim_reg - 2), 0)
    # reg_means <- c(stats::rnorm(num_add_sims, mean = 0, sd = x_sd_among),
    #                reg_means)
    # if (dim_reg > 1) reg_means <- c(0, reg_means); message("Appr. simul. type2")
  } else if (intercept && isTRUE(ZERO_PATTERN)) {
    if (all(zero_dummy_to_use == 0)) {
      reg_means   <- rep(0, times = dim_reg)
      x_sd_within <- 0
    } else {
      stop("Not yet implemented.")
    }
    # if (dim_reg < 2) stop("Can't simulate intercept&policy_dummy for dim < 2!")
    # num_add_sims <- max(min(2, dim_reg - 3), 0)
    # reg_means <- c(stats::rnorm(num_add_sims, mean = 0, sd = x_sd_among),
    #                reg_means)
    # if (dim_reg > 2) reg_means <- c(1, 0, reg_means); message("Approx. type3")
    # if (dim_reg > 1) reg_means <- c(1, reg_means); message("Approx. type4")
  }
  # last_reg_mean <- x_level * (1 - sum(phi_x))
  # last_reg_mean <- last_reg_mean - sum(reg_means * bet_reg[-dim_reg])
  # last_reg_mean <- last_reg_mean/bet_reg[dim_reg]
  # reg_means     <- c(reg_means, last_reg_mean)

  reg_len       <- length(reg_means)
  regs          <- matrix(stats::rnorm(TT*reg_len,
                                       mean = reg_means,
                                       sd = x_sd_within),
                          nrow = TT,
                          ncol = reg_len,
                          byrow = TRUE)

  if (isTRUE(intercept) && isFALSE(ZERO_PATTERN))  return(cbind(1, regs))
  if (isFALSE(intercept) && isFALSE(ZERO_PATTERN)) return(regs)
  if (isFALSE(intercept) && isTRUE(ZERO_PATTERN))  return(regs) # regs[, 1] <- policy_dummy
  if (isTRUE(intercept) && isTRUE(ZERO_PATTERN))   return(regs) # regs[, 1:2] <- c(1, policy_dummy)
}
get_pattern_policy_zeros <- function(zero_pattern, TT) {
  if (is.null(zero_pattern)) return(NULL)
  stopifnot(`Arg. 'zero_pattern' must be logical or numeric` =
              is.logical(zero_pattern) || is.numeric(zero_pattern))
  # Case I: numeric zero_pattern - policy dummies and/or/no structural zeros
  if (is.numeric(zero_pattern)) {
    if (zero_pattern == 1) {
      dummy_to_use <- rep(0, times = TT)
    }
    if (zero_pattern == 2) {
      dummy_to_use <- c(rep(0, times = round(0.5*TT, digits = 0)),
                        rep(1, times = TT - round(0.5*TT, digits = 0)))
    } else if (zero_pattern == 3) {
      dummy_to_use <- c(rep(1, times = TT - round(0.5*TT, digits = 0)),
                        rep(0, times = round(0.5*TT, digits = 0)))
    } else if (zero_pattern == 4) {
      my_third <-  round(1/3*TT, digits = 0)
      dummy_to_use <- c(rep(1, times = my_third),
                        rep(0, times = my_third),
                        rep(1, times = TT - 2*my_third))
    } else if (zero_pattern == 5) {
      my_third <-  round(1/3*TT, digits = 0)
      dummy_to_use <- c(rep(0, times = my_third),
                        rep(1, times = my_third),
                        rep(0, times = TT - 2*my_third))
    } else {
      stop(paste0("Undefined zero patterns: please use 1L to 4L for different",
                  "patterns and check the doc/help for their meaning!"))
    }
  } else if (is.logical(zero_pattern)) {
    # Case I: logical zero_pattern - only true structural zeros or no zeros
    if (isTRUE(zero_pattern)) {
      dummy_to_use <- rep(0, times = TT)
    } else if (isFALSE(zero_pattern)) {
      dummy_to_use <- rep(1, times = TT)
    }
  } else {
    stop("Unknown behavior/case; cannot implement zero pattern.")
  }
  return(dummy_to_use)
}
get_zero_pattern_flag <- function(zero_dummy_to_use) {
  if (is.null(zero_dummy_to_use) || all(zero_dummy_to_use == 1)) {
    return(FALSE)
  } else {
    return(TRUE)
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
