#' Deduces from vector of parameter names, which type of modelling to employ
#'
#' @param get_modelling_reg_types names of parameters list elements; the list
#'   of true parameter values is taken
#'
#' @return logical vector of dimension 4 given \code{TRUE} or \code{FALSE} if
#'   (in this order) modelling of z-regressors, u-regressors, z spline regegressors
#'   or u spline regegressors should be performed
#'
get_modelling_reg_types <- function(get_modelling_reg_types) {
  correct_names <- c("bet_z", "bet_u", "bet_z_spl", "bet_u_spl")
  if (length(intersect(get_modelling_reg_types, correct_names)) == 0) {
    stop("The 'par_true' argument must have correct names: choose from 'bet_z', 'bet_u', 'bet_z_spl' or 'bet_u_spl'! ")
  }
  out <- vector("logical", 4)
  out[1] <- correct_names[1] %in% get_modelling_reg_types
  out[2] <- correct_names[2] %in% get_modelling_reg_types
  out[3] <- correct_names[3] %in% get_modelling_reg_types
  out[4] <- correct_names[5] %in% get_modelling_reg_types
  return(out)
}
#' Generates random effects per cross sectional unit (e.g. US-state)
#'
#' @param DD dimension of the latent state process
#' @param NN number of cross sectional units
#' @param from_IW logical; if \code{TRUE}, then random effects (per dimension of
#'   the state process \code{d=1,...,DD}) are generated from an IW-distribution
#'   with internally specified degrees of freedom and scale matrix s.th. per d,
#'   the NN different random effect vectors are i.i.d
#' @param num_re integer vector of dimension \code{DD} if \code{from_IW == TRUE}
#'   or a scalar integer value if \code{from_IW == FALSE} that specifies the
#'   number of random effects per component \code{d=1,...,DD}
#' @param seed_no integer; random seed to set at the beginning, set \code{NULL}
#'   if not required
#'
#' @return a list of \code{DD} components, each containing a matrix of dimension
#'   \code{num_re[d]} $times$ \code{NN}
#'
#' @export
generate_bet_u <- function(DD, NN, from_IW = FALSE,
                           num_re = rep(2, times = DD),
                           seed_no = 42) {
  if (from_IW) {
    stopifnot(is.numeric(num_re) && (length(num_re) == DD))
    if (!is.null(seed_no))  set.seed(seed_no)
    n0u <- num_re + 5
    D0u <- vector("list", DD)
    true_bet_u <- vector("list", DD)
    for (d in 1:DD) {
      D0u[[d]] <- solve(diag(1:num_re[d]*(10 + 1/num_re[d])))
      D0u[[d]] <- (stats::rWishart(1, n0u[d], D0u[[d]]))[, , 1]
      true_bet_u[[d]] <- matrix(0, nrow = num_re[d], ncol = NN)
      for (n in 1:NN) {
        true_bet_u[[d]][, n] <- MASS::mvrnorm(n = 1,
                                              mu = rep(0, times = num_re[d]),
                                              Sigma = D0u[[d]])
      }
    }
    return(list(true_bet_u, D0u))
  } else {
    stopifnot(is.numeric(num_re) && (length(num_re) == 1))
    true_bet_u <- vector("list", DD)
    vals <- matrix(1:(num_re*DD)*c(-1, 1), nrow = num_re*DD, ncol = NN)
    for (d in 1:DD) {
      true_bet_u[[d]] <- vals[1:num_re + num_re*(d - 1), ]
    }
    return(true_bet_u)
  }
}
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
#' @param include_intercept logical vector of dimension \code{DD}; if \code{TRUE}
#'   include an intercept at the cross sectional unit for component \code{d}
#' @param include_policy logical vector of dimension \code{DD}; if \code{TRUE}
#'   include a policy dummy at the cross sectional unit for component \code{d}
#' @param include_zeros numeric vector of dimension \code{DD} with values 1,
#'   2, 3 or 4:
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
#' @param seed_no integer; random seed to set at the beginning, set \code{NULL}
#'   if not required
#' @param plot_measurements logical; if \code{TRUE}, measurements are plotted
#'   per cross sectional unit \code{n=1,...,N}
#' @param plot_states logical; if \code{TRUE}, latent states are plotted
#'   per cross sectional unit \code{n=1,...,N}
#'
#' @return NN-dimensional list of lists of two: \code{[[1]]} -> regressors and
#'   \code{[[2]]} -> latent states
#' @export
generate_data_t_n <- function(distribution = "dirichlet",
                              NN, TT, DD,
                              par_true,
                              x_levels,
                              x_log_scale,
                              include_intercept,
                              include_policy,
                              include_zeros,
                              plot_measurements = FALSE,
                              plot_states = FALSE,
                              seed_no = NULL) {
  if (x_log_scale) {
    x_levels <- log(x_levels)
  }
  if (!is.null(seed_no)) set.seed(seed_no)

  modelling_reg_types <- names(par_true)
  modelling_reg_types <- get_modelling_reg_types(modelling_reg_types)

  data_part1 <- array(0, c(TT, DD, NN))
  if (distribution == "multinomial" || distribution == "mult-diri") {
    data_part2 <- matrix(0, nrow = TT, ncol = NN)
  }

  states <- array(0, c(TT, DD, NN))

  if (modelling_reg_types[1]) {
    dim_bet_z <- sum(sapply(par_true[["bet_z"]], length))
    z <- array(0, c(TT, dim_bet_z, NN))
  }
  if (modelling_reg_types[2]) {
    dim_bet_u <- sum(sapply(par_true[["bet_u"]], length))
    u <- array(0, c(TT, dim_bet_u, NN))
  }
  for (n in 1:NN) {
    par_true_current <- list(par_true[[1]][, n],
                             par_true[[2]][, n])
    # browser()
    if (modelling_reg_types[1]) par_true_current$bet_z <- par_true[["bet_z"]]
    if (modelling_reg_types[2]) par_true_current$bet_u <- par_true[["bet_u"]]

    out_data_tmp <- generate_data_t(distribution = distribution,
                                    TT = TT, DD = DD, n = n,
                                    par_true = par_true_current,
                                    x_levels = x_levels[, n],
                                    x_log_scale = x_log_scale,
                                    include_intercept = include_intercept,
                                    include_policy = include_policy[, n, drop = TRUE],
                                    include_zeros = include_zeros,
                                    modelling_reg_types = modelling_reg_types)
    data_part1[, , n] <- out_data_tmp[[1]]$yraw
    if (distribution == "multinomial" || distribution == "mult-diri") {
      data_part2[, n]   <- out_data_tmp[[1]]$num_counts
    }
    if (modelling_reg_types[1]) {
      z[, , n] <- out_data_tmp$z

    }
    if (modelling_reg_types[2]) {
      u[, , n] <- out_data_tmp$u
    }
    states[, , n] <- out_data_tmp$x
    plot_data_per_n(DD,
                    yraw = data_part1[, , n],
                    x = states[, , n],
                    x_log_scale = x_log_scale,
                    x_levels = x_levels[, n],
                    plot_measurements = plot_measurements,
                    plot_states       = plot_states)
  }
  out_data <- vector("list", 3)
  names(out_data) <- c("data", "regs", "states")
  if (distribution == "dirichlet") {
    out_data[[1]] <- list(yraw = data_part1)
  }
  if (distribution == "multinomial" || distribution == "mult-diri") {
    out_data[[1]] <- list(yraw = data_part1, num_counts = data_part2)
  }
  if (modelling_reg_types[1]) {
    out_data[[2]]$z <- z
  }
  if (modelling_reg_types[2]) {
    out_data[[2]]$u <- u
  }
  out_data[[3]] <- states

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
#' @param n current cross sectional unit i.e. an integer \code{n=1,...,NN}
#' @param par_true list of true parameters that describe the latent state
#'   process
#' @param x_levels vector of target "means"/"levels" of the states around which
#'   they fluctuate
#' @param x_log_scale logical; if \code{TRUE}, then \code{x_levels[i]} ist taken
#'   as log and the random number generation for the states is performed in logs
#' @param include_intercept logical vector of dimension \code{DD}; if \code{TRUE}
#'   include an intercept at the cross sectional unit for component \code{d}
#' @param include_policy logical vector of dimension \code{DD}; if \code{TRUE}
#'   include a policy dummy at the cross sectional unit for component \code{d}
#' @param include_zeros numeric vector of dimension \code{DD} with values 1,
#'   2, 3 or 4:
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
#' @param modelling_reg_types logical vector of dimension 4 where each component
#'   indicates, if \code{TRUE}, that the corresponding regressor type is to be
#'   generated in the following order:
#'   \itemize{
#'     \item{modelling_reg_types, component 1: }{includes z regressors}
#'     \item{modelling_reg_types, component 2: }{includes u regressors}
#'     \item{modelling_reg_types, component 3: }{includes z spline regressors}
#'     \item{modelling_reg_types, component 4: }{includes u spline regressors}
#'   }
#'
#' @return a list of two: \code{[[1]]} -> regressors and \code{[[2]]} -> latent
#'   states
generate_data_t <- function(distribution = "dirichlet",
                            TT, DD, n,
                            par_true,
                            x_levels,
                            x_log_scale,
                            include_intercept,
                            include_policy,
                            include_zeros,
                            modelling_reg_types) {
  x <- matrix(nrow = TT, ncol = DD, 0)
  sig_sq_x <- par_true[[1]]
  phi_x <- par_true[[2]]
  if (modelling_reg_types[1]) {
    z <- list()
    bet_z <- list()
    for (d in 1:DD) {
      bet_z[[d]] <- par_true[["bet_z"]][[d]]
    }
  } else {
    bet_z <- NULL
  }
  if (modelling_reg_types[2]) {
    u <- list()
    bet_u <- list()
    for (d in 1:DD) {
      bet_u[[d]] <- par_true[["bet_u"]][[d]]
    }
  } else {
    bet_u <- NULL
  }


  if (distribution %in% c("multinomial", "mult-diri", "mult-gen-diri")) {
    num_counts <- sample(x = 80000:120000, size = TT)
  }

  x_sd_level <- c(0.0125, 0.1, 0.025, 0.1, 0.1, 0.1)

  for (d in 1:DD) {
    # browser()
    res <- generate_x_z(TT = TT,
                        phi_x = phi_x[d],
                        sig_sq_x = sig_sq_x[d],
                        bet_z = bet_z[[d]],
                        bet_u = bet_u[[d]][, n],
                        modelling_reg_types = modelling_reg_types,
                        x_level = x_levels[d],
                        x_sd = x_sd_level[d],
                        x_log_scale = x_log_scale,
                        intercept      = include_intercept[d],
                        policy_dummy   = include_policy[d],
                        zero_pattern   = include_zeros[d])
    # browser()
    x[, d] <- res$x
    if (modelling_reg_types[1]) {
      z[[d]] <- res$z
    }
    if (modelling_reg_types[2]) {
      u[[d]] <- res$u
    }
  }
  if (x_log_scale) {
    x <- exp(x)
  }
  alphas <- x

  out_data_all <- list()
  if (distribution == "dirichlet") {
    yraw <- my_rdirichlet(alpha = alphas)
    if (sum(rowSums(yraw)) != TT) {
      stop("Something is wrong with Dirichelet: fractions don't sum up to 1!")
    }
    out_data_all$data <- list(yraw = yraw)
  }
  if (distribution == "multinomial") {
    alphas <- alphas/rowSums(alphas)
    yraw <- my_rmultinomial(probs = alphas, num_counts = num_counts)
    out_data_all$data <- list(yraw = yraw, num_counts = num_counts)
  }
  if (distribution == "mult-diri") {
    yraw <- my_rmult_diri(alpha = alphas, num_counts = num_counts)
    out_data_all$data <- list(yraw = yraw, num_counts = num_counts)
  }
  if (modelling_reg_types[1]) {
    z <- Reduce(cbind, z)
    out_data_all$z <- z
  }
  if (modelling_reg_types[2]) {
    u <- Reduce(cbind, u)
    out_data_all$u <- u
  }
  out_data_all$x <- x
  return(out_data_all)
}
#' Simulates regressors and latent states.
#'
#' The core function that simulates the regressors and latent states given
#' certain parameters.
#'
#' @param TT number of time periods
#' @param phi_x autoregressive parameter of latent states
#' @param sig_sq_x standard deviation of error term in latent state transition
#'   equation
#' @param bet_z regressor coefficients/parameters of latent state process
#'   referring to all cross secitonal units  and are non-linear i.e. splines
#' @param bet_u regressor coefficients/parameters of latent state process;
#'   random effects that refer to each individual cross sectional unit
#' @param bet_z_spl regressor coefficients/parameters of latent state process
#'   referring to all cross secitonal units
#' @param bet_u_spl regressor coefficients/parameters of latent state process;
#'   random effects that refer to each individual cross sectional unit and are
#'   non-linear i.e. splines
#' @param modelling_reg_types logical vector of dimension 4 where each component
#'   indicates, if \code{TRUE}, that the corresponding regressor type is to be
#'   generated in the following order:
#'   \itemize{
#'     \item{modelling_reg_types, component 1: }{includes z regressors}
#'     \item{modelling_reg_types, component 2: }{includes u regressors}
#'     \item{modelling_reg_types, component 3: }{includes z spline regressors}
#'     \item{modelling_reg_types, component 4: }{includes u spline regressors}
#'   }
#' @param x_level target levels of latent states around which they fluctuate (a
#'   tuning parameter that ensures that states at each multivariate component of
#'   the response fluctuate around a particular level)
#' @param x_sd standard deviation allowed for \code{x_level} target (a tuning
#'   parameter that ensures that deviations from the \code{x_level} target are
#'   not too severe)
#' @param x_log_scale logical; if \code{TRUE} process is simulated on the
#'   log-scale
#' @param intercept logical; if \code{TRUE}, includes an intercept term i.e. a
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
#' @return a \code{length(beta_x)}x\code{TT} matrix of regressors and a \code{TT}-dimensional vector of latent states
generate_x_z <- function(TT,
                         phi_x,
                         sig_sq_x,
                         bet_z = NULL,
                         bet_u = NULL,
                         bet_z_spl = NULL,
                         bet_u_spl = NULL,
                         modelling_reg_types,
                         x_level,
                         x_sd,
                         x_log_scale,
                         intercept,
                         policy_dummy = FALSE,
                         zero_pattern = 1,
                         drift = FALSE) {
  dim_z <- length(bet_z)
  dim_u <- length(bet_u)
  dim_reg <- dim_z + dim_u
  bet_reg <- c(bet_z, bet_u)
  x <- rep(0, TT)
  # BEGINNING OF REGRESSOR SIMULATION: --------------------------------------
  if (modelling_reg_types[1] && !modelling_reg_types[2]) {
    z <- generate_reg_vals(TT = TT,
                           bet_reg = bet_reg,
                           dim_reg = dim_reg,
                           phi_x = phi_x,
                           x_level = x_level,
                           x_sd = x_sd,
                           intercept = intercept,
                           policy_dummy = policy_dummy,
                           zero_pattern = zero_pattern)
    u <- NULL
  }
  if (!modelling_reg_types[1] && modelling_reg_types[2]) {
    u <- generate_reg_vals(TT = TT,
                           bet_reg = bet_reg,
                           dim_reg = dim_reg,
                           phi_x = phi_x,
                           x_level = x_level,
                           x_sd = x_sd,
                           intercept = intercept,
                           policy_dummy = policy_dummy,
                           zero_pattern = zero_pattern)
    z <- NULL
  }
  if (modelling_reg_types[1] && modelling_reg_types[2]) {
    # browser()
    x_level_split <- x_level*c(0.5, 0.5)
    z <- generate_reg_vals(TT = TT,
                           bet_reg = bet_z,
                           dim_reg = dim_z,
                           phi_x = phi_x,
                           x_level = x_level_split[1],
                           x_sd = x_sd,
                           intercept = FALSE,
                           policy_dummy = FALSE,
                           zero_pattern = NULL)
    u <- generate_reg_vals(TT = TT,
                           bet_reg = bet_u,
                           dim_reg = dim_u,
                           phi_x = phi_x,
                           x_level = x_level_split[2],
                           x_sd = x_sd,
                           intercept = intercept,
                           policy_dummy = policy_dummy,
                           zero_pattern = zero_pattern)
  }
  reg_all <- cbind(z, u)
  # END OF REGRESSOR SIMULATION: --------------------------------------------
  xinit <- x_level
  x[1] <- f(x_tt = xinit, regs = reg_all[1, ], phi_x = phi_x, bet_reg = bet_reg)
  x[1] <- x[1] + sqrt(sig_sq_x)*stats::rnorm(n = 1)

  for (t in 1:TT) {
    if (t < TT) {
      x[t + 1] <- f(x_tt = x[t], regs = reg_all[t + 1, ],
                    phi_x = phi_x, bet_reg = bet_reg)
      x[t + 1] <- x[t + 1] + sqrt(sig_sq_x)*stats::rnorm(n = 1)
    }
  }
  # browser()

  if (sum(any(x <= 0)) & x_log_scale == FALSE) {
    stop("some state process (x1_t, x2_t, ... or x6_t) not positive!")
  }
  out <- list()
  out$x <- x
  if (modelling_reg_types[1]) {
    out$z <- z
  }
  if (modelling_reg_types[2]) {
    out$u <- u
  }
  return(out)
}
#' Generates the actual regressor values
#'
#' Generates all types of regressors values including \code{z} (standard
#' "overall regressors"), \code{u} i.e. random effect regressors affecting only
#' the cross sectional unit and the correspodning spline version \code{z_spl}
#' and \code{u_spl}
#'
#' @param TT number of time periods
#' @param bet_reg merged betas including thos for \code{z} (\code{bet_z}),
#'   \code{u} (\code{bet_u}), \code{z_spl} (\code{bet_z_spl}) and
#'   \code{u_spl} (\code{bet_u_spl})
#' @param dim_reg overall dimension including all regressors (i.e. length of
#'   \code{bet_reg})
#' @param phi_x autoregressive paramter
#' @param x_level target level of the states to simulate
#' @param x_sd standard deviations for the regressor values
#' @param intercept logical; if \code{TRUE}, then an intercept is included
#' @param policy_dummy logical; if \code{TRUE}, then a policz dummy is included
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
                              x_level, x_sd,
                              intercept,
                              policy_dummy,
                              zero_pattern = NULL) {
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
      stop("Undefined zero patterns: please use 1L to 4L for different patterns and check the documentation/help for their meaning!")
    }
  }
  if (dim_reg == 1) {
    if (intercept) {
      const_level <- x_level * (1 - phi_x)/bet_reg
      regs <- matrix(const_level, nrow = TT, ncol = dim_reg)
    } else if (policy_dummy) {
      regs <- dummy_to_use
    } else {
      const_mean <- x_level * (1 - phi_x)/bet_reg
      regs          <- matrix(stats::rnorm(T*dim_reg, mean = const_mean, sd = x_sd),
                              nrow = TT,
                              ncol = dim_reg,
                              byrow = TRUE)
    }
    return(regs)
  }
  if (dim_reg == 2) {
    if (intercept && policy_dummy) {
      regs <- cbind(rep(1, times = TT),
                    dummy_to_use)
      return(regs)
    }
    if (intercept && !policy_dummy) {
      reg_means <- 1
    } else if (!intercept && policy_dummy) {
      reg_means <- dummy_to_use
    } else if (!intercept && !policy_dummy) {
      reg_means <- stats::rnorm(dim_reg - 1, mean = 0, sd = 3)
    }
    last_reg_mean <- x_level * (1 - phi_x) - sum(reg_means * bet_reg[-dim_reg])
    last_reg_mean <- last_reg_mean/bet_reg[dim_reg]
    reg_means     <- c(reg_means, last_reg_mean)
    regs          <- matrix(stats::rnorm(TT*dim_reg,
                                         mean = reg_means,
                                         sd = x_sd),
                            nrow = TT,
                            ncol = dim_reg,
                            byrow = TRUE)
    return(regs)
  }
  if (dim_reg > 2) {
    if (intercept && !policy_dummy) {
      reg_means <- stats::rnorm(dim_reg - 2, mean = 0, sd = 3)
      reg_means <- c(1, reg_means)

      last_reg_mean <- x_level * (1 - phi_x) - sum(reg_means * bet_reg[-dim_reg])
      last_reg_mean <- last_reg_mean/bet_reg[dim_reg]
      reg_means     <- c(reg_means, last_reg_mean)
      reg_len       <- length(reg_means)
      regs          <- matrix(stats::rnorm(TT*reg_len,
                                           mean = reg_means,
                                           sd = x_sd),
                              nrow = TT,
                              ncol = reg_len,
                              byrow = TRUE)
        regs[, 1] <- 1
    } else if (!intercept && policy_dummy ) {
      reg_means <- stats::rnorm(dim_reg - 2, mean = 0, sd = 3)
      last_reg_mean <- x_level * (1 - phi_x) - sum(reg_means * bet_reg[-c(1, dim_reg)])
      last_reg_mean <- last_reg_mean/bet_reg[dim_reg]
      reg_means     <- c(reg_means, last_reg_mean)
      reg_len       <- length(reg_means)
      regs          <- matrix(stats::rnorm(TT*reg_len,
                                           mean = reg_means,
                                           sd = x_sd),
                              nrow = TT,
                              ncol = reg_len,
                              byrow = TRUE)
      regs <- cbind(dummy_to_use, regs)
    } else if (policy_dummy && intercept) {
      reg_means <- stats::rnorm(dim_reg - 3, mean = 0, sd = 3)
      reg_means <- c(1, reg_means)

      last_reg_mean <- x_level * (1 - phi_x) - sum(reg_means * bet_reg[-c(2, dim_reg)])
      last_reg_mean <- last_reg_mean/bet_reg[dim_reg]
      reg_means     <- c(reg_means, last_reg_mean)
      reg_len       <- length(reg_means)
      regs          <- matrix(stats::rnorm(TT*reg_len,
                                           mean = reg_means,
                                           sd = x_sd),
                              nrow = TT,
                              ncol = reg_len,
                              byrow = TRUE)
      regs <- cbind(1, dummy_to_use, regs[, -c(1)])
    } else {
      reg_means <- stats::rnorm(dim_reg - 1, mean = 0, sd = 3)
      last_reg_mean <- x_level * (1 - phi_x) - sum(reg_means * bet_reg[-dim_reg])
      last_reg_mean <- last_reg_mean/bet_reg[dim_reg]
      reg_means     <- c(reg_means, last_reg_mean)
      reg_len       <- length(reg_means)
      regs          <- matrix(stats::rnorm(TT*reg_len,
                                           mean = reg_means,
                                           sd = x_sd),
                              nrow = TT,
                              ncol = reg_len,
                              byrow = TRUE)
    }
  return(regs)
  }
}
plot_data_per_n <- function(DD,
                            yraw, x,
                            x_log_scale,
                            x_levels,
                            plot_measurements,
                            plot_states) {
  if (plot_measurements && plot_states)  graphics::par(mfrow = c(2,2))
  if (plot_measurements && !plot_states) graphics::par(mfrow = c(2,1))
  if (!plot_measurements && plot_states) graphics::par(mfrow = c(2,1))
  if (plot_measurements) {
    names_title <- "Measurement components"
    names_ylab  <- "measurements: y_t's"
    names_xlab <- paste0("ya1_t (black),", " ya2_t (red),",
                         " ya3_t (green),", " ya4_t (blue)",
                         " ya5_t (turq.),", " and", " ya6_t (pink)")
    all_measurms <- yraw
    graphics::matplot(all_measurms,
                      type = "l",
                      lty = 1,
                      lwd = 1,
                      main = names_title,
                      ylab = names_ylab,
                      xlab = names_xlab)
    graphics::matplot(all_measurms/rowSums(all_measurms),
                      type = "l",
                      lty = 1,
                      lwd = 1,
                      main = names_title,
                      ylab = names_ylab,
                      xlab = names_xlab)
  }
  if (plot_states) {
    names_title <- "True States"
    names_ylab  <- "states: xt's"
    names_xlab  <- c("x1_t (black),", "x2_t (red),",
                     "x3_t (green),",  "x4_t (blue)",
                     "x5_t (turq.)", "and", "x6_t (pink)")
    all_states <- x
    if (x_log_scale) {
      all_states <- log(all_states)
    }
    graphics::matplot(all_states/rowSums(all_states),
                      type = "l",
                      lty = 2,
                      lwd = 1,
                      main = names_title,
                      ylab = names_ylab,
                      xlab = paste(names_xlab))
    graphics::matplot(all_states,
                      type = "l",
                      lty = 2,
                      lwd = 1,
                      main = names_title,
                      ylab = names_ylab,
                      xlab = names_xlab)
    graphics::par(mfrow = c(6,2))
    for (d in 1:DD) {
      graphics::plot(all_states[, d],
                     type = "l",
                     lty = 2,
                     lwd = 1,
                     main = names_title,
                     ylab = names_ylab,
                     xlab = names_xlab,
                     col = d)
      graphics::abline(h = x_levels[d], lty = 1, lwd = 2)
      graphics::abline(h = mean(all_states[, d, drop = TRUE]),
                       lty = 1, lwd = 1, col = d)

      graphics::plot((all_states/rowSums(all_states))[, d],
                     type = "l",
                     lty = 2,
                     lwd = 1,
                     main = names_title,
                     ylab = names_ylab,
                     xlab = names_xlab,
                     col = d)
      graphics::abline(h = (x_levels/sum(x_levels))[d], lty = 1, lwd = 3)
      graphics::abline(h = mean((all_states/rowSums(all_states))[, d, drop = TRUE]),
                       lty = 1, lwd = 1, col = d)
    }
  }
  graphics::par(mfrow = c(1, 1))
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
#' @param regs regressor values e.g. z_{t}, or u_{t} or both (matrix)
#' @param phi_x autoregressive parameter phi
#' @param bet_reg regressor parameters/coefficients e.g. at z_{t}, or u_{t} or
#'   both (matrix)
#'
#' @return \code{T}-dimensional vector of (deterministically computed) state
#'   transitions (conditional means)
f <- function(x_tt, regs, phi_x, bet_reg) {
  # xt <- phi_x*xtt
  x_t <- phi_x*x_tt + regs %*% bet_reg
  # x_t <- phi_x*x_tt + z %*% bet_z + u %*% bet_u

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
