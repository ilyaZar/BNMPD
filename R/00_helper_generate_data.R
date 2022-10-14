#' Generates output container for true parameter values
#'
#' Arguments specify true parameter values, dimension of the model and which
#' regressor types to use. For random effects (beta_u-type regressors) a seed
#' makes sure to produce the same setup. Output container can be used when
#' plotting diagnostics (in a simulation study) as well as generating simulated
#' data sets based on the true parameter values.
#'
#' @param dim_model a vector with three components: \code{NN x TT x DD}
#' @param sig_sq a vector of length \code{DD} and strictly positive elements
#' @param phi a vector of length \code{DD} with elements between 0 and 1
#' @param bet_z a list of length \code{DD} with (possibly) varying number of
#'    regressors
#' @param SIMUL_Z_BETA logical; if \code{TRUE}, then standard continuous
#'   (z-type) covariates are included
#' @param SIMUL_U_BETA logical; if \code{TRUE}, then random effects (u-type
#'    regressors) are simulated based on the argument \code{seed_taken}
#' @param NUM_BETA_U integer giving the number of random effects per cross
#'   section to simulate
#' @param seed_taken the seed used drawing random effects
#'
#' @return a list of four elements each giving the corresponding true parameters
#'    container for the model size/dimension
#' @export
generate_true_params <- function(dim_model,
                                 sig_sq,
                                 phi,
                                 bet_z,
                                 SIMUL_Z_BETA,
                                 SIMUL_U_BETA,
                                 NUM_BETA_U,
                                 seed_taken) {
  check_args <- (missing(sig_sq) || missing(phi) || missing(bet_z) ||
                   missing(SIMUL_U_BETA) || missing(SIMUL_Z_BETA) ||
                   missing(seed_taken) || missing(dim_model))
  if (check_args) stop("Missing arguments")
  if (SIMUL_U_BETA && missing(NUM_BETA_U)) stop("Number of REs not specified!")
  # 1. Data settings: -------------------------------------------------------
  NN <- dim_model[1]   # Cross sectional length
  cat(crayon::green("Setting dimension "), crayon::yellow("NN"),
      crayon::green("to "), crayon::red(NN), crayon::green("!\n"))
  TT <- dim_model[2]  # Time series length
  cat(crayon::green("Setting dimension "), crayon::yellow("TT"),
      crayon::green("to "), crayon::red(TT), crayon::green("!\n"))
  DD <- dim_model[3]
  cat(crayon::green("Setting dimension "), crayon::yellow("DD"),
      crayon::green("to "), crayon::red(DD), crayon::green("!\n"))
  # 2. Set up parameter values: ---------------------------------------------
  check_sig_sq <- all(sig_sq > 0) && all(length(sig_sq) == DD)
  if (!check_sig_sq) stop("'sig_sq' >0 and a vector of length DD ...\n")
  true_sig_sq <- matrix(sig_sq, nrow = DD, ncol = NN)
  #
  check_phi <- all(phi > 0) && all(phi < 1) && all(length(phi) == DD)
  if (!check_phi) stop("phi > 0 or phi < 1 or not a vector of length DD ...\n")
  true_phi <- matrix(phi, nrow = DD, ncol = NN)
  #
  if (SIMUL_Z_BETA) {
    check_bet_z <- is.list(bet_z) && length(bet_z) == DD
    if (!check_bet_z) stop("'bet_z' not a list or not of length = DD ...\n")
    true_bet_z <- bet_z
  } else {
    true_bet_z <- NULL
  }
  #
  if (SIMUL_U_BETA) {
    true_out_u <- BNMPD::generate_bet_u(DD, NN,
                                        TRUE, rep(NUM_BETA_U,
                                                  times = DD),
                                        seed_no = seed_taken)
    true_bet_u <- true_out_u[[1]]
    true_D0u_u <- true_out_u[[2]]
  } else {
    true_bet_u <- NULL
    true_D0u_u <- NULL
  }
  par_trues <- list(sig_sq = true_sig_sq,
                    phi = true_phi,
                    bet_z = true_bet_z,
                    bet_u = true_bet_u,
                    vcm_u = true_D0u_u)
  return(true_params = par_trues)
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
#' @return a named list of two elements:
#'   \itemize{
#'      \item{\code{true_bet_u}:}{a list of length \code{DD} of matrices of
#'                                dimension \code{num_re x NN} with random
#'                                effects}
#'      \item{\code{true_vcm_u}:}{a list of length \code{DD} of diagonal
#'                                matrices of dimension \code{num_re x num_re}
#'                                which are the covariance matrices of random
#'                                effects}
#'   }
#'
#' @export
generate_bet_u <- function(DD, NN,
                           from_IW = FALSE,
                           num_re = rep(2, times = DD),
                           seed_no = 42) {
  true_bet_u <- vector("list", DD)
  if (from_IW) {
    stopifnot(is.numeric(num_re) && (length(num_re) == DD))
    if (!is.null(seed_no))  set.seed(seed_no)
    n0u <- num_re + 1
    D0u <- vector("list", DD)
    for (d in 1:DD) {
      # D0u[[d]] <- solve(diag(1:num_re[d]*(10 + 1/num_re[d]),
      #  nrow = num_re[d]))
      D0u[[d]] <- solve(diag(1:num_re[d]/10*(0.01 + 1/num_re[d]),
                             nrow = num_re[d]))
      D0u[[d]] <- solve((stats::rWishart(1, n0u[d], D0u[[d]]))[, , 1])
      # D0u[[d]] <- diag(1, num_re[d])
      true_bet_u[[d]] <- matrix(0, nrow = num_re[d], ncol = NN)
      # browser()
      for (n in 1:NN) {
        true_bet_u[[d]][, n] <- MASS::mvrnorm(n = 1,
                                              mu = rep(0, times = num_re[d]),
                                              Sigma = D0u[[d]])
      }
    }
    return(list(true_bet_u = true_bet_u,
                true_vcm_u = D0u))
  } else {
    stopifnot(is.numeric(num_re) && (length(num_re) == 1))
    vals <- matrix(1:(num_re*DD)*c(-1, 1), nrow = num_re*DD, ncol = NN)
    for (d in 1:DD) {
      true_bet_u[[d]] <- vals[1:num_re + num_re*(d - 1), , drop = FALSE]
    }
    return(true_bet_u = true_bet_u)
  }
}
#' Deduces from vector of parameter names, which type of modelling to employ
#'
#' @param get_modelling_reg_types a named list of parameter names; the list
#'   of true parameter values is taken usually taken but empty named list
#'   elements are also possible
#'
#' @return a named, logical vector of dimension 4 given \code{TRUE} or
#'   \code{FALSE} if (in this order) modeling of z-regressors, u-regressors, z
#'   spline regressors or u spline regressors should be performed (with names
#'   of return vector set to these variants)
#'
get_modelling_reg_types <- function(get_modelling_reg_types) {
  correct_names <- c("bet_z", "bet_u", "bet_z_spl", "bet_u_spl")
  if (length(intersect(get_modelling_reg_types, correct_names)) == 0) {
    stop(paste0("The 'par_true' argument must have correct names: choose from",
                "'bet_z', 'bet_u', 'bet_z_spl' or 'bet_u_spl'! "))
  }
  out <- vector("logical", 4)
  out[1] <- correct_names[1] %in% get_modelling_reg_types
  out[2] <- correct_names[2] %in% get_modelling_reg_types
  out[3] <- correct_names[3] %in% get_modelling_reg_types
  out[4] <- correct_names[5] %in% get_modelling_reg_types

  names(out) <- c("z-regressors",
                  "u-regressors",
                  "z-spline-regressors",
                  "u-spline-regressors")
  return(out)
}
#' Generate Dirichlet target levels for simulation study
#'
#' A sequence of fractions should not be too erratic over time and the different
#' fractions not too far away from each other. This function achieves the latter
#' by computing target fractions levels that are neither too small, nor too
#' large so that fractions are not too far apart.
#'
#' The tuning parameter list, as given by default values of the argument, works
#' nicely for DD = 3, see the examples section below. The resulting values - 30,
#' 30, 40 - are reasonable target parameters for the Dirichlet taken as
#' `alpha_1=30`, `alpha_2=30`, `and alpha_3=40`.
#'
#'
#' @param DD integer giving the multivariate dimension
#' @param NN number of cross sectional units (repetitions of target values)
#' @param tuning_parameters a set of tuning parameters that generate a
#'   reasonably spaced sequence of target values
#'
#' @return a matrix of dimension \code{NN x DD}, with each row being the target
#'   levels (currently all the same)
#' @export
#'
#' @examples
#' get_dirichlet_levels(DD = 3, NN = 4)
get_dirichlet_levels <- function(DD, NN,
                                 tuning_parameters = list(seq_start = 0.3,
                                                          seq_step = 0.1,
                                                          seq_rep = 2,
                                                          seq_scale = 100)) {
  seq_start <- tuning_parameters$seq_start
  seq_step  <- tuning_parameters$seq_step
  seq_rep   <- tuning_parameters$seq_rep

  tuned_vec <- rep(seq(from = seq_start,
                       to = seq_start + seq_step * DD,
                       by = seq_step),
                   each = seq_rep)
  tuned_vec <- tuned_vec[1:DD]
  tuned_vec <- tuned_vec/sum(tuned_vec)

  matrix(tuned_vec * 100,
         nrow = DD,
         ncol = NN)

}
#' Generates panel data for various models.
#'
#' The data is a panel data from a Dirichlet, Multinomial or
#' Dirichlet-multinomial distribution. The mulitivariate draws can vary along
#' the time and cross-sectional dimensions. This is because the parameter of the
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
#' @param include_intercept logical vector of dimension \code{DD}; if
#'   \code{TRUE} include an intercept at the cross sectional unit for component
#'   \code{d}
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
#'   per cross sectional unit \code{n=1,...,N} with a joint plot of all
#'   components together
#' @param plot_states_each_d logical; if \code{TRUE}, latent states are
#'   plotted per cross sectional unit \code{n=1,...,N} with a separate plot for
#'   each component
#'
#' @return NN-dimensional list of lists of two: \code{[[1]]} -> regressors and
#'   \code{[[2]]} -> latent states
#' @export
generate_data_t_n <- function(distribution,
                              NN, TT, DD,
                              par_true,
                              x_levels,
                              x_log_scale,
                              include_intercept = NULL,
                              include_policy,
                              include_zeros = NULL,
                              plot_measurements = FALSE,
                              plot_states = FALSE,
                              plot_states_each_d = FALSE,
                              seed_no = NULL) {
  if (is.null(include_intercept)) {
    include_intercept <- rep(FALSE, times = DD)
  }
  if (is.null(include_policy)) {
    include_policy <- matrix(FALSE, nrow = DD, ncol = NN)
  # policy_modelling    <- cbind(c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
  #                              c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
  #                              c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
  #                              c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
  #                              c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
  #                              c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE))
  }
  if (is.null(include_zeros)) {
    include_zeros <- NULL
  }
  # zero_modelling      <- c(1:4, 1, 2) # rep(4, times = DD) # c(1, 2, 3, 4, 1, 2)
  densitities_supported <- c("multinomial", "dirichlet-mult",
                             "gen-dirichlet-mult", "gen-dirichlet",
                             "dirichlet")
  if (!(distribution %in% densitities_supported)) {
    stop(paste0("Argument to distribution must be one of: ",
                paste0(densitities_supported, collapse = ", "), "!"))
  }
  if (x_log_scale) {
    x_levels <- log(x_levels)
  }
  if (!is.null(seed_no)) set.seed(seed_no)

  modelling_reg_types <- names(par_true)
  modelling_reg_types <- get_modelling_reg_types(modelling_reg_types)

  data_part1 <- array(0, c(TT, DD, NN))
  if (distribution == "multinomial" || distribution == "dirichlet-mult") {
    data_part2 <- matrix(0, nrow = TT, ncol = NN)
  }

  states <- array(0, c(TT, DD, NN))

  if (modelling_reg_types[1]) {
    dim_bet_z <- sum(sapply(par_true[["bet_z"]], length))
    z <- array(0, c(TT, dim_bet_z, NN))
  }
  if (modelling_reg_types[2]) {
    if (NN == 1) {
      dim_bet_u <- sum(sapply(par_true[["bet_u"]], length))
    } else {
      dim_bet_u <- sum(sapply(par_true[["bet_u"]], nrow))
    }
    u <- array(0, c(TT, dim_bet_u, NN))
  }
  for (n in 1:NN) {
    # browser()
    par_true_current <- list(sig_sq = par_true[["sig_sq"]][, n],
                             phi = par_true[["phi"]][, n])
    if (modelling_reg_types[1]) par_true_current$bet_z <- par_true[["bet_z"]]
    if (modelling_reg_types[2]) {
      # par_true[["bet_u"]]
      par_true_current$bet_u <- lapply(par_true[["bet_u"]], `[`, i = , j = n)
    }

    out_data_tmp <- generate_data_t(distribution = distribution,
                                    TT = TT, DD = DD,
                                    par_true = par_true_current,
                                    x_levels = x_levels[, n],
                                    x_log_scale = x_log_scale,
                                    include_intercept = include_intercept,
                                    include_policy = include_policy[, n,
                                                                   drop = TRUE],
                                    include_zeros = include_zeros,
                                    modelling_reg_types = modelling_reg_types)
    data_part1[, , n] <- out_data_tmp[[1]]$yraw
    if (distribution == "multinomial" || distribution == "dirichlet-mult") {
      data_part2[, n]   <- out_data_tmp[[1]]$num_counts
    }
    if (modelling_reg_types[1]) {
      z[, , n] <- out_data_tmp$z

    }
    if (modelling_reg_types[2]) {
      u[, , n] <- out_data_tmp$u
    }
    states[, , n] <- out_data_tmp$x
    # browser()
    if (isTRUE(any(plot_measurements, plot_states, plot_states_each_d))) {
      plot_data_per_n(DD, yraw = data_part1[, , n], x = states[, , n],
                      x_log_scale = x_log_scale,
                      x_levels = x_levels[, n],
                      plot_measurements = plot_measurements,
                      plot_states       = plot_states,
                      plot_states_each_d  = plot_states_each_d,
                      cs = n)
    }
  }
  out_data <- vector("list", 3)
  names(out_data) <- c("data", "regs", "states")
  if (distribution == "dirichlet") {
    out_data[[1]] <- list(yraw = data_part1)
  }
  if (distribution == "multinomial" || distribution == "dirichlet-mult") {
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
#' Dirichlet, Multinomial or Dirichlet-multinomial distribution. The
#' multivariate draws can vary along the time dimension. This is because the
#' parameter of the distributions that generate the draws are modeled as a
#' function of regressors and latent states where the regressors and latent
#' states can vary over time.
#'
#' @param distribution specifies the distribution; "Dirichlet", "multinomial" or
#'   "dirichlet-multinomial"
#' @param TT number of time periods
#' @param DD number of shares/fractions (for Dirichlet or Dirichlet-multinomial)
#'   or the number of categories for a multinomial distribution
# @param n current cross sectional unit i.e. an integer \code{n=1,...,NN}
#' @param par_true list of true parameters that describe the latent state
#'   process
#' @param x_levels vector of target "means"/"levels" of the states around which
#'   they fluctuate
#' @param x_log_scale logical; if \code{TRUE}, then \code{x_levels[i]} is taken
#'   as log and the random number generation for the states is performed in logs
#' @param include_intercept logical vector of dimension \code{DD}; if
#'   \code{TRUE} include an intercept at the cross sectional unit for component
#'   \code{d}
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
generate_data_t <- function(distribution,
                            TT, DD,
                            par_true,
                            x_levels,
                            x_log_scale,
                            include_intercept,
                            include_policy,
                            include_zeros,
                            modelling_reg_types) {
  # browser()
  x <- matrix(nrow = TT, ncol = DD, 0)
  sig_sq_x <- par_true[["sig_sq"]]
  phi_x    <- par_true[["phi"]]
  if (modelling_reg_types[1]) {
    z <- list()
    bet_z <- par_true[["bet_z"]]
  } else {
    bet_z <- NULL
  }
  if (modelling_reg_types[2]) {
    u <- list()
    bet_u <- par_true[["bet_u"]]
  } else {
    bet_u <- NULL
  }

  if (distribution %in% c("multinomial", "dirichlet-mult",
                          "gen-dirichlet-mul")) {
    num_counts <- sample(x = 80000:120000, size = TT)
  }

  # reg_sd_levels <- c(0.0125, 0.1, 0.025, 0.1, 0.1, 0.1)
  # reg_sd_levels <- rep(0.5, times = DD)
  reg_sd_levels <- rep(0.05, times = DD)
  bet_sd_level  <- 1

  for (d in 1:DD) {
    res <- generate_x_z_u(TT = TT,
                          phi_x = phi_x[d],
                          sig_sq_x = sig_sq_x[d],
                          bet_z = bet_z[[d]],
                          bet_u = bet_u[[d]],
                          modelling_reg_types = modelling_reg_types,
                          x_level = x_levels[d],
                          reg_sd = reg_sd_levels[d],
                          bet_sd = bet_sd_level,
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
    # browser()
    yraw <- my_rdirichlet(alpha = alphas)
    if (sum(rowSums(yraw)) != TT) {
      stop("Something is wrong with Dirichelet: fractions don't sum up to 1!")
    }
    out_data_all$data <- list(yraw = yraw)
  }
  if (distribution == "multinomial") {
    alphas <- alphas / rowSums(alphas)
    yraw <- my_rmultinomial(probs = alphas, num_counts = num_counts)
    out_data_all$data <- list(yraw = yraw, num_counts = num_counts)
  }
  if (distribution == "dirichlet-mult") {
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
#' @param x_level target levels i.e. stationary mean/level of latent states (
#'   ensures that states at each multivariate component of the response
#'   fluctuates around that particular level too)
#' @param reg_sd standard deviations of the regressor values (tunable)
#' @param bet_sd standard deviations of random coefficient value draws (tunable)
#' @param modelling_reg_types named logical vector of dimension 4 where each
#'   component indicates, if \code{TRUE}, that the corresponding regressor type
#'   is to be generated in the following order:
#'   \itemize{
#'     \item{modelling_reg_types, component 1: }{includes z regressors}
#'     \item{modelling_reg_types, component 2: }{includes u regressors}
#'     \item{modelling_reg_types, component 3: }{includes z spline regressors}
#'     \item{modelling_reg_types, component 4: }{includes u spline regressors}
#'   }
#' @param x_log_scale logical; if \code{TRUE} process is simulated on the
#'   log-scale
#' @param intercept logical; if \code{TRUE} includes an intercept term i.e. a
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
                           x_level,
                           reg_sd,
                           bet_sd,
                           modelling_reg_types,
                           x_log_scale,
                           intercept,
                           policy_dummy = FALSE,
                           zero_pattern = NULL,
                           drift = FALSE) {
  dim_z <- length(bet_z)
  dim_u <- length(bet_u)
  dim_reg <- dim_z + dim_u
  bet_reg <- c(bet_z, bet_u)
  x <- rep(0, TT)
  # browser()
  # BEGINNING OF REGRESSOR SIMULATION: --------------------------------------
  if (modelling_reg_types[1] && !modelling_reg_types[2]) {
    z <- generate_reg_vals(TT = TT,
                           bet_reg = bet_reg,
                           dim_reg = dim_reg,
                           phi_x = phi_x,
                           x_level = x_level,
                           reg_sd = reg_sd,
                           bet_sd = bet_sd,
                           intercept = intercept,
                           policy_dummy = policy_dummy,
                           zero_pattern = NULL)
    u <- NULL
  }
  if (!modelling_reg_types[1] && modelling_reg_types[2]) {
    u <- generate_reg_vals(TT = TT,
                           bet_reg = bet_reg,
                           dim_reg = dim_reg,
                           phi_x = phi_x,
                           x_level = x_level,
                           reg_sd = reg_sd,
                           bet_sd = bet_sd,
                           intercept = intercept,
                           policy_dummy = policy_dummy,
                           zero_pattern = NULL)
    z <- NULL
  }
  if (modelling_reg_types[1] && modelling_reg_types[2]) {
    # browser()
    x_level_split <- x_level*c(0.7, 0.3)
    z <- generate_reg_vals(TT = TT,
                           bet_reg = bet_z,
                           dim_reg = dim_z,
                           phi_x = phi_x,
                           x_level = x_level_split[1],
                           reg_sd = reg_sd,
                           bet_sd = bet_sd,
                           intercept = FALSE,
                           policy_dummy = FALSE,
                           zero_pattern = NULL)
    u <- generate_reg_vals(TT = TT,
                           bet_reg = bet_u,
                           dim_reg = dim_u,
                           phi_x = phi_x,
                           x_level = x_level_split[2],
                           reg_sd = reg_sd,
                           bet_sd = bet_sd,
                           intercept = intercept,
                           policy_dummy = policy_dummy,
                           zero_pattern = NULL)
  }
  reg_all <- cbind(z, u)
  # END OF REGRESSOR SIMULATION: --------------------------------------------
  # browser()
  xinit <- x_level
  x[1] <- f(x_tt = xinit, regs = reg_all[1, ], phi_x = phi_x, bet_reg = bet_reg)
  x[1] <- x[1] + sqrt(sig_sq_x)*stats::rnorm(n = 1)

  # browser()
  for (t in 1:TT) {
    if (t < TT) {
      x[t + 1] <- f(x_tt = x[t], regs = reg_all[t + 1, ],
                    phi_x = phi_x, bet_reg = bet_reg)
      x[t + 1] <- x[t + 1] + sqrt(sig_sq_x)*stats::rnorm(n = 1)
    }
  }

  if (sum(any(x <= 0)) & x_log_scale == FALSE) {
    stop("some state process (x1_t, x2_t, ... or xD_t) not positive!")
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
#' @param reg_sd standard deviations for the regressor values
#' @param reg_sd standard deviations for the regressor values
#' @param bet_sd standard deviations for random coefficient value draws
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
                              x_level, reg_sd, bet_sd,
                              intercept,
                              policy_dummy,
                              zero_pattern = NULL) {
  # browser()
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
  }
  if (dim_reg == 1) {
    if (intercept) {
      const_level <- x_level * (1 - phi_x)/bet_reg
      regs <- matrix(const_level, nrow = TT, ncol = dim_reg)
    } else if (policy_dummy) {
      regs <- dummy_to_use
    } else {
      const_mean <- x_level * (1 - phi_x)/bet_reg
      regs          <- matrix(stats::rnorm(TT*dim_reg,
                                           mean = const_mean,
                                           sd = reg_sd),
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
    if (!intercept && !policy_dummy) {
      reg_means <- stats::rnorm(dim_reg - 1, mean = 0, sd = bet_sd)
      last_reg_mean <- x_level * (1 - phi_x)
      last_reg_mean <- last_reg_mean - sum(reg_means * bet_reg[-dim_reg])
      last_reg_mean <- last_reg_mean/bet_reg[dim_reg]
      reg_means     <- c(reg_means, last_reg_mean)
      regs          <- matrix(stats::rnorm(TT*dim_reg,
                                           mean = reg_means,
                                           sd = reg_sd),
                              nrow = TT,
                              ncol = dim_reg,
                              byrow = TRUE)
      return(regs)
    }
    if (intercept && !policy_dummy) {
      reg_means <- 1
    } else if (!intercept && policy_dummy) {
      reg_means <- dummy_to_use[1]
    }
    last_reg_mean <- x_level * (1 - phi_x) - sum(reg_means * bet_reg[-dim_reg])
    last_reg_mean <- last_reg_mean/bet_reg[dim_reg]
    reg_means     <- c(reg_means, last_reg_mean)
    regs          <- matrix(stats::rnorm(TT*dim_reg,
                                         mean = reg_means,
                                         sd = reg_sd),
                            nrow = TT,
                            ncol = dim_reg,
                            byrow = TRUE)
    if (intercept && !policy_dummy) {
      regs[, 1] <- 1
      return(regs)
    } else if (!intercept && policy_dummy) {
      regs[, 1] <- dummy_to_use[1]
      return(regs)
    }
  }
  if (dim_reg > 2) {
    if (intercept && !policy_dummy) {
      reg_means <- stats::rnorm(dim_reg - 2, mean = 0, sd = bet_sd)
      reg_means <- c(1, reg_means)

      last_reg_mean <- x_level * (1 - phi_x)
      last_reg_mean <- last_reg_mean - sum(reg_means * bet_reg[-dim_reg])
      last_reg_mean <- last_reg_mean/bet_reg[dim_reg]
      reg_means     <- c(reg_means, last_reg_mean)
      reg_len       <- length(reg_means)
      regs          <- matrix(stats::rnorm(TT*reg_len,
                                           mean = reg_means,
                                           sd = reg_sd),
                              nrow = TT,
                              ncol = reg_len,
                              byrow = TRUE)
      regs[, 1] <- 1
    } else if (!intercept && policy_dummy ) {
      reg_means <- stats::rnorm(dim_reg - 2, mean = 0, sd = bet_sd)
      last_reg_mean <- x_level * (1 - phi_x)
      last_reg_mean <- last_reg_mean - sum(reg_means * bet_reg[-c(1, dim_reg)])
      last_reg_mean <- last_reg_mean/bet_reg[dim_reg]
      reg_means     <- c(reg_means, last_reg_mean)
      reg_len       <- length(reg_means)
      regs          <- matrix(stats::rnorm(TT*reg_len,
                                           mean = reg_means,
                                           sd = reg_sd),
                              nrow = TT,
                              ncol = reg_len,
                              byrow = TRUE)
      regs <- cbind(dummy_to_use, regs)
    } else if (policy_dummy && intercept) {
      reg_means <- stats::rnorm(dim_reg - 3, mean = 0, sd = bet_sd)
      reg_means <- c(1, reg_means)

      last_reg_mean <- x_level * (1 - phi_x)
      last_reg_mean <- last_reg_mean- sum(reg_means * bet_reg[-c(2, dim_reg)])
      last_reg_mean <- last_reg_mean/bet_reg[dim_reg]
      reg_means     <- c(reg_means, last_reg_mean)
      reg_len       <- length(reg_means)
      regs          <- matrix(stats::rnorm(TT*reg_len,
                                           mean = reg_means,
                                           sd = reg_sd),
                              nrow = TT,
                              ncol = reg_len,
                              byrow = TRUE)
      regs <- cbind(1, dummy_to_use, regs[, -c(1)])
    } else {
      reg_means <- stats::rnorm(dim_reg - 1, mean = 0, sd = bet_sd)
      last_reg_mean <- x_level * (1 - phi_x)
      last_reg_mean <- last_reg_mean- sum(reg_means * bet_reg[-dim_reg])
      last_reg_mean <- last_reg_mean/bet_reg[dim_reg]
      reg_means     <- c(reg_means, last_reg_mean)
      reg_len       <- length(reg_means)
      regs          <- matrix(stats::rnorm(TT*reg_len,
                                           mean = reg_means,
                                           sd = reg_sd),
                              nrow = TT,
                              ncol = reg_len,
                              byrow = TRUE)
    }
    return(regs)
  }
}
#' Function to plot the data
#'
#' @param DD number of shares/fractions (for dirichlet or dirichlet-multinomial)
#'   or the number of categories for a multinomial distribution
#' @param yraw a numeric matrix of dimension \code{TTxDD} giving the simulated
#'   values (not the num_counts if a compoung distributions i.e. just the
#'   dirichlet shares e.g.)
#' @param x a numeric matrix of dimension \code{TTxDD}; simulated latent states
#' @param x_log_scale logical; if \code{TRUE}, simulation of latent state
#'   process is was perfomred on the log-scale
#' @param x_levels a numeric vector of length \code{DD} giving the target
#' dirichlet levels i.e. latent state levels (stationary mean of the latent
#' state process etc.)
#' @param plot_measurements logical; if \code{TRUE}, measurements are plotted
#'   per cross sectional unit \code{n=1,...,N}
#' @param plot_states logical; if \code{TRUE}, latent states are plotted
#'   per cross sectional unit \code{n=1,...,N} with a joint plot of all
#'   components together
#' @param plot_states_each_d logical; if \code{TRUE}, latent states are
#'   plotted per cross sectional unit \code{n=1,...,N} with a seperate plot for
#'   each component
#' @param cs integer giving the cross sectional observation (1,2,...NN)
#'
#' @return invisible return; pure side effects function that plots the data
plot_data_per_n <- function(DD,
                            yraw, x,
                            x_log_scale,
                            x_levels,
                            plot_measurements,
                            plot_states,
                            plot_states_each_d,
                            cs) {
  if (plot_measurements && plot_states)  graphics::par(mfrow = c(2, 1))
  if (plot_measurements && !plot_states) graphics::par(mfrow = c(2, 1))
  if (!plot_measurements && plot_states) graphics::par(mfrow = c(2, 1))
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
                      main = paste0(names_title, " (levels; N = ", cs, ")"),
                      ylab = names_ylab,
                      xlab = names_xlab)
    graphics::matplot(all_measurms/rowSums(all_measurms),
                      type = "l",
                      lty = 1,
                      lwd = 1,
                      main = paste0(names_title, " (fractions; N = ", cs, ")"),
                      ylab = names_ylab,
                      xlab = names_xlab)
  }
  if (plot_states || plot_states_each_d) {
    names_title <- "True States"
    names_ylab  <- "states: xt's"
    names_xlab  <- paste0("x1_t (black),", "x2_t (red),",
                          "x3_t (green),", "x4_t (blue)",
                          "x5_t (turq.)", "and", "x6_t (pink)")
    all_states <- x
    if (x_log_scale) {
      all_states <- log(all_states)
    }
  }
  if (plot_states) {
    graphics::matplot(all_states,
                      type = "l",
                      lty = 2,
                      lwd = 1,
                      main = paste0(names_title, " (levels; N = ", cs, ")"),
                      ylab = names_ylab,
                      xlab = paste(names_xlab))
    graphics::matplot(all_states/rowSums(all_states),
                      type = "l",
                      lty = 2,
                      lwd = 1,
                      main = paste0(names_title, " (fractions; N = ", cs, ")"),
                      ylab = names_ylab,
                      xlab = names_xlab)
  }
  if (plot_states_each_d) {
    graphics::par(mfrow = c(ceiling(DD/2), 2))
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

      # graphics::plot((all_states/rowSums(all_states))[, d],
      #                type = "l",
      #                lty = 2,
      #                lwd = 1,
      #                main = names_title,
      #                ylab = names_ylab,
      #                xlab = names_xlab,
      #                col = d)
      # graphics::abline(h = (x_levels/sum(x_levels))[d], lty = 1, lwd = 3)
      # graphics::abline(h = mean((all_states/rowSums(all_states))[, d,
      #                                                           drop = TRUE]),
      #                  lty = 1, lwd = 1, col = d)
    }
  }
  graphics::par(mfrow = c(1, 1))
  return(invisible(DD))
}
#' Generates random samples from dirichlet distribution
#'
#' Generates random samples from dirichlet distribution; the dimension of the
#' Dirichlet distribution (i.e. the number of shares or fractions) is taken as
#' the number of columns in the alpha-matrix; the number of samples are taken as
#' the number of rows of the alpha matrix; hence, each row of \code{alpha}
#' corresponds to \code{D} alpha parameters for which a \code{D}-dimensional
#' random draw from the Dirichlet is generated.
#'
#' @param alpha a matrix of alpha parameters of a Dirichlet distribution
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
#
#
#
#' Save simulated data and true parameter values used to generate it.
#'
#' @param pth_to_write path to output directory
#' @param fn_data_set file name (\code{.csv}-ending required) for simulated data
#' @param fn_true_states file name for R-containter object that stores true
#'   states
#' @param fn_zero_states file name for R-containter object that stores true
#'   states
#' @param data_sim simulated data set as produced via [generate_data_t_n()]
#' @param true_params true parameter values used to generate data (to determine
#'   containter/data sizes) as produced by output from [generate_true_params()]
#' @param dim_model a vector with three components: \code{NN x TT x DD}
#'
#' @return side effect function; saves simulated data and true params to output
#'   path
#'
#' @export
save_simulated_data <- function(pth_to_write,
                                fn_data_set,
                                fn_true_states,
                                fn_zero_states,
                                data_sim,
                                dim_model,
                                true_params) {
  SIMUL_U_BETA <- !is.null(data_sim$regs$u)
  NN <- dim_model[1]   # Cross sectional length
  cat(crayon::green("Setting dimension "), crayon::yellow("NN"),
      crayon::green("to "), crayon::red(NN), crayon::green("!"))
  TT <- dim_model[2]  # Time series length
  cat(crayon::green("Setting dimension "), crayon::yellow("TT"),
      crayon::green("to "), crayon::red(TT), crayon::green("!"))
  DD <- dim_model[3]
  cat(crayon::green("Setting dimension "), crayon::yellow("DD"),
      crayon::green("to "), crayon::red(DD), crayon::green("!"))

  seq_regs_z <- lapply(true_params$bet_z, length)
  num_regs_z <- sum(unlist(seq_regs_z))
  seq_regs_u <- lapply(true_params$bet_u, nrow)
  num_regs_u <- sum(unlist(seq_regs_u))

  ncol_out <- DD + num_regs_z + num_regs_u
  data_out <- matrix(0, nrow = TT * NN, ncol = ncol_out)
  data_out <- as.data.frame(data_out)

  names_z_reg <- paste0(paste0(paste0("Z_", sapply(lapply(seq_regs_z,
                                                          function(x) {
                                                            seq_len(x)}),
                                                   as.character)), "_"),
                        rep(1:DD, unlist(seq_regs_z)))
  if (SIMUL_U_BETA) {
    names_u_reg <- paste0(paste0(paste0("U_", sapply(lapply(seq_regs_u,
                                                            function(x) {
                                                              seq_len(x)}),
                                                     as.character)), "_"),
                          rep(1:DD, unlist(seq_regs_u)))
  } else {
    names_u_reg <- NULL
  }
  names(data_out) <- c(paste0("Y", 1:DD), names_z_reg, names_u_reg)
  vals_cs  <- as.character(paste0("cs_", rep(seq_len(NN), each = TT)))
  vals_ts  <- rep(1:TT, times = NN)
  cs_ts    <-tibble::tibble(CS = vals_cs, TS = vals_ts)
  data_out <- dplyr::bind_cols(cs_ts, data_out)
  for(n in 1:NN) {
    id_rows <- TT*(n - 1) + (1:TT)
    offset_col <- 2
    id_col_y <- 1:DD + offset_col
    id_col_z <- DD + 1:(num_regs_z) + offset_col
    id_col_u <- DD + num_regs_z + 1:(num_regs_u) + offset_col
    data_out[id_rows, id_col_y] <- data_sim$data$yraw[, , n]
    data_out[id_rows, id_col_z] <- data_sim$regs$z[, , n]
    if (SIMUL_U_BETA) {
      data_out[id_rows, id_col_u] <- data_sim$regs$u[, , n]
    }
  }
  write.csv(data_out, file = file.path(pth_to_write, fn_data_set))
  true_states <- data_sim$states
  save(true_states, file = file.path(pth_to_write, fn_true_states))
  zero_states <- true_states
  zero_states[, , ] <- 0
  save(zero_states,  file = file.path(pth_to_write, fn_zero_states))
}
#' Generates consistent file names for simulation study
#'
#' @param fn_data character giving the filename of the simulated data set
#'   (saved in \code{.csv}-format)
#' @param fn_true_states character giving the filename of the simulated true
#'   states (saved in \code{.RData}-format)
#' @param fn_zero_states  character giving the filename of the states set all to
#'   zero (saved in \code{.RData}-format)
#' @param dim_model numeric vector of 3 elements: \code{NN x TT x DD}
#' @param SIMUL_Z_BETA logical; if \code{TRUE} Z-type regressors are simulated
#'   which is reflected in the naming of simulated data sets
#' @param SIMUL_U_BETA logical; if \code{TRUE} U-type regressors (i.e. random
#'   effects) are simulated which is reflected in the naming of simulated data
#'   sets
#' @param num_z_regs numeric vector giving number of Z-type regressors; parsed
#'   as comma-seperated sequence of numbers i.e. for \code{num_z_regs = 1:3} we
#'   get "withLIN,1,2,3" in the file names
#' @param num_u_regs numeric vector giving number of U-type regressors; parsed
#'   as comma-seperated sequence of numbers i.e. for \code{num_u_regs = 1:3} we
#'   get "withRE,1,2,3" in the file names
#'
#' @return a list of 3:
#'    \itemize{
#'    \item\code{fn_data_set}: file name of the data set
#'    \item\code{fn_data_set}: file name for true state values
#'    \item\code{fn_data_set}: file name for zero state values
#'    }
#' @export
get_file_names_simul_data <- function(fn_data,
                                      fn_true_states,
                                      fn_zero_states,
                                      dim_model,
                                      SIMUL_Z_BETA,
                                      SIMUL_U_BETA,
                                      num_z_regs,
                                      num_u_regs) {
  tmp_fn <- paste0("NN", dim_model[1],
                   "_TT", dim_model[2],
                   "_DD", dim_model[3])
  if (SIMUL_Z_BETA) {
    tmp_fn <- paste0(tmp_fn, "_", "withLIN", paste0(",", num_z_regs,
                                                       collapse = ""))
  } else {
    tmp_fn <- paste0(tmp_fn, "_", "noLIN")
  }
  if (SIMUL_U_BETA) {
    tmp_fn <- paste0(tmp_fn, "_", "withRE", paste0(",", num_u_regs,
                                                      collapse = ""))
  } else {
    tmp_fn <- paste0(tmp_fn, "_", "noRE")
  }
  fn_data_set <- paste0(fn_data, "_", tmp_fn, ".csv")
  fn_true_val <- paste0(fn_true_states, "_", tmp_fn, ".RData")
  fn_zero_val <- paste0(fn_zero_states, "_", tmp_fn, ".RData")
  return(list(fn_data_set = fn_data_set,
              fn_true_val = fn_true_val,
              fn_zero_val = fn_zero_val))
}
