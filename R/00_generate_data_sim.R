#' Generates multivariate (non-linear) panel data for various models.
#'
#' The measurements (responses, dependent variables) have a Dirichlet,
#' Multinomial or Dirichlet-Multinomial distribution. The multivariate draws
#' can vary along the time and cross-sectional dimensions. This is because the
#' parameters of the distributions that generate the draws are modelled as a
#' function of regressors and latent states where the regressors and latent
#' states can vary over time and cross section.
#'
#' @param true_params an object of class "\code{trueParams}" which is a list of
#'   true parameter values and meta information such as distribution,
#'   model dimension, the seed under which the true parameters are generated,
#'   the logical indicators that describe which parameters to generate and their
#'   lengths, as provided by the function [new_trueParams()]
#' @param distribution specifies the distribution: "dirichlet", "gen_dirichlet",
#'    multinomial", "dirichlet_mult", "gen_dirichlet_mult", or "normal" (the
#'    latter generates the latent states without link-function and measurement/
#'    response transformations, which is useful e.g. when testing the pure Gibbs
#'    sampler)
#' @param x_levels target "mean" levels of the states around which they
#'   fluctuate
#' @param X_LOG_SCALE logical; if \code{TRUE}, \code{x_levels} are taken as logs
#'   and the random number generation for the states is performed in logs
#' @param options_include a list of options for various effects:
#'   \itemize{
#'      \item{\code{intercept: }}{logical vector of dimension \code{DD};
#'      if \code{TRUE} include an intercept at the cross sectional unit for
#'      component \code{d}}
#'      \item{\code{policy: }}{logical vector of dimension \code{DD}; if
#'      \code{TRUE} include a policy dummy at the cross sectional unit for
#'      component \code{d}}
#'      \item{\code{zeros: }}{numeric vector of dimension \code{DD} with
#'      values 1, 2, 3 or 4:
#'      \itemize{
#'         \item{1: }{A dummy pattern that starts at the beginning with zeros
#'         and jumps after half of the overall time period}
#'         \item{2: }{A dummy pattern that starts at the beginning with ones and
#'         plummets to zeros after half of the overall time period}
#'         \item{3: }{A dummy pattern that starts at the beginning with one,
#'         then plummets to zeros after a third of the overall time periods, and
#'         then reverts back to ones for the last third of the time}
#'         \item{4: }{A dummy pattern that starts at the beginning with zeros,
#'         then jumps to ones after a third of the overall time periods, and
#'         then reverts back to zeros for the last third of the time}
#'        }
#'        }
#'   }
#' @param seed_no integer; random seed set at the beginning of the simulation
#'   for reproducibility purposes, set \code{NULL} if not required
#' @param options_plot a list of options for plotting data after simulation:
#'   \itemize{
#'   \item{\code{plt_y: }}{logical; if \code{TRUE}, measurements are plotted
#'   per cross sectional unit \code{n=1,...,N}}
#'   \item{\code{plt_x: }}{logical; if \code{TRUE}, latent states are plotted
#'   per cross sectional unit \code{n=1,...,N} with a joint plot of all
#'   components together}
#'   \item{\code{plt_x_per_d: }}{logical; if \code{TRUE}, latent states are
#'   plotted per cross sectional unit \code{n=1,...,N} with a separate plot for
#'   each component}
#'   }
#'
#' @return a named list of two elements, the first one being a list of three:
#'   \itemize{
#'   \item{\code{data: }}{a list of at most two elements (if \code{distribution}
#'   is of type Dirichlet, the second element is \code{NULL}, but otherwise has
#'   the total number of counts e.g. in a Dirichlet-multinomial model)}
#'   \item{\code{regs: }}{list of regressors with two elements: 'z' and 'u'}
#'   \item{\code{states: }}{simulated latent states}
#'   }
#'   and the second being the first argument i.e. an instance of class
#'   \code{trueParams} that is used to generate simulated data (from the
#'   corresponding true parameter values as stored in this object)
#'
#' @export
new_dataSim <- function(true_params,
                        distribution,
                        x_levels,
                        X_LOG_SCALE,
                        options_include = list(intercept = NULL,
                                               policy = NULL,
                                               zeros = NULL),
                        options_plot = list(measurements = FALSE,
                                            states = FALSE,
                                            states_each_d = FALSE),
                        seed_no = NULL) {
  check_class_true_params(true_params)
  check_true_params_distribution(true_params)

  NN <- get_dimension(true_params, "NN")
  TT <- get_dimension(true_params, "TT")
  DD <- get_dimension(true_params, "DD")

  opt1      <- get_opt_include(options_include, NN, DD)
  x_levels  <- get_x_levels(true_params, distribution, x_levels, X_LOG_SCALE)
  seed_no   <- get_seed_no(true_params, seed_no); set.seed(seed_no);
  reg_types <- get_modelling_reg_types(true_params)

  x <- generate_y_x_containter(NN = NN, TT = TT, DD = DD)
  z <- generate_z_u_container(true_params[["beta_z_lin"]],
                              NN = NN, TT = TT, DD = DD,
                              cnt_name = "z",
                              reg_types[["z-linear-regressors"]])
  u <- generate_z_u_container(true_params[["beta_u_lin"]],
                              NN = NN, TT = TT, DD = DD,
                              cnt_name = "u",
                              reg_types[["u-linear-regressors"]])
  for (n in 1:NN) {
    par_true_current <- get_par_true_n(true_params, reg_types, n)
    out_data_tmp <- generate_data_t(TT = TT, DD = DD,
                                    par_true = par_true_current,
                                    x_levels = x_levels[, n],
                                    options_include = opt1[[n]],
                                    modelling_reg_types = reg_types)
    if (reg_types[["z-linear-regressors"]]) {
      z[, , n] <- out_data_tmp$z
    }
    if (reg_types[["u-linear-regressors"]]) {
      u[, , n] <- out_data_tmp$u
    }
    x[, , n] <- out_data_tmp$x_states
  }
  y <- get_measurements(x, X_LOG_SCALE, distribution)
  if (any(sapply(options_plot, isTRUE))) {
    for (n in 1:NN) {
      plot_data_per_n(DD,
                      yraw = y[["part1"]][, , n, drop = FALSE],
                      x = x[, , n, drop = FALSE],
                      x_levels = x_levels[, n],
                      plot_measurements = options_plot$plt_y,
                      plot_states       = options_plot$plt_x,
                      plot_states_each_d  = options_plot$plt_x_per_d,
                      cs = n)
    }
  }
  out_data_sim <- structure(get_output_data_simul(y, x, z, u, reg_types),
                            TRUE_PARAMS = true_params,
                            SEED_NO = seed_no,
                            class = "dataSim")
  return(out_data_sim)
}
#' Returns "trueParams" that comes with the \code{class} "dataSim"
#'
#' The "trueParams" accompanying the simulated data class is stored as an
#' attribute. This helper retrieves this attribute with appropriate class
#' checks.
#'
#' @param data_sim object of \code{class} "dataSim"
#'
#' @return an object of \code{class} "trueParams" accompanying this data class
#' @export
get_true_params_obj <- function(data_sim) {
  check_class_data_sim(data_sim)
  check_class_true_params(attr(data_sim, which = "TRUE_PARAMS"))
  attr(data_sim, which = "TRUE_PARAMS")
}
check_class_data_sim <- function(obj) {
  checker <- inherits(obj, "dataSim")
  stopifnot(`Arg. "obj" must be of class "dataSim".` = checker)
  check_class_true_params(attr(obj, which = "TRUE_PARAMS"))
  return(invisible(obj))
}
#' Helper function to get good values for state levels to simulation
#'
#' As a side effect, checks if argument 'distribution' is valid.
#'
#' @inheritParams new_dataSim
#'
#' @return resulting target levels, that might need adjustment e.g. taking to
#'    the log-scale or standardization with the variance
get_x_levels <- function(true_params, distribution, x_levels, X_LOG_SCALE) {
  densitities_supported <- c("multinomial", "dirichlet_mult",
                             "gen_dirichlet_mult", "gen_dirichlet",
                             "dirichlet", "normal")
  if (!(distribution %in% densitities_supported)) {
    stop(paste0("Argument to distribution must be one of: ",
                paste0(densitities_supported, collapse = ", "), "!"))
  }
  if (X_LOG_SCALE && distribution == "normal") x_levels <- log(x_levels)
  if (X_LOG_SCALE && any(distribution %in% c("dirichlet", "gen_dirichlet",
                                             "multinomial", "dirichlet_mult",
                                             "gen_dirichlet_mult"))) {
    if (distribution %in% c("gen_dirichlet", "gen_dirichlet_mult")) {
      browser()
      sig_tmp  <- get_params(true_params, n = 1,
                             name_par = "sig_sq", DD_TYPE = "A")
    } else {
      sig_tmp  <- get_params(true_params, n = 1, name_par = "sig_sq")
      sig_tmp <- sig_tmp[, 1]
    }
    x_levels <- log(x_levels) - sig_tmp/2
  }
  return(x_levels)
}
#' Helper function to set the seed
#'
#' @inheritParams new_dataSim
#'
#' @return a
get_seed_no <- function(true_params, seed_no) {
  if (is.null(seed_no)) seed_no <- get_seed(true_params)
  seed_no
}
#' Options list of included effects.
#'
#' Refer to intercept specifications (either at z-type regressors or random
#' effects), policy dummies, or zero specifications.
#'
#' @param includes a list of three elements named: 'intercept',
#'   'policy', and 'zeros'
#' @inheritParams new_dataSim
#'
#' @return a list of the same structure as includes but with elements adjusted
#'   for model dimension; \code{includes} is a list of \code{NULL} elements,
#'   then the default specifications are returned (see the function body)
get_opt_include <- function(includes, NN, DD) {
  intercept <- includes$intercept
  policy    <- includes$policy
  zeros     <- includes$zeros
  if (is.null(intercept)) {
    intercept <- list()

    intercept$at_z <- rep(FALSE, times = DD)
    intercept$at_u <- rep(FALSE, times = DD)

    names(intercept$at_z) <- paste0("d_", seq_len(DD))
    names(intercept$at_u) <- paste0("d_", seq_len(DD))
  }
  if (is.null(policy)) {
    policy <- matrix(FALSE, nrow = DD, ncol = NN)
    rownames(policy) <- paste0("d_", seq_len(DD))
    colnames(policy) <- paste0("n_", seq_len(NN))
    # policy_modelling    <- cbind(c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
    #                              c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
    #                              c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
    #                              c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
    #                              c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
    #                              c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE))
  }
  if (is.null(zeros)) {
    zeros <- NULL
    # zeros      <- c(1:4, 1, 2) # rep(4, times = DD) # c(1, 2, 3, 4, 1, 2)
    # names(zeros) <- paste0("d_", seq_len(DD))
  }
  out_opt <- vector("list", length = NN)
  for (n in 1:NN) {
    out_opt[[n]] <- list(intercept = intercept,
                         policy = policy[, n, drop = TRUE],
                         zeros = zeros)
  }
  names(out_opt) <- paste0("N_", 1:NN)
  return(out_opt)
}

#' Generate container for measurements and states of appropriate dimension.
#'
#' For \code{distribution} of appropriate type, where the data-part is compound
#' of two types, a list of two elements (not \code{NULL}) is returned; otherwise
#' a list element can be \code{NULL} indicating that this element is not needed.
#'
#' @inheritParams new_dataSim
#'
#' @return a list of two; \code{data} and \code{states} where the former can
#'   itself be a list of two elements (part1 and part2 if data is compound e.g.
#'   of counts and total counts as for the Multinomial distribution)
generate_y_x_containter <- function(NN, TT, DD) {
  tmp_names <- get_x_y_containter_names(NN = NN, TT = TT, DD = DD)

  out_cnt <- array(0, c(TT, DD, NN))
  dimnames(out_cnt) <- tmp_names

  return(out_cnt)
}
get_x_y_containter_names <- function(NN, TT, DD) {
  names_out <- list(paste0("t_", seq_len(TT)),
                    paste0("d_", seq_len(DD)),
                    paste0("n_", seq_len(NN)))
  return(names_out)
}
#' Generate container for regressors.
#'
#' The multivariate dimension \code{DD} need not be specified as necessary
#' information is directly inferred from container of true parameters
#' (\code{par}).
#'
#' @param pars container of true parameter values as passed to main function
#'   [new_dataSim(true_params = ...)]
#' @inheritParams new_dataSim
#' @param cnt_name a character: either "z" for z-type regressors or "u" for the
#'   random effects container; other imput gives an error
#' @param reg_type output from [get_modelling_reg_types()] specifying the type
#'   of regressors (effects) to generate
#'
#' @return a named list of two elements: \code{z} and \code{u} for z-type or
#'   u-type regressors; elements can be \code{NULL} whenever the corresponding
#'   regressor type is not needed.
generate_z_u_container <- function(pars, NN, TT, DD, cnt_name, reg_type) {
  if (reg_type) {
    if(cnt_name == "z") {
      dim_bet <- sapply(pars, length)
    } else if (cnt_name == "u") {
      if (NN == 1) {
        dim_bet <- sapply(pars, length)
        num_bet <- sum(dim_bet)
      } else {
        dim_bet <- sapply(pars, nrow)
        num_bet <- sum(dim_bet)
      }
    } else {
      stop("Unknown container name.")
    }
    num_bet <- sum(dim_bet)
    names_cnt <- paste0(paste0(cnt_name, unlist(sapply(dim_bet, seq_len))),
                        "_d", rep(1:DD, unlist(dim_bet)))
    tmp_names <- list(paste0("t_", seq_len(TT)),
                      names_cnt,
                      paste0("n_", seq_len(NN)))
    cnt <- array(0, c(TT, num_bet, NN))
    dimnames(cnt) <- tmp_names
  } else {
    cnt <- NULL
  }
  return(cnt)
}
#' Subset container of true parameters by current cross sectional unit \code{n}
#'
#' @inheritParams generate_z_u_container
#' @param n integer giving the cross sectional unit for which to subset
#'
#' @return subset of true parameter container for current cross sectional unit
get_par_true_n <- function(pars, reg_types, n) {
  pars_out <- list(sig_sq = pars[["sig_sq"]][, n],
                   phi = lapply(pars[["phi"]], `[`, i = , j = n))
  if (reg_types[["z-linear-regressors"]]) {
    pars_out$beta_z_lin <- pars[["beta_z_lin"]]
  }
  if (reg_types[["u-linear-regressors"]]) {
    pars_out$beta_u_lin <- lapply(pars[["beta_u_lin"]], `[`, i = , j = n)
  }
  return(pars_out)
}
#' Produce output container of the main simulation function
#'
#' Takes input from intermediate containers; internal function.
#'
#' @param cnt_data output container of measurements
#' @param cnt_states output container of latent states
#' @param cnt_z output container of z-type regressors
#' @param cnt_u output container of u-type-regressors
#' @param reg_types regressor type specification
#'
#' @return a names list of appropriate structure to return from main function
get_output_data_simul <- function(cnt_data,
                                  cnt_states,
                                  cnt_z, cnt_u,
                                  reg_types) {
  dist <- attr(cnt_data, which = "distribution")
  out_data <- vector("list", 3)
  if (dist %in% c("dirichlet", "normal")) {
    out_data[[1]] <- list(yraw = cnt_data[["part1"]])
  } else if (dist == "multinomial" || dist == "dirichlet_mult") {
    out_data[[1]] <- list(yraw = cnt_data[["part1"]],
                          num_counts = cnt_data[["part2"]])
  } else {
    stop("Unknown distribution attribute of data.")
  }
  attr(out_data[[1]], which = "model_type_obs") <- toupper(dist)
  if (reg_types[["z-linear-regressors"]]) {
    out_data[[2]]$z <- cnt_z
  } else{
    out_data[[2]]$z <- NULL
  }
  if (reg_types[["u-linear-regressors"]]) {
    out_data[[2]]$u <- cnt_u
  } else {
    out_data[[2]]$u <- NULL
  }
  out_data[[3]] <- cnt_states

  lat_type_names <- c("auto", "lin", "re", "splZ", "splU")
  lat_type_names <- paste0(lat_type_names[reg_types], collapse = "_")
  attr(out_data[[3]], which = "model_type_lat") <- lat_type_names
  names(out_data) <- c("data", "regs", "states")
  return(out_data)
}
get_measurements <- function(x_states, X_LOG_SCALE, distribution) {
  tmp_dim <- dim(x_states)
  TT <- tmp_dim[1]
  DD <- tmp_dim[2]
  NN <- tmp_dim[3]

  data_part1 <- generate_y_x_containter(NN = NN, TT = TT, DD = DD)

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
#
#
#
