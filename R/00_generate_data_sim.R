#' Generates multivariate panel data for various measurement models
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
#' @inheritParams new_trueParams
#' @param x_levels target "mean" levels of the states around which they
#'   fluctuate
#' @param X_LOG_SCALE logical; if \code{TRUE}, \code{x_levels} are taken as logs
#'   and the random number generation for the states is performed in logs
#' @param options_include a list of options for various effects:
#'   \itemize{
#'      \item{\code{intercept: }}{a list of logical vectors as defined by
#'      [generate_ic_list()] and checked via [check_ic_list()] to represent
#'      intercept generation;
#'      list elements are logical, and if `TRUE`, include an intercept at the
#'      cross sectional unit for component `d`; `FALSE` does not include an
#'      intercept}
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
#'   the total number of counts e.g. in a Dirichlet-Multinomial model)}
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
                                               # policy = NULL,
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

  check_ic_to_dist(distribution, options_include$intercept, DD)

  opt1      <- set_opt_include(distribution, options_include, NN, DD)
  x_levels  <- set_x_levels(true_params, x_levels, X_LOG_SCALE)
  reg_types <- get_modelling_reg_types(true_params)

  x <- generate_y_x_containter(distribution, NN = NN, TT = TT, DD = DD)
  z <- generate_z_u_container(true_params, NN, TT, DD, "z")
  u <- generate_z_u_container(true_params, NN, TT, DD, "u")

  seed_no <- set_seed_no(true_params, seed_no)
  set.seed(seed_no)

  for (n in 1:NN) {
    out_data_tmp <- generate_data_t(nn = n, TT = TT, DD = DD,
                                    true_params = true_params,
                                    x_levels = x_levels[, n],
                                    options_include = opt1[[n]],
                                    modelling_reg_types = reg_types)
    if (reg_types[["z-linear-regressors"]]) {
      z[, , n] <- out_data_tmp$z_regs
    }
    if (reg_types[["u-linear-regressors"]]) {
      u[, , n] <- out_data_tmp$u_regs
    }
    x[, , n] <- out_data_tmp$x_states
  }
  y <- generate_measurements(x, X_LOG_SCALE, distribution,
                             get_dimension(true_params, dim = "all"),
                             options_include = opt1)
  if (any(sapply(options_plot, isTRUE))) {
    for (n in 1:NN) {
      plot_data_per_n(distribution, DD,
                      yraw = y[["part1"]][, , n, drop = FALSE],
                      x = x[, , n, drop = FALSE],
                      x_levels = x_levels[, n],
                      plot_measurements = options_plot$plt_y,
                      plot_states       = options_plot$plt_x,
                      plot_states_each_d  = options_plot$plt_x_per_d,
                      cs = n)
    }
  }
  out_data_sim_core <- get_output_data_simul(y, x, z, u, reg_types)
  out_data_sim <- structure(out_data_sim_core,
                            TRUE_PARAMS = true_params,
                            SEED_NO = seed_no,
                            class = c(class(out_data_sim_core), "dataSim"))
  return(out_data_sim)
}
check_class_data_sim <- function(obj) {
  checker <- inherits(obj, "dataSim")
  stopifnot(`Arg. "obj" must be of class "dataSim".` = checker)
  check_class_true_params(attr(obj, which = "TRUE_PARAMS"))
  return(invisible(obj))
}
#' Returns "trueParams" object that comes with the \code{class} "dataSim"
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
#' Get z-type regressors from an instance of "dataSim"
#'
#' @inheritParams get_true_params_obj
#'
#' @return the data container of the z-type regressor
#' @export
get_regs_z <- function(data_sim) {
  check_class_data_sim(data_sim)
  data_sim[["regs"]][["z"]]
}
#' Get u-type regressors from an instance of "dataSim"
#'
#' @inheritParams get_true_params_obj
#'
#' @return the data container of the u-type regressor
#' @export
get_regs_u <- function(data_sim) {
  check_class_data_sim(data_sim)
  data_sim[["regs"]][["u"]]
}
#' Get measurements/observations i.e. the data from an instance of "dataSim"
#'
#' @inheritParams get_true_params_obj
#' @param type a character giving the "type" of measurement/observational
#'    data to return; either "raw" or "count"
#'
#' @return the measurement/observation data set of appropriate type e.g. the
#'   "raw" fractions if the distribution is a Dirichlet or component counts and
#'    total counts of distribution is of Multinomial type (e.g.
#'    Dirichlet-Multinomial)
#' @export
get_data <- function(data_sim, type) {
  check_class_data_sim(data_sim)
  stopifnot(`Unknown value for arg. 'type'` = type %in% c("raw", "count"))
  if (type == "raw") {
    return(data_sim[["data"]][["yraw"]])
  } else if (type == "count") {
    return(data_sim[["data"]][["num_counts"]])
  }
}
#' Get latent states from an instance of "dataSim"
#'
#' @inheritParams get_true_params_obj
#'
#' @return the latent states which are simulated
#' @export
get_sim_latent_states <- function(data_sim) {
  check_class_data_sim(data_sim)
  data_sim[["states"]]
}
#' Get type of measurements/observations of the "dataSim" object
#'
#' @inheritParams get_true_params_obj
#'
#' @return a character string of either: "NORMAL", "DIRICHLET",
#'    "GEN_DIRICHLET", "MULTINOMIAL", "DIRICHLET_MULT", or "GEN_DIRICHLET_MULT"
#' @export
get_type_obs <- function(data_sim) {
  check_class_data_sim(data_sim)
  attr(data_sim, "model_type_obs")
}
#' Get type of latent states of the "dataSim" object
#'
#' @inheritParams get_true_params_obj
#'
#' @return a character string of either: "auto", "lin", "re", "splZ", "splU" or
#'    a combination thereof
#' @export
get_type_lat <- function(data_sim) {
  check_class_data_sim(data_sim)
  attr(data_sim, "model_type_lat")
}
#' Generic function to model dimension meta info
#'
#' Dispatches on class \code{trueParams}, \code{trueParamsDirichlet},
#' \code{trueParamsGenDirichlet} etc. or \code{dataSim}.
#'
#' Applicable to object of class "dataSim", "trueParams" or derived classed of
#' the latter. The dimensions of the implied model for which the "trueParams"
#' object is meant to be used are returned.
#'
#' @inheritParams get_seed
#' @param dim a character string; either of "NN", "TT", "DD" or "all" which
#'   returns all three
#'
#' @return dimension for object of class \code{class} "trueParams"
#' @export
get_dimension <- function(obj, dim = NULL) {
  UseMethod("get_dimension")
}
#' S3 method for generic 'get_dimension' for class "dataSim"
#'
#' See [get_dimension()] for details.
#'
#' @inheritParams get_dimension
#'
#' @return dimension for object of class (underlying) \code{class} "trueParams"
get_dimension.dataSim <- function(obj, dim = NULL) {
  check_class_data_sim(obj)
  true_params <- get_true_params_obj(obj)
  get_dimension.trueParams(true_params, dim)
}
