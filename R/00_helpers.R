msg_dim_ready <- function(dim, print_name) {
  cat(crayon::green("Setting dimension "), crayon::yellow(print_name),
      crayon::green("to "), crayon::red(dim), crayon::green("!\n"))
}
#' Checks if class or argument is admissible w.r.t. distribution name
#'
#' Used internally only; match testing is performed for allowed class names like
#' "Dirichlet", "GenDirichletMult" etc. OR based on argument values for the
#' '\code{distribution} argument. The latter pops up in many internal functions
#' and admissible values are e.g. "dirichlet", "gen_dirichlet" etc. See function
#' definition for more details.
#'
#' @param x  object to check for underlying class name or argument value
#' @param type character: either 'arg' then argument check is performed or
#'    'class', then class-check is performed
#' @param FORCE_SCALAR logical
#'
#' @return pure side effect-function checking for validity of first argument;
#'   throws error if invalid
check_distribution <- function(x, type = "arg", FORCE_SCALAR = FALSE) {
  stopifnot(`Arg. 'FORCE_SCALAR' must be logical` = is.logical(FORCE_SCALAR))
  stopifnot(`Arg. 'type' must be either 'class' or 'arg'` =
              type %in% c("class", "arg"))
  if (isTRUE(FORCE_SCALAR)) stopifnot(`Arg. 'x' not length 1` = length(x) == 1)
  if (type == "arg") {
    dist_names <- c("dirichlet", "gen_dirichlet", "multinomial",
                    "dirichlet_mult", "gen_dirichlet_mult",
                    "normal")
    stopifnot(`Distrubtion arg. not supported` = x %in% dist_names)
  } else if (type == "class") {
    dist_names_class <- c("Dirichlet", "GenDirichlet", "Multinomial",
                          "DirichletMult", "GenDirichletMult",
                          "NORMAL")
    stopifnot(`Object class not supported` = x %in% dist_names_class)
  } else {
    stop("Unknown value for arg. 'type': use either 'arg' or 'class'.")
  }
  return(invisible(x))
}
#' Generic function to get the seed
#'
#' Dispatches on class \code{trueParams} or \code{dataSim}.
#'
#' @param obj an object of class \code{trueParams} or \code{dataSim}; see
#'    details of corresponding S3 methods
#' @param type a character: either "data_sim", returning the seed used for data
#'   simulation, or "true_params" returning the seed for generation of the
#'   underlying \code{trueParams} object; this is only useful for
#'   \code{get_seed.dataSim} i.e. \code{get_seed()} applied to type
#'   \code{dataSim}; defaults to "data_sim" but will be ignored when first
#'   argument is \code{class=="trueParams"}
#'
#'
#' @return an \code{integer} giving the underlying seed number(s)
#' @export
get_seed <- function(obj, type = "data_sim") {
  stopifnot(`Arg. 'type' must be either 'data_sim' or 'true_params'` =
              (type %in% c("data_sim", "true_params")))
  UseMethod("get_seed")
}
#' Helper function to get the seed
#'
#' Allows access to two types of seeds: data simulation seed and true parameter
#' seeds, even when accessing the `dataSim`-class (because the latter stores
#' both seeds!)
#'
#' @inheritParams get_seed
#'
#' @return an \code{integer} giving the underlying seed number(s)
#' @export
get_seed.dataSim <- function(obj, type = "data_sim") {
  check_class_data_sim(obj)

  if (type == "data_sim") return(attr(obj, "SEED_NO"))
  if (type == "true_params") return(get_seed(get_true_params_obj(obj)))
  return(invisible(NULL))
}
#' Access to meta information for object of class `trueParams`
#'
#' Specifically, getting seed number used during `bet_u_lin` and VCM value
#' simulation.
#'
#' @inheritParams get_seed
#'
#' @return seed number of attribute of \code{class} "trueParams"; the random
#'    seed during beta_u type value construction and corresponding VCMs
#' @export
get_seed.trueParams <- function(obj, type = NULL) {
  check_class_true_params(obj)
  attr(obj, which = "meta_info")[["SEED_NO"]]
}
set_seed_no <- function(true_params, seed_no) {
  if (is.null(seed_no)) seed_no <- get_seed(true_params)
  return(seed_no)
}
progress_any <- function(iter,
                         iter_max,
                         settings = list(digits = 0,
                                         repeat_every = 10),
                         msg = NULL) {
  percentage_progress <- iter / iter_max * 100

  chunck_to_print  <- floor(iter_max / settings$repeat_every)
  repeat_every_seq <- seq(from = chunck_to_print,
                          to = iter_max,
                          by = chunck_to_print)
  check_print_progress <- iter %in% repeat_every_seq
  if (check_print_progress && is.null(msg)) {
    cat(crayon::green("Progress: "),
        crayon::yellow(
          paste0(round(percentage_progress, digits = settings$digits), "%.\n")))
  } else if (check_print_progress && !is.null(msg)) {
    cat(crayon::green("Progress - "),
        crayon::green(paste0(msg, ": ")),
        crayon::yellow(
          paste0(round(percentage_progress, digits = settings$digits), "%.\n")))
  }
}
cry <- function(x) {
  crayon::yellow(x)
}
crg <- function(x) {
  crayon::green(x)
}