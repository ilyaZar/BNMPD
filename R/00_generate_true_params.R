#' Generates output container for true parameter values
#'
#' Arguments specify true parameter values, dimension of the model and which
#' regressor types to use. For random effects (beta_u-type regressors) a seed
#' makes sure to produce the same setup. This output container can be used when
#' plotting diagnostics (in a simulation study) as well as generating simulated
#' data sets based on the true parameter values.
#'
#' @param dim_model a vector with three components: \code{NN x TT x DD}
#' @param sig_sq a vector of length \code{DD} and strictly positive elements; if
#'    \code{NULL}, then defaults are generated (see [get_default_sig_sq])
#' @param phi a vector of length \code{DD} with elements between 0 and 1; if
#'    \code{NULL}, then defaults are generated (see [get_default_phi])
#' @param beta_z_lin a list of length \code{DD} with (possibly) varying number
#'   of regressors; if \code{NULL}, then defaults are generated (see
#'   [get_default_beta_z_lin])
#' @param settings_pars a list with named elements of the form:
#'   \itemize{
#'   \item{SIMUL_PHI: }{logical; if \code{TRUE}, then autoregressive process of
#'   order \code{order_p_vec} is generated}
#'   \item{SIMUL_Z_BETA: }{logical; if \code{TRUE}, then standard continuous
#'   (z-type) covariates are included}
#'    \item{SIMUL_U_BETA: }{logical; if \code{TRUE}, then random effects (u-type
#'    regressors) are simulated based on the argument \code{seed_taken}}
#'   \item{num_z_regs: }{integer giving the number of z-type regressors (either
#'   the same number for all multivariate components, or a vector of length
#'   equal to the number of multivariate components)}
#'   \item{num_u_regs: }{integer giving the number of random effects per
#'   cross section to simulate (either the same number for all multivariate
#'   components, or a vector of length equal to the number of multivariate
#'   components; currently not varying along the cross sectional units but
#'   extension is possible)}
#'   \item{order_p_vec: }{a vector of length DD with integers >= 0 giving the
#'   order of auto-regression per component d=1,...,DD}
#'   }
#' @param options a list with settings for true parameter generation with first
#'   component being a \code{dwn_scl} value for the variance parameter
#'   \code{sigma} and the second an optional list with two logical vectors of
#'   length DD indicating whether the corresponding component (either at z or u)
#'   should have an intercept
#' @param seed_taken the seed used for drawing random effects
#' @inheritParams generate_data_t_n
#'
#' @return an object of class "\code{trueParams}" which is a list of three: true
#'   parameter values as a list, meta information such as model dimension, the
#'   logical indicators that describe which parameters to generate and their
#'   number:
#'   \itemize{
#'   \item{\code{SIMUL_PHI:}}{ logical; if \code{TRUE} autoregressive latent
#'   laten states have been simulated}
#'   \item{\code{SIMUL_Z_BETA:}}{ logical; if \code{TRUE} z-regressors have been
#'   simulated}
#'   \item{\code{SIMUL_U_BETA:}}{ logical; if \code{TRUE} u-regressors have been
#'   simulated}
#'   \item{\code{num_z_regs:}}{ number of z-type regressors or \code{NULL}}
#'   \item{\code{num_z_regs:}}{ number of u-type regressors or \code{NULL}}
#'   \item{\code{order_p_vec:}}{ an integer vector giving the autoregressive
#'   lag order for each \code{d} in \code{1,...,DD}}
#'   }
#' @export
new_trueParams <- function(distribution,
                           model_dim,
                           sig_sq = NULL,
                           phi = NULL,
                           beta_z_lin = NULL,
                           settings_pars = list(SIMUL_PHI = TRUE,
                                                SIMUL_Z_BETA = TRUE,
                                                SIMUL_U_BETA = TRUE,
                                                num_z_regs = 3,
                                                num_u_regs = 2,
                                                order_p_vec = 1),
                           options = list(dwn_scl = 10,
                                          intercepts = NULL),
                           seed_taken) {
  SIMUL_PHI    <- settings_pars$SIMUL_PHI
  SIMUL_Z_BETA <- settings_pars$SIMUL_Z_BETA
  SIMUL_U_BETA <- settings_pars$SIMUL_U_BETA
  num_z_regs   <- settings_pars$num_z_regs
  num_u_regs   <- settings_pars$num_u_regs
  order_p_vec  <- settings_pars$order_p_vec
  check_args <- (missing(SIMUL_U_BETA) || missing(SIMUL_Z_BETA) ||
                   missing(SIMUL_PHI) || missing(seed_taken) ||
                   missing(model_dim))
  if (check_args) stop("Missing arguments")
  if (SIMUL_U_BETA && missing(num_u_regs)) stop("Number of REs not specified!")
  # 1. Data settings: -------------------------------------------------------
  NN <- model_dim[1]   # Cross sectional length
  cat(crayon::green("Setting dimension "), crayon::yellow("NN"),
      crayon::green("to "), crayon::red(NN), crayon::green("!\n"))
  TT <- model_dim[2]  # Time series length
  cat(crayon::green("Setting dimension "), crayon::yellow("TT"),
      crayon::green("to "), crayon::red(TT), crayon::green("!\n"))
  DD <- model_dim[3]
  cat(crayon::green("Setting dimension "), crayon::yellow("DD"),
      crayon::green("to "), crayon::red(DD), crayon::green("!\n"))
  if (distribution %in% c("gen_dirichlet", "gen_dirichlet_mult")) {
    DD2 <- get_DD2(distribution = distribution, DD = model_dim[[3]])
    cat(crayon::green("Setting !internal! dimension "), crayon::yellow("DD2"),
        crayon::green("to "), crayon::red(DD2),
        crayon::green(paste0("for ", distribution," type parameters!\n")))
    model_dim <- c(model_dim, DD2 = DD2)
  }
  # 2. Set up parameter values: ---------------------------------------------
  true_sig_sq <- new_sig_sq_x(
    distribution,
    sig_sq,
    DD,
    NN,
    options$dwn_scl)
  true_phi <- new_phi(
    SIMUL_PHI,
    distribution,
    phi,
    DD,
    NN,
    order_p_vec)
  true_bet_z  <- new_bet_z(
    SIMUL_Z_BETA,
    distribution,
    beta_z_lin,
    DD, num_z_regs,
    options$intercepts$at_z)
  tmp_u  <- new_bet_vcm_u(
    SIMUL_U_BETA,
    distribution,
    DD, NN,
    num_u_regs,
    seed_taken,
    options$intercepts$at_u)

  true_bet_u <- tmp_u[["true_bet_u"]]
  true_D0u_u <- tmp_u[["true_D0u_u"]]

  par_trues <- list(sig_sq = true_sig_sq,
                    phi = true_phi,
                    beta_z_lin = true_bet_z,
                    beta_u_lin = true_bet_u,
                    vcm_u_lin = true_D0u_u)

  structure(par_trues,
            meta_info = list(MODEL_TYPE = distribution,
                             MODEL_DIM = model_dim,
                             PAR_SETTINGS = settings_pars,
                             SEED_NO = seed_taken),
            class = get_class_true_param(distribution))
}
#' @param trueParam object of `class` "trueParam" and one of its subclasses i.e.
#'    "trueParamDirichlet", "trueParamGenDirichlet" etc.
get_meta_info_all <- function(trueParam) {

}
#' Sets true values (default or user supplied) for parameter phi
#'
#' @inheritParams new_trueParams
#' @param DD number of multivariate components
#' @param NN number of cross sectional units
#'
#' @return a list of length \code{DD} with each element being a matrix of
#'   dimension \code{order_p_vec[d] x NN} that stores the phi's per number of
#'   autoregressions, number of cross section and per component \code{d}.
#' @export
new_phi <- function(SIMUL_PHI, distribution, phi, DD, NN, order_p_vec) {
  if (SIMUL_PHI) {
    order_p_vec <- get_order_p_vec(order_p_vec, DD)
    if (!is.null(phi)) {
      out_phi <- get_manual_phi(phi, DD, order_p_vec, NN)
    } else {
      out_phi <- get_default_phi(distribution, DD, NN, order_p_vec)
    }
  } else {
    out_phi <- NULL
  }
  structure(out_phi,
            class = c("true_phi"))
}
#' Sets true values (default or user supplied) for parameter sig_sq_x
#'
#' @inheritParams new_trueParams
#' @inheritParams new_phi
#'
#' @return a matrix of dimension \code{DD x NN} of true parameter values for
#'   sig_sq_x
new_sig_sq_x <- function(distribution, sig_sq, DD, NN, dwn_scl) {
  if (distribution %in% c("dirichlet", "multinomial", "dirichlet_mult")) {
    if (!is.null(sig_sq)) {
      check_sig_sq <- all(sig_sq > 0) && all(length(sig_sq) == DD)
      if (!check_sig_sq) stop("'sig_sq' > 0 and a vector of length DD ...\n")
      out_sig_sq <- matrix(sig_sq, nrow = DD, ncol = NN)
    } else {
      tmp_sig_sq  <- get_default_sig_sq(distribution, DD, dwn_scl)
      out_sig_sq <- matrix(tmp_sig_sq, nrow = DD, ncol = NN)
    }
    rownames(out_sig_sq) <- paste0("D", 1:DD)
    colnames(out_sig_sq) <- paste0("N", 1:NN)
    return(structure(out_sig_sq,
                     class = c("true_sig_sq", "matrix", "array")))
  } else if (distribution %in% c("gen_dirichlet")) {
    if (!is.null(sig_sq)) {
      check_sig_sq <- all(sig_sq > 0) && all(length(sig_sq) == DD)
      if (!check_sig_sq) stop("'sig_sq' > 0 and a vector of length DD ...\n")
      out_sig_sq <- array(sig_sq, c(DD, NN, 2))
    } else {
      tmp_sig_sq <- get_default_sig_sq(distribution, DD, dwn_scl)
      out_sig_sq <- array(tmp_sig_sq, c(DD, NN, 2))
    }
    dimnames(out_sig_sq) <- list(paste0("D", 1:DD),
                                 paste0("N", 1:NN),
                                 c("A", "B"))
    return(structure(out_sig_sq,
                     class = c("true_sig_sq", "matrix", "array")))
  } else if (distribution %in% c("gen_dirichlet_mult")) {
    if (!is.null(sig_sq)) {
      check_sig_sq <- all(sig_sq > 0) && all(length(sig_sq) == DD)
      if (!check_sig_sq) stop("'sig_sq' > 0 and a vector of length DD ...\n")
      out_sig_sq <- array(sig_sq, c(DD - 1, NN, 2))
    } else {
      tmp_sig_sq <- get_default_sig_sq(distribution, DD, dwn_scl)
      out_sig_sq <- array(tmp_sig_sq, c(DD - 1, NN, 2))
    }
    dimnames(out_sig_sq) <- list(paste0("D", 1:(DD - 1)),
                                 paste0("N", 1:NN),
                                 c("A", "B"))
    return(structure(out_sig_sq,
                     class = c("true_sig_sq", "matrix", "array")))
  } else {
    stop("Uknown distribution argument.")
  }
  return(invisible(distribution))
}
#' Sets true values (default or user supplied) for parameter bet_z
#'
#' @inheritParams new_trueParams
#' @inheritParams new_phi
#'
#' @return a list of length \code{DD} each element being a vector of the
#'   corresponding number of regressors
new_bet_z <- function(SIMUL_Z_BETA,
                      distribution,
                      beta_z_lin,
                      DD, num_z_regs,
                      intercepts) {
  if (SIMUL_Z_BETA) {
    if (!is.null(beta_z_lin)) {
      out_bet_z <- get_manual_bet_z(beta_z_lin, DD)
    } else {
      if (is.null(num_z_regs)) stop("Num. of 'beta_z_lin' regressors required.")
      out_bet_z <- get_default_beta_z_lin(
        distribution,
        DD,
        num = num_z_regs,
        intercepts)
    }
  } else {
    out_bet_z <- NULL
  }
  structure(out_bet_z,
            class = "true_bet_z")
}
#' Sets true values (default or user supplied) for parameter bet_u
#'
#' @inheritParams new_trueParams
#' @inheritParams new_phi
#'
#' @return a list of two elements for bet_u and vcm_bet_u each again a list of
#'   dimension \code{DD} with first element being a matrix of dimension
#'   \code{num_u_regs x NN} (matrix of true random effects) and the second
#'   matrix \code{num_u_regs x num_u_regs} (covariance matrix for random
#'   effects)
#' @export
new_bet_vcm_u <- function(SIMUL_U_BETA, distribution,
                          DD, NN,
                          num_u_regs,
                          seed_taken,
                          intercepts) {
  if (SIMUL_U_BETA) {
    DD <- get_DD(distribution, DD)
    num_reg_seq <- get_num_reg_seq(num_u_regs, DD)
    true_out_u  <- generate_bet_u(
      distribution,
      DD, NN, TRUE, num_reg_seq,
      seed_no = seed_taken)
    true_bet_u <- true_out_u[[1]]
    true_D0u_u <- true_out_u[[2]] # nolint: object_name_linter.
  } else {
    true_bet_u <- NULL
    true_D0u_u <- NULL # nolint: object_name_linter.
  }
  structure(list(true_bet_u = true_bet_u, true_D0u_u = true_D0u_u),
            class = "true_bet_vcm_u")
}
#' Generic methods to extract parameters
#'
#' Dispatches on attribute-class of \code{trueParam} access appropriate
#' parameters for cross sectional unit \code{n}n
#'
#' @param true_params object of `class` "trueParams" and one of its subclasses
#'   i.e. "trueParamDirichlet", "trueParamGenDirichlet" etc.
#' @param n cross sectional unit; an integer or vector of integers of
#'    appropriate values matching (sub-)indices of the cross sectional dimension
#'    if not then R-type out of bounds error is thrown
#' @param DD multivariate dimension; an integer or vector of integers of
#'    appropriate values matching (sub-)indices of the multivariate dimension;
#'    if not then R-type out of bounds error is thrown
#' @param name_par optional string giving the parameter name; if \code{NULL}
#'    then all parameters are returned; must be of either 'sig_sq', 'phi',
#'    'beta_z_lin', 'beta_u_lin', or 'vcm_u_lin'
#' @param DD_TYPE optional; if not \code{NULL}, then it must be "A", "B" or "AB"
#'    to use DD-A, DD-B or both slice(s) of multivariate component; only
#'    applicable for special distributions such as generalized Dirichlet or
#'    generalized Multinomial Dirichlet
#'
#' @return sliced parameters for some cross sectional unit \code{n}
#' @export
get_params <- function(true_params, n = NULL, DD = NULL,
                       name_par = NULL, DD_TYPE = NULL) {
  stopifnot(`Arg. 'n' is either NULL or a single number.` = length(n) <= 1)
  stopifnot(`Arg. 'DD' is either NULL or a single number.` = length(DD) <= 1)
  stopifnot(`Arg. 'name_par' is either NULL or a single number.` = length(name_par) <= 1)
  stopifnot(`Arg. 'DD_TYPE' is either NULL or a single number.` = length(DD_TYPE) <= 1)
  UseMethod("get_params")
}
#' Method for class 'trueParamsDirichlet' derived from 'trueParams'
#'
#' Used for the Dirichlet distribution.
#'
#' @inheritParams get_params
#'
#' @return sliced parameters for some cross sectional unit \code{n}
#' @export
get_params.trueParamsDirichlet <- function(true_params,
                                           n = NULL, DD = NULL,
                                           name_par = NULL,
                                           DD_TYPE = NULL) {
  if (is.null(n)) n <- seq_len(nrow(true_params$sig_sq))
  reg_types <- get_modelling_reg_types(true_params)
  pars_out  <- true_params %>%
    get_default_pars(n, DD, reg_types) %>%
    get_par_name(name_par)
  return(pars_out)
}
#' Method for class 'trueParamsGenDirichlet' derived from 'trueParams'
#'
#' Used for the generalized Dirichlet distribution.
#'
#' @inheritParams get_params
#'
#' @return sliced parameters for some cross sectional unit \code{n}
#' @export
get_params.trueParamsGenDirichlet <- function(true_params,
                                              n = NULL, DD = NULL,
                                              name_par = NULL,
                                              DD_TYPE = NULL) {
  if (is.null(n)) n <- seq_len(dim(true_params$sig_sq)[1])
  if (missing(DD_TYPE)) stop("Must set arg. 'DD_TYPE' for gen. Dirichlet!" )
  reg_types <- get_modelling_reg_types(true_params)
  pars_out  <- true_params %>%
    get_special_pars(n, DD, reg_types, special_type = DD_TYPE) %>%
    get_par_name(name_par)
  return(pars_out)
}
#' Method for class 'trueParamsDirichletMult' derived from 'trueParams'
#'
#' Used for the Dirichlet Multinomial distribution.
#'
#' @inheritParams get_params
#'
#' @return sliced parameters for some cross sectional unit \code{n}
#' @export
get_params.trueParamsDirichletMult <- function(true_params,
                                               n = NULL, DD = NULL,
                                               name_par = NULL,
                                               DD_TYPE = NULL) {
  if (is.null(n)) n <- seq_len(nrow(true_params$sig_sq))
  reg_types <- get_modelling_reg_types(true_params)
  pars_out  <- true_params %>%
    get_default_pars(n, DD, reg_types) %>%
    get_par_name(name_par)
  return(pars_out)
}
get_default_pars <- function(true_params, n, DD = NULL, reg_types) {
  pars_out <- list(sig_sq = true_params[["sig_sq"]][, n, drop = FALSE],
                   phi = lapply(true_params[["phi"]], `[`, i = , j = n))
  if (reg_types[["z-linear-regressors"]]) {
    pars_out$beta_z_lin <- true_params[["beta_z_lin"]]
  }
  if (reg_types[["u-linear-regressors"]]) {
    pars_out$beta_u_lin <- lapply(true_params[["beta_u_lin"]], `[`, i = , j = n)
    pars_out$vcm_u_lin  <- true_params[["vcm_u_lin"]]
  }
  if (!is.null(DD)) pars_out <- get_dd_slice(pars_out, DD)
  return(pars_out)
}
get_special_pars <- function(true_params, n, DD = NULL,
                             reg_types, special_type) {
  stopifnot(`Wrong arg. to 'special_type':` = special_type %in% c("A", "B"))
  # stopifnot(`Wrong arg. to 'special_type':` = special_type %in% c("A", "B", "AB"))
  # if (special_type == "AB") special_type <- c("A", "B")
  pars_out <- list(sig_sq = true_params[["sig_sq"]][, n, special_type, drop = FALSE],
                   phi = lapply(true_params[["phi"]][[special_type]],
                                `[`, i = , j = n))
  if (reg_types[["z-linear-regressors"]]) {
    pars_out$beta_z_lin <- true_params[["beta_z_lin"]][[special_type]]
  }
  if (reg_types[["u-linear-regressors"]]) {
    pars_out$beta_u_lin <- lapply(true_params[["beta_u_lin"]][[special_type]],
                                  `[`, i = , j = n)
    pars_out$vcm_u_lin  <- true_params[["vcm_u_lin"]][[special_type]]
  }
  if (!is.null(DD)) pars_out <- get_dd_slice(pars_out, DD, special_type)
  return(pars_out)
}
get_dd_slice <- function(tmp_pars_out, DD, special_type = NULL) {
  pars_out <- tmp_pars_out
  if (is.null(DD)) return(pars_out)
  if (!is.null(special_type)) {
    pars_out$sig_sq <- pars_out$sig_sq[DD, , special_type]
  } else {
    pars_out$sig_sq <- pars_out$sig_sq[DD, ]
  }
  pars_out$phi        <- pars_out$phi[[DD]]
  pars_out$vcm_u_lin  <- pars_out$vcm_u_lin[[DD]]
  pars_out$beta_z_lin <- pars_out$beta_z_lin[[DD]]
  pars_out$beta_u_lin <- pars_out$beta_u_lin[[DD]]
  return(pars_out)
}
get_par_name <- function(tmp_pars_out, name) {
  if (is.null(name)) return(tmp_pars_out)
  tmp_pars_out[[name]]
}
