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
new_trueParams <- function(dim_model,
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
                   missing(dim_model))
  if (check_args) stop("Missing arguments")
  if (SIMUL_U_BETA && missing(num_u_regs)) stop("Number of REs not specified!")
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
  true_sig_sq <- new_sig_sq_x(sig_sq, DD, NN, options$dwn_scl)
  true_phi    <- new_phi(SIMUL_PHI, phi, DD, NN, order_p_vec)
  true_bet_z  <- new_bet_z(SIMUL_Z_BETA, beta_z_lin,
                                      DD, num_z_regs,
                                      options$intercepts$at_z)
  tmp_u       <- new_bet_vcm_u(SIMUL_U_BETA, DD, NN,
                               num_u_regs, seed_taken,
                               options$intercepts$at_u)
  true_bet_u <- tmp_u[["true_bet_u"]]
  true_D0u_u <- tmp_u[["true_D0u_u"]]

  par_trues <- list(sig_sq = true_sig_sq,
                    phi = true_phi,
                    beta_z_lin = true_bet_z,
                    beta_u_lin = true_bet_u,
                    vcm_u_lin = true_D0u_u)
  class(par_trues) <- "trueParams"
  attr(par_trues, which = "meta_info") <- list(MODEL_DIM = dim_model,
                                               PAR_SETTINGS = settings_pars,
                                               SEED_NO = seed_taken)
  return(true_params = par_trues)
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
new_phi <- function(SIMUL_PHI, phi, DD, NN, order_p_vec) {
  if (SIMUL_PHI) {
    if (length(order_p_vec) == 1) {
      order_p_vec <- rep(order_p_vec, times = DD)
    } else if (length(order_p_vec) != DD) {
      stop("Autoregressive order_p must be either length 1 or DD!")
    }
    if (!is.null(phi)) {
      if (length(phi) != DD) {
        stop("Phi params must be passed as a list of length DD!")
      }
      out_phi <- vector("list", DD)
      for (d in 1:DD) {
        check_phi <- all(phi[[d]] >= 0) && all(phi[[d]] < 1)
        if (!check_phi) {
          stop("phi > 0 or phi < 1 or not a vector of length DD ...\n")
        }
        out_phi[[d]] <- matrix(phi[[d]], nrow = order_p_vec[d], ncol = NN)
      }
    } else {
      out_phi  <- get_default_phi(DD, NN, order_p_vec)
    }
  } else {
    out_phi <- NULL
  }
  structure(out_phi,
            class = c("true_phi", list))
}
#' Sets true values (default or user supplied) for parameter sig_sq_x
#'
#' @inheritParams new_trueParams
#' @inheritParams new_phi
#'
#' @return a matrix of dimension \code{DD x NN} of true parameter values for
#'   sig_sq_x
new_sig_sq_x <- function(sig_sq, DD, NN, dwn_scl) {
  if (!is.null(sig_sq)) {
    check_sig_sq <- all(sig_sq > 0) && all(length(sig_sq) == DD)
    if (!check_sig_sq) stop("'sig_sq' > 0 and a vector of length DD ...\n")
    out_sig_sq <- matrix(sig_sq, nrow = DD, ncol = NN)
  } else {
    tmp_sig_sq  <- get_default_sig_sq(DD, dwn_scl)
    out_sig_sq <- matrix(tmp_sig_sq, nrow = DD, ncol = NN)
  }
  rownames(out_sig_sq) <- paste0("D", 1:DD)
  colnames(out_sig_sq) <- paste0("N", 1:NN)
  structure(out_sig_sq,
            class = c("true_sig_sq", "matrix", "array"))
}
#' Sets true values (default or user supplied) for parameter bet_z
#'
#' @inheritParams new_trueParams
#' @inheritParams new_phi
#'
#' @return a list of length \code{DD} each element being a vector of the
#'   corresponding number of regressors
new_bet_z <- function(SIMUL_Z_BETA, beta_z_lin, DD, num_z_regs,
                                 intercepts) {
  if (SIMUL_Z_BETA) {
    if (!is.null(beta_z_lin)) {
      check_bet_z <- is.list(beta_z_lin) && length(beta_z_lin) == DD
      if (!check_bet_z) stop("'beta_z_lin' not a list or not of length = DD...")
      out_bet_z <- beta_z_lin
    } else {
      if (is.null(num_z_regs)) stop("Num. of 'beta_z_lin' regressors required.")
      out_bet_z <- get_default_beta_z_lin(DD, num = num_z_regs, intercepts)
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
new_bet_vcm_u <- function(SIMUL_U_BETA, DD, NN,
                          num_u_regs, seed_taken,
                          intercepts) {
  if (SIMUL_U_BETA) {
    if (is.null(num_u_regs)) stop("Specify number of U-type regressors.")
    num_reg_seq <- if (length(num_u_regs) == 1) {
      num_reg_seq <- rep(num_u_regs, times = DD)
    } else if (length(num_u_regs) == DD) {
      num_reg_seq <- num_u_regs
    } else {
      stop(paste0("Number of random effects must be length, in which case it",
                  "gets recycled, or DD!"))
    }
    true_out_u <- BNMPD::generate_bet_u(DD, NN, TRUE, num_reg_seq,
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
#' Get default values for true \code{sig_sq} parameters in simulation
#'
#' @param DD number of multivariate components
#' @param dwn_scl control parameter that scales variance upwards/downwards
#'
#' @return a vector of length \code{DD}
get_default_sig_sq <- function(DD, dwn_scl) { # nolint: object_name_linter.
  str_scl <- 3.1 / dwn_scl
  add_scl <- 0.1 / dwn_scl

  # str_scl <- 2.1 / dwn_scl # nolint: commented_code_linter.
  # add_scl <- 1.1 / dwn_scl # nolint: object_name_linter.
  return(str_scl + add_scl * 0:(DD - 1))
}
#' Get default values for true \code{phi} parameters in simulation
#'
#' @inheritParams new_trueParams
#' @inheritParams new_phi
#'
#' @return a list of length \code{DD} with elements being matrices of dimension
#'   \code{order_p_vec[d] x NN} containing the true parameter values for phi
get_default_phi <- function(DD, NN, order_p_vec) {
  if (any(order_p_vec > 4)) stop("Need more default phis ...")
  # possible_phis <- matrix(c(0.15, 0.45, 0.25, 0.35,
  #                           0.25, 0.20, 0.35, 0.25,
  #                           0.35, 0.20, 0.20, 0.15,
  #                           0.15, 0.10, 0.10, 0.15),
  #                         nrow = 4, ncol = 4)
  possible_phis <- matrix(c(0.50, 0.50, 0.50, 0.50,
                            0.50, 0.50, 0.50, 0.50,
                            0.50, 0.50, 0.50, 0.50,
                            0.50, 0.50, 0.50, 0.50),
                          nrow = 4, ncol = 4)
  possible_phis <- rbind(possible_phis, possible_phis, possible_phis)
  out_phi <- vector("list", DD)
  for (d in seq_len(DD)) {
    tmp_phis <- possible_phis[d, 1:order_p_vec[d]]
    tmp_phis <- matrix(tmp_phis, nrow = order_p_vec[d], ncol = NN)
    rownames(tmp_phis) <- paste0("p", 1:order_p_vec[d])
    colnames(tmp_phis) <- paste0("NN", 1:NN)
    out_phi[[d]] <- tmp_phis
  }
  return(out_phi)
}
#' Get default values for true \code{beta_z_lin} parameters in simulation
#'
#' @inheritParams get_default_sig_sq
#' @param num inter giving number of regressor components (if
#'   \code{length(num) == 1}, then each code{beta_z_lin} gets the same number;
#'   else must be a vector satisfying \code{length(num) == DD})
#'
#' @return a list of length \code{DD}, with number of elements = \code{num}, or,
#'   if \code{num} is a vector of length \code{DD}, \code{num[d]}, for
#'   \code{d=1,...,DD}
get_default_beta_z_lin <- function(DD, num, intercepts) {
  num_reg_len <- length(num)
  tmp1 <- c(0.325, -0.44)  # tmp values large, first component negative
  tmp2 <- c(0.327, -0.435)    # tmp values large, first component positive
  tmp3 <- c(0.33, -0.43)   # tmp values small, first component negative
  tmp4 <- c(0.332, -0.415) # tmp values small, first component positive
  if (num_reg_len == 1) {
    tmp_neg_pos_large <- rep(tmp1, length.out = num)
    tmp_pos_neg_large <- rep(tmp2, length.out = num)
    tmp_neg_pos_small <- rep(tmp3, length.out = num)
    tmp_pos_neg_small <- rep(tmp4, length.out = num)

  } else if (num_reg_len == DD) {
    tmp_neg_pos_large <- rep(tmp1, length.out = num[1])
    tmp_pos_neg_large <- rep(tmp2, length.out = num[2])
    tmp_neg_pos_small <- rep(tmp3, length.out = num[3])
    tmp_pos_neg_small <- rep(tmp4, length.out = num[4])
  } else {
    stop("If 'num_z_regs' is a vector, it must be of length 'DD'...")
  }
  list_vals <- list(tmp_neg_pos_large,
                    tmp_pos_neg_large,
                    tmp_neg_pos_small,
                    tmp_pos_neg_small)
  list_vals <- rep(list_vals, length.out = DD)
  list_vals <- scale_up_intercept(list_vals, 100, intercepts)
  return(list_vals)
}
scale_up_intercept <- function(vals_list, scl_up, intercept_ids) {
  DD <- length(intercept_ids)
  scl <- rep(1, times = DD)
  scl[intercept_ids] <- scl_up
  for (d in 1:DD) {
    vals_list[[d]][1] <- vals_list[[d]][1] * scl[d]
  }
  return(vals_list)
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
                           from_IW,
                           num_re,
                           vcm_u_scl = 0.03,
                           rel_var_to_cov = 5,
                           n0u = 50, # n0u <- num_re + 1
                           seed_no = 42) {
  true_bet_u <- vector("list", DD)
  if (from_IW) {
    stopifnot(is.numeric(num_re) && (length(num_re) == DD))
    if (!is.null(seed_no))  set.seed(seed_no)
    D0u <- vector("list", DD)
    for (d in 1:DD) {
      D0u[[d]] <- get_hyper_prior_vcm(vcm_u_scl, rel_var_to_cov, num_re[d])
      D0u[[d]] <- solveme((stats::rWishart(1, n0u, D0u[[d]]))[, , 1])
      colnames(D0u[[d]]) <- paste0("U", 1:num_re[d])
      rownames(D0u[[d]]) <- paste0("U", 1:num_re[d])
      true_bet_u[[d]] <- matrix(0, nrow = num_re[d], ncol = NN)
      colnames(true_bet_u[[d]]) <- paste0("N", 1:NN)
      rownames(true_bet_u[[d]]) <- paste0("U", 1:num_re[d])
      for (n in 1:NN) {
        true_bet_u[[d]][, n] <- MASS::mvrnorm(n = 1,
                                              # mu = c(2, rep(0.5, times = num_re[d] - 1)),
                                              mu = c(1.25, rep(0.75, times = num_re[d] - 1)),
                                              Sigma = D0u[[d]])
      }
    }
    return(list(true_bet_u = true_bet_u,
                true_vcm_u = D0u))
  } else {
    warning("Not simulating RE with VCM distributed as ~IW().")
    stopifnot(is.numeric(num_re) && (length(num_re) == 1))
    vals <- matrix(1:(num_re*DD)*c(-1, 1), nrow = num_re*DD, ncol = NN)
    for (d in 1:DD) {
      true_bet_u[[d]] <- vals[1:num_re + num_re*(d - 1), , drop = FALSE]
    }
    return(true_bet_u = true_bet_u)
  }
}
get_hyper_prior_vcm <- function(vcm_u_scl, rel_var_to_cov, num_re) {
  # vcm_u_scl_adj <- vcm_u_scl * (1 + 1/num_re)
  cov_entry <- -1.0
  var_entry <- abs(cov_entry) * rel_var_to_cov
  out_mat <- matrix(cov_entry, nrow = num_re, ncol = num_re)
  diag(out_mat) <- rep(var_entry, times = num_re)
  out_mat <- out_mat * vcm_u_scl
  # solveme(out_mat)
  out_mat
}
#' Write a yaml config file for a BNMPD-model class
#'
#' @inheritParams get_file_name_main
#' @param model_type a named or unnamed vector of two characters specifying the
#'   model type observations (e.g. 'DIRICHLET") and the model type latent states
#'   i.e. "lin_re" giving the specification of linear and random effect
#'   covariates in the latent state process
#' @param pth_to_yaml path to yaml file to write to (including file-name and
#'   '.yaml' extension)
#'
#' @return side effect function that writes a 'yaml'-config file to path
#' @export
generate_yaml_model_defintion <- function(model_type,
                                          dim_model,
                                          par_settings,
                                          pth_to_yaml) {
  SIMUL_PHI    <-par_settings[["SIMUL_PHI"]]
  SIMUL_Z_BETA <-par_settings[["SIMUL_Z_BETA"]]
  SIMUL_U_BETA <-par_settings[["SIMUL_U_BETA"]]
  num_z_regs   <-par_settings[["num_z_regs"]]
  num_u_regs   <-par_settings[["num_u_regs"]]
  order_p      <-par_settings[["order_p"]]
  check_args_to_generate_yaml_json(dim_model, par_settings)
  tmp_nms <- names(dim_model)
  dim_model <- as.integer(dim_model)
  names(dim_model) <- tmp_nms
  out_list <- list(model_type_obs = model_type_to_list(model_type)[1],
                   model_type_lat = model_type_to_list(model_type)[2],
                   dimension = dimensions_to_list(dim_model),
                   cross_section_used = cross_section_to_list(dim_model),
                   time_series_used = time_series_to_list(dim_model))
  for (d in seq_len(dim_model[["DD"]])) {
    # browser()
    name_list_elem <- paste0("D", ifelse(d < 10, paste0(0, d), d))
    out_list[[name_list_elem]] <- DD_to_list(par_settings, d)
  }

  yaml::write_yaml(out_list, file = pth_to_yaml)
}
check_args_to_generate_yaml_json <- function(dim = NULL, regspc = NULL) {
  if (!is.null(dim)) {
    check_dim <- !(is.null(names(dim)))
    stopifnot(`Arg. 'dim' is unnamed; use NN, TT, DD as names.` = check_dim)
    check_dim <- all(names(dim) %in% c("NN", "TT", "DD"))
    stopifnot(`Unknown names of arg. 'dimensions' ... ` = check_dim)
    check_dim <- length(dim) == 3
    stopifnot(`Arg. 'dimensions' not of length 3.` = check_dim)
  }
  if (!is.null(regspc)) {
    stopifnot(`Argument 'par_settings' is not a list ` = is.list(regspc))
    check_names <- all(unique(names(regspc)) %in% c("SIMUL_PHI",
                                                    "SIMUL_Z_BETA",
                                                    "SIMUL_U_BETA",
                                                    "num_z_regs",
                                                    "num_u_regs",
                                                    "order_p_vec"))
    stopifnot(`Arg. 'par_settings' has wrong names` = check_names)
    check_types <- all(unlist(lapply(regspc, typeof)) %in% c("logical",
                                                             "logical",
                                                             "logical",
                                                             "double",
                                                             "double",
                                                             "double"))
    stopifnot(`Arg. 'par_settings' has wrong element types.` = check_types)
  }
}
model_type_to_list <- function(mod) {
  if (is.null(names(mod))) {
    msg <- paste0("Arg. model_type is unnamed: ",
                  "using ", mod[1], " as 'model_type_obs' and ",
                  mod[2], " as 'model_type_lat'.")
    warning(msg)
    return(c(mod[1], mod[2]))
  } else {
    stopifnot(`Unknown model names` = all(names(mod) %in% c("model_type_obs",
                                                            "model_type_lat")))
  }
  c(mod["model_type_obs"], mod["model_type_lat"])
}
dimensions_to_list <- function(dim) {
  list(num_cross_section = dim[["NN"]],
       num_time_periods = dim[["TT"]],
       num_mult_comp = dim[["DD"]])
}
cross_section_to_list <- function(dim) {
  tmp_cs_var_val <- paste0("paste0('cs_', 1:", dim[["NN"]], ")")
  attr(tmp_cs_var_val, which = "quoted") <- TRUE
  tmp_cs_var_lab <- paste0("paste0('cs_', 1:", dim[["NN"]], ")")
  attr(tmp_cs_var_lab, which = "quoted") <- TRUE
  list(cs_name_var = "CS",
       cs_name_lab = "cross_section",
       cs_var_val =  tmp_cs_var_val,
       cs_var_lab =  tmp_cs_var_lab)
}
time_series_to_list <- function(dim) {
  tmp_ts_var_val <- paste0("seq(from = 1L, to = ", dim[["TT"]], "L, by = 1L)")
  attr(tmp_ts_var_val, which = "quoted") <- TRUE
  tmp_ts_var_lab <- paste0("as.character(seq(from = 1L, to = ",
                           dim[["TT"]], "L, by = 1L))")
  attr(tmp_ts_var_lab, which = "quoted") <- TRUE
  list(ts_name_var = "TS",
       ts_name_lab = "time_series",
       ts_var_val =  tmp_ts_var_val,
       ts_var_lab =  tmp_ts_var_lab)
}
DD_to_list <- function(regspc, d_tmp) {
  out_list <- list(y_var = paste0("Y", d_tmp),
                   y_lab = paste0("Y", d_tmp))
  if (regspc$SIMUL_Z_BETA) {
    out_list$z_reg <-  get_reg_to_sublist(regspc$num_z_regs, d_tmp, "Z")
  }
  if (regspc$SIMUL_U_BETA) {
    out_list$u_reg <-  get_reg_to_sublist(regspc$num_u_regs, d_tmp, "U")
  }
  return(out_list)
}
get_reg_to_sublist <- function(num, d_tmp, name) {
  list(lab = paste0(name, "_", seq_len(num), "_", d_tmp),
       var = paste0(name, "_", seq_len(num), "_", d_tmp))
}
#' Generate config file \code{setup_inits.json} from simulated data
#'
#' Parses parameter object used to generate simulated data (generated by
#' [new_trueParams()]).
#'
#' @param params_used output generated by [new_trueParams()] or any other
#'   list type object that follows the structure of instances from this class
#' @param pth_to_json character giving the path to write to the resulting
#'   \code{setup_inits.json}-file
#'
#' @return pure side-effect function that generates the \code{json}-file
#' @export
generate_setup_init_json <- function(params_used, pth_to_json) {
  params_used$sig_sq <- params_used$sig_sq[, 1]
  DD <- length(params_used$sig_sq)

  list_json <- list()
  for(d in seq_len(DD)) {
    name_list_elem <- paste0("D", ifelse(d < 10, paste0(0, d), d))

    par_avail <- get_par_avail(params_used, num_d = d)
    par_avail_names <- names(par_avail)

    par_to_list <- list()
    for (j in par_avail_names) {
      par_to_list[[j]] <- par_to_list(par_avail[[j]], j, d)
    }
    list_json[[name_list_elem]] <- par_to_list
  }
  jsonlite::write_json(list_json, pth_to_json, digits = 8, pretty = TRUE)
}
par_to_list <- function(par_val, par_name, num_d) {
  if (par_name == "phi") {
    num_comp <- seq_len(length(par_val))
    lab <- paste0("phi_{", num_comp, ",", num_d,"}")
    var <- paste0(par_name, "_", num_comp, "_", num_d)
    val <- par_val
  }
  if (par_name == "beta_z_lin") {
    num_comp <- seq_along(par_val)
    lab <- paste0("beta_{z", num_comp, ",", num_d, "}^{lin}")
    var <- paste0("Z_", num_comp, "_", num_d)
    val <- par_val
  }
  if (par_name == "beta_u_lin") {
    if (is.null(nrow(par_val))) {
      num_comp <- 1
    } else {
      num_comp <- seq_len(nrow(par_val))
    }
    lab <- paste0("beta_{u", num_comp, ",", num_d, "}^{lin}")
    var <- paste0("U_", num_comp, "_", num_d)
    val <- par_val
  }
  if (par_name == "vcm_u_lin") {
    lab <- paste0("vcm_beta_u{", num_d,"}")
    var <- paste0("buVCM", num_d)
    val <- par_val
  }
  if (par_name == "sig_sq") {
    lab <- paste0("sigma_{", num_d,"}")
    var <- paste0(par_name, num_d)
    val <- par_val
  }
  list(lab =  lab,
       var =  var,
       val =  val)
}
get_par <- function(par_to_check, par_name, num_d) {
  if (par_name == "sig_sq") {
    if (!is.null(par_to_check)) return(par_to_check[num_d])
  } else if (par_name == "phi") {
    if (!is.null(par_to_check)) return(par_to_check[[num_d]][, 1])
  } else if (par_name %in% c("beta_z_lin", "beta_u_lin", "vcm_u_lin")) {
    if (!is.null(par_to_check)) return(par_to_check[[num_d]])
  } else {
    stop("Unknown parameter name.")
  }
}
get_par_avail <- function(params_used, num_d) {
  check_pars <- c("phi", "beta_z_lin", "beta_u_lin", "sig_sq", "vcm_u_lin")
  num_pars <- length(check_pars)
  id_avail  <- NULL
  par_avail <- vector("list", num_pars)
  for (j in seq_len(num_pars)) {
    if(isTRUE(!is.null((params_used[[check_pars[j]]])))) {
      par_avail[[j]] <- get_par(params_used[[check_pars[j]]],
                                check_pars[j], num_d)
      id_avail <- c(id_avail, j)
    }
  }
  par_avail <- par_avail[id_avail]
  names(par_avail) <- check_pars[id_avail]
  return(par_avail)
}
update_settings_project_yaml <- function(pth, proj_no, proj_name, notes) {
  tmp_proj_yaml <- yaml::read_yaml(pth)
  if (!is.null(proj_no)) {
    tmp_proj_yaml$project_name <- proj_no
  }
  if (!is.null(proj_name)) {
    tmp_proj_yaml$project_name <- proj_name
  }
  if (!is.null(notes)) {
    tmp_proj_yaml$project_notes <- notes
  }
  yaml::write_yaml(tmp_proj_yaml, file = pth)
  return(invisible(pth))
}
