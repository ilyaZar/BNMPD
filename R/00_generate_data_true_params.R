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
#' @param beta_z_lin a list of length \code{DD} with (possibly) varying number of
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
                                 beta_z_lin,
                                 SIMUL_Z_BETA,
                                 SIMUL_U_BETA,
                                 NUM_BETA_U,
                                 seed_taken) {
  check_args <- (missing(sig_sq) || missing(phi) || missing(beta_z_lin) ||
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
  check_phi <- all(phi >= 0) && all(phi < 1) && all(length(phi) == DD)
  if (!check_phi) stop("phi > 0 or phi < 1 or not a vector of length DD ...\n")
  true_phi <- matrix(phi, nrow = DD, ncol = NN)
  #
  if (SIMUL_Z_BETA) {
    check_bet_z <- is.list(beta_z_lin) && length(beta_z_lin) == DD
    if (!check_bet_z) stop("'beta_z_lin' not a list or not of length = DD ...\n")
    true_bet_z <- beta_z_lin
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
                    beta_z_lin = true_bet_z,
                    beta_u_lin = true_bet_u,
                    vcm_u_lin = true_D0u_u)
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
                           num_re,
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
      D0u[[d]] <- solve(diag(1:num_re[d]*100/10*(0.01 + 1/num_re[d]),
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
#' Write a yaml config file for a BNMPD-model class
#'
#' @param model_type a named or unnamed vector of two character specifying the
#'   model type observations (e.g. 'DIRICHLET") and the model type latent states
#'   i.e. "lin_reg" giving the specification of linear and random effect
#'   covariates in the latent state process
#' @param dimensions a named numeric vector of three components: 'NN', 'TT',
#'   'DD'
#' @param regressor_specs a list of four:
#'   \itemize{
#'   \item{\code{SIMUL_Z_BETA: }}{logical; if \code{TRUE}, then Z-type regressors
#'    are written to yaml-config}
#'   \item{\code{SIMUL_U_BETA: }}{logical; if \code{TRUE}, then U-type regressors
#'    are written to yaml-config}
#'   \item{\code{num_z_regs: }}{numeric giving the number of z-regressors}
#'   \item{\code{num_u_regs: }}{numeric giving the number of u-regressors}
#'    }
#' @param pth_to_yaml path to yaml file to write to (including file-name and
#'   '.yaml' extension)
#'
#' @return side effect function that writes a 'yaml'-config file to path
#' @export
generate_yaml_model_defintion <- function(model_type,
                                          dimensions,
                                          regressor_specs,
                                          pth_to_yaml) {
  check_args_to_generate_yaml_json(dimensions, regressor_specs)
  tmp_nms <- names(dimensions)
  dimensions <- as.integer(dimensions)
  names(dimensions) <- tmp_nms
  out_list <- list(model_type_obs = model_type_to_list(model_type)[1],
                   model_type_lat = model_type_to_list(model_type)[2],
                   dimension = dimensions_to_list(dimensions),
                   cross_section_used = cross_section_to_list(dimensions),
                   time_series_used = time_series_to_list(dimensions))
  for (d in seq_len(dimensions[["DD"]])) {
    # browser()
    name_list_elem <- paste0("D", ifelse(d<10, paste0(0, d), d))
    out_list[[name_list_elem]] <- DD_to_list(regressor_specs, d)
  }

  yaml::write_yaml(out_list, file = pth_to_yaml)
}
check_args_to_generate_yaml_json <- function(dim = NULL, regspc = NULL) {
  if (!is.null(dim)) {
    check_dim <- !(is.null(names(dim)))
    stopifnot(`Arg. 'dimensions' is unnamed; use NN, TT, DD as names.` = check_dim)
    check_dim <- all(names(dim) %in% c("NN", "TT", "DD"))
    stopifnot(`Unknown names of arg. 'dimensions' ... ` = check_dim)
    check_dim <- length(dim) == 3
    stopifnot(`Arg. 'dimensions' not of length 3.` = check_dim)
  }
  if (!is.null(regspc)) {
    stopifnot(`Argument 'regressor_specs' is not a list ` = is.list(regspc))
    check_names <- all(unique(names(regspc)) %in% c("SIMUL_Z_BETA",
                                                    "SIMUL_U_BETA",
                                                    "num_z_regs",
                                                    "num_u_regs"))
    stopifnot(`Arg. 'regressor_specs' has wrong names` = check_names)
    check_types <- all(unlist(lapply(regspc, typeof)) %in% c("logical",
                                                             "logical",
                                                             "double",
                                                             "double" ))
    stopifnot(`Arg. 'regressor_specs' has wrong element types.` = check_types)
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
#' [generate_true_params()]).
#'
#' @inheritParams generate_yaml_model_defintion
#' @param true_params output generated by [generate_true_params()]
#' @param pth_to_json character giving the path to write to the resulting
#'   \code{setup_inits.json}-file
#'
#' @return pure side-effect function that generates the \code{json}-file
#' @export
generate_setup_init_json <- function(dimension, true_params, pth_to_json) {
  check_args_to_generate_yaml_json(dimension)
  DD <- dimension[["DD"]]

  true_params$phi    <- true_params$phi[, 1]
  true_params$sig_sq <- true_params$sig_sq[, 1]
  list_json <- list()
  for(d in seq_len(DD)) {
    name_list_elem <- paste0("D", ifelse(d < 10, paste0(0, d), d))

    par_avail <- get_par_avail(true_params, num_d = d)
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
    lab <- paste0("phi_{", num_d,"}")
    var <- paste0(par_name, num_d)
    val <- par_val
  }
  if (par_name == "beta_z_lin") {
    num_comp <- seq_len(length(par_val))
    lab <-  paste0(paste0("beta_{z", num_comp, ","), num_d,"}^{lin}")
    var <- paste0(paste0("Z_", num_comp), "_", num_d)
    val <- par_val
  }
  if (par_name == "beta_u_lin") {
    if (is.matrix(par_val)) {
      num_comp <- seq_len(nrow(par_val))
    } else {
      num_comp <- seq_len(length(par_val))
    }
    lab <- paste0(paste0("beta_{u", num_comp, ","), num_d,"}^{lin}")
    var <- paste0(paste0("U_", num_comp), "_", num_d)
    val <- par_val
  }
  if (par_name == "sig_sq") {
    lab <- paste0("sigma_{", num_d,"}")
    var <- paste0(par_name, num_d)
    val <- par_val
  }
  if (par_name == "vcm_u_lin") {
    lab <- paste0("beta_{u-vcm", num_d,"}")
    var <- paste0("buVCM", num_d)
    val <- as.vector(par_val)
  }
  list(lab =  lab,
       var =  var,
       val =  val)
}
get_par <- function(par_to_check, par_name, num_d) {
  if (par_name %in% c("sig_sq", "phi")) {
    tmp_val <- par_to_check[num_d]
    return(tmp_val)
    # if(is.na(tmp_val)) return(NA_real_) or NLL; else return(tmp_val)
  } else if (par_name %in% c("beta_z_lin", "beta_u_lin", "vcm_u_lin")) {
    tmp_val <- par_to_check[[num_d]]
    # if(is.null(tmp_val)) return(tmp_val); else return(tmp_val)
    return(tmp_val)
  } else {
    stop("Unknown parameter name.")
  }
}
get_par_avail <- function(true_params, num_d) {
  check_pars <- c("phi", "beta_z_lin", "beta_u_lin", "vcm_u_lin", "sig_sq")
  id_avail  <- NULL
  par_avail <- list()
  for (j in seq_along(check_pars)) {
    par_avail[[j]] <- get_par(true_params[[check_pars[j]]],
                              check_pars[j],
                              num_d)
    if(!is.null(par_avail[[j]])) id_avail <- c(id_avail, j)
  }
  names(par_avail) <- check_pars[id_avail]
  return(par_avail)
}
