#' Write a yaml config file for a BNMPD-model class
#'
#' @inheritParams new_dataSim
#' @inheritParams set_new_project_name
#' @param pth_project path to top-level project directory
#'
#' @return side effect function that writes a 'yaml'-config file to path
#' @export
generate_yaml_model_defintion <- function(true_params,
                                          pth_project,
                                          model_type) {

  dim_model    <- get_dimension(true_params, "all")
  par_settings <- get_par_settings(true_params)

  pth_to_yaml <- file.path(pth_project, "model",
                           "model-definition", "model_definition.yaml")
  SIMUL_PHI    <- par_settings[["SIMUL_PHI"]]
  SIMUL_Z_BETA <- par_settings[["SIMUL_Z_BETA"]]
  SIMUL_U_BETA <- par_settings[["SIMUL_U_BETA"]]
  num_z_regs   <- par_settings[["num_z_regs"]]
  num_u_regs   <- par_settings[["num_u_regs"]]
  order_p      <- par_settings[["order_p"]]
  check_args_to_generate_yaml_json(dim_model, par_settings)
  dim_model <- set_dim_model_writeable(dim_model)
  out_list <- list(model_type_obs = model_type_to_list(model_type)[1],
                   model_type_lat = model_type_to_list(model_type)[2],
                   dimension = dimensions_to_list(dim_model),
                   cross_section_used = cross_section_to_list(dim_model),
                   time_series_used = time_series_to_list(dim_model))
  for (d in seq_len(dim_model[["DD"]])) {
    name_list_elem <- paste0("D", ifelse(d < 10, paste0(0, d), d))
    out_list[[name_list_elem]] <- DD_to_list(par_settings, d)
  }

  yaml::write_yaml(out_list, file = pth_to_yaml)
}
set_dim_model_writeable <- function(dims) {
  tmp  <- names(dims)
  dims <- as.integer(dims)
  names(dims) <- tmp
  return(dims)
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
#'   list type object that follows the structure of instances from this class;
#'   e.g. 'zero_params' as generated by [get_zero_or_defaults()], see body of
#'   [generate_simulation_study()]
#' @inheritParams generate_yaml_model_defintion
#'
#' @return pure side-effect function that generates the \code{json}-file
#' @export
generate_setup_init_json <- function(params_used, pth_project) {
  pth_to_json <- file.path(pth_project, "model",
                           "model-definition",
                           "setup_inits.json")
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
update_settings_project_yaml <- function(pth_top_level,
                                         proj_name,
                                         proj_no = NULL,
                                         notes = NULL) {
  pth <- file.path(pth_top_level, "model", "settings", "settings_project.yaml")
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
  return(invisible(pth_top_level))
}
