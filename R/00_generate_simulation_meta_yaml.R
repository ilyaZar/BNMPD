#' Write a yaml config file for a BNMPD-model class
#'
#' @inheritParams new_dataSim
#' @inheritParams set_new_project_name
#' @param pth_project path to top-level project directory
#' @param model_type concatenated from [get_type_obs()], [get_type_lat()]:
#'    `c(model_type_obs = get_type_obs(dataSim),
#'       model_type_lat = get_type_lat(dataSim))`
#'
#' @return side effect function that writes a 'yaml'-config file to path
#' @export
generate_yaml_model_defintion <- function(true_params,
                                          pth_project,
                                          model_type) {

  dim_model    <- get_dimension(true_params, "all") %>%
    set_dim_model_writeable()
  par_settings <- get_par_settings(true_params)
  distribution <- get_distribution(true_params)
  check_args_to_generate_yaml_json(dim_model, par_settings)

  pth_to_yaml <- file.path(pth_project, "model",
                           "model-definition", "model_definition.yaml")
  SIMUL_PHI    <- par_settings[["SIMUL_PHI"]]
  SIMUL_Z_BETA <- par_settings[["SIMUL_Z_BETA"]]
  SIMUL_U_BETA <- par_settings[["SIMUL_U_BETA"]]
  num_z_regs   <- par_settings[["num_z_regs"]]
  num_u_regs   <- par_settings[["num_u_regs"]]
  order_p      <- par_settings[["order_p"]]

  out_list <- list(model_type_obs = model_type_to_list(model_type)[1],
                   model_type_lat = model_type_to_list(model_type)[2],
                   dimension = dimensions_to_list(dim_model),
                   cross_section_used = cross_section_to_list(dim_model),
                   time_series_used = time_series_to_list(dim_model))

  out_list <- multi_DD_to_list(out_list, distribution, dim_model, par_settings)

  yaml::write_yaml(out_list, file = pth_to_yaml)

  return(invisible(true_params))
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
    check_dim <- all(names(dim) %in% c("NN", "TT", "DD", "DD2"))
    stopifnot(`Unknown names of arg. 'dimensions' ... ` = check_dim)
    check_dim <- length(dim) %in% c(3, 4)
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
multi_DD_to_list <- function(list_yaml, dist, dim_model, par_settings) {
  DD <- dim_model[["DD"]]
  if (dist %in% c("gen_dirichlet", "gen_dirichlet_mult")) {
    list_yaml <- get_DD_to_list(list_yaml, par_settings, DD, "special")
  } else {
    list_yaml <- get_DD_to_list(list_yaml, par_settings, DD, "default")
  }
  return(list_yaml)
}
get_DD_to_list <- function(list_out, par_settings, DD, type) {
  for (d in seq_len(DD)) {
    dd_elem             <- get_name_list_elem(d)
    if (type == "special") {
      list_out[[dd_elem]] <- get_y_to_sublist(d)
      if (d == DD) break;
      ad <- ifelse(d < 10, paste0(0, d), d)
      eA <- paste0("DA_", ad)
      eB <- paste0("DB_", ad)
      list_out[[dd_elem]][[eA]] <- DD_to_list(par_settings, "A", d)
      list_out[[dd_elem]][[eB]] <- DD_to_list(par_settings, "B", d)
    } else if (type == "default") {
      list_out[[dd_elem]] <- c(get_y_to_sublist(d),
                               DD_to_list(par_settings, "", d))
    } else {
      stop("Unknown value for argument 'type'; either special or default")
    }
  }
  return(list_out)
}
DD_to_list <- function(regspc, type, d_tmp) {
  out_list <- list()
  if (regspc$SIMUL_Z_BETA) {
    out_list$z_reg <-  get_reg_to_sublist(regspc$num_z_regs, type, d_tmp, "Z")
  }
  if (regspc$SIMUL_U_BETA) {
    out_list$u_reg <-  get_reg_to_sublist(regspc$num_u_regs, type, d_tmp, "U")
  }
  return(out_list)
}
get_name_list_elem <- function(d_tmp) {
  paste0("D", ifelse(d_tmp < 10, paste0(0, d_tmp), d_tmp))
}
get_reg_to_sublist <- function(num, type, d_tmp, name) {
  list(lab = paste0(name, "_", seq_len(num), "_", type, d_tmp),
       var = paste0(name, "_", seq_len(num), "_", type, d_tmp))
}
get_y_to_sublist <- function(d_tmp) {
   list(y_var = paste0("Y", d_tmp),
        y_lab = paste0("Y", d_tmp))
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
