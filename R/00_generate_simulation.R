#' Generate directories and files for simulation study
#'
#' @param data_simulation an object of \code{class} "dataSim" as generated by
#'    [new_dataSim()]
#' @param INIT_AT a character: either "trues" or "default" to write
#'   initialization values at true or default values (an identity matrix for the
#'   VCM of bet_u_lin e.g. but the beta_z/beta_u set to 0).
#' @param pth_top_level character giving path to where create the project (top
#'   level directory and not full path with project name, see argument
#'   \code{project_name})
#' @param project_name a named list with pre-pending and appending meta-infos
#'   characterizing the simulation study (e.g. set \code{append = "MCMC_only"})
#' @param overwrite logical; if \code{TRUE} then implied project directory is
#'   overwritten when it exists, otherwise an error is thrown
#' @param TESTING local; if `TRUE` then meta files for model-dir setup are
#'   copied with different relative paths that allows to use this functions
#'   inside tests; `{testthat}` type tests start relative to the `test/testthat`
#'    directory so relative paths need adjustments and `testthat::test()` or
#'    `testthat::test_file()` run smoothly
#'
#' @return side effect function generating directories and files for simulation
#'   study and copying data/template files into corresponding directories
#' @export
generate_simulation_study <- function(data_simulation,
                                      INIT_AT = "trues",
                                      pth_top_level,
                                      project_name = list(prepend = NULL,
                                                          append = NULL),
                                      overwrite = FALSE,
                                      TESTING = FALSE) {
  check_class_data_sim(data_simulation)
  stopifnot(`INIT_AT must be eihter 'true' or 'default'.` =
              INIT_AT %in% c("trues", "default"))
  stopifnot(`'pth_top_level' must be of type character` =
              is.character(pth_top_level))
  stopifnot(`Arg. 'project_name' must be named list` =
              names(project_name) %in% c("prepend", "append"))
  stopifnot(`Arg. 'overwrite' must be logical` = is.logical(overwrite))

  dataSim    <- data_simulation
  trueParams <- get_true_params_obj(dataSim)

  zeroParams <- get_zero_or_defaults(trueParams)
  usedParams <- get_used_or_zeros(trueParams, zeroParams, INIT_AT = INIT_AT)

  model_type <- c(model_type_obs = get_type_obs(dataSim),
                  model_type_lat = get_type_lat(dataSim))

  base_name      <- set_base_name(dataSim)
  project_name   <- set_new_project_name(project_name, base_name,
                                         model_type[["model_type_obs"]])
  pth_to_project <- file.path(pth_top_level, project_name)

  dir_proj_top_level_update(pth_to_project, overwrite = overwrite)
  set_project_dir_structure(pth_to_project)

  generate_yaml_model_defintion(trueParams, pth_to_project, model_type)
  generate_setup_init_json(usedParams, pth_to_project)

  save_simulated_data(pth_to_project, base_name, dataSim,
                      trueParams, zeroParams)

  copy_meta_files(pth_to_project, project_name, TESTING)
  update_settings_project_yaml(pth_to_project, project_name)
  return(invisible(data_simulation))
}
#' Helper that generates zero/default paramters
#'
#' Parameters are of class `trueParams`
#'
#' @inheritParams new_dataSim
#'
#' @return an object of class `trueParams` with entries all zero or some default
#'   values for its data containers
get_zero_or_defaults <- function(true_params) {
  check_class_true_params(true_params)
  zero_params <- true_params
  if (!is.null(true_params[["sig_sq"]])) {
    zero_params <- set_params(zero_params,
                              name_par = "sig_sq",
                              values = 1)
  }
  if (!is.null(true_params[["phi"]])) {
    zero_params <- set_params(zero_params,
                              name_par = "phi",
                              values = 0)
  }
  if (!is.null(true_params[["beta_z_lin"]])) {
    zero_params <- set_params(zero_params,
                              name_par = "beta_z_lin",
                              values = 0)
  }
  if (!is.null(true_params[["beta_u_lin"]])) {
    zero_params <- set_params(zero_params,
                              name_par = "beta_u_lin",
                              values = 0)
  }
  if (!is.null(true_params[["vcm_u_lin"]])) {
    zero_params <- set_params(zero_params,
                              name_par = "vcm_u_lin",
                              values = c(0, 1))
  }
  return(zero_params)
}
#' Generate main part of a file name
#'
#' Which is the 'basename' internally. Takes into consideration the type of
#' simulation, regressor numbers and types, order of the autoregression, seeds
#' etc.
#'
#' @param data_sim an object of class `dataSim`, see [new_dataSim()]
#'
#' @return main part for file names
set_base_name <- function(data_sim) {
  true_params  <- get_true_params_obj(data_sim)
  dim_model    <- get_dimension(true_params, "all")
  par_settings <- get_par_settings(true_params)
  seed_nos_01  <- get_seed(true_params)
  seed_nos_02  <- get_seed(data_sim)

  SIMUL_PHI    <- par_settings[["SIMUL_PHI"]]
  SIMUL_Z_BETA <- par_settings[["SIMUL_Z_BETA"]]
  SIMUL_U_BETA <- par_settings[["SIMUL_U_BETA"]]
  num_z_regs   <- par_settings[["num_z_regs"]]
  num_u_regs   <- par_settings[["num_u_regs"]]
  order_p      <- par_settings[["order_p_vec"]]

  tmp_fn <- paste0("NN", dim_model[1], "_TT", dim_model[2], "_DD", dim_model[3])

  if (SIMUL_PHI) {
    tmp_fn <- paste0(tmp_fn, "_", "withAUTO", paste0(",", unique(order_p),
                                                     collapse = ""))
  } else {
    tmp_fn <- paste0(tmp_fn, "_", "noAUTO")
  }
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
  tmp_fn <- paste0(tmp_fn,
                   "_", "parSEED", seed_nos_01,
                   "_", "simSEED", seed_nos_02)
  return(tmp_fn)
}
dir_proj_top_level_update <- function(pth_top_lvl, overwrite) {
  stopifnot(`Missing args. not permitted` =
              !(missing(pth_top_lvl) && missing(overwrite)))
  if (isFALSE(overwrite)) {
    stopifnot(`Project directory already exists!` = !dir.exists(pth_top_lvl))
  } else if (isTRUE(overwrite)) {
    if (dir.exists(pth_top_lvl)) {
      unlink(pth_top_lvl, recursive = TRUE)
      cat(crayon::magenta(paste0("\n", "Overwriting project structure in: ",
                                 pth_top_lvl, "\n")))
    }
  }
  dir.create(pth_top_lvl, recursive = TRUE)
  return(invisible(pth_top_lvl))
}
#' Generates consistent file names for simulation study
#'
#' @param fn_main_part character giving base part of the name, see
#'   [set_base_name]
#' @param fn_data character giving the filename of the simulated data set
#'   (saved in \code{.csv}-format); default to "sim_data"
#' @param fn_true_states character giving the filename of the simulated true
#'   states (saved in \code{.rds}-format); default to "states_true"
#' @param fn_zero_states  character giving the filename of the states set all to
#'   zero (saved in \code{.rds}-format); default to "states_zero"
#'
#' @return a list of 3:
#'    \itemize{
#'    \item\code{fn_data_set}: file name of the data set
#'    \item\code{fn_data_set}: file name for true state values
#'    \item\code{fn_data_set}: file name for zero state values
#'    }
get_file_names_simul_data <- function(fn_main_part,
                                      fn_data = "sim_data",
                                      fn_true_states = "states_true",
                                      fn_zero_states = "states_zero") {

  fn_data_set <- paste0(fn_data, "_", fn_main_part, ".csv")
  fn_true_val <- paste0(fn_true_states, "_", fn_main_part, ".rds")
  fn_zero_val <- paste0(fn_zero_states, "_", fn_main_part, ".rds")
  return(list(fn_data_set = fn_data_set,
              fn_true_val = fn_true_val,
              fn_zero_val = fn_zero_val))
}
#' Save simulated data and true parameter values used to generate it.
#'
#' @inheritParams generate_yaml_model_defintion
#' @inheritParams set_new_project_name
#' @inheritParams generate_simulation_study
#' @inheritParams generate_yaml_model_defintion
#' @inheritParams new_dataSim
#' @param data_sim an object of class `dataSim`, see [new_dataSim()]
#' @param defl_params default parameter values used to generate data (to
#'   determine container/data sizes); usually set to an instance of 'zeroParams'
#'   as generated by [get_zero_or_defaults()]
#'
#' @return side effect function; saves simulated data and true parameter values
#'   to output path
save_simulated_data <- function(pth_project,
                                base_name,
                                data_sim,
                                true_params,
                                defl_params) {
  pth_to_write <- file.path(pth_project, "model", "input")

  fn_all         <- get_file_names_simul_data(fn_main_part = base_name)
  fn_true_states <- fn_all[["fn_true_val"]]
  fn_zero_states <- fn_all[["fn_zero_val"]]
  fn_data_set    <- file.path("datasets", fn_all[["fn_data_set"]])

  pth_data        <- file.path(pth_to_write, fn_data_set)
  pth_true_states <- file.path(pth_to_write, fn_true_states)
  pth_zero_states <- file.path(pth_to_write, fn_zero_states)
  pth_true_params <- file.path(pth_to_write, "true_params.rds")
  pth_defl_params <- file.path(pth_to_write, "defl_params.rds")
  SIMUL_U_BETA <- !is.null(get_regs_u(data_sim))
  SIMUL_Z_BETA <- !is.null(get_regs_z(data_sim))
  true_states <- get_sim_latent_states(data_sim)
  zero_states <- true_states
  zero_states[, , ] <- 0


  dist_type <- get_dist_type_simul(data_sim)

  NN <- get_dimension(data_sim, "NN") # Cross sectional length
  TT <- get_dimension(data_sim, "TT") # Time series length
  DD <- get_dimension(data_sim, "DD") # mult. comp. length
  msg_dim_ready(NN, "NN")
  msg_dim_ready(TT, "TT")
  msg_dim_ready(DD, "DD")

  tmp_list <- get_names_num_regs_sim(true_params, SIMUL_Z_BETA, SIMUL_U_BETA)
  num_regs_z <- tmp_list$num$num_regs_z
  num_regs_u <- tmp_list$num$num_regs_u
  names_z_reg <- tmp_list$names$names_z_reg
  names_u_reg <- tmp_list$names$names_u_reg

  data_out <- generate_data_sim_csv(data_sim, dist_type,
                                    SIMUL_Z_BETA, SIMUL_U_BETA,
                                    names_z_reg, names_u_reg,
                                    num_regs_z, num_regs_u,
                                    NN, TT, DD)
  utils::write.csv(data_out, file = pth_data, row.names = FALSE)
  saveRDS(true_states, file = pth_true_states)
  saveRDS(true_params, file = pth_true_params)
  saveRDS(zero_states, file = pth_zero_states)
  saveRDS(defl_params, file = pth_defl_params)
}
get_dist_type_simul <- function(data_sim) {
  DISTRIBUTION <- get_type_obs(data_sim)
  if (any(DISTRIBUTION %in% c("DIRICHLET", "GEN_DIRICHLET", "NORMAL"))) {
    dist_type <- "type1"
  } else if (any(DISTRIBUTION %in% c("DIRICHLET_MULT", "GEN_DIRICHLET_MULT",
                                     "MULTINOMIAL"))) {
    dist_type <- "type2"
  }
  return(dist_type)
}
generate_data_sim_csv <- function(data_sim, dist_type,
                                  SIMUL_Z_BETA, SIMUL_U_BETA,
                                  names_z_reg, names_u_reg,
                                  num_regs_z, num_regs_u,
                                  NN, TT, DD) {
  if (dist_type == "type1") {
    ncol_out <- DD + num_regs_z + num_regs_u
    data_out_colnames <- c(paste0("Y", 1:DD), names_z_reg, names_u_reg)

    id_col_y <- 1:DD
    if (SIMUL_Z_BETA) id_col_z <- DD + 1:(num_regs_z)
    if (SIMUL_U_BETA) id_col_u <- DD + num_regs_z + 1:(num_regs_u)
  } else if (dist_type == "type2") {
    ncol_out <- DD + 1 + num_regs_z + num_regs_u
    data_out_colnames <- c(paste0("Y", 1:DD), "num_counts",
                           names_z_reg, names_u_reg)

    id_col_y <- 1:(DD + 1)
    if (SIMUL_Z_BETA) id_col_z <- (DD + 1) + 1:(num_regs_z)
    if (SIMUL_U_BETA) {
      id_col_u <- (DD + 1) + num_regs_z + 1:(num_regs_u)
    }
  }
  data_out <- matrix(0, nrow = TT * NN, ncol = ncol_out)
  data_out <- as.data.frame(data_out)
  names(data_out) <- data_out_colnames
  for (n in 1:NN) {
    id_rows <- TT * (n - 1) + (1:TT)
    if (dist_type == "type1") {
      tmp_data <- get_data(data_sim, type = "raw")[, , n]
    } else if (dist_type == "type2") {
      tmp_data <- cbind(get_data(data_sim, type = "raw")[, , n],
                        get_data(data_sim, type = "count")[, n])
    }
    data_out[id_rows, id_col_y] <- tmp_data
    if (SIMUL_Z_BETA) {
      data_out[id_rows, id_col_z] <- get_regs_z(data_sim)[, , n]
    }
    if (SIMUL_U_BETA) {
      data_out[id_rows, id_col_u] <- get_regs_u(data_sim)[, , n]
    }
  }
  vals_cs  <- as.character(paste0("cs_", rep(seq_len(NN), each = TT)))
  vals_ts  <- rep(1:TT, times = NN)
  cs_ts    <- tibble::tibble(CS = vals_cs, TS = vals_ts)
  data_out <- dplyr::bind_cols(cs_ts, data_out)
  return(data_out)
}
copy_meta_files <- function(pth_to_project, project_name, TESTING = FALSE) {
  tmp_root <- "."
  if (isTRUE(TESTING)) tmp_root <- "./../.."
  if (grepl("NORMAL", project_name)) {
    file.copy(from = file.path(tmp_root,
                               "inst/meta-sources/main_run_T01.R"),
              to = file.path(pth_to_project,
                             paste0("main_run_", project_name, ".R")))
    file.copy(from = file.path(tmp_root,
                               "inst/meta-sources/main_diagnostics_T01.R"),
              to = file.path(pth_to_project,
                             paste0("main_diagnostics_", project_name, ".R")))
  } else {
    file.copy(from = file.path(tmp_root,
                               "inst/meta-sources/main_run_T02.R"),
              to = file.path(pth_to_project,
                             paste0("main_run_", project_name, ".R")))
    file.copy(from = file.path(tmp_root,
                               "inst/meta-sources/main_diagnostics_T02.R"),
              to = file.path(pth_to_project,
                             paste0("main_diagnostics_", project_name, ".R")))
  }
  file.copy(from = file.path(tmp_root,
                             "inst/meta-sources/setup_priors.json"),
            to = file.path(pth_to_project, "model", "model-definition"))
  file.copy(from = file.path(tmp_root,
                             "inst/meta-sources/settings_plattform&sampler.yaml"),
            to = file.path(pth_to_project, "model", "settings"))
  file.copy(from = file.path(tmp_root,
                             "inst/meta-sources/settings_project.yaml"),
            to = file.path(pth_to_project, "model", "settings"))
  copy_bash_script_slurm(tmp_root, pth_to_project, project_name)
}
copy_bash_script_slurm <- function(pth_root, pth_project, nm_project) {
  devel_sh <- readLines(
  file.path(pth_root, "inst/meta-sources/runma_DEVEL.sh"))
  devel_sh <- replace_sh(devel_sh, nm_project)
  cheops_sh <- readLines(
    file.path(pth_root, "inst/meta-sources/runma_CHEOPS.sh"))
  cheops_sh <- replace_sh(cheops_sh, nm_project)
  writeLines(devel_sh, file.path(pth_project, "runma_DEVEL.sh"))
  writeLines(cheops_sh, file.path(pth_project, "runma_CHEOPS.sh"))
}
replace_sh <- function(file, nm_project) {
  out_file <- gsub("<PATH_TO_MODEL>", nm_project, file)
  gsub("<PATH_TO_SCRIPT.R>",
       paste0("main_run_", nm_project, ".R"),
       out_file)
}