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
#' @inheritParams new_dataSim
#'
#' @return a matrix of dimension \code{NN x DD}, with each row being the target
#'   levels (currently all the same)
#' @export
#'
#' @examples
#' get_target_dist_levels(distribution = "dirichlet", DD = 3, NN = 4)
get_target_dist_levels <- function(distribution,
                                   DD, NN,
                                   target_val_fixed = NULL,
                                   tuning_parameters = list(seq_start = 0.3,
                                                            seq_step = 0.025,
                                                            seq_rep = 2,
                                                            seq_scale = 1e4)) {
  DD  <- get_DD(distribution, DD)
  DD2 <- get_DD2(distribution,DD)
  if (is.null(target_val_fixed)) {
    seq_start <- tuning_parameters$seq_start
    seq_step  <- tuning_parameters$seq_step
    seq_rep   <- tuning_parameters$seq_rep

    tuned_vec <- rep(seq(from = seq_start,
                         to = seq_start + seq_step * DD2,
                         by = seq_step),
                     each = seq_rep)
    tuned_vec <- tuned_vec[1:DD2]
    tuned_vec <- tuned_vec/sum(tuned_vec)

    out <- matrix(tuned_vec * tuning_parameters$seq_scale,
                  nrow = DD2,
                  ncol = NN)
  } else {
    out <- matrix(rep(target_val_fixed, times = NN * DD2),
                  nrow = DD2,
                  ncol = NN)
  }
  out <- set_name_target_vals(out, DD, DD2, NN)
  return(out)
}
set_name_target_vals <- function(out, DD, DD2, NN) {
  if (DD == DD2) {
    rownames(out) <- paste0("D", seq_len(DD))
    colnames(out) <- paste0("N", seq_len(NN))
  } else if (2 * DD == DD2) {
    rownames(out) <- c(paste0("A_D", seq_len(DD)),
                       paste0("B_D", seq_len(DD)))
    colnames(out) <- paste0("N", seq_len(NN))
  } else {
    stop("Something went wrong. DD and DD2 dims do not match properly.")
  }
  return(out)
}
get_DD <- function(distribution, DD) {
  switch(distribution,
         "dirichlet" = DD,
         "dirichlet_mult" = DD,
         "gen_dirichlet" = DD,
         "gen_dirichlet_mult" = DD - 1)
}
get_DD2 <- function(distribution, DD) {
  DD <- unname(DD)
  switch(distribution,
         "dirichlet" = DD,
         "dirichlet_mult" = DD,
         "gen_dirichlet" = DD * 2,
         "gen_dirichlet_mult" = DD * 2 - 2)
}
#' Save simulated data and true parameter values used to generate it.
#'
#' @param pth_to_write path to output directory
#' @param fn_data_set file name (\code{.csv}-ending required) for simulated data
#' @param fn_true_states file name for R-container object that stores true
#'   states
#' @param fn_zero_states file name for R-container object that stores true
#'   states
#' @param data_sim simulated data set as produced via [new_dataSim()]
#' @param true_params true parameter values used to generate data (to determine
#'   container/data sizes) as produced by output from [new_trueParams()]
#' @param dim_model a vector with three components: \code{NN x TT x DD}
#'
#' @return side effect function; saves simulated data and true parameter values
#'   to output path
#'
#' @export
save_simulated_data <- function(pth_to_write,
                                fn_data_set,
                                fn_true_states,
                                fn_zero_states,
                                data_sim,
                                dim_model,
                                true_params,
                                defl_params) {
  pth_data        <- file.path(pth_to_write, fn_data_set)
  pth_true_states <- file.path(pth_to_write, fn_true_states)
  pth_zero_states <- file.path(pth_to_write, fn_zero_states)
  pth_true_params <- file.path(pth_to_write, "true_params.RData")
  pth_defl_params <- file.path(pth_to_write, "defl_params.RData")
  SIMUL_U_BETA <- !is.null(data_sim$regs$u)
  SIMUL_Z_BETA <- !is.null(data_sim$regs$z)
  true_states <- data_sim$states
  zero_states <- true_states
  zero_states[, , ] <- 0

  DISTRIBUTION <- attr(data_sim$data, "model_type_obs")
  if (any(DISTRIBUTION %in% c("DIRICHLET", "GEN-DIRICHLET", "NORMAL"))) {
    dist_type <- "type1"
  } else if(any(DISTRIBUTION %in% c("DIRICHLET-MULT", "GEN-DIRICHLET-MULT",
                                    "MULTINOMIAL"))) {
    dist_type <- "type2"
  }

  NN <- dim_model[1] # Cross sectional length
  cat(crayon::green("Setting dimension "), crayon::yellow("NN"),
      crayon::green("to "), crayon::red(NN), crayon::green("!"))
  TT <- dim_model[2]  # Time series length
  cat(crayon::green("Setting dimension "), crayon::yellow("TT"),
      crayon::green("to "), crayon::red(TT), crayon::green("!"))
  DD <- dim_model[3] # mult. comp. length
  cat(crayon::green("Setting dimension "), crayon::yellow("DD"),
      crayon::green("to "), crayon::red(DD), crayon::green("!"))

  tmp_list <- get_names_num_simulated(true_params, DD,
                                      SIMUL_Z_BETA, SIMUL_U_BETA)
  num_regs_z <- tmp_list$num$num_regs_z
  num_regs_u <- tmp_list$num$num_regs_u
  names_z_reg <- tmp_list$names$names_z_reg
  names_u_reg <- tmp_list$names$names_u_reg

  offset_col <- 2
  if (dist_type == "type1") {
    ncol_out <- DD + num_regs_z + num_regs_u
    data_out_colnames <- c(paste0("Y", 1:DD), names_z_reg, names_u_reg)

    id_col_y <- 1:DD + offset_col
    if (SIMUL_Z_BETA) id_col_z <- DD + 1:(num_regs_z) + offset_col
    if (SIMUL_U_BETA) id_col_u <- DD + num_regs_z + 1:(num_regs_u) + offset_col
  } else if (dist_type == "type2") {
    ncol_out <- DD + 1 + num_regs_z + num_regs_u
    data_out_colnames <- c(paste0("Y", 1:DD), "num_counts",
                           names_z_reg, names_u_reg)

    id_col_y <- 1:(DD + 1) + offset_col
    if (SIMUL_Z_BETA) id_col_z <- (DD + 1) + 1:(num_regs_z) + offset_col
    if (SIMUL_U_BETA) id_col_u <- (DD + 1) + num_regs_z + 1:(num_regs_u) + offset_col
  }
  data_out <- matrix(0, nrow = TT * NN, ncol = ncol_out)
  data_out <- as.data.frame(data_out)

  names(data_out) <- data_out_colnames
  vals_cs  <- as.character(paste0("cs_", rep(seq_len(NN), each = TT)))
  vals_ts  <- rep(1:TT, times = NN)
  cs_ts    <-tibble::tibble(CS = vals_cs, TS = vals_ts)
  data_out <- dplyr::bind_cols(cs_ts, data_out)

  for(n in 1:NN) {
    id_rows <- TT*(n - 1) + (1:TT)
    if (dist_type == "type1") {
      tmp_data <- data_sim$data$yraw[, , n]
    } else if(dist_type == "type2") {
      tmp_data <- cbind(data_sim$data$yraw[, , n],
                        data_sim$data$num_counts[, n])
    }
    data_out[id_rows, id_col_y] <- tmp_data
    if (SIMUL_Z_BETA) {
      data_out[id_rows, id_col_z] <- data_sim$regs$z[, , n]
    }
    if (SIMUL_U_BETA) {
      data_out[id_rows, id_col_u] <- data_sim$regs$u[, , n]
    }
  }
  write.csv(data_out, file = pth_data, row.names = FALSE)
  save(true_states, file = pth_true_states)
  save(true_params, file = pth_true_params)
  save(zero_states, file = pth_zero_states)
  save(defl_params, file = pth_defl_params)
}
get_names_num_simulated <- function(true_params, DD, SIMUL_Z, SIMUL_U) {
  if (SIMUL_Z) {
    tmp_list    <- get_name_num(type = "reg_z", true_params$beta_z_lin, DD)
    num_regs_z  <- tmp_list$num
    names_z_reg <- tmp_list$names
  } else {
    num_regs_z  <- 0
    names_z_reg <- NULL
  }
  if (SIMUL_U) {
    tmp_list    <- get_name_num(type = "reg_u", true_params$beta_u_lin, DD)
    num_regs_u  <- tmp_list$num
    names_u_reg <- tmp_list$names
  } else {
    num_regs_u  <- 0
    names_u_reg <- NULL
  }
  out <- list(names = list(names_z_reg = names_z_reg,
                           names_u_reg = names_u_reg),
              num = list(num_regs_z = num_regs_z,
                         num_regs_u = num_regs_u))
  return(out)
}
get_name_num <- function(type, par, DD) {
  if (type == "reg_z") {
    seq_regs <- lapply(par, length)
    tmp_name <- paste0(paste0(paste0("Z_", sapply(lapply(seq_regs,
                                                         function(x) {
                                                           seq_len(x)}),
                                                  as.character)), "_"),
                       rep(1:DD, unlist(seq_regs)))
  }
  if (type == "reg_u") {
    seq_regs <- lapply(par, nrow)
    tmp_name <- paste0(paste0(paste0("U_", sapply(lapply(seq_regs,
                                                         function(x) {
                                                           seq_len(x)}),
                                                  as.character)), "_"),
                       rep(1:DD, unlist(seq_regs)))
  }
  list(num = sum(unlist(seq_regs)),
       names = tmp_name)
}
#' Generates consistent file names for simulation study
#'
#' @param fn_main_part character giving base part of the name, see
#'   [get_file_name_main]
#' @param fn_data character giving the filename of the simulated data set
#'   (saved in \code{.csv}-format)
#' @param fn_true_states character giving the filename of the simulated true
#'   states (saved in \code{.RData}-format)
#' @param fn_zero_states  character giving the filename of the states set all to
#'   zero (saved in \code{.RData}-format)
#'
#' @return a list of 3:
#'    \itemize{
#'    \item\code{fn_data_set}: file name of the data set
#'    \item\code{fn_data_set}: file name for true state values
#'    \item\code{fn_data_set}: file name for zero state values
#'    }
#' @export
get_file_names_simul_data <- function(fn_main_part,
                                      fn_data,
                                      fn_true_states,
                                      fn_zero_states) {

  fn_data_set <- paste0(fn_data, "_", fn_main_part, ".csv")
  fn_true_val <- paste0(fn_true_states, "_", fn_main_part, ".RData")
  fn_zero_val <- paste0(fn_zero_states, "_", fn_main_part, ".RData")
  return(list(fn_data_set = fn_data_set,
              fn_true_val = fn_true_val,
              fn_zero_val = fn_zero_val))
}
#' Generate main part of a file name
#'
#' @param dim_model numeric vector of 3 elements: \code{NN x TT x DD} (can be
#'   taken from attributes of any object from class \code{trueParams},
#'   see [new_trueParams()] and below)
#' @param par_settings list of parameter settings as returned via attributes of
#'   any object from class \code{trueParams} (see [new_trueParams()] for
#'   details)
#' @param seed_nos integer vector with two components giving the seed number
#'   under which the trueParams object (first entry) and simulated data (second
#'   entry) is obtained
#'
#' @return main part for file names
get_file_name_main <- function(dim_model,
                               par_settings,
                               seed_nos) {
  SIMUL_PHI    <-par_settings[["SIMUL_PHI"]]
  SIMUL_Z_BETA <-par_settings[["SIMUL_Z_BETA"]]
  SIMUL_U_BETA <-par_settings[["SIMUL_U_BETA"]]
  num_z_regs   <-par_settings[["num_z_regs"]]
  num_u_regs   <-par_settings[["num_u_regs"]]
  order_p      <-par_settings[["order_p"]]
  tmp_fn <- paste0("NN", dim_model[1],
                   "_TT", dim_model[2],
                   "_DD", dim_model[3])
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
                   "_", "parSEED", seed_nos[1],
                   "_", "simSEED", seed_nos[2])
  return(tmp_fn)
}
#' Subset simulated data from [new_dataSim()]
#'
#' The output from [new_dataSim()] is a list with all data components
#' simulated, so this function returns the same structure but subsetted
#' accordingly.
#'
#' @param data the data set as given by the output from [new_dataSim()]
#' @param names_dim character with names of the dimension to retrieve; any of
#'   "TT", "DD", or "NN", or a vector containing either of these elements
#' @param id_seq for each element in \code{names_dim}, an integer or character
#'   sequence that retrieves the subset of the dimension
#'
#' @return a subset of \code{data}
#' @export
subset_data_simul <- function(data, names_dim, id_seq) {
  out_data_subset <-  data
  dim_data <- dim(data$states)

  dim_tk <- c(t = "TT", d = "DD", n = "NN")
  dim_ns <- names(dim_tk)
  dim_id <- which(dim_tk == names_dim)
  dim_nm <- length(dim_id)

  id_subset <- vector("list", 3)
  iter_id_seq <- 1
  for (i in 1:3) {
    if (i %in% dim_id) {
      id_subset[[i]] <- paste0(dim_ns[i], "_", id_seq[[iter_id_seq]])
      iter_id_seq <- iter_id_seq + 1
    } else {
      id_subset[[i]] <- paste0(dim_ns[i], "_", seq_len(dim_data[[i]]))
    }
  }
  out_data_subset$states <- data$states[id_subset[[1]],
                                        id_subset[[2]],
                                        id_subset[[3]],
                                        drop = FALSE]
  out_data_subset$data[[1]] <- data$data[[1]][id_subset[[1]],
                                              id_subset[[2]],
                                              id_subset[[3]],
                                              drop = FALSE]

  tmp_d_id <- paste0(substr(id_subset[[2]], 1, 1),
                     substr(id_subset[[2]], 3, 3))
  regexp   <- paste0("(", paste0(".*_", tmp_d_id, collapse = "|"), ")")

  tmp_z_id <- grepl(regexp,dimnames(data$regs$z)[[2]])
  out_data_subset$regs$z <- data$regs$z[id_subset[[1]],
                                        tmp_z_id,
                                        id_subset[[3]],
                                        drop = FALSE]
  tmp_u_id <- grepl(regexp, dimnames(data$regs$u)[[2]])
  out_data_subset$regs$u <- data$regs$u[id_subset[[1]],
                                        tmp_u_id,
                                        id_subset[[3]],
                                        drop = FALSE]

  return(out_data_subset)
}
#' Generate directories and files for simulation study
#'
#' @param data_simulation an object of \code{class} "dataSim" as generated by
#'    [new_dataSim()]}
#' @param INIT_AT a character: either "trues" or "default" to write
#'   initialization values at true or default values (an identity matrix for the
#'   VCM of bet_u_lin e.g. but the beta_z/beta_u set to 0).
#' @param pth_to_project character giving path to where the project (top level
#'   directory is to be created)
#' @param project_name character; name of the simulation study used as top-level
#'   directory and in project information
#'
#' @return side effect function generating directories and files for simulation
#'   study and copying data/template files into corresponding directories
#' @export
generate_simulation_study <- function(data_simulation,
                                      INIT_AT = "true",
                                      pth_to_project,
                                      project_name = list(prepend = NULL,
                                                          append = NULL)) {
  browser()
  check_class_data_sim(data_simulation)
  stopifnot(`INIT_AT must be eihter 'true' or 'default'.` =
              INIT_AT %in% c("true", "default"))
  stopifnot(`'pth_to_proj' must be of type character` =
              is.character(pth_to_project))
  stopifnot(`Arg. 'project_name' must be named list` =
              names(project_name) %in% c("prepend", "append"))

  dataSim       <- data_simulation
  trueParams    <- get_true_params_obj(dataSim)
  meta_info_tmp <- attr(trueParams, "meta_info")

  seeds_both <- c(par_seed = meta_info_tmp$SEED_NO,
                  sim_seed = attr(dataSim, which = "SEED_NO"))

  zeroParams   <- get_zero_or_defaults(trueParams)
  usedParams   <- if (INIT_AT == "true") {
    usedParams <- trueParams
  } else {
    usedParams <-zeroParams
  }

  model_type <- c(model_type_obs = attr(dataSim[["data"]],
                                        which = "model_type_obs"),
                  model_type_lat = attr(dataSim[["states"]],
                                        which = "model_type_lat"))
  base_name <- get_file_name_main(dim_model =  meta_info_tmp$MODEL_DIM,
                                  par_settings = meta_info_tmp$PAR_SETTINGS,
                                  seed_nos = seeds_both)
  tmp_name <- paste0(model_type[["model_type_obs"]],
                     "_", base_name)
  project_name <- paste0(c(project_name$prepend, tmp_name, project_name$append),
                         collapse = "_")

  pth_top_lvl <- file.path(pth_to_project, project_name)
  stopifnot(`Project directory already exists!` = !dir.exists(pth_top_lvl))
  dir.create(pth_top_lvl)

  pth_tmp <- c(file.path(pth_top_lvl, "model"),
               file.path(pth_top_lvl, "results"),
               file.path(pth_top_lvl, "model",
                         c("history",
                           "input",
                           "history/log",
                           "model-definition",
                           "output",
                           "settings")),
               file.path(pth_top_lvl, "model", "input", "datasets"),
               file.path(pth_top_lvl, "results",
                         c("diagnostics",
                           "inference",
                           "interpretation",
                           "summary")))
  lapply(pth_tmp, dir.create)
  fn_all <- get_file_names_simul_data(fn_main_part = base_name,
                                      fn_data = "sim_data",
                                      fn_true_states = "states_true",
                                      fn_zero_states = "states_zero")
  generate_yaml_model_defintion(model_type,
                                dim_model = meta_info_tmp$MODEL_DIM,
                                par_settings = meta_info_tmp$PAR_SETTINGS,
                                file.path(pth_top_lvl, "model",
                                          "model-definition",
                                          "model_definition.yaml"))
  generate_setup_init_json(usedParams, file.path(pth_top_lvl, "model",
                                                 "model-definition",
                                                 "setup_inits.json"))
  save_simulated_data(file.path(pth_top_lvl, "model", "input"),
                      file.path("datasets", fn_all[["fn_data_set"]]),
                      fn_all[["fn_true_val"]],
                      fn_all[["fn_zero_val"]],
                      data_sim = dataSim,
                      dim_model = meta_info_tmp$MODEL_DIM,
                      true_params =  trueParams,
                      defl_params = zeroParams)
  copy_meta_files(pth_top_lvl, project_name)
  update_settings_project_yaml(file.path(pth_to_project, project_name, "model",
                                         "settings", "settings_project.yaml"),
                               proj_no = NULL, proj_name = project_name,
                               notes = NULL)
  return(invisible(data_simulation))
}
get_zero_or_defaults <- function(true_params) {
  check_class_true_params(true_params)
  zero_params <- true_params
  if (!is.null(true_params[["sig_sq"]])) {
    zero_params[["sig_sq"]][] <- 1
  }
  if (!is.null(true_params[["phi"]])) {
    zero_params[["phi"]][] <- lapply(zero_params[["phi"]],
                                     function(x) {y <- x; y[] <- 0; y})
  }
  if (!is.null(true_params[["beta_z_lin"]])) {
    zero_params[["beta_z_lin"]][] <- lapply(zero_params[["beta_z_lin"]],
                                            function(x) {y <- x; y[] <- 0; y})
  }
  if (!is.null(true_params[["beta_u_lin"]])) {
    zero_params[["beta_u_lin"]][] <- lapply(zero_params[["beta_u_lin"]],
                                            function(x) {y <- x; y[] <- 0; y})
  }
  if (!is.null(true_params[["vcm_u_lin"]])) {
    zero_params[["vcm_u_lin"]][] <- lapply(zero_params[["vcm_u_lin"]],
                                           function(x) {
                                             y <- x;
                                             y[] <- 0;
                                             diag(y) <- 1;
                                             y}
    )
  }

  return(zero_params)
}
copy_meta_files <- function(pth_to_project, project_name) {
  if (grepl("NORMAL", project_name)) {
    file.copy(from = "./inst/meta-sources/main_run_T01.R",
              to = file.path(pth_to_project,
                             paste0("main_run_", project_name, ".R")))
  } else {
    file.copy(from = "./inst/meta-sources/main_run_T02.R",
              to = file.path(pth_to_project,
                             paste0("main_run_", project_name, ".R")))
  }
  file.copy(from = "./inst/meta-sources/setup_priors.json",
            to = file.path(pth_to_project, "model", "model-definition"))
  file.copy(from = "./inst/meta-sources/settings_plattform&sampler.yaml",
            to = file.path(pth_to_project, "model", "settings"))
  file.copy(from = "./inst/meta-sources/settings_project.yaml",
            to = file.path(pth_to_project, "model", "settings"))
}
