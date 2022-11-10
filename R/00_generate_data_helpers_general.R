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
#'
#' @return a matrix of dimension \code{NN x DD}, with each row being the target
#'   levels (currently all the same)
#' @export
#'
#' @examples
#' get_dirichlet_levels(DD = 3, NN = 4)
get_dirichlet_levels <- function(DD, NN,
                                 tuning_parameters = list(seq_start = 0.3,
                                                          seq_step = 0.1,
                                                          seq_rep = 2,
                                                          seq_scale = 100)) {
  seq_start <- tuning_parameters$seq_start
  seq_step  <- tuning_parameters$seq_step
  seq_rep   <- tuning_parameters$seq_rep

  tuned_vec <- rep(seq(from = seq_start,
                       to = seq_start + seq_step * DD,
                       by = seq_step),
                   each = seq_rep)
  tuned_vec <- tuned_vec[1:DD]
  tuned_vec <- tuned_vec/sum(tuned_vec)

  matrix(tuned_vec * 100,
         nrow = DD,
         ncol = NN)

}
#' Save simulated data and true parameter values used to generate it.
#'
#' @param pth_to_write path to output directory
#' @param fn_data_set file name (\code{.csv}-ending required) for simulated data
#' @param fn_true_states file name for R-container object that stores true
#'   states
#' @param fn_zero_states file name for R-container object that stores true
#'   states
#' @param data_sim simulated data set as produced via [generate_data_t_n()]
#' @param true_params true parameter values used to generate data (to determine
#'   container/data sizes) as produced by output from [generate_true_params()]
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
                                true_params) {
  SIMUL_U_BETA <- !is.null(data_sim$regs$u)
  NN <- dim_model[1] # Cross sectional length
  cat(crayon::green("Setting dimension "), crayon::yellow("NN"),
      crayon::green("to "), crayon::red(NN), crayon::green("!"))
  TT <- dim_model[2]  # Time series length
  cat(crayon::green("Setting dimension "), crayon::yellow("TT"),
      crayon::green("to "), crayon::red(TT), crayon::green("!"))
  DD <- dim_model[3] # mult. comp. length
  cat(crayon::green("Setting dimension "), crayon::yellow("DD"),
      crayon::green("to "), crayon::red(DD), crayon::green("!"))

  seq_regs_z <- lapply(true_params$beta_z_lin, length)
  num_regs_z <- sum(unlist(seq_regs_z))
  seq_regs_u <- lapply(true_params$beta_u_lin, nrow)
  num_regs_u <- sum(unlist(seq_regs_u))

  ncol_out <- DD + num_regs_z + num_regs_u
  data_out <- matrix(0, nrow = TT * NN, ncol = ncol_out)
  data_out <- as.data.frame(data_out)

  names_z_reg <- paste0(paste0(paste0("Z_", sapply(lapply(seq_regs_z,
                                                          function(x) {
                                                            seq_len(x)}),
                                                   as.character)), "_"),
                        rep(1:DD, unlist(seq_regs_z)))
  if (SIMUL_U_BETA) {
    names_u_reg <- paste0(paste0(paste0("U_", sapply(lapply(seq_regs_u,
                                                            function(x) {
                                                              seq_len(x)}),
                                                     as.character)), "_"),
                          rep(1:DD, unlist(seq_regs_u)))
  } else {
    names_u_reg <- NULL
  }
  names(data_out) <- c(paste0("Y", 1:DD), names_z_reg, names_u_reg)
  vals_cs  <- as.character(paste0("cs_", rep(seq_len(NN), each = TT)))
  vals_ts  <- rep(1:TT, times = NN)
  cs_ts    <-tibble::tibble(CS = vals_cs, TS = vals_ts)
  data_out <- dplyr::bind_cols(cs_ts, data_out)
  for(n in 1:NN) {
    id_rows <- TT*(n - 1) + (1:TT)
    offset_col <- 2
    id_col_y <- 1:DD + offset_col
    id_col_z <- DD + 1:(num_regs_z) + offset_col
    id_col_u <- DD + num_regs_z + 1:(num_regs_u) + offset_col
    data_out[id_rows, id_col_y] <- data_sim$data$yraw[, , n]
    data_out[id_rows, id_col_z] <- data_sim$regs$z[, , n]
    if (SIMUL_U_BETA) {
      data_out[id_rows, id_col_u] <- data_sim$regs$u[, , n]
    }
  }
  write.csv(data_out, file = file.path(pth_to_write, fn_data_set),
            row.names = FALSE)
  true_states <- data_sim$states
  save(true_states, file = file.path(pth_to_write, fn_true_states))
  zero_states <- true_states
  zero_states[, , ] <- 0
  save(zero_states,  file = file.path(pth_to_write, fn_zero_states))
}
#' Generates consistent file names for simulation study
#'
#' @param fn_data character giving the filename of the simulated data set
#'   (saved in \code{.csv}-format)
#' @param fn_true_states character giving the filename of the simulated true
#'   states (saved in \code{.RData}-format)
#' @param fn_zero_states  character giving the filename of the states set all to
#'   zero (saved in \code{.RData}-format)
#' @param dim_model numeric vector of 3 elements: \code{NN x TT x DD}
#' @param SIMUL_Z_BETA logical; if \code{TRUE} Z-type regressors are simulated
#'   which is reflected in the naming of simulated data sets
#' @param SIMUL_U_BETA logical; if \code{TRUE} U-type regressors (i.e. random
#'   effects) are simulated which is reflected in the naming of simulated data
#'   sets
#' @param num_z_regs numeric vector giving number of Z-type regressors; parsed
#'   as comma-seperated sequence of numbers i.e. for \code{num_z_regs = 1:3} we
#'   get "withLIN,1,2,3" in the file names
#' @param num_u_regs numeric vector giving number of U-type regressors; parsed
#'   as comma-seperated sequence of numbers i.e. for \code{num_u_regs = 1:3} we
#'   get "withRE,1,2,3" in the file names
#'
#' @return a list of 3:
#'    \itemize{
#'    \item\code{fn_data_set}: file name of the data set
#'    \item\code{fn_data_set}: file name for true state values
#'    \item\code{fn_data_set}: file name for zero state values
#'    }
#' @export
get_file_names_simul_data <- function(fn_data,
                                      fn_true_states,
                                      fn_zero_states,
                                      dim_model,
                                      SIMUL_Z_BETA,
                                      SIMUL_U_BETA,
                                      num_z_regs,
                                      num_u_regs) {
  tmp_fn <- paste0("NN", dim_model[1],
                   "_TT", dim_model[2],
                   "_DD", dim_model[3])
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
  fn_data_set <- paste0(fn_data, "_", tmp_fn, ".csv")
  fn_true_val <- paste0(fn_true_states, "_", tmp_fn, ".RData")
  fn_zero_val <- paste0(fn_zero_states, "_", tmp_fn, ".RData")
  return(list(fn_data_set = fn_data_set,
              fn_true_val = fn_true_val,
              fn_zero_val = fn_zero_val))
}
#' Subset simulated data from [generate_data_t_n()]
#'
#' The output from [generate_data_t_n()] is a list with all data components
#' simulated, so this function returns the same structure but subseted
#' accordingly.
#'
#' @param data the data set as given by the output from [generate_data_t_n()]
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
#' @param pth_to_project character giving path to where the project (top level
#'   directory is to be created)
#'
#' @param project_name character; name of the simulation study used as top-level
#'   directory and in project information
#' @param meta_info either \code{prepend} or \code{append}, defaults to the
#'   former and appends \code{project_name} by measurement type and model
#'   dimensions
#' @param data_simulation a named list of three elements:
#'   \itemize{
#'   \item{\code{data:}}{ the data set as generated by [generate_data_t_n()]}
#'   \item{\code{trueParams:}}{ true parameter values as generated by
#'                              [generate_true_params()]}
#'   \item{\code{dim:}}{ a named numeric vector ("NN", "TT", "DD")
#'   giving the dimensions of the model}
#'   }
#' @param data_simulation_meta a named list of four elements:
#'   \itemize{
#'   \item{\code{SIMUL_Z_BETA:}}{ logical; if \code{TRUE} z-regressors have been
#'   simulated}
#'   \item{\code{SIMUL_U_BETA:}}{ logical; if \code{TRUE} u-regressors have been
#'   simulated}
#'   \item{\code{num_z_regs:}}{ number of z-type regressors or \code{NULL}}
#'   \item{\code{num_z_regs:}}{ number of u-type regressors or \code{NULL}}
#'   }
#'
#' @return side effect function generating directories and files for simulation
#'   study and copying data/template files into corresponding directories
#' @export
generate_simulation_study <- function(pth_to_project,
                                      project_name,
                                      meta_info = "prepend",
                                      data_simulation,
                                      data_simulation_meta) {
  model_type <- c(model_type_obs = attr(data_simulation$dataSim[["data"]],
                                        which = "model_type_obs"),
                  model_type_lat = attr(data_simulation$dataSim[["states"]],
                                        which = "model_type_lat"))
  tmp_name <- paste0(paste0(names(data_simulation$dim),
                            data_simulation$dim), collapse = "_")
  tmp_name <- paste0(model_type[["model_type_obs"]],
                     "_", tmp_name, "_",
                     model_type[["model_type_lat"]])
  if (meta_info == "prepend") {
    project_name <- paste0(tmp_name, "_", project_name)
  } else if  (meta_info == "append") {
    project_name <- paste0(project_name, "_", tmp_name)
  } else {
    stop("Unknown value for arg 'meta_info'; use either 'prepend' or 'append'.")
  }
  pth_top_lvl <- file.path(pth_to_project, project_name)
  stopifnot(`Project directory already exists!` = !dir.exists(pth_top_lvl))
  dir.create(pth_top_lvl)

  pth_tmp <- c(file.path(pth_top_lvl, "model"),
               file.path(pth_top_lvl, "results"),
               file.path(pth_top_lvl, "model",
                         c("history",
                           "input",
                           "log",
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

  fn_all <- get_file_names_simul_data(fn_data = "sim_data",
                                      fn_true_states = "states_true",
                                      fn_zero_states = "states_zero",
                                      dim_model = data_simulation[["dim"]],
                                      data_simulation_meta[["SIMUL_Z_BETA"]],
                                      data_simulation_meta[["SIMUL_U_BETA"]],
                                      data_simulation_meta[["num_z_regs"]],
                                      data_simulation_meta[["num_u_regs"]])
  generate_yaml_model_defintion(model_type,
                                dimensions = data_simulation[["dim"]],
                                data_simulation_meta,
                                file.path(pth_top_lvl, "model",
                                          "model-definition",
                                          "model_definition.yaml"))
  generate_setup_init_json(data_simulation[["dim"]],
                           data_simulation[["trueParams"]],
                           file.path(pth_top_lvl, "model",
                                     "model-definition",
                                     "setup_inits.json"))
  save_simulated_data(file.path(pth_top_lvl, "model", "input"),
                      file.path("datasets", fn_all[["fn_data_set"]]),
                      fn_all[["fn_true_val"]],
                      fn_all[["fn_zero_val"]],
                      data_sim = data_simulation[["dataSim"]],
                      dim_model = data_simulation[["dim"]],
                      true_params = data_simulation[["trueParams"]])
  copy_meta_files(pth_top_lvl)
  return(invisible(pth_to_project))
}
copy_meta_files <- function(pth_to_project) {
  file.copy(from = "./inst/meta-sources/setup_priors.json",
            to = file.path(pth_to_project, "model", "model-definition"))
  file.copy(from = "./inst/meta-sources/settings_plattform&sampler.yaml",
            to = file.path(pth_to_project, "model", "settings"))
  file.copy(from = "./inst/meta-sources/settings_project.yaml",
            to = file.path(pth_to_project, "model", "settings"))
}
