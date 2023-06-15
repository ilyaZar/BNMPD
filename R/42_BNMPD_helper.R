copy_env = function(to_env = NULL, from_env) {
  if (is.null(to_env)) to_env <- new.env()
  for(n in ls(from_env, all.names = TRUE)) {
    assign(n, get(n, from_env), to_env)
  }
  invisible(to_env)
}
copy_env_to_parent = function(from_env) {
  to_env <- parent.frame()
  for(n in ls(from_env, all.names = TRUE)) {
    assign(n, get(n, from_env), to_env)
  }
  invisible(to_env)
}
#' Helper to ease path settings for simulation studies.
#'
#' Automatically retrieves dir under BNMPD/inst/simulation studies for model
#' input settings. This includes paths to true/zero states and true/default
#' parameter settings for (P)MCMC.
#'
#' @param pth_main path to main estimation/model-class folder
#'
#' @return returns a set of 4 paths as a list
#' @export
get_paths_modelBNMPD_input <- function(pth_main) {

  pth_project <- pth_main

  pth_states <- file.path(pth_main, "model/input")
  fn_all     <- list.files(pth_states)

  fn_id_states_true <- grep("states_true", fn_all)
  pth_states_true   <- file.path(pth_states, fn_all[fn_id_states_true])

  fn_id_states_zero <- grep("states_zero", fn_all)
  pth_states_zero   <- file.path(pth_states, fn_all[fn_id_states_zero])

  pth_params_true   <- file.path(pth_main, "model/input", "true_params.rds")
  pth_params_defl   <- file.path(pth_main, "model/input", "defl_params.rds")

  out <- list()
  out$pth_project     <- pth_project
  out$pth_states_zero <- pth_states_zero
  out$pth_states_true <- pth_states_true
  out$pth_params_defl <- pth_params_defl
  out$pth_params_true <- pth_params_true
  return(out)
}
#' Helper to ease path settings and file naming for simulation studies.
#'
#' Automatically retrieves dirs under BNMPD/inst/simulation studies for model
#' results settings. This includes paths to true/zero states and true/default
#' parameter settings for (P)MCMC and plot/table file names.
#'
#' @param pth_main path to main estimation/model-class folder
#'
#' @return returns a set of 2 paths and 2 file names as characters of a list
#' @export
get_paths_modelBNMPD_results <- function(pth_main) {

  pth_table <- file.path(pth_main, "results", "inference")
  pth_plots <- file.path(pth_main, "results", "diagnostics")
  fnm_table <- basename(pth_main)
  fnm_plots <- paste0(basename(pth_main), "_all_plots")

  out <- list()
  out$pth_plots <- pth_plots
  out$pth_table <- pth_table
  out$fnm_table <- fnm_table
  out$fnm_plots <- fnm_plots
  return(out)
}
read_rds <- function(pth_from) {
  if (!is.null(pth_from)) {
    return(readRDS(normalizePath(pth_from)))
  }
}