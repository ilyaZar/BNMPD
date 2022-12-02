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
#' Helper to ease path settings especially for simulation studies.
#'
#' Automatically retrieves dir under BNMPD/inst/simulation studies.
#'
# @param type type of study; works well with 'simulation-study' so
#   so far.
#' @param pth_main path to main estimation/model-class folder
#'
#' @return returns a set of 4 paths as a list
#' @export
get_paths_modelBNMPD <- function(pth_main) {

  pth_project <- pth_main

  pth_states <- file.path(pth_main, "model/input")
  fn_all     <- list.files(pth_states)

  fn_id_states_true <- grep("states_true", fn_all)
  pth_states_true   <- file.path(pth_states, fn_all[fn_id_states_true])

  fn_id_states_zero <- grep("states_zero", fn_all)
  pth_states_zero   <- file.path(pth_states, fn_all[fn_id_states_zero])

  pth_params_true   <- file.path(pth_main, "model/input", "true_params.RData")

  out <- list()
  out$pth_project     <- pth_project
  out$pth_states_zero <- pth_states_zero
  out$pth_states_true <- pth_states_true
  out$pth_params      <- pth_params_true
  return(out)
}
